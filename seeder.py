from __future__ import annotations
import copy
import inspect
from . import utilities as util
from . import harvester as harvest
from pathlib import Path
import subprocess as sp
import shutil  # for copy--will be obviated by python 3.14
import re
import json


default_run_script = """
from mdfarmer.simulate import omm_basic_sim_block_json as runner
runner('config.json')
"""


class ConfigError(ValueError):
    """A mistake in the Farmer's configuration rather than in one clone's disk
    state.

    Clone setup otherwise swallows failures so one corrupt generation directory
    cannot stop an orchestrator minding hundreds of others. A config error is
    not per-clone: it applies to every clone it touches, so it propagates
    instead, rather than booting a campaign quietly short of the seeds asked
    for.
    """


# Config entries from_disk works out per clone. A seed override naming one would
# be silently overwritten, so it is refused instead.
CLONE_DERIVED_CONFIG_KEYS = frozenset((
    'seed_index', 'clone_index', 'gen_index', 'steps', 'steps_per_gen',
    'append', 'new_velocities', 'seed_fn', 'top_fn', 'system_fn'))


def _gen_sort_key(path, sep):
    """Sort key ordering gen directories numerically.

    Anything that is not '<prefix><sep><digits>' sorts to the end, keyed by
    name so the order is still deterministic.
    """
    tail = path.name.rsplit(sep, 1)[-1]
    if tail.isdigit():
        return (0, int(tail), '')
    return (1, 0, path.name)


def _align_tandem_trajs(gen_path: Path, prev_config: dict, target_nset: int):
    """Trim the velocity and force trajectories to target_nset frames.

    True once every tandem trajectory the config names is aligned. False when
    one cannot be: a format that cannot be trimmed, a file that cannot be read,
    or a tandem behind positions, whose missing frames cannot be fabricated
    without a matching checkpoint. The caller cascades on False and redoes the
    generation, which is dead either way.
    """
    for name_key, suffix_key, default_name in (
            ('velocity_name', 'velocity_traj_suffix', 'velocities'),
            ('force_name', 'force_traj_suffix', 'forces')):
        tandem_suffix = prev_config.get(suffix_key)
        if not tandem_suffix:
            continue
        tandem_p = (gen_path / prev_config.get(name_key, default_name)
                    ).with_suffix(tandem_suffix)
        if (not tandem_p.is_file()) or tandem_p.stat().st_size == 0:
            continue
        # Trimmed from what the bytes hold, not from a header that can be
        # ahead of them, so a torn last frame is removed rather than counted.
        try:
            tandem_actual = util.truncate_traj_to_nframes(tandem_p, target_nset)
        except Exception as exc:
            print(f'_align_tandem_trajs: could not trim {tandem_p}: {exc}; '
                  f'cascading.')
            return False
        if tandem_actual != target_nset:
            print(f'_align_tandem_trajs: {tandem_p} reached {tandem_actual} '
                  f'frames, not the {target_nset} positions has; cascading.')
            return False
    return True


def _try_recover_gen(gen_path: Path, *,
                     append_mode: bool,
                     restart_name: str,
                     traj_name: str,
                     traj_suffix: str,
                     write_interval: int,
                     total_steps: int,
                     top_fn: str):
    """Decide whether one gen directory can seed the next launch.

    Returns (gen_index, seed_fn, steps_to_run, append): the gen the launch is
    for, which is one higher than this directory's if this gen is complete;
    the state.xml to load; the steps to take; and whether the trajectory
    reporters open in append mode. Returns None if the gen is unrecoverable,
    and the caller cascades to an older gen, then to the initial seed.

    Unrecoverable means prune and redo: a state.xml that is missing or will
    not parse, or a stepCount that does not sit on a trajectory frame
    boundary, where header and checkpoint disagree and both are suspect. A
    state behind the trajectory is the expected drift from a kill between the
    two reporters, so the trajectory is trimmed to match. A state ahead of it
    loses the trajectory tail but keeps an intact integrator state, so it
    advances to the next gen and seeds from there.

    The step at gen start is counted from the earlier generations' own configs,
    not read from the DCD header: OpenMM's DCDReporter hard-codes istart to
    reportInterval, while state.xml's stepCount accumulates across gens.
    """
    config_p = gen_path / 'config.json'
    if not config_p.is_file():
        return None
    try:
        prev_config = json.loads(config_p.read_text())
    except json.JSONDecodeError:
        print(f'_try_recover_gen: malformed config at {config_p}; skipping gen.')
        return None
    gen_index = prev_config['gen_index']

    restart_p = gen_path / restart_name
    if not util.is_state_xml_usable(restart_p):
        return None
    seed_fn = str(restart_p.resolve())
    try:
        state_step = util.state_xml_step_count(restart_p)
    except ValueError as exc:
        print(f'_try_recover_gen: {exc}; skipping gen.')
        return None

    traj_p = (gen_path / traj_name).with_suffix(traj_suffix)
    if (not traj_p.is_file()) or traj_p.stat().st_size == 0:
        # No traj yet for this gen. Start it now from this state.
        return gen_index, seed_fn, total_steps, False

    # Only DCD has a header we can read and rewrite; an XTC has to be counted
    # instead, and one that ran past its checkpoint is redone rather than cut.
    try:
        nset = util.get_traj_len(str(traj_p), top_fn)
    except Exception as exc:
        print(f'_try_recover_gen: cannot count frames in {traj_p}: '
              f'{type(exc).__name__}: {exc}; skipping gen.')
        return None
    nsavc = write_interval

    if nset == 0:
        return gen_index, seed_fn, total_steps, False

    # Counted from the campaign's own zero: the seed may have arrived already
    # carrying the step count of whatever equilibrated it.
    gen_start_step = prev_config.get('step_origin', 0) + util.steps_before(
        prev_config['traj_dir_top_level'], prev_config['seed_index'],
        prev_config['clone_index'], gen_index, prev_config['dirname_pad'],
        sep=prev_config['sep'])
    state_offset = state_step - gen_start_step
    if state_offset <= 0 or state_offset % nsavc != 0:
        print(f'_try_recover_gen: state stepCount={state_step} not aligned '
              f'to gen start {gen_start_step} at nsavc={nsavc} ({gen_path}); '
              f'cascading.')
        return None
    target_nset = state_offset // nsavc

    if target_nset > nset:
        # State ahead of the trajectory: the integrator state is intact, advance.
        return gen_index + 1, seed_fn, total_steps, False
    if target_nset < nset:
        # Kill between the trajectory reporter and the checkpoint writer.
        print(f'_try_recover_gen: trimming {traj_p} from {nset} to '
              f'{target_nset} frames to match state.xml stepCount.')
        try:
            actual = util.truncate_traj_to_nframes(traj_p, target_nset)
        except Exception as exc:
            print(f'_try_recover_gen: could not truncate {traj_p}: {exc}; '
                  f'cascading.')
            return None
        if actual != target_nset:
            print(f'_try_recover_gen: truncate returned {actual} != '
                  f'target {target_nset}; cascading.')
            return None
        nset = actual
        # Reached only once positions were trimmed; a healthy in-flight gen
        # takes the target_nset == nset path, so its tandems are left alone.
        if not _align_tandem_trajs(gen_path, prev_config, target_nset):
            return None

    remaining = total_steps - nset * nsavc
    if remaining > 0:
        if append_mode:
            return gen_index, seed_fn, remaining, True
        # Non-append: preserve the partial traj and try an older gen.
        old_traj_p = traj_p.parent / ('old_' + traj_p.name)
        print(f'_try_recover_gen: non-append mode, renaming partial '
              f'{traj_p} -> {old_traj_p}')
        traj_p.rename(old_traj_p)
        return None
    # Gen complete; advance.
    return gen_index + 1, seed_fn, total_steps, False


class Clone:
    """One trajectory, grown a generation at a time by repeated scheduler jobs.

    Built by a Farmer at boot, and directly by adaptive sampling scripts that
    steer a clone themselves.
    """

    __slots__ = (
        'config', 'dry_run', 'job_number', 'job_number_re', 'job_name_fstring', 'current_seed',
        'current_gen_dir', 'config_p', 'scheduler_script_p', 'compare_keys',
        'scheduler_fstring', 'scheduler', 'traj_list', 'sep', 'dirname_pad',
        'scheduler_kws', 'restarts_per_gen', 'restart_attempts', 'run_script',
        'harvester', 'remaining_steps', 'run_script_name', 'total_steps',
        'preemption_checker', 'node_blocklist', 'progress_fn',
        'scheduler_log_dir', 'last_gen_index', 'reaped_gen')

    def __init__(self,
                 # keys must match the runner function's parameter names.
                 config: dict,

                 # Scheduler to call. Will also be used to name files later.
                 scheduler: str,
                 # BSUB/SBATCH/Scheduler script with anchors for str.format().
                 scheduler_fstring: str,
                 # keys match the fstring's anchors, values the substitutions.
                 scheduler_kws: dict,
                 # path to the serialized xml this clone will grow from.
                 seed_fn: str,
                 # if true, run inspect.cleandoc on scheduler_fstring first.
                 cleandoc_sched_fstring=True,
                 restarts_per_gen=3,
                 # the scheduler's assigned job number for the current job.
                 job_number=None,
                 dirname_pad=2,
                 sep='-',
                 run_script=default_run_script,
                 # regex to extract job_number from submission call output.
                 job_number_re='[1-9][0-9]*',
                 # anchors for config keys, rendering the job name each gen.
                 job_name_fstring=None,
                 # Joined into job_name_fstring when that is None. Seed, clone
                 # and gen must appear in that order or reassociation fails.
                 job_name_elements=(
                     '{title}', '{seed_index}', '{clone_index}',
                     '{gen_index}'),
                 # config keys __eq__ and __hash__ compare two clones on.
                 compare_keys=('seed_index', 'clone_index'),
                 harvester=None,
                 # Callable jid -> bool; a preemption it reports costs no
                 # restart_attempt.
                 preemption_checker=None,
                 # Shared BadNodeRegistry. A 0-step abort whose scheduler log
                 # names a node-local failure bars that node from later jobs.
                 node_blocklist=None,
                 # Full per-gen step count. On a resume config['steps'] holds
                 # only the remainder, which would otherwise shorten later gens.
                 steps_per_gen=None,
                 # Callable(gen_dir, **context) -> steps a gen still owes. None
                 # counts frames (OpenMM); GROMACS passes gmx_gen_progress.
                 progress_fn=None,
                 # 0-based index of the last generation; None means no limit.
                 last_gen_index=None,
                 dry_run=False
                 ):
        # REQUIRED ARGS below here
        self.config = config  # dict keys and values must be json serializable.
        # Before set_seed, whose error paths report self.get_tag().
        self.compare_keys = compare_keys
        if steps_per_gen is None:
            self.total_steps = config['steps']
        else:
            self.total_steps = steps_per_gen
        self.remaining_steps = config['steps']
        seed_p = Path(seed_fn)
        # An equilibrated seed's state.xml already carries a step count, and
        # every later one counts on from it. Recorded once so recovery can tell
        # this campaign's own steps from the seed's history.
        config.setdefault('step_origin', util.state_xml_origin(seed_p))
        self.set_seed(seed_p)
        self.scheduler = scheduler
        if cleandoc_sched_fstring:
            self.scheduler_fstring = inspect.cleandoc(scheduler_fstring)
        else:
            self.scheduler_fstring = scheduler_fstring
        self.scheduler_kws = scheduler_kws
        try:
            self.run_script_name = self.scheduler_kws['run_script_name']
        except KeyError:
            self.run_script_name = 'run.py'
            scheduler_kws['run_script_name'] = self.run_script_name

        # Args with defaults below here
        self.dry_run = dry_run
        self.restarts_per_gen = restarts_per_gen
        # this should always start at zero, since it's incremented below.
        self.restart_attempts = 0
        self.current_gen_dir = util.dir_seeds_clones_gens(
            Path(self.config['traj_dir_top_level']),
            self.config['seed_index'],
            self.config['clone_index'],
            self.config['gen_index'],
            self.config['dirname_pad'],
            sep=self.config['sep'],
            mkdir=True
        )
        self.job_number = job_number
        self.sep = sep
        if job_name_fstring:
            self.job_name_fstring = job_name_fstring
        else:
            self.job_name_fstring = self.sep.join(job_name_elements)
        self.dirname_pad = dirname_pad
        self.job_number_re = re.compile(job_number_re)
        self.harvester = harvester
        # The harvest reads the full generation length from config.json, and
        # config['steps'] holds only the steps still owed on a resume.
        if config.get('steps_per_gen') is None:
            config['steps_per_gen'] = self.total_steps
        elif config['steps_per_gen'] != self.total_steps:
            raise ConfigError(
                f"config['steps_per_gen']={config['steps_per_gen']} disagrees "
                f'with steps_per_gen={self.total_steps}; the harvest would '
                'place this clone\'s frames at the wrong global index.')
        # Catch an incommensurable generation length now, while it is free to
        # fix: it breaks harvested trajectories at every seam, and never loudly.
        if harvester is not None:
            downsample_frq = (getattr(harvester, 'run_config', None)
                              or {}).get('downsample_frq')
            if downsample_frq:
                try:
                    harvest.check_commensurability(
                        self.total_steps, config['write_interval'],
                        downsample_frq)
                except ValueError as exc:
                    raise ConfigError(str(exc)) from exc
        self.preemption_checker = preemption_checker
        self.node_blocklist = node_blocklist
        # The generation already harvested, so it never gets a second harvest.
        self.reaped_gen = None
        self.progress_fn = progress_fn
        self.last_gen_index = last_gen_index
        self.run_script = run_script
        # Where this clone's scheduler log lands. None means its own gen dir; a
        # ClonePack points every member at the one log a packed job writes.
        self.scheduler_log_dir = None
        self.scheduler_script_p = None  # always redefined each run

    @classmethod
    def from_disk(cls,
                  tdir: Path,
                  seed_index: int,
                  clone_index: int,
                  *,
                  # This seed's starting structure, used when no gen recovers.
                  initial_seed_fn: str,
                  # Resolved, validated top and system paths for this seed.
                  top_fn: str,
                  system_fn: str,
                  # The .gro or .pdb this seed's generation 0 starts from. Per
                  # seed, since one shared value can only suit one of them.
                  structure_fn: str = None,
                  # Config entries laid over the template for this seed. May
                  # not name a key from_disk works out per clone.
                  config_overrides: dict = None,
                  # The Farmer's config_template. Deepcopied before mutation.
                  config_template: dict,
                  # Scheduler context (passed straight through to __init__).
                  scheduler: str,
                  scheduler_fstring: str,
                  scheduler_kws: dict,
                  # Discovery / Clone wiring.
                  dirname_pad: int,
                  sep: str,
                  job_number_re: str,
                  job_name_fstring: str,
                  harvester=None,
                  preemption_checker=None,
                  node_blocklist=None,
                  restarts_per_gen=3,
                  # Passed through to __init__; None means no generation limit.
                  last_gen_index=None,
                  # GROMACS hooks. None gives the OpenMM defaults, which are
                  # _try_recover_gen and default_run_script.
                  recover_fn=None,
                  run_script=None,
                  progress_fn=None,
                  # (seed, clone, gen) -> jid for jobs currently queued, so they
                  # can be re-associated after an orchestrator restart.
                  rep_dict: dict = None,
                  dry_run: bool = False):
        """Build a Clone from the on-disk state for one seed and clone index.

        Decides which gen to run next, what to seed it from, and whether to
        append or start fresh, by falling through the extant gens newest to
        oldest and then back to initial_seed_fn if none is recoverable.
        """
        if rep_dict is None:
            rep_dict = {}
        append_mode = config_template['append']
        steps_per_gen = config_template['steps']

        clone_dir = util.dir_seeds_clones(
            tdir, seed_index, clone_index, dirname_pad, sep=sep, mkdir=False)
        if clone_dir.is_dir():
            # By generation number, not by name: names only sort right at equal
            # width, so gen-999 would come after gen-1000. Unparseable go last.
            gen_paths = sorted(clone_dir.iterdir(),
                               key=lambda p: _gen_sort_key(p, sep))
        else:
            gen_paths = []

        # A gen whose job is still running is never recovered: recovery trims
        # and renames files that job holds open. Bind to the live job instead.
        live_gens = sorted(gen for six, cix, gen in rep_dict
                           if (six, cix) == (seed_index, clone_index))
        _recover = recover_fn if recover_fn is not None else _try_recover_gen
        recovered = None
        if not live_gens:
            for gen_path in reversed(gen_paths):
                recovered = _recover(
                    gen_path,
                    append_mode=append_mode,
                    restart_name=config_template['restart_name'],
                    traj_name=config_template['traj_name'],
                    traj_suffix=config_template['traj_suffix'],
                    write_interval=config_template['write_interval'],
                    total_steps=steps_per_gen,
                    top_fn=top_fn,
                )
                if recovered is not None:
                    break

        if live_gens:
            gen_index = live_gens[-1]
            live_dir = util.dir_seeds_clones_gens(
                tdir, seed_index, clone_index, gen_index, dirname_pad,
                sep=sep, mkdir=False)
            restart_p = live_dir / config_template['restart_name']
            seed_fn = (str(restart_p.resolve())
                       if restart_p.is_file() and restart_p.stat().st_size
                       else initial_seed_fn)
            steps_this_launch = steps_per_gen
            append_now = True
            is_internal_restart = True
            print(f'seed {seed_index} clone {clone_index} gen {gen_index} is '
                  'still running; leaving its files alone.')
        elif recovered is not None:
            gen_index, seed_fn, steps_this_launch, append_now = recovered
            is_internal_restart = True
        else:
            gen_index = 0
            seed_fn = initial_seed_fn
            steps_this_launch = steps_per_gen
            append_now = False
            is_internal_restart = False

        config = copy.deepcopy(config_template)
        config['seed_index'] = seed_index
        config['clone_index'] = clone_index
        config['gen_index'] = gen_index
        config['steps'] = steps_this_launch
        config['append'] = append_now
        if is_internal_restart:
            config['new_velocities'] = False
        else:
            config['new_velocities'] = True
        config['system_fn'] = system_fn
        config['top_fn'] = top_fn
        if structure_fn is not None:
            config['structure_fn'] = structure_fn
        if config_overrides:
            clashing = sorted(set(config_overrides) & CLONE_DERIVED_CONFIG_KEYS)
            if clashing:
                raise ConfigError(
                    f'config_overrides may not set {clashing}: from_disk '
                    'derives those per clone and would overwrite them.')
            config.update(config_overrides)
        # Keeps the StateDataReporter's % complete honest on a resume.
        if isinstance(config.get('state_data_kwargs'), dict):
            config['state_data_kwargs'] = {
                **config['state_data_kwargs'],
                'totalSteps': steps_this_launch,
            }

        jid = rep_dict.get((seed_index, clone_index, gen_index))

        # The full generation length: the GROMACS runner needs it for an
        # absolute step target even when config['steps'] is only a remainder.
        config['steps_per_gen'] = steps_per_gen

        # OpenMM callers pass nothing and get the OpenMM runner.
        if run_script is None:
            run_script = default_run_script
        return cls(
            config,
            scheduler,
            scheduler_fstring,
            scheduler_kws,
            seed_fn,
            steps_per_gen=steps_per_gen,
            job_number=jid,
            job_number_re=job_number_re,
            job_name_fstring=job_name_fstring,
            restarts_per_gen=restarts_per_gen,
            dirname_pad=dirname_pad,
            sep=sep,
            harvester=harvester,
            preemption_checker=preemption_checker,
            node_blocklist=node_blocklist,
            progress_fn=progress_fn,
            last_gen_index=last_gen_index,
            dry_run=dry_run,
            run_script=run_script,
        )

    def __hash__(self):
        """Hash the config entries named by compare_keys."""
        return hash(tuple(self.config[k] for k in self.compare_keys))

    def __eq__(self, other: Clone) -> bool:
        """Equal when the hashes match.

        Two clones compared on different compare_keys will nearly always
        differ, whatever else they have in common.
        """
        return hash(self) == hash(other)

    def get_tag(self):
        """Return a string of the config entries that identify this clone."""
        return ' '.join((f'{k}: {self.config[k]}'
                         for k in self.compare_keys))

    def set_seed(self, seed_fp):
        if seed_fp.is_file():
            if seed_fp.stat().st_size == 0:
                print('ERROR: CLONE', self.get_tag(),
                      'found seed, but it is an empty file.')
                raise IOError(seed_fp)
            self.current_seed = seed_fp
            self.config['seed_fn'] = str(seed_fp)
        else:
            print('ERROR: CLONE', self.get_tag(), 'could not find seed.')
            raise FileNotFoundError(seed_fp)

    def holds_own_progress(self, restart_p):
        """True if this restart file is newer than the seed we would copy in."""
        return (restart_p.is_file() and restart_p.stat().st_size
                and restart_p.stat().st_mtime
                > self.current_seed.stat().st_mtime)

    def target_step(self):
        """The absolute step this generation ends at, counted from the chain."""
        return util.steps_before(
            self.config['traj_dir_top_level'], self.config['seed_index'],
            self.config['clone_index'], self.config['gen_index'],
            self.config['dirname_pad'], sep=self.config['sep']
        ) + self.total_steps

    def check_copy_set_restart_seed(self):
        """Put a copy of this clone's seed in the current gen directory.

        A restart file newer than the seed is this generation's own progress,
        so it is kept: copying an older seed over it rewinds the generation.
        Either way the config ends up naming the file beside the job.
        """
        if self.current_seed.parent == self.current_gen_dir:
            return
        cg_seed_p = self.current_gen_dir / self.config['restart_name']
        if self.holds_own_progress(cg_seed_p):
            print(f'{cg_seed_p} is newer than its seed; keeping the progress '
                  'this generation has already made.')
        else:
            shutil.copy(self.current_seed, cg_seed_p)
        self.set_seed(cg_seed_p)

    def was_preempted(self):
        """True if the scheduler reports this clone's last job was preempted.

        The restart that follows one does not count against restarts_per_gen.
        """
        if self.preemption_checker is None:
            return False
        if self.job_number is None:
            return False
        return self.preemption_checker(self.job_number)

    def spend_restart(self):
        """Charge one restart of this generation, and say whether to go on.

        The one place a clone is written off for good: every failure to advance
        funnels through here, and False means restarts_per_gen is spent. Any
        advance that gets somewhere clears the counter, so the budget is per
        generation and not per campaign.
        """
        if self.restart_attempts >= self.restarts_per_gen:
            print(self.current_gen_dir, 'has been restarted',
                  self.restart_attempts, 'times. Aborting this clone.')
            return False
        self.restart_attempts += 1
        return True

    def note_failed_submission(self, count_as_restart):
        """Report a submission that did not take, as a strike not a give-up.

        A refused sbatch or an unreadable job id is a failed attempt at this
        generation, so it costs a restart and not the rest of the campaign.
        plow_harrow_plant already charged one unless the caller waived that.
        """
        return count_as_restart or self.spend_restart()

    def plow_harrow_plant(self, overwrite=False, count_as_restart=True):
        """Make this generation's directory and write the files a launch needs.

        The generation comes from config['gen_index'], so a caller starting a
        new one has to increment that first. Returns False once this generation
        has used up restarts_per_gen. count_as_restart False, as after a
        preemption, skips the restart_attempts increment.
        """
        # Subsequent gens restart from the previous positions and velocities.
        if self.config['gen_index'] != 0:
            self.config['new_velocities'] = False

        self.current_gen_dir = util.dir_seeds_clones_gens(
            Path(self.config['traj_dir_top_level']),
            self.config['seed_index'],
            self.config['clone_index'],
            self.config['gen_index'],
            self.config['dirname_pad'],
            sep=self.config['sep'],
            mkdir=True
        )
        # A fresh directory needs the previous seed copied into it.
        self.check_copy_set_restart_seed()
        # The absolute step this generation must reach. Counted here, not in a
        # runner: the Clone is what knows where this generation sits in the
        # chain, and the runners are handed the number.
        self.config['target_step'] = self.target_step()
        # Rewritten even when overwrite is False: a stale config on disk would
        # relaunch a half-finished gen from scratch.
        config_p = self.current_gen_dir / 'config.json'
        util.write_json_atomic(config_p, self.config, indent=4)

        job_name = self.job_name_fstring.format(**self.config)
        scheduler_script = self.scheduler_fstring.format(job_name=job_name,
                                                         **self.scheduler_kws)
        self.scheduler_script_p = (
            self.current_gen_dir / self.scheduler).with_suffix('.sh')
        if overwrite or not self.scheduler_script_p.is_file():
            self.scheduler_script_p.write_text(scheduler_script)
        run_script_p = self.current_gen_dir / self.run_script_name
        if overwrite or not run_script_p.is_file():
            run_script_p.write_text(self.run_script)
        if not count_as_restart:
            return True
        return self.spend_restart()

    def start_current(self, overwrite=False, count_as_restart=True,
                      submit=True):
        """Prepare this generation's launch, and (unless submit is False)
        submit it.

        False only once this generation's restart budget is spent: a submission
        that is refused costs a restart and still returns True.

        submit=False is what lets a ClonePack do every member's preparation,
        gen directory, seed copy, config.json, run script, and then send a
        single sbatch for the whole pack.
        """
        should_launch = self.plow_harrow_plant(
            overwrite=overwrite, count_as_restart=count_as_restart)
        print(self.get_tag(), 'should_launch',
              should_launch, 'dry_run', self.dry_run)
        if not submit:
            return should_launch
        if should_launch and not self.dry_run:
            print('launching', self.get_tag())
            # SIMULATION RUNS HERE. OUTPUT SCANNED FOR JOB NUMBER.
            # capture_output, not check_output, so stderr can be surfaced.
            with self.scheduler_script_p.open() as f:
                result = sp.run(
                    self.scheduler, stdin=f, cwd=self.current_gen_dir,
                    text=True, capture_output=True)
            if result.returncode != 0:
                # Submission failed: bad QOS, account, scheduler hiccup or
                # script. Worth another tick, so it costs a restart.
                print(f'{self.scheduler} call for {self.get_tag()} returned '
                      f'exit code {result.returncode}')
                print('  stdout:', result.stdout)
                print('  stderr:', result.stderr)
                return self.note_failed_submission(count_as_restart)
            # Assumes the scheduler prints something on submission, whose first
            # job_number_re match is the job number.
            match = self.job_number_re.search(result.stdout)
            if match is None:
                # The job is running; only its id is lost. Say so, so it can be
                # cancelled rather than left writing into an abandoned dir.
                print(f'{self.scheduler} SUBMITTED a job for {self.get_tag()} '
                      f'in {self.current_gen_dir}, but no job number could be '
                      f'read from its output. That job is RUNNING AND '
                      f'UNTRACKED: find and cancel it by hand.')
                print('  stdout:', result.stdout)
                print('  stderr:', result.stderr)
                return self.note_failed_submission(count_as_restart)
            self.job_number = int(match.group(0))
            print('Started:', self.get_tag())
        return should_launch

    def start_next(self, overwrite=False, submit=True):
        """Advance this clone one generation and launch it."""
        # An absolute path, so the restart file is found from the next gen dir.
        new_seed = self.current_gen_dir/self.config['restart_name']
        self.set_seed(new_seed.resolve())
        self.restart_attempts = 0
        # The full step count again, for the clone and for config.json.
        self.remaining_steps = self.total_steps
        self.config['steps'] = self.total_steps
        # A new generation is never a resume; an append carried over from the
        # last one starts it as though continuing a run it never began.
        self.config['append'] = False
        # Increment before building, since the point is to start the next one.
        self.config['gen_index'] += 1
        attempted_launch = self.start_current(overwrite=overwrite,
                                              submit=submit)
        return attempted_launch

    @property
    def current_gen(self):
        """The generation this clone would run next.

        Read from the config rather than tracked alongside it, so an adaptive
        sampling script that redirects a clone by moving config['gen_index']
        cannot leave the two disagreeing.
        """
        return self.config['gen_index']

    @property
    def is_done(self):
        """True once this clone has finished its last configured generation.

        Always False when last_gen_index is None: there is then no limit.
        """
        return (self.last_gen_index is not None
                and self.current_gen > self.last_gen_index)

    def gen_remaining_steps(self):
        """How many steps this generation still owes.

        progress_fn answers if one was given, otherwise frames are counted.
        """
        if self.progress_fn is not None:
            return self.progress_fn(
                self.current_gen_dir,
                total_steps=self.total_steps,
                gen_index=self.config['gen_index'],
                restart_name=self.config['restart_name'],
                traj_name=self.config['traj_name'],
                traj_suffix=self.config['traj_suffix'],
                write_interval=self.config['write_interval'],
                top_fn=self.config['top_fn'],
            )
        traj_p = (self.current_gen_dir / self.config['traj_name']
                  ).with_suffix(self.config['traj_suffix'])
        if not traj_p.is_file():
            return self.total_steps
        return util.calx_remaining_steps(
            str(traj_p), self.config['top_fn'], self.total_steps,
            self.config['write_interval'])

    def check_start_gen(self, scheduler_report: set, overwrite=False,
                        submit=True):
        """Look at where this clone is and launch whatever it needs next.

        False when no launch was attempted because this generation has spent
        its restart budget.
        """
        if self.job_number in scheduler_report:
            print('Job', self.job_number, 'still running',
                  self.job_name_fstring.format(**self.config))
            return True

        # A preemption should not burn a restart_attempt; genuine failures
        # (segfault, OOM, GPU error) still count.
        count_as_restart = not self.was_preempted()

        previous_remaining = self.remaining_steps
        self.remaining_steps = self.gen_remaining_steps()
        # A launch that got somewhere is not a restart: that is how a gen longer
        # than one allocation finishes, and progress clears the budget too.
        if self.remaining_steps < previous_remaining:
            count_as_restart = False
            self.restart_attempts = 0

        if self.remaining_steps <= 0:
            # Generation finished. start_next resets this clone's counters.
            print('Preparing to move to next generation!')
            # do any automated traj postprocessing encoded by harvester
            if self.harvester and self.reaped_gen != self.config['gen_index']:
                print('running harvester!')
                # Recorded before the attempt: one harvest is all this gen gets,
                # even if a failed start_next brings us back here next tick.
                self.reaped_gen = self.config['gen_index']
                try:
                    self.harvester.reap(
                        self.current_gen_dir, dry_run=self.dry_run)
                except Exception as exc:
                    # A harvest is post-processing; losing it must not stop the
                    # simulation campaign from advancing.
                    print(f'harvester failed for {self.get_tag()}: '
                          f'{type(exc).__name__}: {exc}; continuing.')
            if (self.last_gen_index is not None
                    and self.config['gen_index'] >= self.last_gen_index):
                # The last generation asked for; count it done instead of
                # starting one more. is_done stops anything building its dir.
                self.config['gen_index'] += 1
                return True
            return self.start_next(overwrite=overwrite, submit=submit)

        if self.remaining_steps >= self.total_steps:
            # Nothing ran. Read the scheduler log for a node-local cause before
            # the next submission overwrites it; a match steers later jobs off.
            if self.node_blocklist is not None:
                self.node_blocklist.scan_and_record(
                    self.scheduler_log_dir or self.current_gen_dir,
                    self.get_tag())
            # A trajectory of an unstarted sim can be a zero-frame file, which
            # breaks many appenders; clear it so the next launch starts clean.
            traj_p = (self.current_gen_dir / self.config['traj_name']
                      ).with_suffix(self.config['traj_suffix'])
            if traj_p.is_file():
                print(f'{traj_p} found, but zero steps. Removing and '
                      f'attempting restart number {self.restart_attempts}.')
                traj_p.unlink()
            else:
                print(f'{traj_p} not found, attempting start number '
                      f'{self.restart_attempts} for this gen.')
            self.config['steps'] = self.total_steps
            self.config['append'] = False
            return self.start_current(
                overwrite=overwrite, count_as_restart=count_as_restart,
                submit=submit)

        # Partially run: continue it.
        self.config['steps'] = self.remaining_steps
        self.config['append'] = True
        # A continuation always reloads the checkpoint's momenta; redrawing
        # them here would splice a thermal discontinuity into the trajectory.
        self.config['new_velocities'] = False
        self.check_copy_set_restart_seed()
        return self.start_current(
            overwrite=overwrite, count_as_restart=count_as_restart,
            submit=submit)


class ClonePack:
    """K Clones that share one GPU, one sbatch job, and one generation step.

    The cluster's Slurm exposes only a gpu gres, with no mps and no shard, so
    it cannot co-schedule two independent jobs onto one card. Packing therefore
    has to happen inside a single job, which is a level above Clone: each
    member does everything check_start_gen does *except* submit, and then the
    pack submits once for all of them.

    Clone is deliberately untouched by this. The pack drives members through
    the same code path a solo clone uses, check_start_gen with submit False,
    so a packed generation and a solo generation prepare identically. Only the
    submission is shared.

    Members must come from ONE condition and system, and preferably run
    identical steps: the job holds the card until its slowest member finishes,
    so mismatched per-step costs waste GPU time.

    That rule is about straggler cost, not data safety. Failure is per-member
    (gmx_pack collects K outcomes and the tender fails exactly one clone) and
    so is recovery, from mdrun -cpi off that member's own checkpoint, so a bad
    job costs a lost block that gets redone rather than a damaged dataset.
    Packing across conditions is a throughput decision, not an unsafe one.
    """

    def __init__(self, clones, pack_dir, scheduler, scheduler_fstring,
                 scheduler_kws, run_script, cpus_per_task,
                 job_name_fstring=None, job_number_re='[1-9][0-9]*',
                 job_number=None, dry_run=False,
                 pack_manifest_name='pack.json',
                 run_script_name='run.py',
                 member_cores=None,
                 # If True, unequal steps per generation are not even noted:
                 # you are saying every launch ends at -maxh, not at a target.
                 wallclock_matched=False,
                 sep=None,
                 job_name_elements=('{title}', '{seed_index}',
                                    '{clone_index}', '{gen_index}')):
        if not clones:
            raise ValueError('a ClonePack needs at least one Clone')
        steps = {c.total_steps for c in clones}
        if len(steps) != 1 and not wallclock_matched:
            print(f'NOTE: pack members run {sorted(steps)} steps per '
                  'generation. The job holds the card until its slowest member '
                  'finishes, and one Harvester cannot serve two step lengths: '
                  'a clone whose steps_per_gen disagrees with the harvester '
                  'config loses that harvest silently. Give each step length '
                  'its own Harvester.')
        self.clones = list(clones)
        self.pack_dir = Path(pack_dir)
        self.pack_dir.mkdir(parents=True, exist_ok=True)
        self.scheduler = scheduler
        self.scheduler_fstring = inspect.cleandoc(scheduler_fstring)
        # Held by reference, not copied, so nodes blocked later still reach this
        # pack. The pack's own keys are laid over it at submit time.
        self.scheduler_kws = scheduler_kws
        self.pack_scheduler_kws = {'cpus': int(cpus_per_task),
                                   'run_script_name': run_script_name}
        self.run_script = run_script
        self.run_script_name = run_script_name
        self.cpus_per_task = int(cpus_per_task)
        # One core count per member, in member order; None splits evenly.
        if member_cores is not None:
            member_cores = [int(c) for c in member_cores]
            if len(member_cores) != len(self.clones):
                raise ValueError(
                    f'member_cores has {len(member_cores)} entries for '
                    f'{len(self.clones)} members')
        self.member_cores = member_cores
        self.pack_manifest_name = pack_manifest_name
        self.job_number_re = re.compile(job_number_re)
        # The name must parse as seed, clone and gen, the way Farmer
        # re-associates at boot; one that does not gets a second job here.
        self.sep = clones[0].sep if sep is None else sep
        self.job_name_fstring = (job_name_fstring
                                 or self.sep.join(job_name_elements))
        # Member 0's config renders that name, so its job number is the pack's.
        if job_number is None:
            job_number = next((c.job_number for c in self.clones
                               if c.job_number is not None), None)
        self.job_number = job_number
        for clone in self.clones:
            clone.job_number = job_number
            # A packed job writes one scheduler log, here. No member has one
            # in its own gen dir, so a node scan rooted there finds nothing.
            clone.scheduler_log_dir = self.pack_dir
        self.dry_run = dry_run
        # Members that exhausted their restart budget. Kept out of self.clones,
        # which member_cores and the pack's job name index positionally.
        self.retired = set()

    @property
    def current_gen(self):
        """Where the pack is: the laggard among the members still live.

        A retired member's frozen current_gen would otherwise pin it forever.
        """
        live = [c.current_gen for i, c in enumerate(self.clones)
                if i not in self.retired]
        return min(live) if live else min(c.current_gen for c in self.clones)

    def get_tag(self):
        # Retired members are marked, so a short pack reads differently.
        return 'pack[' + ' | '.join(
            c.get_tag() + (' RETIRED' if i in self.retired else '')
            for i, c in enumerate(self.clones)) + ']'

    def __hash__(self):
        return hash(tuple(hash(c) for c in self.clones))

    def __eq__(self, other):
        return isinstance(other, ClonePack) and hash(self) == hash(other)

    def _job_name(self):
        head = dict(self.clones[0].config)
        return self.job_name_fstring.format(**head)

    def live_indexes(self):
        """Positions of the members that have not been retired, in pack order."""
        return [i for i in range(len(self.clones)) if i not in self.retired]

    def retire(self, index):
        """Drop a member that has spent its restart budget, keeping its slot.

        It stays in self.clones, which member_cores and the pack's job name
        index positionally; only live_indexes stops naming it.
        """
        self.retired.add(index)
        print(f'{self.clones[index].get_tag()}: exhausted its restart budget; '
              'retiring it from the pack.')

    def spend_restart(self):
        """Charge every live member one restart; False once none are left.

        A member whose budget runs out retires instead of failing the pack, so
        this says what Clone.spend_restart says: whether there is anything left
        to try.
        """
        for index in self.live_indexes():
            if not self.clones[index].spend_restart():
                self.retire(index)
        return bool(self.live_indexes())

    def check_start_gen(self, scheduler_report: set, overwrite=False):
        """Advance every member, then submit one job for the pack.

        False only once every member has spent its restart budget; a pack whose
        submission is refused charges the budget and returns True.
        """
        if self.job_number in scheduler_report:
            print('Pack job', self.job_number, 'still running', self._job_name())
            return True

        live_indexes = self.live_indexes()
        if not live_indexes:
            print(f'{self.get_tag()}: every member retired; failing pack.')
            return False

        member_configs, member_indexes = [], []
        tried = newly_finished = 0
        for index in live_indexes:
            clone = self.clones[index]
            if clone.is_done:
                continue
            tried += 1
            try:
                ok = clone.check_start_gen(scheduler_report,
                                           overwrite=overwrite, submit=False)
            except Exception as exc:
                print(f'ERROR preparing pack member {clone.get_tag()}: '
                      f'{type(exc).__name__}: {exc}')
                continue
            if not ok:
                # Its restart budget is spent, and a member that never runs can
                # never earn it back; retire it so it stops pinning the pack.
                self.retire(index)
                continue
            if clone.is_done:
                # That call finished this member's last gen; nothing to submit.
                newly_finished += 1
                continue
            member_configs.append(clone.current_gen_dir / 'config.json')
            member_indexes.append(index)

        still_live = self.live_indexes()
        if not still_live:
            print(f'{self.get_tag()}: every member retired; failing pack.')
            return False
        if all(self.clones[i].is_done for i in still_live):
            print(f'{self.get_tag()}: every live member finished its '
                  'generations.')
            return True
        if not member_configs:
            # Every live member raised on the way in, so none was charged for
            # the attempt; charge them here rather than retrying forever.
            print(f'{self.get_tag()}: no member could be prepared.')
            return self.spend_restart()
        unfinished = tried - newly_finished
        if len(member_configs) < unfinished:
            # Launch the healthy members rather than shrinking the pack for
            # good, which would leave the card underpacked from here on.
            print(f'{self.get_tag()}: {unfinished - len(member_configs)} of '
                  f'{unfinished} members could not be prepared; launching '
                  'the rest.')

        from . import gmx_pack
        # Subset the core split to the members that actually prepared.
        member_cores = (None if self.member_cores is None
                        else [self.member_cores[i] for i in member_indexes])
        gmx_pack.write_pack_manifest(
            self.pack_dir, member_configs, cpus_per_task=self.cpus_per_task,
            member_cores=member_cores,
            pack_manifest_name=self.pack_manifest_name)
        (self.pack_dir / self.run_script_name).write_text(self.run_script)
        script_p = (self.pack_dir / self.scheduler).with_suffix('.sh')
        script_p.write_text(self.scheduler_fstring.format(
            job_name=self._job_name(),
            **{**self.scheduler_kws, **self.pack_scheduler_kws}))

        if self.dry_run:
            print(f'{self.get_tag()}: dry run, wrote {script_p} and manifest.')
            return True

        with script_p.open() as f:
            result = sp.run(self.scheduler, stdin=f, cwd=self.pack_dir,
                            text=True, capture_output=True)
        if result.returncode != 0:
            # Preparing the members already charged their budgets for this
            # attempt, so the pack retries until those run out.
            print(f'{self.scheduler} call for {self.get_tag()} returned '
                  f'exit code {result.returncode}')
            print('  stdout:', result.stdout)
            print('  stderr:', result.stderr)
            return True
        match = self.job_number_re.search(result.stdout)
        if match is None:
            # As for a solo clone: the pack's job is running, untracked.
            print(f'{self.scheduler} SUBMITTED the job for {self.get_tag()} in '
                  f'{self.pack_dir}, but no job number could be read from its '
                  f'output. That job is RUNNING AND UNTRACKED: find and cancel '
                  f'it by hand.')
            print('  stdout:', result.stdout)
            print('  stderr:', result.stderr)
            return True
        self.job_number = int(match.group(0))
        # Every member answers to the pack's job id, so per-member checks work.
        for clone in self.clones:
            clone.job_number = self.job_number
        print('Started pack:', self.get_tag(), 'as job', self.job_number)
        return True
