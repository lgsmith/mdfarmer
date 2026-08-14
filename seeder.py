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


# Sort key that orders gen directories numerically. Anything that isn't
# '<prefix><sep><digits>' sorts to the end, keyed by name so the order is still
# deterministic.
def _gen_sort_key(path, sep):
    tail = path.name.rsplit(sep, 1)[-1]
    if tail.isdigit():
        return (0, int(tail), '')
    return (1, 0, path.name)


# Examine one gen directory and decide whether it can seed the next launch.
#
# Returns (gen_index, seed_fn, steps_to_run, append) on success:
#   - gen_index: which gen the next launch is for (may equal this dir's
#     gen, or be one higher if this gen is complete and we advance).
#   - seed_fn: absolute path to the state.xml to load from.
#   - steps_to_run: how many sim steps the next launch should take.
#   - append: whether trajectory reporters should open in append mode.
#
# Returns None if this gen is unrecoverable. Caller cascades to an older
# gen, then falls back to the initial seed.
#
# Two failures collapse into "unrecoverable" (β policy — prune & redo):
#   1. state.xml unparseable or missing — torn-mid-write or never landed.
#   2. state.xml stepCount doesn't sit on a DCD frame boundary — header
#      doesn't match the checkpoint, both probably bad.
#
# For state-behind-DCD (kill between DCDReporter and CheckpointReporter,
# expected single-frame drift), we truncate the DCD to align and resume.
# For state-ahead-of-DCD (buffered DCD frames lost at kill before
# flush-per-frame existed), we accept state.xml as authoritative and
# advance to the next gen — the lost DCD tail is unrecoverable but the
# integrator state is intact and the next gen can seed from it.
#
# Note on the gen-relative step math: OpenMM's DCDReporter hard-codes the
# DCD header's `istart` to `reportInterval`, regardless of cumulative
# simulation step. state.xml's stepCount, by contrast, accumulates across
# gens. So we derive the absolute step at gen start from gen_index and
# total_steps rather than trusting the DCD header's istart.
def _try_recover_gen(gen_path: Path, *,
                     append_mode: bool,
                     restart_name: str,
                     traj_name: str,
                     traj_suffix: str,
                     write_interval: int,
                     total_steps: int,
                     top_fn: str):
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

    # Frame accounting differs by format. Only DCD exposes a header we can read
    # (and rewrite) directly; for XTC -- which is the config template's DEFAULT
    # traj_suffix -- we count frames and cannot truncate, so a trajectory that
    # has run ahead of its checkpoint has to be redone rather than trimmed.
    # Demanding a DCD header unconditionally made every .xtc generation look
    # unrecoverable, cascading each clone all the way back to gen 0 and
    # overwriting trajectories that were perfectly good.
    is_dcd = traj_suffix == '.dcd'
    if is_dcd:
        try:
            info = util.dcd_header_info(traj_p)
        except Exception as exc:
            print(f'_try_recover_gen: bad DCD header at {traj_p}: {exc}; '
                  f'skipping gen.')
            return None
        nset, nsavc = info['nset'], info['nsavc']
    else:
        try:
            nset = util.get_traj_len(str(traj_p), top_fn)
        except Exception as exc:
            print(f'_try_recover_gen: cannot count frames in {traj_p}: '
                  f'{type(exc).__name__}: {exc}; skipping gen.')
            return None
        nsavc = write_interval
    if nset == 0:
        return gen_index, seed_fn, total_steps, False

    gen_start_step = gen_index * total_steps
    state_offset = state_step - gen_start_step
    if state_offset <= 0 or state_offset % nsavc != 0:
        print(f'_try_recover_gen: state stepCount={state_step} not aligned '
              f'to gen start {gen_start_step} at nsavc={nsavc} ({gen_path}); '
              f'cascading.')
        return None
    target_nset = state_offset // nsavc

    if target_nset > nset:
        # State ahead of DCD — legacy preempt-buffer-drift from before
        # FlushingDCDReporter existed. state.xml's integrator state is
        # intact; advance and let the next gen seed from it.
        return gen_index + 1, seed_fn, total_steps, False
    if target_nset < nset and not is_dcd:
        # Trajectory ahead of the checkpoint, and this format cannot be trimmed
        # in place. Advancing anyway would leave frames past the checkpoint that
        # the next generation re-simulates from an earlier point -- a backward
        # jump in the concatenated trajectory. Redo the generation instead.
        print(f'_try_recover_gen: {traj_p} has {nset} frames but state.xml is '
              f'at frame {target_nset}, and {traj_suffix} cannot be truncated; '
              f'cascading so this gen is redone rather than left discontiguous.')
        return None
    if target_nset < nset:
        # Kill between DCDReporter and CheckpointReporter writes.
        print(f'_try_recover_gen: trimming {traj_p} from {nset} to '
              f'{target_nset} frames to match state.xml stepCount.')
        actual = util.truncate_dcd_to_nframes(traj_p, target_nset)
        if actual != target_nset:
            print(f'_try_recover_gen: truncate returned {actual} != '
                  f'target {target_nset}; cascading.')
            return None
        nset = actual
        # Keep the parallel velocity/force DCDs frame-aligned with positions.
        # This lives INSIDE the position-truncation branch on purpose: it only
        # runs when the position DCD is ahead of state (an unclean kill), which
        # never happens for a healthy in-flight gen (positions and state advance
        # together, so that gen takes the target_nset == nset path and skips
        # this entirely). So it inherits the position truncation's in-flight
        # safety and won't touch a live clone's tandem files during graceful
        # re-association. Trim any tandem that is ahead; if one is behind we
        # can't fabricate the missing frame without a matching checkpoint, so
        # cascade — this is a dead gen, so redoing it is safe.
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
            if tandem_suffix != '.dcd':
                print(f'_try_recover_gen: cannot frame-align {tandem_p} '
                      f'(only .dcd truncation is supported); cascading.')
                return None
            try:
                tandem_nset = util.dcd_header_info(tandem_p)['nset']
            except Exception as exc:
                print(f'_try_recover_gen: bad tandem DCD header at {tandem_p}: '
                      f'{exc}; cascading.')
                return None
            if tandem_nset == target_nset:
                continue
            if tandem_nset < target_nset:
                print(f'_try_recover_gen: {tandem_p} has {tandem_nset} frames, '
                      f'behind positions/state ({target_nset}); cannot realign '
                      f'without a matching checkpoint, cascading.')
                return None
            print(f'_try_recover_gen: trimming {tandem_p} from {tandem_nset} '
                  f'to {target_nset} frames to match positions.')
            tandem_actual = util.truncate_dcd_to_nframes(tandem_p, target_nset)
            if tandem_actual != target_nset:
                print(f'_try_recover_gen: tandem truncate returned '
                      f'{tandem_actual} != target {target_nset}; cascading.')
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
    __slots__ = (
        'config', 'dry_run', 'job_number', 'job_number_re', 'job_name_fstring', 'current_seed',
        'current_gen_dir', 'current_gen', 'config_p', 'scheduler_script_p', 'compare_keys',
        'scheduler_fstring', 'scheduler', 'traj_list', 'sep', 'dirname_pad',
        'scheduler_kws', 'restarts_per_gen', 'restart_attempts', 'run_script',
        'harvester', 'remaining_steps', 'run_script_name', 'total_steps',
        'preemption_checker', 'node_blocklist', 'progress_fn')

    # This should mostly be used by the init function, and by adaptive sampling scripts.

    def __init__(self,
                 # keys should match argument parameter names from runner
                 # function, 4x:omm_basic_sim_block
                 config: dict,

                 # Scheduler to call. Will also be used to name files later.
                 scheduler: str,
                 # BSUB/SBATCH/Scheduler script with anchors for str.format().
                 scheduler_fstring: str,
                 # Keys should match fstring anchors. Values should be desired substitution.
                 scheduler_kws: dict,
                 # should be path to the serialized xml this clone will grow from.
                 seed_fn: str,
                 # If true, call inspect.cleandoc on sched. fstring prior to binding it to self.
                 cleandoc_sched_fstring=True,
                 restarts_per_gen=3,
                 # will store the scheduler's assigned job  number for current job.
                 job_number=None,
                 dirname_pad=2,
                 sep='-',
                 run_script=default_run_script,
                 # regex to extract job_number from submission call output.
                 job_number_re='[1-9][0-9]*',
                 # A string with anchors for keys from the config to fill in  job_name at
                 # each generation. This is arg to -J flag in bsub.sh
                 job_name_fstring=None,
                 # if job_name_fstring is none, join iterable of job_name_elements
                 # and save as sef.job_name_fstring. Default vals are recommended min.
                 # If each of these are not in job name, with seed, clone, gen in that order
                 # reassociation from a killed orchestrator will fail.
                 job_name_elements=(
                     '{title}', '{seed_index}', '{clone_index}',
                     '{gen_index}'),
                 # Compare keys are used by eq and hash to determine whether two clones are equal.
                 compare_keys=('seed_index', 'clone_index'),
                 harvester=None,
                 # Callable jid -> bool, returns True if the named job was
                 # preempted by the scheduler. If provided, preemption restarts
                 # don't count against restarts_per_gen.
                 preemption_checker=None,
                 # Shared BadNodeRegistry. When a 0-step abort is detected
                 # and the gen's scheduler log matches a node-local
                 # failure pattern, scan_and_record harvests the node and
                 # excludes it from subsequent submissions. May be None
                 # (older callers and tests).
                 node_blocklist=None,
                 # The full per-gen step count from the template. On a resume
                 # config['steps'] is the steps-remaining-this-gen, not the
                 # template total, so we have to track the total separately
                 # or start_next will reset gens to the shortened count. If
                 # None, fall back to config['steps'] for backwards compat.
                 steps_per_gen=None,
                 # Callable(gen_dir, **context) -> steps still owed by that
                 # generation. None uses the frame-count inference, which is
                 # right for the OpenMM reporters. GROMACS passes
                 # gmx_simulate.gmx_gen_progress, which reads the step the
                 # runner recorded -- frame counting is off by one write
                 # interval there, because gmx writes a frame at step 0.
                 progress_fn=None,
                 dry_run=False
                 ):
        # REQUIRED ARGS below here
        self.config = config  # dict keys and values must be json serializable.
        if steps_per_gen is None:
            self.total_steps = config['steps']
        else:
            self.total_steps = steps_per_gen
        self.remaining_steps = config['steps']
        seed_p = Path(seed_fn)
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
        self.current_gen = self.config['gen_index']
        self.job_number = job_number
        self.sep = sep
        if job_name_fstring:
            self.job_name_fstring = job_name_fstring
        else:
            self.job_name_fstring = self.sep.join(job_name_elements)
        self.dirname_pad = dirname_pad
        self.job_number_re = re.compile(job_number_re)
        self.compare_keys = compare_keys
        self.harvester = harvester
        # The harvest reads the full generation length out of config.json, and
        # config['steps'] is not it -- on a resume that has been narrowed to the
        # steps still owed. Pin the untouched value here so every construction
        # path records it, not just from_disk.
        if config.get('steps_per_gen') is None:
            config['steps_per_gen'] = self.total_steps
        elif config['steps_per_gen'] != self.total_steps:
            raise ValueError(
                f"config['steps_per_gen']={config['steps_per_gen']} disagrees "
                f'with steps_per_gen={self.total_steps}; the harvest would '
                'place this clone\'s frames at the wrong global index.')
        # Config time is the only free moment to catch a generation length that
        # is not a whole number of write intervals, or a frame count that is not
        # a whole number of downsample periods. Both break the spacing of the
        # harvested streams at every seam, and neither ever fails loudly.
        if harvester is not None:
            downsample_frq = (getattr(harvester, 'run_config', None)
                              or {}).get('downsample_frq')
            if downsample_frq:
                harvest.check_commensurability(
                    self.total_steps, config['write_interval'], downsample_frq)
        self.preemption_checker = preemption_checker
        self.node_blocklist = node_blocklist
        self.progress_fn = progress_fn
        self.run_script = run_script
        self.scheduler_script_p = None  # always redefined each run

    # Construct a Clone by walking the on-disk state for (seed_index,
    # clone_index) under tdir. Decides which gen to run next, what to seed
    # it from, and whether to append or start fresh. Falls through extant
    # gens newest -> oldest, then falls back to initial_seed_fn if nothing
    # is recoverable.
    @classmethod
    def from_disk(cls,
                  tdir: Path,
                  seed_index: int,
                  clone_index: int,
                  *,
                  # The user-provided initial structure for this seed. Used
                  # when no on-disk gen is recoverable.
                  initial_seed_fn: str,
                  # Resolved, validated top and system paths for this seed.
                  top_fn: str,
                  system_fn: str,
                  # The Farmer's full config_template. Read-only here; we
                  # deepcopy before mutating.
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
                  # GROMACS support hooks. None -> OpenMM defaults:
                  #   recover_fn -> _try_recover_gen (state.xml/DCD recovery)
                  #   run_script -> Clone's default_run_script (omm runner)
                  recover_fn=None,
                  run_script=None,
                  progress_fn=None,
                  # (seed, clone, gen) -> jid for jobs currently in the
                  # scheduler queue, so we can re-associate after an
                  # orchestrator restart.
                  rep_dict: dict = None,
                  dry_run: bool = False):
        if rep_dict is None:
            rep_dict = {}
        append_mode = config_template['append']
        steps_per_gen = config_template['steps']

        clone_dir = util.dir_seeds_clones(
            tdir, seed_index, clone_index, dirname_pad, sep=sep, mkdir=False)
        if clone_dir.is_dir():
            # Sorted by the gen NUMBER, not the directory name: a lexicographic
            # sort agrees with numeric order only while every index has the same
            # width, so at dirname_pad=3 'gen-999' sorts after 'gen-1000' and
            # recovery walks back from the wrong generation. Names that don't
            # parse (stray files, harvest output) sort last and are skipped by
            # the recover function's own config.json check.
            gen_paths = sorted(clone_dir.iterdir(),
                               key=lambda p: _gen_sort_key(p, sep))
        else:
            gen_paths = []

        _recover = recover_fn if recover_fn is not None else _try_recover_gen
        recovered = None
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

        if recovered is not None:
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
        # Keep the StateDataReporter progress display consistent with the
        # actual run length so % complete isn't misleading on a resume.
        if isinstance(config.get('state_data_kwargs'), dict):
            config['state_data_kwargs'] = {
                **config['state_data_kwargs'],
                'totalSteps': steps_this_launch,
            }

        jid = rep_dict.get((seed_index, clone_index, gen_index))

        # steps_per_gen is the untouched full generation length; the GROMACS
        # runner needs it to compute an absolute cumulative step target even
        # when config['steps'] has been narrowed to a remainder.
        config['steps_per_gen'] = steps_per_gen

        # Only override Clone's default run_script when one was supplied, so
        # OpenMM callers keep the default and GROMACS callers get their runner.
        run_script_kw = {} if run_script is None else {'run_script': run_script}
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
            dirname_pad=dirname_pad,
            sep=sep,
            harvester=harvester,
            preemption_checker=preemption_checker,
            node_blocklist=node_blocklist,
            progress_fn=progress_fn,
            dry_run=dry_run,
            **run_script_kw,
        )

    # Compare a value across self config and other config in other Clone.
    def conf_value_eq(self, other: Clone, conf_key: str) -> bool:
        return self.config[conf_key] == other.config[conf_key]

    # Two clones should be the same if their config has the same seed, clone, and title in it.
    # Customize by providing a set of keys to compare.
    def __hash__(self):
        return hash(tuple(self.config[k] for k in self.compare_keys))

    # Use has to define equality;
    # Note if the two clones are using different compare keys this will nearly always be false
    def __eq__(self, other: Clone) -> bool:
        return hash(self) == hash(other)

    # Return a string representing key features of this clone
    def get_tag(self):
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

    def check_copy_set_restart_seed(self):
        # Check if there's a state save matching current state here
        seed_p = self.current_seed
        seed_dir = seed_p.parent
        if seed_dir != self.current_gen_dir:
            cg_seed_p = self.current_gen_dir/self.config['restart_name']
            # Critical! Copy the old seed file into the new dir!
            shutil.copy(seed_p, cg_seed_p)
            self.set_seed(cg_seed_p)

    # Returns True if the scheduler reports this clone's last job was
    # preempted (so the upcoming restart shouldn't count against the
    # per-gen restart_attempts budget).
    def was_preempted(self):
        if self.preemption_checker is None:
            return False
        if self.job_number is None:
            return False
        return self.preemption_checker(self.job_number)

    # note this gets the gen index from config then builds the dir for that gen
    # So, if you want to start a new generation, you have to increment/change
    # self.config['gen_index'] before calling this.
    # Tries to set up directory, and checks if we've gone over the number of restart limits
    # Returns a bool based on success (True) or failure (False) of these efforts.
    # If count_as_restart is False (e.g. last job was preempted), skip the
    # restart_attempts increment so the budget isn't burned by preemption.
    def plow_harrow_plant(self, overwrite=False, count_as_restart=True):
        # If we're running subsequent generations, we want to restart from prev.
        # positions and velocities.
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
        # If we've made a fresh directory this should copy the
        # previous seed into the new directory.
        self.check_copy_set_restart_seed()
        # config.json is ALWAYS rewritten from the in-memory config, which is the
        # authority on what this launch should do. Keeping a stale file when
        # overwrite=False (the Farmer default) meant a resume computed the right
        # `steps`, `append` and `seed_fn`, wrote none of them, and the job read
        # the first attempt's config instead -- so a partially-run generation
        # relaunched as if from scratch. Written via a temp file so a reader (or
        # a job starting concurrently) never sees a half-written config.
        config_p = self.current_gen_dir / 'config.json'
        tmp_p = config_p.with_name(config_p.name + '.tmp')
        with tmp_p.open('w') as f:
            json.dump(self.config, f, indent=4)
        tmp_p.replace(config_p)

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
        if self.restart_attempts < self.restarts_per_gen:
            self.restart_attempts += 1
            return True
        else:
            print(self.current_gen_dir, 'has been restarted',
                  self.restart_attempts, 'times. Aborting this clone.')
            return False

    def start_current(self, overwrite=False, count_as_restart=True,
                      submit=True):
        """Prepare this generation's launch, and (unless `submit` is False)
        submit it.

        `submit=False` is what lets a ClonePack do every member's preparation --
        gen directory, seed copy, config.json, run script -- and then send a
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
            # sp.run with capture_output=True (not check_output) so stderr is
            # captured and surfaced — check_output leaves err.stderr=None,
            # which made the previous failure print uninformative.
            with self.scheduler_script_p.open() as f:
                result = sp.run(
                    self.scheduler, stdin=f, cwd=self.current_gen_dir,
                    text=True, capture_output=True)
            if result.returncode != 0:
                # Submission failed (bad QOS, account, scheduler hiccup,
                # malformed script). Return False so the Farmer marks this
                # clone failed via mark_clone_failed and keeps tending the
                # rest, rather than letting one bad sbatch kill the
                # orchestrator.
                print(f'{self.scheduler} call for {self.get_tag()} returned '
                      f'exit code {result.returncode}')
                print('  stdout:', result.stdout)
                print('  stderr:', result.stderr)
                return False
            scheduler_output = result.stdout
            # NOTE: this assumes that some text is printed when a job is started,
            # and that within that text the first number matching job_number_re
            # is the Job number.
            self.job_number = int(
                self.job_number_re.search(scheduler_output).group(0))
            print('Started:', self.get_tag())
        return should_launch

    def start_next(self, overwrite=False, submit=True):
        # Set seed to be current restart file, but full path so it will be found in next gen dir.
        new_seed = self.current_gen_dir/self.config['restart_name']
        self.set_seed(new_seed.resolve())
        # reset number of restart attempts
        self.restart_attempts = 0
        # reset number of steps to take
        # first for accounting inside the clone
        self.remaining_steps = self.total_steps
        # then with respect to the number of steps to write to the config.json
        self.config['steps'] = self.total_steps
        # A new generation is never a resume. append was set True the first time
        # THIS generation was continued mid-flight, and leaving it set leaked
        # into the next generation, which then launched as if it were resuming a
        # partial run it had never started.
        self.config['append'] = False
        # because we want to start next, increment the gen before building
        self.config['gen_index'] += 1
        self.current_gen += 1
        attempted_launch = self.start_current(overwrite=overwrite,
                                              submit=submit)
        return attempted_launch

    # How many steps this generation still owes. Engine-specific when a
    # progress_fn was supplied (GROMACS reads the step the runner recorded);
    # otherwise inferred from the trajectory's frame count, which is what the
    # OpenMM reporters make true.
    def gen_remaining_steps(self):
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

    # Returns false if launch not attempted because of too many restart attempts
    def check_start_gen(self, scheduler_report: set, overwrite=False,
                        submit=True):
        if self.job_number in scheduler_report:
            print('Job', self.job_number, 'still running',
                  self.job_name_fstring.format(**self.config))
            return True

        # If the scheduler reports this job as preempted, the upcoming restart
        # shouldn't burn a restart_attempt. Genuine failures (segfault, OOM, GPU
        # error) still count.
        count_as_restart = not self.was_preempted()

        previous_remaining = self.remaining_steps
        self.remaining_steps = self.gen_remaining_steps()
        # A launch that moved the generation forward is not a "restart" in the
        # sense the budget is meant to police -- it is the normal way a
        # generation longer than one walltime allocation gets finished. Charging
        # it meant an ultralong generation exhausted restarts_per_gen and had
        # its clone abandoned while it was working perfectly.
        if self.remaining_steps < previous_remaining:
            count_as_restart = False

        if self.remaining_steps <= 0:
            # Generation finished.
            self.restart_attempts = 0
            self.config['steps'] = self.total_steps
            print('Preparing to move to next generation!')
            # do any automated traj postprocessing encoded by harvester
            if self.harvester:
                print('running harvester!')
                try:
                    self.harvester.reap(
                        self.current_gen_dir, dry_run=self.dry_run)
                except Exception as exc:
                    # A harvest is post-processing; losing it must not stop the
                    # simulation campaign from advancing.
                    print(f'harvester failed for {self.get_tag()}: '
                          f'{type(exc).__name__}: {exc}; continuing.')
            return self.start_next(overwrite=overwrite, submit=submit)

        if self.remaining_steps >= self.total_steps:
            # Nothing ran. Last chance to scan the failing job's scheduler log
            # for a node-local cause before the next submission overwrites
            # slurm.out / lsf.out. If a fatal-on-node pattern matched, the
            # registry adds the node to scheduler_kws['exclude_nodes'] so
            # plow_harrow_plant renders a directive that steers off it.
            if self.node_blocklist is not None:
                self.node_blocklist.scan_and_record(
                    self.current_gen_dir, self.get_tag())
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
        self.check_copy_set_restart_seed()
        return self.start_current(
            overwrite=overwrite, count_as_restart=count_as_restart,
            submit=submit)


class ClonePack:
    """K Clones that share one GPU, one sbatch job, and one generation step.

    The cluster's Slurm exposes only a `gpu` gres -- no `mps`, no `shard` -- so
    it cannot co-schedule two independent jobs onto one card. Packing therefore
    has to happen inside a single job, which is a level above `Clone`: each
    member does everything `check_start_gen` does *except* submit, and then the
    pack submits once for all of them.

    `Clone` is deliberately untouched by this. The pack drives members through
    the same code path a solo clone uses (`check_start_gen(..., submit=False)`),
    so a packed generation and a solo generation prepare identically -- only the
    submission is shared.

    Members must come from ONE condition and system, with identical `steps`:
    the job holds the card until its slowest member finishes, so mismatched
    per-step costs waste GPU time, and packing across conditions would let one
    bad job damage two datasets at once.
    """

    def __init__(self, clones, pack_dir, scheduler, scheduler_fstring,
                 scheduler_kws, run_script, cpus_per_task,
                 job_name_fstring=None, job_number_re='[1-9][0-9]*',
                 job_number=None, dry_run=False,
                 pack_manifest_name='pack.json',
                 run_script_name='run.py'):
        if not clones:
            raise ValueError('a ClonePack needs at least one Clone')
        steps = {c.total_steps for c in clones}
        if len(steps) != 1:
            raise ValueError(
                f'pack members must all run the same number of steps per '
                f'generation (got {sorted(steps)}); a shorter member would '
                'leave the card idle waiting for the longer one, and variable '
                'generation lengths complicate the contiguity bookkeeping.')
        self.clones = list(clones)
        self.pack_dir = Path(pack_dir)
        self.pack_dir.mkdir(parents=True, exist_ok=True)
        self.scheduler = scheduler
        self.scheduler_fstring = inspect.cleandoc(scheduler_fstring)
        self.scheduler_kws = dict(scheduler_kws)
        self.scheduler_kws.setdefault('run_script_name', run_script_name)
        self.scheduler_kws['cpus'] = cpus_per_task
        self.run_script = run_script
        self.run_script_name = run_script_name
        self.cpus_per_task = int(cpus_per_task)
        self.pack_manifest_name = pack_manifest_name
        self.job_number = job_number
        self.job_number_re = re.compile(job_number_re)
        self.job_name_fstring = job_name_fstring or '{title}-pack-{seed_index}-{clone_index}'
        self.dry_run = dry_run

    @property
    def current_gen(self):
        # The pack advances together, so the laggard defines where it is.
        return min(c.current_gen for c in self.clones)

    def get_tag(self):
        return 'pack[' + ' | '.join(c.get_tag() for c in self.clones) + ']'

    def __hash__(self):
        return hash(tuple(hash(c) for c in self.clones))

    def __eq__(self, other):
        return isinstance(other, ClonePack) and hash(self) == hash(other)

    def _job_name(self):
        head = dict(self.clones[0].config)
        return self.job_name_fstring.format(**head)

    def check_start_gen(self, scheduler_report: set, overwrite=False):
        """Advance every member, then submit one job for the pack."""
        if self.job_number in scheduler_report:
            print('Pack job', self.job_number, 'still running', self._job_name())
            return True

        prepared, member_configs = [], []
        for clone in self.clones:
            try:
                ok = clone.check_start_gen(scheduler_report,
                                           overwrite=overwrite, submit=False)
            except Exception as exc:
                print(f'ERROR preparing pack member {clone.get_tag()}: '
                      f'{type(exc).__name__}: {exc}')
                ok = False
            prepared.append(ok)
            if ok:
                member_configs.append(clone.current_gen_dir / 'config.json')
        if not any(prepared):
            print(f'{self.get_tag()}: no member could be prepared; failing pack.')
            return False
        if not all(prepared):
            # Relaunch the pack with the members that are still healthy rather
            # than shrinking it permanently: a shrunk pack leaves the card
            # underpacked for the rest of the campaign.
            print(f'{self.get_tag()}: {prepared.count(False)} of '
                  f'{len(prepared)} members could not be prepared; launching '
                  'the rest.')

        from . import gmx_pack
        gmx_pack.write_pack_manifest(
            self.pack_dir, member_configs, cpus_per_task=self.cpus_per_task,
            reps_per_card=len(member_configs),
            pack_manifest_name=self.pack_manifest_name)
        (self.pack_dir / self.run_script_name).write_text(self.run_script)
        script_p = (self.pack_dir / self.scheduler).with_suffix('.sh')
        script_p.write_text(self.scheduler_fstring.format(
            job_name=self._job_name(), **self.scheduler_kws))

        if self.dry_run:
            print(f'{self.get_tag()}: dry run, wrote {script_p} and manifest.')
            return True

        with script_p.open() as f:
            result = sp.run(self.scheduler, stdin=f, cwd=self.pack_dir,
                            text=True, capture_output=True)
        if result.returncode != 0:
            print(f'{self.scheduler} call for {self.get_tag()} returned '
                  f'exit code {result.returncode}')
            print('  stdout:', result.stdout)
            print('  stderr:', result.stderr)
            return False
        match = self.job_number_re.search(result.stdout)
        if match is None:
            print(f'could not parse a job number from {result.stdout!r}')
            return False
        self.job_number = int(match.group(0))
        # Every member answers to the pack's job id, so the Farmer's
        # still-running check works per member as well as per pack.
        for clone in self.clones:
            clone.job_number = self.job_number
        print('Started pack:', self.get_tag(), 'as job', self.job_number)
        return True
