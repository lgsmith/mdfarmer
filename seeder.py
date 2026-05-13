from __future__ import annotations
import copy
import inspect
from . import utilities as util
from pathlib import Path
import subprocess as sp
import shutil  # for copy--will be obviated by python 3.14
import re
import json


default_run_script = """
from mdfarmer.simulate import omm_basic_sim_block_json as runner
runner('config.json')
"""


# A "usable" seed is a state file that exists and is non-empty. Zero-byte
# state.xml files are a common preempt failure mode (job died before the
# CheckpointReporter completed a write) and must not be picked as seeds.
def _is_usable_seed(p: Path) -> bool:
    if not p.is_file():
        return False
    if p.stat().st_size == 0:
        return False
    return True


# Examine one gen directory and decide whether it can seed the next launch.
#
# Returns (gen_index, seed_fn, steps_to_run, append) on success:
#   - gen_index: which gen the next launch is for (may equal this dir's
#     gen, or be one higher if this gen is complete and we advance).
#   - seed_fn: absolute path to the state.xml to load from.
#   - steps_to_run: how many sim steps the next launch should take.
#   - append: whether trajectory reporters should open in append mode.
#
# Returns None if this gen is unrecoverable (missing/corrupt config, no
# usable state.xml, partial traj in non-append mode). Caller should try
# an older gen, then fall back to the initial seed.
#
# Side effect (non-append mode only): a partial traj in this gen gets
# renamed to old_<name> so the next launch starts cleanly. This matches
# the documented "cut the ends off the tree" reset behavior.
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
    if not _is_usable_seed(restart_p):
        return None
    seed_fn = str(restart_p.resolve())

    traj_p = (gen_path / traj_name).with_suffix(traj_suffix)
    if (not traj_p.is_file()) or traj_p.stat().st_size == 0:
        # No (or zero-size) traj yet; start this gen fresh from the in-dir state.
        return gen_index, seed_fn, total_steps, False

    remaining = util.calx_remaining_steps(
        str(traj_p), top_fn, total_steps, write_interval)
    if remaining > 0:
        if append_mode:
            return gen_index, seed_fn, remaining, True
        # Non-append mode: preserve the partial traj for later inspection
        # and fall through to an older gen.
        old_traj_p = traj_p.parent / ('old_' + traj_p.name)
        print(f'_try_recover_gen: non-append mode, renaming partial '
              f'{traj_p} -> {old_traj_p}')
        traj_p.rename(old_traj_p)
        return None
    # Gen is complete; advance to gen+1 seeded from this gen's terminal state.
    return gen_index + 1, seed_fn, total_steps, False


class Clone:
    __slots__ = (
        'config', 'dry_run', 'job_number', 'job_number_re', 'job_name_fstring', 'current_seed',
        'current_gen_dir', 'current_gen', 'config_p', 'scheduler_script_p', 'compare_keys',
        'scheduler_fstring', 'scheduler', 'traj_list', 'sep', 'dirname_pad',
        'scheduler_kws', 'restarts_per_gen', 'restart_attempts', 'run_script',
        'harvester', 'remaining_steps', 'run_script_name', 'total_steps',
        'preemption_checker')

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
                 # The full per-gen step count from the template. On a resume
                 # config['steps'] is the steps-remaining-this-gen, not the
                 # template total, so we have to track the total separately
                 # or start_next will reset gens to the shortened count. If
                 # None, fall back to config['steps'] for backwards compat.
                 steps_per_gen=None,
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
        self.preemption_checker = preemption_checker
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
            gen_paths = sorted(clone_dir.iterdir())
        else:
            gen_paths = []

        recovered = None
        for gen_path in reversed(gen_paths):
            recovered = _try_recover_gen(
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
            dry_run=dry_run,
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
        config_p = self.current_gen_dir / 'config.json'
        if overwrite or not config_p.is_file():
            with config_p.open('w') as f:
                json.dump(self.config, f, indent=4)
        else:
            print("There is already a config at:", str(config_p.resolve()),
                  'Proceeding with old config; did not ask for overwrite.')

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

    def start_current(self, overwrite=False, count_as_restart=True):
        should_launch = self.plow_harrow_plant(
            overwrite=overwrite, count_as_restart=count_as_restart)
        print(self.get_tag(), 'should_launch',
              should_launch, 'dry_run', self.dry_run)
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

    def start_next(self, overwrite=False):
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
        # because we want to start next, increment the gen before building
        self.config['gen_index'] += 1
        self.current_gen += 1
        attempted_launch = self.start_current(overwrite=overwrite)
        return attempted_launch

    # Returns false if launch not attempted because of too many restart attempts
    def check_start_gen(self, scheduler_report: set, overwrite=False):
        if self.job_number in scheduler_report:
            print('Job', self.job_number, 'still running',
                  self.job_name_fstring.format(**self.config))
            no_failure = True
        else:
            # If the scheduler reports this job as preempted, the upcoming
            # restart shouldn't burn a restart_attempt. Genuine failures
            # (segfault, OOM, GPU error) still count.
            if self.was_preempted():
                count_as_restart = False
            else:
                count_as_restart = True
            traj_p = (self.current_gen_dir / self.config['traj_name']
                      ).with_suffix(self.config['traj_suffix'])
            if traj_p.is_file():
                self.remaining_steps = util.calx_remaining_steps(
                    str(traj_p),
                    self.config['top_fn'],
                    self.total_steps,
                    self.config['write_interval']
                )
                # Traj of an unstarted sim can be an empty file, which breaks many appenders.
                if self.remaining_steps == self.total_steps:
                    print(
                        f'{traj_p} found, but zero steps. ',
                        f'Removing and attempting restart number {self.restart_attempts}.',
                    )
                    traj_p.unlink()
                    no_failure = self.start_current(
                        overwrite=overwrite, count_as_restart=count_as_restart)
                # if this, the trajectory has steps remaining before it is a full gen.
                # Run those.
                elif self.remaining_steps > 0:
                    self.config['steps'] = self.remaining_steps
                    # set the reporters to append
                    self.config['append'] = True
                    self.check_copy_set_restart_seed()
                    no_failure = self.start_current(
                        overwrite=overwrite, count_as_restart=count_as_restart)
                else:  # trajectory was created, and finished running.
                    self.restart_attempts = 0
                    self.config['steps'] = self.total_steps
                    print('Preparing to move to next generation!')
                    # do any automated traj postprocessing encoded by harvester
                    if self.harvester:
                        print('running harvester!')
                        self.harvester.reap(
                            self.current_gen_dir, dry_run=self.dry_run)
                    no_failure = self.start_next(overwrite=overwrite)
            else:  # no trajectory file, start from here just as if we'd found an empty file.
                print(f'{traj_p} not found, attempting start number '
                      f'{self.restart_attempts} for this gen.')
                no_failure = self.start_current(
                    overwrite=overwrite, count_as_restart=count_as_restart)
        # else:  # this triggers if no JN bound to clone ==> we should start one
            # no_failure = self.start_current(overwrite=overwrite)
        return no_failure
