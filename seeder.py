from __future__ import annotations
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


class Clone:
    __slots__ = (
        'config', 'dry_run', 'job_number', 'job_number_re', 'job_name_fstring', 'current_seed',
        'current_gen_dir', 'current_gen', 'config_p', 'scheduler_script_p', 'compare_keys',
        'scheduler_fstring', 'scheduler', 'traj_list', 'sep', 'dirname_pad',
        'scheduler_kws', 'restarts_per_gen', 'restart_attempts', 'run_script',
        'harvester', 'remaining_steps', 'run_script_name', 'total_steps')

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
                 dry_run=False
                 ):
        # REQUIRED ARGS below here
        self.config = config  # dict keys and values must be json serializable.
        self.total_steps = config['steps']
        self.remaining_steps = self.total_steps
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
        self.run_script = run_script
        self.scheduler_script_p = None  # always redefined each run

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
            self.config['seed_fn'] = str(seed_fp.name)
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

    # note this gets the gen index from config then builds the dir for that gen
    # So, if you want to start a new generation, you have to increment/change
    # self.config['gen_index'] before calling this.
    # Tries to set up directory, and checks if we've gone over the number of restart limits
    # Returns a bool based on success (True) or failure (False) of these efforts.
    def plow_harrow_plant(self, overwrite=False):
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
        conf_dict = self.config.copy()
        # get rid of inputs the basic generation cant ingest:
        del conf_dict['traj_dir_top_level']
        del conf_dict['sep']
        del conf_dict['dirname_pad']
        # set seed to
        if overwrite or not config_p.is_file():
            with config_p.open('w') as f:
                json.dump(conf_dict, f, indent=4)
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
        if self.restart_attempts < self.restarts_per_gen:
            self.restart_attempts += 1
            return True
        else:
            print(self.current_gen_dir, 'has been restarted',
                  self.restart_attempts, 'times. Aborting this clone.')
            return False

    def start_current(self, overwrite=False):
        should_launch = self.plow_harrow_plant(overwrite=overwrite)
        print(self.get_tag(), 'should_launch',
              should_launch, 'dry_run', self.dry_run)
        if should_launch and not self.dry_run:
            print('launching', self.get_tag())
            # SIMULATION RUNS HERE. OUTPUT SCANNED FOR JOB NUMBER.
            try:
                with self.scheduler_script_p.open() as f:
                    scheduler_output = sp.check_output(
                        self.scheduler, stdin=f, cwd=self.current_gen_dir, text=True)
            except sp.CalledProcessError as err:
                print(f'{self.scheduler} call threw error',
                      err.stdout, err.stderr)
                raise
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
                    no_failure = self.start_current(overwrite=overwrite)
                # if this, the trajectory has steps remaining before it is a full gen.
                # Run those.
                elif self.remaining_steps > 0:
                    self.config['steps'] = self.remaining_steps
                    # set the reporters to append
                    self.config['append'] = True
                    self.check_copy_set_restart_seed()
                    no_failure = self.start_current(overwrite=overwrite)
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
                no_failure = self.start_current(overwrite=overwrite)
        return no_failure
