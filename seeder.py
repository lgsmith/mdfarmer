import inspect
import utilities as util
from pathlib import Path
import subprocess as sp
import re
import json

from __future__ import annotations

default_run_script = """
from mdfarmer.simulate import omm_basic_sim_block_json as runner
runner('config.json')
"""


class Clone:
    __slots__ = (
        'config', 'job_number', 'job_number_re', 'job_name_fstring', 'current_seed',
        'current_gen_dir', 'config_p', 'scheduler_script_p', 'compare_keys',
        'scheduler_fstring', 'scheduler', 'traj_list', 'sep', 'dirname_pad',
        'scheduler_kws', 'restarts_per_gen', 'restart_attempts', 'run_script')

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
                 seed: str,
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
                 # If each of these are not in job name, with run, clone, gen in that order
                 # reassociation from a killed orchestrator will fail.
                 job_name_elements=(
                     '{title}', '{run_index}', '{clone_index}',
                     '{gen_index}'),
                 # Compare keys are used by eq and hash to determine whether two clones are equal.
                 compare_keys=('title', 'run_index', 'clone_index')
                 ):
        # REQUIRED ARGS below here
        self.config = config  # dict keys and values must be json serializable.
        seed_p = Path(seed)
        if seed_p.is_file():
            self.current_seed = seed_p
            self.config['seed'] = seed
        else:
            raise FileNotFoundError(seed_p)
        self.scheduler = scheduler
        if cleandoc_sched_fstring:
            self.scheduler_fstring = inspect.cleandoc(scheduler_fstring)
        else:
            self.scheduler_fstring = scheduler_fstring
        self.scheduler_kws = scheduler_kws
        self.run_script = run_script

        # Args with defaults below here
        self.restarts_per_gen = restarts_per_gen
        # this should always start at zero, since it's incremented below.
        self.restart_attempts = 0
        self.current_gen_dir = util.dir_runs_clones_gens(
            Path(self.config['traj_dir_top_level']),
            self.config['run_index'],
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
        self.compare_keys = compare_keys

    # Compare a value across self config and other config in other Clone.
    def conf_value_eq(self, other: Clone, conf_key: str) -> bool:
        return self.config[conf_key] == other.config[conf_key]

    # Two clones should be the same if their config has the same run, clone, and title in it.
    # Customize by providing a set of keys to compare.
    def __hash__(self):
        return hash(tuple(self.config[k] for k in self.compare_keys))

    # Use has to define equality;
    # Note if the two clones are using different compare keys this will nearly always be false
    def __eq__(self, other: Clone) -> bool:
        return hash(self) == hash(other)

    def set_seed(self, seed):
        if seed.is_file():
            self.current_seed = seed
            self.config['seed'] = str(seed)
        else:
            raise FileNotFoundError(seed)
    # note this gets the gen index from config then builds the dir for that gen
    # So, if you want to start a new generation, you have to increment/change
    # self.config['gen_index'] before calling this.

    def plow_harrow_plant(self, overwrite=False):
        # If we're running subsequent generations, we want to restart from prev.
        # positions and velocities.
        if self.config['gen_index'] != 0:
            self.config['initial'] = None
            self.config['new_velocities'] = False

        self.current_gen_dir = util.dir_runs_clones_gens(
            Path(self.config['traj_dir_top_level']),
            self.config['run_index'],
            self.config['clone_index'],
            self.config['gen_index'],
            self.config['dirname_pad'],
            sep=self.config['sep'],
            mkdir=True
        )

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
        run_script_p = self.current_gen_dir / 'run.py'
        if overwrite or not run_script_p.is_file():
            run_script_p.write_text(self.run_script)

    def check_seed(self):
        # Check if there's a state save matching current state here
        if self.current_seed.is_file():
            seed_dir = self.current_seed.parent
            if seed_dir != self.current_gen_dir:
                cg_seed_p = self.current_gen_dir/self.config['restart_name']
                if cg_seed_p.is_file():
                    self.current_seed = cg_seed_p
            self.config['seed'] = str(self.current_seed)
        else:
            raise FileNotFoundError(self.current_seed)

    def start_current(self, overwrite=False):
        self.plow_harrow_plant(overwrite=overwrite)
        # SIMULATION RUNS HERE. OUTPUT SCANNED FOR JOB NUMBER.
        try:
            with self.scheduler_script_p.open() as f:
                scheduler_output = sp.check_output(
                    self.scheduler, stdin=f, cwd=self.current_gen_dir, text=True)
        except sp.CalledProcessError as err:
            print(f'{self.scheduler} call threw error', err.stdout, err.stderr)
            raise
        # NOTE: this assumes that some text is printed when a job is started,
        # and that within that text the first number matching job_number_re
        # is the Job number.
        self.job_number = int(
            self.job_number_re.search(scheduler_output).group(0))
        print('started job', self.job_number, 'run', self.config['run_index'],
              'clone', self.config['clone_index'], 'gen',
              self.config['gen_index'])

    def start_next(self, overwrite=False, seed=None):
        # because we want to start next, increment the gen before building
        self.config['gen_index'] += 1
        if seed:
            self.current_seed = seed
        self.start_current(overwrite=overwrite)

    def check_start_gen(self, scheduler_report: set, overwrite=False):
        if self.job_number:
            if self.job_number in scheduler_report:
                print('Job', self.job_number, 'still running',
                      self.job_name_fstring.format(**self.config))
            else:
                traj_p = (self.current_gen_dir / self.config['traj_name']
                          ).with_suffix(self.config['traj_suffix'])
                if traj_p.is_file():
                    total_steps = self.config['steps']
                    remaining_steps = util.calx_remaining_steps(
                        str(traj_p),
                        self.config['prmtop_fn'],
                        total_steps,
                        self.config['write_interval']
                    )
                    # Traj of an unstarted sim can be an empty file, which breaks many appenders.
                    if remaining_steps == total_steps:
                        print(
                            '{} found, but zero steps. Removing and attempting restart {} '
                            'for this gen.'.format(
                                traj_p, self.restart_attempts)
                        )
                        traj_p.unlink()
                        self.restart_attempts += 1
                        self.start_current(overwrite=overwrite)
                    # if this, the trajectory has steps remaining before it is a full gen. Run those.
                    elif remaining_steps > 0:
                        if self.restarts_per_gen > self.restart_attempts:
                            self.restart_attempts += 1
                            self.config['steps'] = remaining_steps
                            self.check_seed()
                            self.start_current(overwrite=overwrite)
                            # need to reset total for next round
                            self.config['steps'] = total_steps
                        else:
                            print(
                                'Tried to restart from '
                                '{} too many times.'.format(
                                    self.current_gen_dir)
                            )

                    else:  # trajectory was created, and finished running.
                        self.restart_attempts = 0
                        self.start_next(overwrite=overwrite)
                else:  # no trajectory file, start from here just as if we'd found an empty file.
                    print(
                        '{} not found, problem here. Attempting restart {} '
                        'for this gen.'.format(traj_p, self.restart_attempts)
                    )
                    self.restart_attempts += 1
                    self.start_current(overwrite=overwrite)
        else:  # this triggers if no JN bound to clone ==> we should start one
            self.start_current(overwrite=overwrite)
        return self.config['gen_index']
