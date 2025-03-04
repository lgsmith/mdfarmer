import subprocess as sp
import utilities as util
from pathlib import Path
import json
import seeder
import simulate as sims
import time




class Farmer:
    __slots__ = ('runs_clones_list', 'n_runs', 'n_clones', 'n_gens', 'runner',
                 'config_template', 'jn_regex', 'current_jids', 'overwrite',
                 'active_clone_threshold', 'active_clone_set', 'runs_first',
                 'job_name_fstring', 'job_number_re', 'finished_clones',
                 'quiet', 'sep', 'dirname_pad', 'seed_state_fns', 'scheduler',
                 'scheduler_report_cmd', 'scheduler_fstring', 'scheduler_kws',
                 'scheduler_assoc_rep_cmd')

    def update_jids(self):
        jids_string = sp.check_output(
            self.scheduler_report_cmd,
            shell=True, text=True,
            executable='/bin/bash').strip()
        print('jids_string:')
        print(jids_string)
        self.current_jids = set(map(int, jids_string.split()))
        print(self.current_jids)

    def check_path_config(self, key):
        p = Path(self.config_template[key])
        if p.is_file():
            return p
        else:
            raise FileNotFoundError(p)

    def make_clone_check_running(self, tdir, run_index, clone_index, rep_dict):
        # First, determine which (if any) generations for each run and clone have finished.
        clone_dir = util.dir_runs_clones(tdir, run_index,
                                      clone_index,
                                      self.dirname_pad,
                                      mkdir=False)
        config = util.merge_args_defaults_dict(
            self.runner, **self.config_template)
        if clone_dir.is_dir():
            # relies on padding to cause lex sort to yield highest
            # gen dir as last element
            try:
                highest_gen = sorted(clone_dir.iterdir())[-1]
            except IndexError:
                highest_gen = None
        else:
            print('No clone-dir found.')
            highest_gen = None
        if highest_gen:
            config_p = highest_gen / 'config.json'
            if config_p.is_file():
                with config_p.open() as f:
                    prev_config = json.load(f)
                gen_index = prev_config['gen_index']
                try:
                    traj_name = prev_config['traj_name']
                    traj_suff = prev_config['traj_suffix']
                except KeyError:
                    traj_name = config['traj_name']
                    traj_suff = config['traj_suffix']
                traj_p = (highest_gen / traj_name).with_suffix(
                    traj_suff)
                if traj_p.is_file():
                    if self.config_template['append']:
                        config['seed'] = traj_p.parent/self.config_template['state_fn']
                        remaining_steps = util.calx_remaining_steps(
                            str(traj_p),
                            config['prmtop_fn'],
                            prev_config['steps'],
                            prev_config['write_interval']
                        )
                        if remaining_steps > 0:
                            config['steps'] = remaining_steps
                            # change the config's seed to look at the checkpoint/state file in the traj_dir
                        else:  # this generation is done: increment gen counter.
                            gen_index += 1
                    else:
                        old_traj_p = traj_p.parent / 'old_' + traj_p.name
                        print(
                            'Traj found, but not operating in append mode.')
                        print('Moving', str(traj_p), 'to',
                                str(old_traj_p))
                        traj_p.rename(old_traj_p)
                else:
                    print('No traj found from looking in {}.'.format(
                        config_p))
                    print('Starting from fresh dir for run {run_index},'
                          ' clone {clone_index},'
                          ' gen {gen_index}.'.format(**prev_config))
            else:
                print(
                    'Could not find config, parsing gen_index from path.')
                gen_index = int(highest_gen.name.split(self.sep)[-1])
        else:
            # if we get here in control flow,  then we start fresh
            gen_index = 0
        # Now try to see if there's a currently running job with this
        # run_index, clone_index, and gen_index
        try:
            jid = rep_dict[(run_index, clone_index, gen_index)]
        except KeyError:
            jid = None

        # Finally, build the clone with whichver setting needed to be started.
        config['run_index'] = run_index
        config['clone_index'] = clone_index
        config['gen_index'] = gen_index
        if gen_index == 0:
            config['initial'] = self.seed_state_fns[run_index]
            config['new_velocities'] = True
        else:
            config['initial'] = None
            config['new_velocities'] = False
        clone = seeder.Clone(
            config,
            self.scheduler,
            self.scheduler_fstring,
            self.scheduler_kws,
            job_number=jid,
            job_number_re=self.job_number_re,
            job_name_fstring=self.job_name_fstring,
            dirname_pad=self.dirname_pad,
            sep=self.sep,
        )
        if jid:
            self.active_clone_set.add()

        return clone

    def __init__(self, n_runs: int, n_clones: int, n_gens: int,
                 config_template: dict,
                 seed_structure_fns: list,  # len(seed_structure_fns) == n_runs
                 scheduler: str,
                 scheduler_fstring: str,
                 scheduler_kws: dict,
                 scheduler_report_cmd: str,
                 scheduler_assoc_rep_cmd: str,
                 traj_list=None,
                 quiet=False,
                 active_clone_threshold=50,
                 dirname_pad=3,
                 job_number_re='[1-9][0-9]*',
                 sep='-',
                 runs_first=True,
                 job_name_elements=(
                     '{title}', '{run_index}', '{clone_index}',
                     '{gen_index}'),
                 overwrite=False,
                 runner=sims.omm_generation
                 ):
        self.n_runs = n_runs
        self.n_clones = n_clones
        self.n_gens = n_gens
        self.overwrite = overwrite
        self.config_template = config_template
        self.runner = runner
        self.runs_first = runs_first
        self.scheduler = scheduler
        self.scheduler_kws = scheduler_kws
        self.scheduler_fstring = scheduler_fstring
        self.scheduler_report_cmd = scheduler_report_cmd
        self.scheduler_assoc_rep_cmd = scheduler_assoc_rep_cmd
        self.job_number_re = job_number_re
        # Ensure that all file-names are in the config as full paths
        self.config_template['traj_dir_top_level'] = str(
            Path(self.config_template['traj_dir_top_level']).resolve()
        )
        self.config_template['prmtop_fn'] = str(
            self.check_path_config('prmtop_fn').resolve()
        )
        seed_structure_fps = [Path(s).resolve()
                              for s in seed_structure_fns]
        self.seed_state_fns = []
        for p in seed_structure_fps:
            if p.is_file():
                self.seed_state_fns.append(str(p))
            else:
                print(p)
                raise FileNotFoundError

        self.sep = sep
        self.config_template['sep'] = self.sep
        self.dirname_pad = dirname_pad
        self.config_template['dirname_pad'] = self.dirname_pad
        self.quiet = quiet
        # important to pass this down through the clones
        self.job_name_fstring = self.sep.join(job_name_elements)
        self.current_jids = set()
        self.finished_clones = dict()
        self.active_clone_threshold = active_clone_threshold
        self.active_clone_set = set()
        try:
            traj_list = self.config_template['traj_list']
            self.config_template['traj_list'] = str(Path(traj_list).resolve())
        except KeyError:
            if traj_list:
                self.config_template['traj_list'] = traj_list
            else:  # needs to be a fullpath to traj_list to append to, as string
                self.config_template['traj_list'] = str(Path('traj_list.txt')
                                                        .resolve())

        # Update in case we're restarting an orchestrator that had been running
        self.update_jids()
        rep_dict = {}
        if self.current_jids:
            # try to re-associate jobs if possible.
            report = sp.check_output(self.scheduler_assoc_rep_cmd,
                                     shell=True,
                                     executable='/bin/bash',
                                     text=True).strip().split('\n')
            print('re-association scheduler report:')
            print(report)
            for line in report:
                ls = line.split()
                jid = int(ls[0])
                rix, cix, gix = map(int, ls[1].split(self.sep)[1:])
                rep_dict[(rix, cix, gix)] = jid

        tdir = Path(self.config_template['traj_dir_top_level'])
        self.runs_clones_list = []
        # try to figure out how 'complete' all extant clones are
        if self.runs_first:
            # order list so that we try to start an unfinished clone in each run before moving to the next clone.
            for clone_index in range(self.n_clones):
                clone_queue = []
                for run_index in range(self.n_runs):
                    clone_queue.append(self.make_clone_check_running(
                        tdir, run_index, clone_index, rep_dict))
                self.runs_clones_list.append(clone_queue)
        else:
            # else order list so to start each unfinished clone in a given run before moving to the next run.
            for run_index in range(self.n_runs):
                clone_queue = []
                for clone_index in range(self.n_clones):
                    clone_queue.append(self.make_clone_check_running(
                        tdir, run_index, clone_index, rep_dict))
                self.runs_clones_list.append(clone_queue)

    def launch(self, sleep=None, update_jids=True):
        retlist = []  # note, this will be flat
        if update_jids:
            self.update_jids()
        for run_index, clone_list in enumerate(self.runs_clones_list):
            for clone_index, clone in enumerate(clone_list):
                if sleep:
                    time.sleep(sleep)

                next_up_gen = clone.config['gen_index']
                # this is only a check to see if the clone has run enough gens.
                enough_gens = next_up_gen > self.n_gens
                retlist.append(enough_gens)
                if enough_gens:
                    self.finished_clones[
                        (run_index, clone_index, next_up_gen)
                    ] = clone
                    clone_list.remove(clone)
                    try:
                        self.active_clone_set.remove(clone)
                    except KeyError:
                        print(
                            f'done_before_launch: {run_index},{clone_index},{next_up_gen}')
                elif clone in self.active_clone_set:
                    clone.check_start_gen(
                        self.current_jids, overwrite=self.overwrite
                    )
                elif len(self.active_clone_set) < self.active_clone_threshold:
                    self.active_clone_set.add(clone)
                    clone.check_start_gen(
                        self.current_jids, overwrite=self.overwrite
                    )

        return retlist

    # whether you're starting or restarting, this is probably what you want
    # if you'd just like a 'minder' process to start all your sims and keep them
    # running until they've gotten through all the generations.
    def start_seed_minder(self, update_interval=120):
        completed_list = self.launch(sleep=None)
        brake_file_p = Path('stop')
        # this needs to be while all(list of T/F for completed runs/clones)
        while not all(completed_list):
            if brake_file_p.is_file():
                print(
                    f'Brake file detected: {brake_file_p} Stopping submission loop.')
                return False
            completed_list = self.launch(sleep=None)
            if not self.quiet:
                print('The following (run clone gen) are complete:', *(
                    ','.join(map(str, self.finished_clones))))
            time.sleep(update_interval)
        # If we get here, that means the minder thinks all clones have finished.
        return True


class Adaptive(Farmer):
    __slots__ = ('current_gen', 'ranker', 'seeds', 'rank_jobscript')

    def __init__(self, n_runs, n_clones, n_gens, config_template, seed_structure_fns, scheduler, scheduler_fstring,
                 scheduler_kws, scheduler_report_cmd, scheduler_assoc_rep_cmd, traj_list=None, quiet=False,
                 active_clone_threshold=50, dirname_pad=3, job_number_re='[1-9][0-9]*', sep='-', runs_first=True,
                 job_name_elements=('{title}', '{run_index}', '{clone_index}', '{gen_index}'), overwrite=False, 
                 runner=sims.omm_generation, current_gen=None):

        super().__init__(n_runs, n_clones, n_gens, config_template, seed_structure_fns, scheduler, scheduler_fstring, 
                         scheduler_kws, scheduler_report_cmd, scheduler_assoc_rep_cmd, traj_list, quiet, 
                         active_clone_threshold, dirname_pad, job_number_re, sep, runs_first, 
                         job_name_elements, overwrite, runner):  # type: ignore

        if current_gen:
            self.current_gen = current_gen
        else:
            current_youngest_gen = self.n_gens
            for clone_list in self.runs_clones_list:
                for clone in clone_list:
                    clone_gen = clone.config['gen_index']
                    if current_youngest_gen > clone_gen:
                        current_youngest_gen = clone_gen
            self.current_gen = current_youngest_gen

    def launch(self, sleep=None, update_jids=True):
        retlist = []  # note, this will be flat
        if update_jids:
            self.update_jids()
        for run_index, clone_list in enumerate(self.runs_clones_list):
            for clone_index, clone in enumerate(clone_list):
                if sleep:
                    time.sleep(sleep)

                next_up_gen = clone.config['gen_index']
                # Check if the clone is up to date with the current gen of adaptive sampling.
                enough_gens = next_up_gen > self.current_gen
                retlist.append(enough_gens)
                if enough_gens:
                    self.finished_clones[
                        (run_index, clone_index, next_up_gen)
                    ] = clone
                    clone_list.remove(clone)
                    try:
                        self.active_clone_set.remove(clone)
                    except KeyError:
                        print(
                            f'done_before_launch: {run_index},{clone_index},{next_up_gen}')
                elif clone in self.active_clone_set:
                    clone.check_start_gen(
                        self.current_jids, overwrite=self.overwrite
                    )
                elif len(self.active_clone_set) < self.active_clone_threshold:
                    self.active_clone_set.add(clone)
                    clone.check_start_gen(
                        self.current_jids, overwrite=self.overwrite
                    )
        # if there are no clones in active clone list, after adding some and removing others, the gen is complete
        if len(self.active_clone_set) == 0:
            self.current_gen += 1
            self.rank_and_seed()

        # Add the number of gens to be run as an additional condition that must be true for clone-minder to return.
        retlist.append(self.current_gen > self.n_gens + 1)

        return retlist
