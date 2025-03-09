import subprocess as sp
from . import utilities as util
from pathlib import Path
import json
from . import seeder
from . import simulate as sims
import time
import copy


class Farmer:
    __slots__ = ('priority_ordered_clones', 'n_seeds', 'n_clones', 'n_gens', 'runner',
                 'config_template', 'jn_regex', 'current_jids', 'overwrite',
                 'active_clone_threshold', 'active_clone_set', 'failed_clone_set', 'seeds_first',
                 'job_name_fstring', 'job_number_re', 'finished_clones', 'harvester',
                 'quiet', 'sep', 'dirname_pad', 'seed_state_fns', 'scheduler',
                 'scheduler_report_cmd', 'scheduler_fstring', 'scheduler_kws',
                 'scheduler_assoc_rep_cmd', 'system_fns', 'top_fns')

    def update_jids(self):
        jids_string = sp.check_output(
            self.scheduler_report_cmd,
            shell=True, text=True,
            executable='/bin/bash').strip()
        print('jids_string:')
        print(jids_string)
        self.current_jids = set(map(int, jids_string.split()))
        print(self.current_jids)

    def check_path(self, p: Path):
        if p.is_file():
            return p
        else:
            raise FileNotFoundError(p)

    def check_path_config(self, key):
        p = Path(self.config_template[key])
        self.check_path(p)

    # move clone off the active list
    def mark_clone_failed(self, clone):
        self.failed_clone_set.add(clone)
        try:
            self.active_clone_set.remove(clone)
        except KeyError:  # if clone isn't in active set that's OK.
            pass
        print('FAILED CLONE:', clone.get_tag())

    # Check if clone has finished all its generations. 
    # Remove from clone_list, and active set, and add to finished set. 
    # return True if finished, False if not.
    def check_mark_clone_finished(self, clone):
        next_up_gen = clone.current_gen
        enough_gens = next_up_gen > self.n_gens
        if enough_gens:
            print('Finished:', clone.get_tag())
            self.finished_clones.add(clone)
            try:
                self.active_clone_set.remove(clone)
            except KeyError:
                print('done_before_launch', clone.get_tag())
        return enough_gens

    def make_clone_check_running(self, tdir, seed_index, clone_index, rep_dict):
        # First, determine which (if any) generations for each run and clone have finished.
        clone_dir = util.dir_seeds_clones(tdir, seed_index,
                                          clone_index,
                                          self.dirname_pad,
                                          mkdir=False)
        config = copy.deepcopy(self.config_template)
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
            try:
                prev_config_raw = config_p.read_text()
                prev_config = json.loads(prev_config_raw)
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
                        # Because jobs will be launched from traj_p.parent file should be there
                        config['seed_fn'] = self.config_template['restart_fn']
                        remaining_steps = util.calx_remaining_steps(
                            str(traj_p),
                            config['top_fn'],
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
                    print('Starting from fresh dir for seed {seed_index},'
                          ' clone {clone_index},'
                          ' gen {gen_index}.'.format(**prev_config))
            except FileNotFoundError:
                print(
                    'Could not find config, parsing gen_index from path.')
                gen_index = int(highest_gen.name.split(self.sep)[-1])
            except json.decoder.JSONDecodeError:
                print('Config at', config_p,
                      'appears to contain malformed json; here is the raw string:')
                print(prev_config_raw, '\nProceeding to obtain gen index from path.')
                gen_index = int(highest_gen.name.split(self.sep)[-1])
        else:
            # if we get here in control flow,  then we start fresh
            gen_index = 0
        # Now try to see if there's a currently running job with this
        # seed_index, clone_index, and gen_index
        try:
            jid = rep_dict[(seed_index, clone_index, gen_index)]
        except KeyError:
            jid = None

        # Finally, build the clone with whichver setting needed to be started.
        config['seed_index'] = seed_index
        config['clone_index'] = clone_index
        config['gen_index'] = gen_index
        # Pick top and sys based on seed index:
        sys_fn = self.system_fns[seed_index]
        config['system_fn'] = str(self.check_path(Path(sys_fn)).resolve())

        top_fn = self.top_fns[seed_index]
        config['top_fn'] = str(self.check_path(Path(top_fn)).resolve())

        if gen_index == 0:
            config['new_velocities'] = True
        else:
            config['new_velocities'] = False
        clone = seeder.Clone(
            config,
            self.scheduler,
            self.scheduler_fstring,
            self.scheduler_kws,
            self.seed_state_fns[seed_index],
            job_number=jid,
            job_number_re=self.job_number_re,
            job_name_fstring=self.job_name_fstring,
            dirname_pad=self.dirname_pad,
            sep=self.sep
        )
        if jid:
            self.active_clone_set.add(clone)

        return clone

    def __init__(self, n_seeds: int, n_clones: int, n_gens: int,
                 config_template: dict,
                 # len(seed_structure_fns) == n_seeds
                 seed_structure_fns: list,
                 system_fns: list,
                 top_fns: list,
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
                 seeds_first=True,
                 job_name_elements=(
                     '{title}', '{seed_index}', '{clone_index}',
                     '{gen_index}'),
                 overwrite=False,
                 harvester=None,
                 runner=sims.omm_generation
                 ):
        self.n_seeds = n_seeds
        self.n_clones = n_clones
        self.n_gens = n_gens
        self.overwrite = overwrite
        self.config_template = config_template
        self.system_fns = system_fns
        self.top_fns = top_fns
        self.runner = runner
        self.seeds_first = seeds_first
        self.scheduler = scheduler
        self.scheduler_kws = scheduler_kws
        self.scheduler_fstring = scheduler_fstring
        self.scheduler_report_cmd = scheduler_report_cmd
        self.scheduler_assoc_rep_cmd = scheduler_assoc_rep_cmd
        self.job_number_re = job_number_re
        self.harvester = harvester
        # Ensure that all file-names are in the config as full paths
        self.config_template['traj_dir_top_level'] = str(
            Path(self.config_template['traj_dir_top_level']).resolve()
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
        self.finished_clones = set()
        self.active_clone_threshold = active_clone_threshold
        self.active_clone_set = set()
        self.failed_clone_set = set()
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
                six, cix, gix = map(int, ls[1].split(self.sep)[1:])
                rep_dict[(six, cix, gix)] = jid

        tdir = Path(self.config_template['traj_dir_top_level'])
        self.priority_ordered_clones = []
        # try to figure out how 'complete' all extant clones are
        if self.seeds_first:
            # order list so that we try to start an unfinished clone in each
            # seed before moving to the next clone.
            for clone_index in range(self.n_clones):
                clone_queue = []
                for seed_index in range(self.n_seeds):
                    clone_queue.append(self.make_clone_check_running(
                        tdir, seed_index, clone_index, rep_dict))
                self.priority_ordered_clones.append(clone_queue)
        else:
            # else order list so to start each unfinished clone in
            # a given seed before moving to the next seed.
            for seed_index in range(self.n_seeds):
                clone_queue = []
                for clone_index in range(self.n_clones):
                    clone_queue.append(self.make_clone_check_running(
                        tdir, seed_index, clone_index, rep_dict))
                self.priority_ordered_clones.append(clone_queue)

    def launch(self, sleep=None, update_jids=True):
        still_running = []  # note, this will be flat
        if update_jids:
            self.update_jids()
        for clone_list in self.priority_ordered_clones:
            # Record one False for a fully emptied clone-list
            if not clone_list:
                still_running.append(False)
            else:
                clone_indexes_to_remove = []
                for i, clone in enumerate(clone_list):
                    if sleep:
                        time.sleep(sleep)
                    # This probably shouldn't happen, but it's worth checking for
                    if clone in self.finished_clones or \
                          clone in self.failed_clone_set:
                        still_running.append(False)
                        clone_indexes_to_remove.append(i)
                    elif self.check_mark_clone_finished(clone):
                        still_running.append(False)
                        clone_indexes_to_remove.append(i)
                    # If clone is in active set, it may have just finished a generation. 
                    elif clone in self.active_clone_set:
                        # Try to start another.
                        if clone.check_start_gen( self.current_jids, overwrite=self.overwrite):
                            still_running.append(True)
                        else:
                            self.mark_clone_failed(clone)
                            still_running.append(False)
                            clone_indexes_to_remove.append(i)

                    # This condition arises when there are few enough active clones 
                    # that we could launch more.
                    elif len(self.active_clone_set) < self.active_clone_threshold:
                        #  So we try to launch another.
                        if clone.check_start_gen(self.current_jids, overwrite=self.overwrite):
                            self.active_clone_set.add(clone)
                            still_running.append(True)
                        else:
                            self.mark_clone_failed(clone, clone_list)
                            still_running.append(False)
                            clone_indexes_to_remove.append(i)
                    else:
                        print('WARNING:', clone.get_tag(), 'is not accounted for by launch logic.')
                # Because we are changing the length of the list, this must be done in reverse order
                clone_indexes_to_remove.reverse()
                for i in clone_indexes_to_remove:
                    del clone_list[i]
        return still_running

    # whether you're starting or restarting, this is probably what you want
    # if you'd just like a 'minder' process to start all your sims and keep them
    # running until they've gotten through all the generations.
    def start_tending_fields(self, update_interval=120):
        completed_list = self.launch(sleep=None)
        print('completed_list', *completed_list)
        brake_file_p = Path('stop')
        # this needs to be while all(list of T/F for completed seeds/clones)
        while not all(completed_list):
            if brake_file_p.is_file():
                print(
                    f'Brake file detected: {brake_file_p.resolve()} Stopping submission loop.')
                return False
            completed_list = self.launch(sleep=None)
            if not self.quiet:
                print('The following (seed clone gen) are complete:', *(
                    ','.join(map(str, self.finished_clones))))
            time.sleep(update_interval)
        # If we get here, that means the minder thinks all clones have finished.
        return True


class Adaptive(Farmer):
    __slots__ = ('current_gen', 'ranker', 'seeds', 'rank_jobscript')

    def __init__(self, n_seeds, n_clones, n_gens, config_template, seed_structure_fns, scheduler, scheduler_fstring,
                 scheduler_kws, scheduler_report_cmd, scheduler_assoc_rep_cmd, traj_list=None, quiet=False,
                 active_clone_threshold=50, dirname_pad=3, job_number_re='[1-9][0-9]*', sep='-', seeds_first=True,
                 job_name_elements=('{title}', '{seed_index}', '{clone_index}', '{gen_index}'), overwrite=False,
                 runner=sims.omm_generation, current_gen=None):

        super().__init__(n_seeds, n_clones, n_gens, config_template, seed_structure_fns, scheduler, scheduler_fstring,
                         scheduler_kws, scheduler_report_cmd, scheduler_assoc_rep_cmd, traj_list, quiet,
                         active_clone_threshold, dirname_pad, job_number_re, sep, seeds_first,
                         job_name_elements, overwrite, runner)

        if current_gen:
            self.current_gen = current_gen
        else:
            current_youngest_gen = self.n_gens
            for clone_list in self.priority_ordered_clones:
                for clone in clone_list:
                    clone_gen = clone.config['gen_index']
                    if current_youngest_gen > clone_gen:
                        current_youngest_gen = clone_gen
            self.current_gen = current_youngest_gen

    def launch(self, sleep=None, update_jids=True):
        retlist = []  # note, this will be flat
        if update_jids:
            self.update_jids()
        for seed_index, clone_list in enumerate(self.priority_ordered_clones):
            for clone_index, clone in enumerate(clone_list):
                if sleep:
                    time.sleep(sleep)

                next_up_gen = clone.config['gen_index']
                # Check if the clone is up to date with the current gen of adaptive sampling.
                enough_gens = next_up_gen > self.current_gen
                retlist.append(enough_gens)
                if enough_gens:
                    self.finished_clones[
                        (seed_index, clone_index, next_up_gen)
                    ] = clone
                    clone_list.remove(clone)
                    try:
                        self.active_clone_set.remove(clone)
                    except KeyError:
                        print(
                            f'done_before_launch: {seed_index},{clone_index},{next_up_gen}')
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
