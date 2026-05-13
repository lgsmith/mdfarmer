import subprocess as sp
from . import utilities as util
from pathlib import Path
from . import seeder
from . import simulate as sims
import time


class Farmer:
    __slots__ = ('priority_ordered_clones', 'n_seeds', 'n_clones', 'n_gens', 'runner', 'jids_file',
                 'config_template', 'jn_regex', 'current_jids', 'dry_run', 'overwrite',
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
        print(f'jids_string:\n{jids_string}')
        self.current_jids = set(map(int, jids_string.split()))
        self.jids_file.write_text(jids_string)

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

    # Build one Clone via disk-state discovery, isolating failures so one
    # corrupt clone dir can't kill orchestrator boot.
    def _setup_one_clone(self, tdir, seed_index, clone_index, rep_dict):
        try:
            clone = seeder.Clone.from_disk(
                tdir, seed_index, clone_index,
                initial_seed_fn=self.seed_state_fns[seed_index],
                top_fn=self.top_fns[seed_index],
                system_fn=self.system_fns[seed_index],
                config_template=self.config_template,
                scheduler=self.scheduler,
                scheduler_fstring=self.scheduler_fstring,
                scheduler_kws=self.scheduler_kws,
                dirname_pad=self.dirname_pad,
                sep=self.sep,
                job_number_re=self.job_number_re,
                job_name_fstring=self.job_name_fstring,
                harvester=self.harvester,
                preemption_checker=util.preemption_checkers.get(self.scheduler),
                rep_dict=rep_dict,
                dry_run=self.dry_run,
            )
        except Exception as exc:
            print(f'Skipping clone seed={seed_index} clone={clone_index} '
                  f'during setup: {type(exc).__name__}: {exc}')
            return None
        if clone.job_number is not None:
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
                 runner=sims.omm_generation,
                 dry_run=False
                 ):
        self.n_seeds = n_seeds
        self.n_clones = n_clones
        self.n_gens = n_gens
        self.overwrite = overwrite
        self.config_template = config_template
        # Resolve and validate top / system paths once up front so per-clone
        # setup doesn't repeat the work and a bad file fails loudly at boot.
        self.system_fns = [str(self.check_path(Path(p)).resolve())
                           for p in system_fns]
        self.top_fns = [str(self.check_path(Path(p)).resolve())
                        for p in top_fns]
        self.runner = runner
        self.seeds_first = seeds_first
        self.scheduler = scheduler
        self.scheduler_kws = scheduler_kws
        self.scheduler_fstring = scheduler_fstring
        self.scheduler_report_cmd = scheduler_report_cmd
        self.scheduler_assoc_rep_cmd = scheduler_assoc_rep_cmd
        self.job_number_re = job_number_re
        self.harvester = harvester
        self.dry_run = dry_run
        # Ensure that all file-names are in the config as full paths
        self.config_template['traj_dir_top_level'] = str(
            Path(self.config_template['traj_dir_top_level']).resolve()
        )
        self.config_template['integrator_xml'] = str(
            self.check_path(Path(self.config_template['integrator_xml'])).resolve()
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
        self.jids_file = Path(f'{config_template["title"]}-jids.txt')
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

        # One scheduler query at boot, used to populate both current_jids and
        # rep_dict. Calling scheduler_report_cmd and scheduler_assoc_rep_cmd
        # separately let a job appear / disappear between the two, producing a
        # stale jid that would be bound to a Clone and double-launched on the
        # first tick.
        self.current_jids = set()
        rep_dict = {}
        assoc_raw = sp.check_output(self.scheduler_assoc_rep_cmd,
                                    shell=True,
                                    executable='/bin/bash',
                                    text=True).strip()
        print('boot re-association scheduler report:')
        print(assoc_raw)
        for line in assoc_raw.split('\n'):
            if not line.strip():
                continue
            ls = line.split()
            try:
                jid = int(ls[0])
            except (ValueError, IndexError):
                continue
            self.current_jids.add(jid)
            try:
                six, cix, gix = map(int, ls[1].split(self.sep)[1:])
                rep_dict[(six, cix, gix)] = jid
            except (ValueError, IndexError):
                # Job name doesn't fit our seed-clone-gen suffix scheme;
                # leave it in current_jids so we don't relaunch over it
                # but don't try to bind it to a Clone.
                continue
        self.jids_file.write_text(' '.join(map(str, sorted(self.current_jids))))

        tdir = Path(self.config_template['traj_dir_top_level'])
        self.priority_ordered_clones = []
        if self.seeds_first:
            outer_range = range(self.n_clones)
            inner_range = range(self.n_seeds)
            indexed = lambda outer, inner: (inner, outer)  # (seed, clone)
        else:
            outer_range = range(self.n_seeds)
            inner_range = range(self.n_clones)
            indexed = lambda outer, inner: (outer, inner)  # (seed, clone)
        for outer in outer_range:
            clone_queue = []
            for inner in inner_range:
                seed_index, clone_index = indexed(outer, inner)
                clone = self._setup_one_clone(
                    tdir, seed_index, clone_index, rep_dict)
                if clone is not None:
                    clone_queue.append(clone)
            self.priority_ordered_clones.append(clone_queue)

    def launch(self, sleep=None, update_jids=True):
        still_running = []  # note, this will be flat
        if update_jids:
            self.update_jids()
        for clone_list in self.priority_ordered_clones:
            # Record one False for a fully emptied clone-list
            if not clone_list:
                print('not clonelist-triggered')
                still_running.append(False)
            else:
                clone_indexes_to_remove = []
                for i, clone in enumerate(clone_list):
                    print('starting into clone loop for clone index',
                          i, clone.get_tag())
                    if sleep:
                        time.sleep(sleep)
                    # This probably shouldn't happen, but it's worth checking for
                    if clone in self.finished_clones or \
                            clone in self.failed_clone_set:
                        print('clone is finished clones or failed clones')
                        still_running.append(False)
                        clone_indexes_to_remove.append(i)
                    elif self.check_mark_clone_finished(clone):
                        print('clone was just marked finished')
                        still_running.append(False)
                        clone_indexes_to_remove.append(i)
                    # If clone is in active set, it may have just finished a generation.
                    elif clone in self.active_clone_set:
                        print('clone is in active clone list')
                        # Try to start another.
                        did_start = clone.check_start_gen(
                            self.current_jids, overwrite=self.overwrite)
                        if did_start:
                            still_running.append(True)
                        else:
                            self.mark_clone_failed(clone)
                            still_running.append(False)
                            clone_indexes_to_remove.append(i)

                    # This condition arises when there are few enough active clones
                    # that we could launch more.
                    elif len(self.active_clone_set) < self.active_clone_threshold:
                        print(
                            'there are some more active clones, let us launch', clone.get_tag())
                        #  So we try to launch another.
                        if clone.check_start_gen(self.current_jids, overwrite=self.overwrite):
                            print('started clone, adding to active_clone_set')
                            self.active_clone_set.add(clone)
                            still_running.append(True)
                        else:
                            self.mark_clone_failed(clone)
                            still_running.append(False)
                            clone_indexes_to_remove.append(i)
                    else:
                        print('WARNING:', clone.get_tag(),
                              'is not accounted for by launch logic.')
                # Because we are changing the length of the list, this must be done in reverse order
                clone_indexes_to_remove.reverse()
                for i in clone_indexes_to_remove:
                    del clone_list[i]
        return still_running

    # whether you're starting or restarting, this is probably what you want
    # if you'd just like a 'minder' process to start all your sims and keep them
    # running until they've gotten through all the generations.
    def start_tending_fields(self, update_interval=120):
        still_running = self.launch(sleep=None, update_jids=False)
        brake_file_p = Path('stop')
        # this needs to be while all(list of T/F for completed seeds/clones)
        print('still_running:', *still_running, flush=True)
        # If dry run, short circuit the tending loop.
        if self.dry_run:
            still_running = [False]
        # check on how the jobs are doing, see if you can launch more!
        # AKA the tending loop.
        while any(still_running):
            if brake_file_p.is_file():
                print(
                    f'Brake file detected: {brake_file_p.resolve()} Stopping submission loop.')
                return False
            time.sleep(update_interval)
            still_running = self.launch(sleep=None)
            if not self.quiet:
                print('The following (seed clone gen) are complete:', ', '.join((
                    map(lambda c: c.get_tag(), self.finished_clones))))
            print('STILL RUNNING:', *still_running, flush=True)
        # If we get here, that means the minder thinks all clones have finished.
        return True


# class Adaptive(Farmer):
#     __slots__ = ('current_gen', 'ranker', 'seeds', 'rank_jobscript')

#     def __init__(self, n_seeds, n_clones, n_gens, config_template, seed_structure_fns, scheduler, scheduler_fstring,
#                  scheduler_kws, scheduler_report_cmd, scheduler_assoc_rep_cmd, traj_list=None, quiet=False,
#                  active_clone_threshold=50, dirname_pad=3, job_number_re='[1-9][0-9]*', sep='-', seeds_first=True,
#                  job_name_elements=('{title}', '{seed_index}', '{clone_index}', '{gen_index}'), overwrite=False,
#                  runner=sims.omm_generation, current_gen=None):

#         super().__init__(n_seeds, n_clones, n_gens, config_template, seed_structure_fns, scheduler, scheduler_fstring,
#                          scheduler_kws, scheduler_report_cmd, scheduler_assoc_rep_cmd, traj_list, quiet,
#                          active_clone_threshold, dirname_pad, job_number_re, sep, seeds_first,
#                          job_name_elements, overwrite, runner)

#         if current_gen:
#             self.current_gen = current_gen
#         else:
#             current_youngest_gen = self.n_gens
#             for clone_list in self.priority_ordered_clones:
#                 for clone in clone_list:
#                     clone_gen = clone.config['gen_index']
#                     if current_youngest_gen > clone_gen:
#                         current_youngest_gen = clone_gen
#             self.current_gen = current_youngest_gen

#     def launch(self, sleep=None, update_jids=True):
#         retlist = []  # note, this will be flat
#         if update_jids:
#             self.update_jids()
#         for seed_index, clone_list in enumerate(self.priority_ordered_clones):
#             for clone_index, clone in enumerate(clone_list):
#                 if sleep:
#                     time.sleep(sleep)

#                 next_up_gen = clone.config['gen_index']
#                 # Check if the clone is up to date with the current gen of adaptive sampling.
#                 enough_gens = next_up_gen > self.current_gen
#                 retlist.append(enough_gens)
#                 if enough_gens:
#                     self.finished_clones[
#                         (seed_index, clone_index, next_up_gen)
#                     ] = clone
#                     clone_list.remove(clone)
#                     try:
#                         self.active_clone_set.remove(clone)
#                     except KeyError:
#                         print(
#                             f'done_before_launch: {seed_index},{clone_index},{next_up_gen}')
#                 elif clone in self.active_clone_set:
#                     clone.check_start_gen(
#                         self.current_jids, overwrite=self.overwrite
#                     )
#                 elif len(self.active_clone_set) < self.active_clone_threshold:
#                     self.active_clone_set.add(clone)
#                     clone.check_start_gen(
#                         self.current_jids, overwrite=self.overwrite
#                     )
#         # if there are no clones in active clone list, after adding some and removing others, the gen is complete
#         if len(self.active_clone_set) == 0:
#             self.current_gen += 1
#             self.rank_and_seed()

#         # Add the number of gens to be run as an additional condition that must be true for clone-minder to return.
#         retlist.append(self.current_gen > self.n_gens + 1)

#         return retlist
