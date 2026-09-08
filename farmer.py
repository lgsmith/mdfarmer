import subprocess as sp
import traceback
from . import utilities as util
from pathlib import Path
from . import seeder
from . import simulate as sims
from . import gmx_simulate as gmx
from . import gmx_pack
import time

# Stands in for a seed index that runs off the end of one of the input lists.
MISSING_ENTRY = '<no entry>'

# Seconds between one submission and the next. 0 hands the scheduler the whole
# tick at once, which some of them mind and some do not.
LAUNCH_INTERVAL = 0

# Where every finished generation's trajectory path is appended, if the
# campaign does not name the file itself.
TRAJ_LIST_NAME = 'traj_list.txt'


def missing_seed_inputs(seed_structure_fns, system_fns, top_fns,
                        missing_entry=MISSING_ENTRY):
    """Which of each seed's three input files are not on disk yet.

    Maps seed index to the paths missing for it and omits a seed whose inputs
    are all there, so an empty mapping means the whole campaign can start. An
    index past the end of one of the lists has no path to name and is reported
    as missing_entry.
    """
    fn_lists = (list(seed_structure_fns), list(system_fns), list(top_fns))
    missing = {}
    for seed_index in range(max((len(fns) for fns in fn_lists), default=0)):
        absent = [missing_entry if seed_index >= len(fns)
                  else str(fns[seed_index])
                  for fns in fn_lists
                  if seed_index >= len(fns)
                  or not Path(fns[seed_index]).is_file()]
        if absent:
            missing[seed_index] = absent
    return missing


def ready_seed_count(seed_structure_fns, system_fns, top_fns):
    """How many seeds, counting from index 0, have every input in place.

    The leading run and not the total: a seed index is an on-disk identity, so
    a campaign grows by appending and a gap cannot be closed by renumbering
    the seeds past it. Pass this as Farmer(n_seeds=...) to run what is ready.
    """
    missing = missing_seed_inputs(seed_structure_fns, system_fns, top_fns)
    n_given = min(len(seed_structure_fns), len(system_fns), len(top_fns))
    return next((i for i in range(n_given) if i in missing), n_given)


def campaign_path(value, traj_dir_top_level, default_name):
    """An absolute path for value, anchored to the campaign directory.

    A relative value hangs off traj_dir_top_level rather than off whatever
    directory the driver script was started in, so the campaign names its own
    files and the whole tree can be moved to another cluster. An absolute one
    is left where it was put. The result is absolute either way: every
    generation runs from its own directory, so a relative path would give each
    of them a private file.
    """
    return str((Path(traj_dir_top_level) / (value or default_name)).resolve())


def campaign_and_recorded_jobs(report, title, sep, known_jids,
                               index_count=util.JOB_NAME_INDEX_COUNT):
    """A queue report split into this campaign's jobs and its orphaned ids.

    Returns (jobs, orphans). jobs is what campaign_jobs bound by name, as
    (job_id, indices) pairs. orphans are the queued ids in known_jids whose
    name does not split into a title and indices at all, so nothing can be
    bound to them: ours by record rather than by name. A name that does split
    belongs to whichever campaign its title names and is never orphaned here,
    so a reused id running somebody else's job is not claimed back.
    """
    jobs = util.campaign_jobs(report, title, sep=sep, index_count=index_count)
    unbindable = {jid for jid, name in util.parse_scheduler_report(report)
                  if util.split_job_name(name, sep=sep,
                                         index_count=index_count) is None}
    return jobs, (unbindable & set(known_jids)) - {jid for jid, _ in jobs}


def seed_slice(name, values, n_seeds):
    """The first n_seeds entries of a seed-indexed list, or None for None.

    A list longer than n_seeds is a campaign running a ready prefix; one
    shorter has no entry for a seed that will be built, which raises here
    rather than as an IndexError once clone setup reaches that seed.
    """
    if values is None:
        return None
    values = list(values)
    if len(values) < n_seeds:
        raise ValueError(
            f'{name} has {len(values)} entries for {n_seeds} seeds')
    return values[:n_seeds]


class Farmer:
    __slots__ = ('priority_ordered_clones', 'n_seeds', 'n_clones', 'n_gens', 'runner', 'jids_file',
                 'config_template', 'jn_regex', 'current_jids', 'dry_run', 'overwrite',
                 'active_clone_threshold', 'active_set', 'failed_clone_set', 'seeds_first',
                 'job_name_fstring', 'job_number_re', 'finished_clones', 'harvester',
                 'quiet', 'sep', 'dirname_pad', 'seed_state_fns', 'scheduler',
                 'scheduler_report_cmd', 'scheduler_fstring', 'scheduler_kws',
                 'scheduler_assoc_rep_cmd', 'system_fns', 'top_fns',
                 'node_blocklist', 'run_script', 'recover_fn', 'progress_fn',
                 'launch_interval', 'restarts_per_gen', 'pack_size',
                 'pack_grouping', 'pack_cpus_per_task', 'pack_scheduler_fstring',
                 'pack_run_script', 'pack_member_cores',
                 'seed_config_overrides', 'seed_labels')

    def write_jids(self):
        """Mirror the live job ids to jids_file, for the next tick to read."""
        self.jids_file.write_text(' '.join(map(str, sorted(self.current_jids))))

    def recorded_jids(self):
        """The job ids an earlier tick wrote to jids_file, or an empty set.

        A missing or unreadable file is not an error: a campaign's first boot
        has none, and the record is only ever a safety net over the job names.
        """
        try:
            text = self.jids_file.read_text()
        except OSError:
            return set()
        return {int(f) for f in text.split() if f.isascii() and f.isdigit()}

    def known_jids(self):
        """Every job id this campaign has reason to believe is its own.

        The job numbers its clones or packs are bound to, which a tender learns
        as it submits them, plus what an earlier tick recorded, which is all a
        fresh boot has to go on.
        """
        bound = {unit.job_number
                 for queue in self.priority_ordered_clones for unit in queue
                 if unit.job_number is not None}
        return bound | self.recorded_jids()

    def update_jids(self):
        """Refresh the set of job ids the scheduler says are ours and alive.

        False means the answer could not be trusted and current_jids was left
        alone. A query that failed prints nothing, and reading that as "no jobs
        running" would relaunch every live clone on top of itself.
        """
        trusted, report = util.scheduler_query(self.scheduler_report_cmd)
        if not trusted:
            print(f'WARNING: keeping the previous {len(self.current_jids)} job '
                  'ids and skipping this tick rather than relaunching live jobs.')
            return False
        jobs, orphans = campaign_and_recorded_jobs(
            report, self.config_template['title'], self.sep, self.known_jids())
        if orphans:
            print(f'NOTE: queued jobs {sorted(orphans)} carry names this '
                  'campaign cannot parse, but it recorded them as its own; '
                  'counting them as live rather than launching over them.')
        self.current_jids = {jid for jid, _ in jobs} | orphans
        self.write_jids()
        return True

    def check_path(self, p: Path):
        if p.is_file():
            return p
        else:
            raise FileNotFoundError(p)

    def check_seed_map(self, tdir, seed_map_name=util.SEED_MAP_NAME):
        """Check this boot's seed labels against the campaign's record of them.

        Without seed_labels nothing is bound, but a record left by an earlier
        boot means the guard was in use and is now off, which is worth saying.
        """
        seed_map_p = Path(tdir) / seed_map_name
        if self.seed_labels is None:
            if seed_map_p.is_file():
                print(f'WARNING: {seed_map_p} records what each seed index '
                      'means, but this boot passed no seed_labels, so a '
                      'reordered seed_structure_fns will not be caught.')
            return None
        return util.check_seed_map(seed_map_p, self.seed_labels)

    def mark_clone_failed(self, clone):
        """Move a clone out of active_set and onto the failed one.

        active_set holds whichever unit is being scheduled, so the argument is
        a Clone in an unpacked campaign and a ClonePack in a packed one.
        """
        self.failed_clone_set.add(clone)
        try:
            self.active_set.remove(clone)
        except KeyError:  # if clone isn't in active set that's OK.
            pass
        print('FAILED CLONE:', clone.get_tag())

    def check_mark_clone_finished(self, clone):
        """True when the clone has run all n_gens generations.

        A finished clone is moved out of active_set and into
        finished_clones. Either holds whichever unit is being scheduled: a
        Clone, or a ClonePack once packing.
        """
        next_up_gen = clone.current_gen
        enough_gens = next_up_gen >= self.n_gens
        if enough_gens:
            print('Finished:', clone.get_tag())
            self.finished_clones.add(clone)
            try:
                self.active_set.remove(clone)
            except KeyError:
                print('done_before_launch', clone.get_tag())
        return enough_gens

    def reassociate_running_jobs(self):
        """Map our live job ids back to the (seed, clone, gen) they belong to.

        One query answers both which jobs are alive and whose they are: two
        would let a job come or go in between, binding a dead id to a Clone.
        Raises if the scheduler cannot be reached, since booting blind would
        submit a second job into every live generation directory. Only a name
        carries the indices a Clone is bound by, so a recorded id queued under
        an unreadable one is counted as live and warned about, not bound.
        """
        self.current_jids = set()
        rep_dict = {}
        trusted, assoc_raw = util.scheduler_query(self.scheduler_assoc_rep_cmd)
        if not trusted:
            raise RuntimeError(
                'the scheduler could not be queried at boot, so which of this '
                "campaign's jobs are still running is unknown. Booting anyway "
                'would submit a second job into every live generation '
                'directory. Fix the query and start again.')
        print('boot re-association scheduler report:')
        print(assoc_raw)
        jobs, orphans = campaign_and_recorded_jobs(
            assoc_raw, self.config_template['title'], self.sep,
            self.known_jids())
        for jid, key in jobs:
            self.current_jids.add(jid)
            # The tender cannot cancel either job, so all it can do is say so.
            if key in rep_dict:
                print(f'WARNING: jobs {rep_dict[key]} and {jid} are both '
                      f'queued for seed/clone/gen {key}. Two jobs in one '
                      'generation directory will corrupt it -- cancel one by '
                      'hand.')
            rep_dict[key] = jid
        if orphans:
            print(f'WARNING: an earlier tick recorded jobs {sorted(orphans)} '
                  'as this campaign\'s and they are still queued, but their '
                  'names carry no seed/clone/gen, so no clone can be bound to '
                  'them. A second job may launch into their generation '
                  'directories -- cancel them by hand.')
        self.current_jids |= orphans
        self.write_jids()
        return rep_dict

    def pack_template(self):
        """The template a pack submits: the given one, else the MPS default."""
        family = util.scheduler_families.get(self.scheduler, self.scheduler)
        return (self.pack_scheduler_fstring
                or util.basic_scheduler_fstrings_mps[family])

    def check_preempt_template(self):
        """Check that the template about to be submitted traps SIGTERM.

        A preempted job is killed outright unless its submit script traps
        SIGTERM and touches PREEMPT_SIGTERM, which is what the simulation
        watches for so it can shut down on a whole frame. Packing submits the
        pack template and never the solo one, and a packed member watches the
        sentinel whether or not handle_preempt is set. A solo template without
        the trap raises; a pack template without it only warns.
        """
        packing = bool(self.pack_size or self.pack_grouping)
        if not packing and not self.config_template.get('handle_preempt'):
            return
        fstring = self.pack_template() if packing else self.scheduler_fstring
        if 'PREEMPT_SIGTERM' in fstring and 'trap' in fstring:
            return
        if packing:
            print('WARNING: the pack template has no SIGTERM trap that '
                  'touches PREEMPT_SIGTERM, so a preempted pack loses the '
                  'block every member is running. Use '
                  'basic_scheduler_fstrings_mps[<scheduler>] or add an '
                  'equivalent trap+background+wait pattern.')
            return
        raise ValueError(
            'handle_preempt is set (via Farmer(handle_preempt=True) or '
            'config_template["handle_preempt"]) but scheduler_fstring '
            'lacks a SIGTERM trap that touches PREEMPT_SIGTERM. Use '
            'basic_scheduler_fstrings_preempt[<scheduler>] or include '
            'an equivalent trap+background+wait pattern in your custom '
            'template.')

    def select_engine(self, runner, run_script, recover_fn, progress_fn):
        """Fill in the run_script, recover_fn and progress_fn matching runner.

        A hand-supplied set that disagrees with runner is refused rather than
        half applied: the engine that runs is the one named in the run script.
        """
        self.runner = runner
        self.run_script = run_script
        self.recover_fn = recover_fn
        self.progress_fn = progress_fn
        if self.runner is gmx.gmx_generation:
            gmx_defaults = (('run_script', gmx.default_gmx_run_script),
                            ('recover_fn', gmx.gmx_try_recover_gen),
                            ('progress_fn', gmx.gmx_gen_progress))
            for name, default in gmx_defaults:
                if getattr(self, name) is None:
                    setattr(self, name, default)
            custom = [name for name, default in gmx_defaults
                      if getattr(self, name) is not default]
            if custom:
                print(f'NOTE: runner=gmx_generation with custom {custom}; '
                      'make sure they implement the GROMACS contract.')
        elif self.runner is sims.omm_generation:
            gmx_pieces = (gmx.default_gmx_run_script, gmx.gmx_try_recover_gen,
                          gmx.gmx_gen_progress)
            for name in ('run_script', 'recover_fn', 'progress_fn'):
                if getattr(self, name) in gmx_pieces:
                    raise ValueError(
                        f'runner is the OpenMM omm_generation but {name} is the '
                        'GROMACS one. Pass runner=gmx_generation for a GROMACS '
                        'campaign; mixing them runs one engine with the '
                        "other's recovery logic.")

    def _setup_one_clone(self, tdir, seed_index, clone_index, rep_dict):
        """Build one Clone from its directory, or None if that clone alone fails.

        Isolating the failure keeps one corrupt clone directory from killing
        boot. A ConfigError is re-raised instead: it is the configuration, not
        this clone, and every clone it touches will hit it.
        """
        try:
            clone = seeder.Clone.from_disk(
                tdir, seed_index, clone_index,
                initial_seed_fn=self.seed_state_fns[seed_index],
                top_fn=self.top_fns[seed_index],
                system_fn=self.system_fns[seed_index],
                structure_fn=self.seed_state_fns[seed_index],
                config_overrides=(None if self.seed_config_overrides is None
                                  else self.seed_config_overrides[seed_index]),
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
                node_blocklist=self.node_blocklist,
                restarts_per_gen=self.restarts_per_gen,
                last_gen_index=self.n_gens - 1,
                rep_dict=rep_dict,
                run_script=self.run_script,
                recover_fn=self.recover_fn,
                progress_fn=self.progress_fn,
                dry_run=self.dry_run,
            )
        except seeder.ConfigError:
            raise
        except Exception as exc:
            print(f'Skipping clone seed={seed_index} clone={clone_index} '
                  f'during setup: {type(exc).__name__}: {exc}')
            return None
        if clone.job_number is not None:
            self.active_set.add(clone)
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
                 # Where finished trajectories are listed. Relative to
                 # traj_dir_top_level; None -> its 'traj_list.txt'.
                 traj_list=None,
                 quiet=False,
                 # Slots for running clones at once, or for packs once packing.
                 active_clone_threshold=50,
                 # Seconds to wait between one submission and the next.
                 launch_interval=LAUNCH_INTERVAL,
                 dirname_pad=3,
                 job_number_re='[1-9][0-9]*',
                 # Where live job ids are mirrored; None -> '<title>-jids.txt'.
                 jids_file=None,
                 # Dead launches a gen may take before its clone is dropped.
                 restarts_per_gen=3,
                 # MPS packing: this many clones share one job and one GPU.
                 pack_size=None,
                 # callable(clones) -> groups, choosing who packs with whom.
                 pack_grouping=None,
                 # A number, or callable(group) -> number for uneven packs.
                 pack_cpus_per_task=None,
                 pack_scheduler_fstring=None,
                 pack_run_script=None,
                 # A list, or callable(group) -> list, of cores per member.
                 pack_member_cores=None,
                 # One dict per seed, laid over config_template. n_seeds long.
                 seed_config_overrides=None,
                 # One label per seed, checked against the campaign's record of
                 # what each seed index means. n_seeds long; None checks nothing.
                 seed_labels=None,
                 sep='-',
                 seeds_first=True,
                 job_name_elements=(
                     '{title}', '{seed_index}', '{clone_index}',
                     '{gen_index}'),
                 overwrite=False,
                 harvester=None,
                 runner=sims.omm_generation,
                 # run.py body for each gen dir. None -> runner's own default.
                 run_script=None,
                 # Disk-recovery classifier. None -> runner's own default.
                 recover_fn=None,
                 # How a gen's progress is measured. None -> its frame count.
                 progress_fn=None,
                 dry_run=False,
                 # Shut down cleanly on preemption. Template must trap SIGTERM.
                 handle_preempt=False,
                 # Where learned bad nodes are written, and reloaded from.
                 bad_node_persist='bad_nodes.txt',
                 # Log substrings marking a node-local failure, or the default.
                 bad_node_patterns=None,
                 ):
        self.n_seeds = n_seeds
        # More inputs than seeds asked for is a campaign running a ready prefix.
        if len(seed_structure_fns) > n_seeds:
            print(f'NOTE: {len(seed_structure_fns)} seed structures given for '
                  f'n_seeds={n_seeds}; running the first {n_seeds} and leaving '
                  'the rest for a later boot with a larger n_seeds.')
        self.n_clones = n_clones
        self.restarts_per_gen = restarts_per_gen
        self.pack_size = pack_size
        self.pack_grouping = pack_grouping
        self.pack_cpus_per_task = pack_cpus_per_task
        self.pack_scheduler_fstring = pack_scheduler_fstring
        self.pack_run_script = pack_run_script
        self.pack_member_cores = pack_member_cores
        self.seed_config_overrides = seed_slice(
            'seed_config_overrides', seed_config_overrides, n_seeds)
        self.seed_labels = seed_slice('seed_labels', seed_labels, n_seeds)
        # gen-seed is base + stride * seed_index + clone_index.
        gen_seed_stride = config_template.get('gen_seed_stride',
                                              gmx.GEN_SEED_STRIDE)
        if n_seeds > 1 and n_clones > gen_seed_stride:
            raise ValueError(
                f'gen_seed_stride={gen_seed_stride} is not larger than '
                f'n_clones={n_clones}; seeds would share velocity seeds.')
        self.n_gens = n_gens
        self.overwrite = overwrite
        self.config_template = config_template
        # Resolve and check once here, so a bad file fails loudly at boot.
        self.system_fns = [str(self.check_path(Path(p)).resolve())
                           for p in seed_slice('system_fns', system_fns,
                                               n_seeds)]
        self.top_fns = [str(self.check_path(Path(p)).resolve())
                        for p in seed_slice('top_fns', top_fns, n_seeds)]
        self.select_engine(runner, run_script, recover_fn, progress_fn)
        self.seeds_first = seeds_first
        self.scheduler = scheduler
        self.scheduler_kws = scheduler_kws
        self.scheduler_fstring = scheduler_fstring
        # Before anything formats scheduler_fstring, so exclude_nodes is set.
        self.node_blocklist = util.BadNodeRegistry(
            bad_node_persist, scheduler, self.scheduler_kws,
            patterns=bad_node_patterns)
        # Either route to handle_preempt must reach the check below.
        if handle_preempt or self.config_template.get('handle_preempt'):
            self.config_template['handle_preempt'] = True
        self.check_preempt_template()
        self.scheduler_report_cmd = scheduler_report_cmd
        self.scheduler_assoc_rep_cmd = scheduler_assoc_rep_cmd
        self.job_number_re = job_number_re
        self.harvester = harvester
        self.dry_run = dry_run
        # Ensure that all file-names are in the config as full paths
        self.config_template['traj_dir_top_level'] = str(
            Path(self.config_template['traj_dir_top_level']).resolve()
        )
        # integrator_xml is OpenMM-only; GROMACS drivers omit it.
        if self.config_template.get('integrator_xml'):
            self.config_template['integrator_xml'] = str(
                self.check_path(Path(self.config_template['integrator_xml'])).resolve()
            )

        # Both keys are optional; not every engine's template carries them.
        util.check_whole_frames(self.config_template.get('steps'),
                                self.config_template.get('write_interval'),
                                source='config_template')

        # Whoever set this meant "flush per frame", already the default.
        if self.config_template.get('buffering') == 0:
            print('WARNING: config_template["buffering"] == 0 will make each '
                  'DCDFile.writeModel struct.pack a separate syscall and slow '
                  'trajectory writes dramatically. FlushingDCDReporter already '
                  'flushes the kernel buffer after every frame; remove the '
                  'buffering=0 entry.')

        self.seed_state_fns = [
            str(self.check_path(Path(s).resolve()))
            for s in seed_slice('seed_structure_fns', seed_structure_fns,
                                n_seeds)]
        # Before any directory is made, so a re-index refuses rather than runs.
        self.check_seed_map(self.config_template['traj_dir_top_level'])

        self.sep = sep
        self.config_template['sep'] = self.sep
        self.dirname_pad = dirname_pad
        self.config_template['dirname_pad'] = self.dirname_pad
        self.quiet = quiet
        self.jids_file = Path(
            jids_file if jids_file is not None
            else f'{config_template["title"]}-jids.txt')
        # important to pass this down through the clones
        self.job_name_fstring = self.sep.join(job_name_elements)
        self.current_jids = set()
        self.finished_clones = set()
        if active_clone_threshold < 1:
            raise ValueError(
                f'active_clone_threshold={active_clone_threshold} leaves every '
                'clone waiting for a slot that never opens. It must be at '
                'least 1.')
        self.active_clone_threshold = active_clone_threshold
        self.launch_interval = launch_interval
        # Whichever unit is being scheduled: Clones, or packs once packing.
        self.active_set = set()
        self.failed_clone_set = set()
        # Anchored to the campaign directory, so the tree can be moved whole.
        self.config_template['traj_list'] = campaign_path(
            self.config_template.get('traj_list') or traj_list,
            self.config_template['traj_dir_top_level'], TRAJ_LIST_NAME)

        # Empty before boot re-association, which asks it what ids are bound.
        self.priority_ordered_clones = []
        rep_dict = self.reassociate_running_jobs()

        tdir = Path(self.config_template['traj_dir_top_level'])
        if self.seeds_first:
            # One queue per clone index, holding every seed of it.
            queue_keys = [[(seed_index, clone_index)
                           for seed_index in range(self.n_seeds)]
                          for clone_index in range(self.n_clones)]
        else:
            # One queue per seed, holding every clone of it.
            queue_keys = [[(seed_index, clone_index)
                           for clone_index in range(self.n_clones)]
                          for seed_index in range(self.n_seeds)]
        for keys in queue_keys:
            clone_queue = []
            for seed_index, clone_index in keys:
                clone = self._setup_one_clone(
                    tdir, seed_index, clone_index, rep_dict)
                if clone is not None:
                    clone_queue.append(clone)
            self.priority_ordered_clones.append(clone_queue)
        # Count now: packing replaces the clones in these queues with packs.
        asked_for = self.n_seeds * self.n_clones
        built = sum(len(queue) for queue in self.priority_ordered_clones)
        if built < asked_for:
            print(f'WARNING: {asked_for - built} of {asked_for} clones could '
                  'not be set up; this campaign will be that much smaller '
                  'than asked for.')
        if self.pack_size or self.pack_grouping:
            self.build_packs(tdir)

    def group_clones(self, clones):
        """Partition every built Clone into pack-sized groups.

        pack_grouping supplies the policy; without one the flat priority
        order is cut into consecutive runs of pack_size. Either way every
        clone must land in exactly one group. A clone left out of
        the plan would never be submitted, and one in two packs would get two
        jobs in its generation directory.
        """
        if self.pack_grouping is not None:
            groups = [list(g) for g in self.pack_grouping(clones)]
        else:
            groups = [clones[i:i + self.pack_size]
                      for i in range(0, len(clones), self.pack_size)]

        def key(clone):
            return (clone.config['seed_index'], clone.config['clone_index'])

        packed = [key(c) for group in groups for c in group]
        duplicated = sorted({k for k in packed if packed.count(k) > 1})
        if duplicated:
            raise ValueError(f'seed/clone {duplicated} appear in two packs')
        missing = sorted(set(map(key, clones)) - set(packed))
        if missing:
            raise ValueError(f'clones built but never packed: {missing}')
        if self.pack_size:
            wrong = [len(g) for g in groups if len(g) != self.pack_size]
            if wrong:
                raise ValueError(
                    f'pack_size={self.pack_size} but groups of size {wrong} '
                    'were produced')
        return groups

    def build_packs(self, tdir):
        """Replace the clone queues with ClonePacks, one queue per pack.

        A pack answers every call launch makes on a Clone, so the tending
        loop is unchanged. active_set is rebuilt over the packs because
        _setup_one_clone populated it with the individual Clones.
        """
        fstring = self.pack_template()
        run_script = self.pack_run_script or gmx_pack.default_gmx_pack_run_script
        cpus = self.pack_cpus_per_task or self.scheduler_kws.get('cpus')
        if not cpus:
            raise ValueError(
                'packing needs pack_cpus_per_task, or a "cpus" entry in '
                'scheduler_kws for the pack template to fill in')
        flat = [c for queue in self.priority_ordered_clones for c in queue]
        packs = []
        for group in self.group_clones(flat):
            # A callable lets unlike packs ask for unlike amounts of core.
            group_cpus = cpus(group) if callable(cpus) else cpus
            tag = self.sep.join(
                f's{c.config["seed_index"]:0{self.dirname_pad}d}'
                f'c{c.config["clone_index"]:0{self.dirname_pad}d}'
                for c in group)
            member_cores = self.pack_member_cores
            if callable(member_cores):
                member_cores = member_cores(group)
            packs.append(seeder.ClonePack(
                group, Path(tdir) / 'packs' / f'pack{self.sep}{tag}',
                self.scheduler, fstring, self.scheduler_kws,
                run_script=run_script, cpus_per_task=group_cpus,
                member_cores=member_cores, sep=self.sep,
                job_number_re=self.job_number_re,
                dry_run=self.dry_run))
        self.priority_ordered_clones = [[pack] for pack in packs]
        self.active_set = {pack for pack in packs
                           if pack.job_number is not None}
        if packs:
            biggest = max(len(pack.clones) for pack in packs)
            print(f'NOTE: {len(packs)} packs of up to {biggest} clones. '
                  f'active_clone_threshold={self.active_clone_threshold} '
                  'counts what active_set holds, which is now packs rather '
                  f'than clones, so up to {self.active_clone_threshold * biggest} '
                  'clones will run at once.')
        return packs

    def _safe_check_start_gen(self, clone):
        """check_start_gen, with any exception charged to the restart budget.

        It touches the filesystem, the scheduler and the engine binary, any of
        which can raise, and one raise must not kill a tender that has been
        minding a campaign for weeks. A raise leaves the attempt uncounted, so
        it is counted here.
        """
        try:
            return clone.check_start_gen(
                self.current_jids, overwrite=self.overwrite)
        except Exception as exc:
            print(f'ERROR advancing clone {clone.get_tag()}: '
                  f'{type(exc).__name__}: {exc}')
            traceback.print_exc()
            return clone.spend_restart()

    def advance_clone(self, clone):
        """Try to move a clone on by one tick; False once it has been failed.

        Every failure to advance is charged to the clone's own restart budget,
        so a False here means that budget is spent and there is nothing left
        to retry.
        """
        if self._safe_check_start_gen(clone):
            return True
        self.mark_clone_failed(clone)
        return False

    def resolve_launch_interval(self, sleep):
        """The campaign's launch_interval, or a deprecated sleep= overriding it.

        sleep was the per-call name for the same wait; it still works so that
        existing driver scripts keep running, and says so once.
        """
        if sleep is None:
            return self.launch_interval
        print(f'WARNING: sleep={sleep} is deprecated. Pass '
              f'Farmer(launch_interval={sleep}) instead; honouring sleep for '
              'now.')
        return sleep

    def launch(self, sleep=None, update_jids=True):
        """One tick: advance, finish or fail every clone the queues still hold.

        Returns one flat True/False per clone looked at, True meaning it is
        still part of the campaign. A clone that is waiting for a free slot
        counts as running, not as a failure.
        """
        launch_interval = self.resolve_launch_interval(sleep)
        still_running = []  # note, this will be flat
        if update_jids and not self.update_jids():
            # Scheduler unreachable; presume all running, retry next tick.
            return [True] * max(1, sum(len(cl) for cl in
                                       self.priority_ordered_clones))
        for queue_index, clone_list in enumerate(self.priority_ordered_clones):
            # Record one False for a fully emptied clone-list
            if not clone_list:
                print('not clonelist-triggered')
                still_running.append(False)
                continue
            # What this queue keeps for the next tick; the rest are dropped.
            survivors = []
            for clone in clone_list:
                print('starting into clone loop for', clone.get_tag())
                # Spaces the submissions out, for a scheduler that minds a
                # whole campaign landing at once.
                if launch_interval:
                    time.sleep(launch_interval)
                # This probably shouldn't happen, but it's worth checking for
                if clone in self.finished_clones or \
                        clone in self.failed_clone_set:
                    print('clone is finished clones or failed clones')
                    still_running.append(False)
                elif self.check_mark_clone_finished(clone):
                    print('clone was just marked finished')
                    still_running.append(False)
                # In the active set, so it may have just finished a gen.
                elif clone in self.active_set:
                    print('clone is in active clone list')
                    # Try to start another.
                    if self.advance_clone(clone):
                        survivors.append(clone)
                        still_running.append(True)
                    else:
                        still_running.append(False)

                # Few enough active clones that we could launch another.
                elif len(self.active_set) < self.active_clone_threshold:
                    print(
                        'there are some more active clones, let us launch', clone.get_tag())
                    #  So we try to launch another.
                    if self.advance_clone(clone):
                        print('started clone, adding to active_set')
                        self.active_set.add(clone)
                        survivors.append(clone)
                        still_running.append(True)
                    else:
                        still_running.append(False)
                # Every slot is taken, so this clone waits its turn.
                else:
                    if not self.quiet:
                        print(clone.get_tag(),
                              'is waiting for a free slot.')
                    survivors.append(clone)
                    still_running.append(True)
            self.priority_ordered_clones[queue_index] = survivors
        return still_running

    def start_tending_fields(self, update_interval=120, sleep=None):
        """Mind the whole campaign: launch every clone and keep it going.

        Whether you are starting or restarting, this is probably what you
        want. Returns True once every clone has finished, and False if a
        'stop' brake file halted the loop or any clone was given up on.
        Raises if no clone could be built at all. sleep is the deprecated
        spelling of Farmer(launch_interval=...) and overrides it for this run.
        """
        self.launch_interval = self.resolve_launch_interval(sleep)
        if not self.priority_ordered_clones or not any(
                self.priority_ordered_clones):
            # All of them failing is not the campaign finishing.
            raise RuntimeError(
                'No clones could be set up; nothing to tend. Check the '
                'per-clone setup errors printed above (missing structure, '
                'topology, .mdp, or an unreadable checkpoint).')
        still_running = self.launch(update_jids=False)
        brake_file_p = Path('stop')
        print('still_running:', *still_running, flush=True)
        # If dry run, short circuit the tending loop.
        if self.dry_run:
            still_running = [False]
        # The tending loop: see how the jobs do, launch where there is room.
        while any(still_running):
            if brake_file_p.is_file():
                print(
                    f'Brake file detected: {brake_file_p.resolve()} Stopping submission loop.')
                return False
            time.sleep(update_interval)
            # Losing the tender leaves every running job unminded.
            try:
                still_running = self.launch()
            except Exception as exc:
                print(f'ERROR in tending loop: {type(exc).__name__}: {exc}')
                traceback.print_exc()
                print('Continuing; will retry next tick.', flush=True)
                still_running = [True]
                continue
            if not self.quiet:
                print('The following (seed clone gen) are complete:', ', '.join((
                    map(lambda c: c.get_tag(), self.finished_clones))))
            print('STILL RUNNING:', *still_running, flush=True)
        # They have all stopped; only good news if they stopped by finishing.
        if self.failed_clone_set:
            print(f'WARNING: {len(self.failed_clone_set)} clone(s) failed and '
                  f'{len(self.finished_clones)} finished: '
                  + ', '.join(c.get_tag() for c in self.failed_clone_set))
        return not self.failed_clone_set
