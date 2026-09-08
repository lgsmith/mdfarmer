"""What the tender does at boot and on a tick, in the corners that bite.

Every Farmer here is a dry run over placeholder files, so no GROMACS and no
scheduler are needed: the decisions under test are the orchestrator's own.
"""
import contextlib
import io
import sys
from pathlib import Path

import harness
from harness import Suite

from mdfarmer import utilities as util
from mdfarmer import gmx_simulate as gs
from mdfarmer import farmer as fm

N_SEEDS = 1
N_CLONES = 2
N_GENS = 2
STEPS_PER_GEN = 1000
WRITE_INTERVAL = 100
CPUS = 8
PACK_SIZE = 2
# Clones per pack-threshold check, and how many clones a short campaign asks
# for when only half of them build.
PACK_THRESHOLD = 3
SHORT_CLONES = 4
FAILURE_LIMIT = 3


def make_template(work, steps_per_gen=STEPS_PER_GEN,
                  write_interval=WRITE_INTERVAL):
    return gs.gmx_config_template(
        traj_dir_top_level=str(work / 'farm'), title='fm',
        structure_fn=str(work / 'a.gro'), mdp_fn=str(work / 'base.mdp'),
        dirname_pad=2, sep='_', traj_name='prod', traj_suffix='.xtc',
        restart_name='state.cpt', steps=steps_per_gen,
        steps_per_gen=steps_per_gen, write_interval=write_interval,
        temperature=300, gen_seed_base=1, mdrun_args=['-nb', 'gpu'],
        traj_list=str(work / 'tl.txt'))


def make_farmer(work, template=None, n_seeds=N_SEEDS, n_clones=N_CLONES,
                n_gens=N_GENS, cpus=CPUS, cls=None, scheduler_fstring=None,
                **kwargs):
    return (fm.Farmer if cls is None else cls)(
        n_seeds=n_seeds, n_clones=n_clones, n_gens=n_gens,
        config_template=make_template(work) if template is None else template,
        seed_structure_fns=[str(work / 'a.gro')] * n_seeds,
        system_fns=[str(work / 'base.mdp')] * n_seeds,
        top_fns=[str(work / 'topol.top')] * n_seeds,
        scheduler='sbatch',
        scheduler_fstring=(util.basic_scheduler_fstrings['slurm']
                           if scheduler_fstring is None else scheduler_fstring),
        scheduler_kws=dict(gpu_line='', queue_name='gpu', exclude_nodes='',
                           cpus=cpus, run_script_name='run.py'),
        scheduler_report_cmd='true', scheduler_assoc_rep_cmd='true',
        sep='_', dirname_pad=2, runner=gs.gmx_generation, dry_run=True,
        overwrite=True, jids_file=work / 'jids.txt', **kwargs)


class StubClone:
    """Answers the handful of calls launch makes, and reports what it did.

    `succeeds` is what its check_start_gen returns, so a test can stage a run
    of failures and then a recovery. It carries the same restart budget a real
    Clone does, since that is what the tender now consults.
    """

    def __init__(self, tag='stub', current_gen=0, succeeds=False,
                 restarts_per_gen=FAILURE_LIMIT, raises=False):
        self.tag = tag
        self.current_gen = current_gen
        self.succeeds = succeeds
        self.raises = raises
        self.restarts_per_gen = restarts_per_gen
        self.restart_attempts = 0
        self.starts = 0

    def get_tag(self):
        return self.tag

    def spend_restart(self):
        if self.restart_attempts >= self.restarts_per_gen:
            return False
        self.restart_attempts += 1
        return True

    def check_start_gen(self, scheduler_report, overwrite=False):
        self.starts += 1
        if self.raises:
            raise RuntimeError('the filesystem stalled')
        if self.succeeds:
            # A real advance clears the generation's budget; so does this one.
            self.restart_attempts = 0
        return self.succeeds


def tend(farmer, clone):
    """Put one stub clone in the tender's care, as an already-running clone."""
    farmer.priority_ordered_clones = [[clone]]
    farmer.active_set = {clone}
    farmer.finished_clones = set()
    farmer.failed_clone_set = set()
    return farmer


class HalfBuildingFarmer(fm.Farmer):
    """A Farmer whose odd-numbered clones all fail to build."""

    __slots__ = ()

    def _setup_one_clone(self, tdir, seed_index, clone_index, rep_dict):
        if clone_index % 2:
            print(f'Skipping clone seed={seed_index} clone={clone_index}')
            return None
        return super()._setup_one_clone(tdir, seed_index, clone_index,
                                        rep_dict)


def captured(call):
    """Run call(), returning (result, everything it printed)."""
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        result = call()
    return result, out.getvalue()


def main(n_clones=N_CLONES, pack_size=PACK_SIZE,
         pack_threshold=PACK_THRESHOLD, short_clones=SHORT_CLONES,
         failure_limit=FAILURE_LIMIT):
    suite = Suite('farmer_findings')
    work = harness.workdir('farmer_findings')
    for name in ('a.gro', 'topol.top', 'base.mdp'):
        (work / name).write_text('placeholder\n')

    suite.section('the preempt trap is checked in the submitted template')
    try:
        captured(lambda: make_farmer(work, handle_preempt=True))
        suite.check('a solo template with no trap is refused', False,
                    '-> no exception')
    except ValueError as exc:
        suite.check('a solo template with no trap is refused',
                    'PREEMPT_SIGTERM' in str(exc), f'-> {str(exc)[:60]}')
    farmer, _ = captured(lambda: make_farmer(
        work, handle_preempt=True,
        scheduler_fstring=util.basic_scheduler_fstrings_preempt['slurm']))
    suite.check('a solo template with a trap boots and arms the reporter',
                farmer.config_template['handle_preempt'] is True)
    trapless = '\n'.join(
        line for line in util.basic_scheduler_fstrings_mps['slurm'].splitlines()
        if 'PREEMPT_SIGTERM' not in line and 'trap ' not in line)
    _, log = captured(lambda: make_farmer(
        work, pack_size=pack_size, pack_scheduler_fstring=trapless))
    suite.check('a pack template with no trap warns rather than refusing',
                'WARNING' in log and 'pack template' in log)
    _, log = captured(lambda: make_farmer(work, pack_size=pack_size))
    suite.check('the stock pack template passes the same check',
                'pack template' not in log)
    label = 'packing checks the pack template, not the unsubmitted one'
    try:
        _, log = captured(lambda: make_farmer(
            work, handle_preempt=True, pack_size=pack_size,
            scheduler_fstring=util.basic_scheduler_fstrings['slurm']))
        suite.check(label, 'pack template' not in log)
    except ValueError as exc:
        suite.check(label, False, f'-> refused over the solo one: {exc}'[:60])

    suite.section('packing changes what active_clone_threshold counts')
    _, log = captured(lambda: make_farmer(
        work, pack_size=pack_size, active_clone_threshold=pack_threshold))
    suite.check('boot says how many clones the threshold now allows',
                f'{pack_threshold * pack_size} clones will run at once' in log)

    suite.section('a generation that is not a whole number of write intervals')
    ragged = make_template(work, steps_per_gen=STEPS_PER_GEN + 1)
    try:
        captured(lambda: make_farmer(work, template=ragged))
        suite.check('a ragged generation length is refused at boot', False,
                    '-> no exception')
    except ValueError as exc:
        suite.check('a ragged generation length is refused at boot',
                    'write_interval' in str(exc), f'-> {str(exc)[:60]}')
    engineless = make_template(work)
    engineless.pop('write_interval')
    farmer, _ = captured(lambda: make_farmer(work, template=engineless))
    suite.check('a template without write_interval still boots',
                'write_interval' not in farmer.config_template)

    suite.section('traj_list is resolved however it arrives')
    relative = 'given-by-hand.txt'
    template = make_template(work)
    template.pop('traj_list')
    farmer, _ = captured(lambda: make_farmer(work, template=template,
                                             traj_list=relative))
    suite.check('a relative traj_list argument is resolved against the cwd',
                farmer.config_template['traj_list']
                == str(Path(relative).resolve()),
                f"-> {farmer.config_template['traj_list']}")
    template = make_template(work)
    template['traj_list'] = ''
    farmer, _ = captured(lambda: make_farmer(work, template=template))
    suite.check('an empty entry in the template falls back to traj_list.txt',
                farmer.config_template['traj_list']
                == str(Path('traj_list.txt').resolve()),
                f"-> {farmer.config_template['traj_list']}")

    suite.section('clones that could not be set up are counted, not just listed')
    farmer, log = captured(lambda: make_farmer(
        work, n_clones=short_clones, cls=HalfBuildingFarmer))
    built = sum(len(queue) for queue in farmer.priority_ordered_clones)
    suite.check('only half the clones were built',
                built == short_clones // 2, f'-> {built}')
    suite.check('boot says how many of how many are missing',
                f'{short_clones - built} of {short_clones} clones' in log)
    suite.check('the campaign still boots on what it has',
                any(farmer.priority_ordered_clones))
    _, log = captured(lambda: make_farmer(work))
    suite.check('a campaign that builds every clone says nothing',
                'could not be set up' not in log)

    suite.section('another campaign\'s job in the same queue')
    farmer, _ = captured(lambda: make_farmer(work))
    farmer.scheduler_assoc_rep_cmd = "printf '77 someone-elses-job\\n'"
    rep_dict, log = captured(farmer.reassociate_running_jobs)
    suite.check('it is not taken for one of ours',
                farmer.current_jids == set(), f'-> {farmer.current_jids}')
    suite.check('no clone is bound to it', rep_dict == {}, f'-> {rep_dict}')
    suite.check('and it is passed over quietly, not warned about every tick',
                'WARNING' not in log, f'-> {log.strip()[:60]}')

    suite.section('our own job with a name that does not parse')
    farmer, _ = captured(lambda: make_farmer(work))
    farmer.scheduler_assoc_rep_cmd = "printf '88 fm_0_0\\n'"
    rep_dict, log = captured(farmer.reassociate_running_jobs)
    suite.check('no clone is bound to it', rep_dict == {}, f'-> {rep_dict}')
    suite.check('boot names the job id and the name it could not parse',
                'WARNING' in log and '88' in log and 'fm_0_0' in log,
                f'-> {log.strip()[:70]}')

    suite.section('two queued jobs for one generation')
    farmer, _ = captured(lambda: make_farmer(work))
    farmer.scheduler_assoc_rep_cmd = "printf '11 fm_0_0_0\\n22 fm_0_0_0\\n'"
    rep_dict, log = captured(farmer.reassociate_running_jobs)
    suite.check('boot names both job ids and the generation they share',
                'WARNING' in log and '11' in log and '22' in log
                and '(0, 0, 0)' in log)
    suite.check('the tender keeps tending, bound to the newer job',
                rep_dict == {(0, 0, 0): 22}, f'-> {rep_dict}')
    suite.check('both jobs are still counted as ours',
                farmer.current_jids == {11, 22}, f'-> {farmer.current_jids}')

    suite.section('a clone waiting for a free slot')
    farmer, _ = captured(lambda: make_farmer(
        work, seeds_first=False, active_clone_threshold=1))
    still_running, log = captured(lambda: farmer.launch(update_jids=False))
    queue = farmer.priority_ordered_clones[0]
    suite.check('both clones are in one queue, so one has to wait',
                len(queue) == n_clones, f'-> {len(queue)}')
    suite.check('the waiting clone counts as still running',
                still_running == [True] * n_clones, f'-> {still_running}')
    suite.check('the waiting clone is not dropped from its queue',
                len(farmer.priority_ordered_clones[0]) == n_clones,
                f'-> {len(farmer.priority_ordered_clones[0])}')
    suite.check('the waiting clone is not marked failed',
                not farmer.failed_clone_set)
    suite.check('waiting is not reported as a launch-logic failure',
                'not accounted for' not in log)

    suite.section('a clone that says its budget is spent is failed at once')
    limit = failure_limit
    farmer, _ = captured(lambda: make_farmer(work))
    stub = StubClone('budget-spent')
    tend(farmer, stub)
    still_running, _ = captured(lambda: farmer.launch(update_jids=False))
    suite.check('the tender keeps no second budget of its own',
                still_running == [False] and stub in farmer.failed_clone_set,
                f'-> {still_running}')
    suite.check('it asked the clone exactly once', stub.starts == 1,
                f'-> {stub.starts}')
    suite.check('it did not charge the clone again on the way out',
                stub.restart_attempts == 0, f'-> {stub.restart_attempts}')

    farmer, _ = captured(lambda: make_farmer(work))
    stub = StubClone('never-launches')
    tend(farmer, stub)
    farmer.active_set = set()
    still_running, _ = captured(lambda: farmer.launch(update_jids=False))
    suite.check('a clone waiting for its first launch is failed the same way',
                still_running == [False] and stub in farmer.failed_clone_set,
                f'-> {still_running}')

    suite.section('a raise while advancing is charged to the restart budget')
    farmer, _ = captured(lambda: make_farmer(work))
    stub = StubClone('raises', restarts_per_gen=limit, raises=True)
    tend(farmer, stub)
    for attempt in range(1, limit + 1):
        still_running, _ = captured(lambda: farmer.launch(update_jids=False))
        suite.check(f'raise {attempt} of {limit} costs a restart, not the clone',
                    still_running == [True]
                    and farmer.priority_ordered_clones == [[stub]]
                    and stub.restart_attempts == attempt
                    and not farmer.failed_clone_set,
                    f'-> {still_running}, {stub.restart_attempts}')
    still_running, _ = captured(lambda: farmer.launch(update_jids=False))
    suite.check('the raise past the budget fails the clone',
                still_running == [False] and stub in farmer.failed_clone_set,
                f'-> {still_running}')

    suite.section('an advance that works clears the restart count')
    farmer, _ = captured(lambda: make_farmer(work))
    stub = StubClone('recovers', restarts_per_gen=limit, raises=True)
    tend(farmer, stub)
    captured(lambda: farmer.launch(update_jids=False))
    stub.raises = False
    stub.succeeds = True
    captured(lambda: farmer.launch(update_jids=False))
    suite.check('the count is cleared by the advance that worked',
                stub.restart_attempts == 0, f'-> {stub.restart_attempts}')
    stub.raises = True
    still_running, _ = captured(lambda: farmer.launch(update_jids=False))
    suite.check('so the next failure starts the count over',
                still_running == [True] and stub.restart_attempts == 1
                and not farmer.failed_clone_set,
                f'-> {still_running}, {stub.restart_attempts}')

    suite.section('a threshold that would leave every clone waiting')
    try:
        captured(lambda: make_farmer(work, active_clone_threshold=0))
        suite.check('active_clone_threshold=0 is refused', False,
                    '-> no exception')
    except ValueError as exc:
        suite.check('active_clone_threshold=0 is refused', True,
                    f'-> {str(exc)[:60]}')

    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
