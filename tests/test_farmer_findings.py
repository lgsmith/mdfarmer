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
                n_gens=N_GENS, cpus=CPUS, **kwargs):
    return fm.Farmer(
        n_seeds=n_seeds, n_clones=n_clones, n_gens=n_gens,
        config_template=make_template(work) if template is None else template,
        seed_structure_fns=[str(work / 'a.gro')] * n_seeds,
        system_fns=[str(work / 'base.mdp')] * n_seeds,
        top_fns=[str(work / 'topol.top')] * n_seeds,
        scheduler='sbatch',
        scheduler_fstring=util.basic_scheduler_fstrings['slurm'],
        scheduler_kws=dict(gpu_line='', queue_name='gpu', exclude_nodes='',
                           cpus=cpus, run_script_name='run.py'),
        scheduler_report_cmd='true', scheduler_assoc_rep_cmd='true',
        sep='_', dirname_pad=2, runner=gs.gmx_generation, dry_run=True,
        overwrite=True, jids_file=work / 'jids.txt', **kwargs)


def captured(call):
    """Run call(), returning (result, everything it printed)."""
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        result = call()
    return result, out.getvalue()


def main(n_clones=N_CLONES):
    suite = Suite('farmer_findings')
    work = harness.workdir('farmer_findings')
    for name in ('a.gro', 'topol.top', 'base.mdp'):
        (work / name).write_text('placeholder\n')

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
