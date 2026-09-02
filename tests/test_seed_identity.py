"""What a seed index means, and what happens when that changes.

`seed_labels` binds each index to a name and records it, so a reordered
seed_structure_fns refuses to boot rather than running one replica into
another's directory. Every Farmer here is a dry run over placeholder files:
no GROMACS, no scheduler.
"""
import contextlib
import io
import json
import sys

import harness
from harness import Suite

from mdfarmer import utilities as util
from mdfarmer import gmx_simulate as gs
from mdfarmer import farmer as fm

N_SEEDS = 3
N_CLONES = 1
N_GENS = 2
STEPS_PER_GEN = 1000
WRITE_INTERVAL = 100
CPUS = 8
LABELS = ['native', 'unfolded', 'misfolded']


def captured(call):
    """Run call(), returning (result, everything it printed)."""
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        result = call()
    return result, out.getvalue()


def raised(call):
    """Run call(), returning the exception it raised or None."""
    try:
        captured(call)
    except Exception as exc:
        return exc
    return None


def seed_inputs(work, n_seeds=N_SEEDS):
    """Three parallel n_seeds-long lists of placeholder structure/mdp/top."""
    work.mkdir(parents=True, exist_ok=True)
    structures, systems, tops = [], [], []
    for seed_index in range(n_seeds):
        names = [f'seed{seed_index}.gro', f'seed{seed_index}.mdp',
                 f'seed{seed_index}.top']
        for name, fns in zip(names, (structures, systems, tops)):
            path = work / name
            path.write_text('placeholder\n')
            fns.append(str(path))
    return structures, systems, tops


def make_template(work, farm):
    return gs.gmx_config_template(
        traj_dir_top_level=str(work / farm), title='fm',
        structure_fn=str(work / 'seed0.gro'), mdp_fn=str(work / 'seed0.mdp'),
        dirname_pad=2, sep='_', traj_name='prod', traj_suffix='.xtc',
        restart_name='state.cpt', steps=STEPS_PER_GEN,
        steps_per_gen=STEPS_PER_GEN, write_interval=WRITE_INTERVAL,
        temperature=300, gen_seed_base=1, mdrun_args=['-nb', 'gpu'],
        traj_list=str(work / 'tl.txt'))


def make_farmer(work, farm, fns, n_seeds=N_SEEDS, n_clones=N_CLONES, **kwargs):
    structures, systems, tops = fns
    return fm.Farmer(
        n_seeds=n_seeds, n_clones=n_clones, n_gens=N_GENS,
        config_template=make_template(work, farm),
        seed_structure_fns=structures, system_fns=systems, top_fns=tops,
        scheduler='sbatch',
        scheduler_fstring=util.basic_scheduler_fstrings['slurm'],
        scheduler_kws=dict(gpu_line='', queue_name='gpu', exclude_nodes='',
                           cpus=CPUS, run_script_name='run.py'),
        scheduler_report_cmd='true', scheduler_assoc_rep_cmd='true',
        sep='_', dirname_pad=2, runner=gs.gmx_generation, dry_run=True,
        overwrite=True, jids_file=work / 'jids.txt', **kwargs)


def seed_map(work, farm, seed_map_name=util.SEED_MAP_NAME):
    return work / farm / seed_map_name


def check_seed_labels(suite, work, labels=LABELS, n_seeds=N_SEEDS):
    suite.section('a seed label is recorded the first time it is given')
    fns = seed_inputs(work)
    farm = 'labelled'
    farmer, _ = captured(lambda: make_farmer(work, farm, fns,
                                             seed_labels=labels))
    suite.check('the record lands under its own name',
                seed_map(work, farm).is_file(),
                f'-> {sorted(p.name for p in (work / farm).iterdir())}')
    recorded = json.loads(seed_map(work, farm).read_text())
    suite.check('every index is written against its label',
                recorded == {str(i): l for i, l in enumerate(labels)},
                f'-> {recorded}')
    suite.check('the farmer keeps the labels it was given',
                farmer.seed_labels == labels)
    suite.check('nothing is left behind under the temp name',
                not list((work / farm).glob('seed_map.json.tmp')))

    suite.section('the same labels boot again')
    exc = raised(lambda: make_farmer(work, farm, fns, seed_labels=labels))
    suite.check('an unchanged campaign is not refused', exc is None,
                f'-> {exc}')

    suite.section('a reordered list is refused')
    swapped = [labels[1], labels[0]] + list(labels[2:])
    exc = raised(lambda: make_farmer(work, farm, fns, seed_labels=swapped))
    suite.check('a swap between two seeds refuses the boot',
                isinstance(exc, ValueError), f'-> {exc!r}'[:70])
    message = str(exc)
    suite.check('the message names both changed indices',
                'seed 0 was' in message and 'seed 1 was' in message,
                f'-> {message[:70]}')
    suite.check('the message names the old and the new label',
                repr(labels[0]) in message and repr(labels[1]) in message)
    suite.check('the message names the file to delete to force the re-index',
                str(seed_map(work, farm)) in message)
    unchanged = json.loads(seed_map(work, farm).read_text())
    suite.check('the record is left as it was',
                unchanged == {str(i): l for i, l in enumerate(labels)},
                f'-> {unchanged}')

    suite.section('a new seed appended past the end is accepted')
    grown = list(labels) + ['extra']
    grown_fns = seed_inputs(work, n_seeds=n_seeds + 1)
    exc = raised(lambda: make_farmer(work, farm, grown_fns,
                                     n_seeds=n_seeds + 1, seed_labels=grown))
    suite.check('growing the campaign is not refused', exc is None,
                f'-> {exc}')
    recorded = json.loads(seed_map(work, farm).read_text())
    suite.check('the new index joins the record',
                recorded == {str(i): l for i, l in enumerate(grown)},
                f'-> {recorded}')

    suite.section('a boot that runs fewer seeds than are recorded')
    exc = raised(lambda: make_farmer(work, farm, grown_fns, n_seeds=1,
                                     seed_labels=grown))
    suite.check('a shorter run is not refused', exc is None, f'-> {exc}')
    recorded = json.loads(seed_map(work, farm).read_text())
    suite.check('the seeds it did not run keep their labels',
                recorded == {str(i): l for i, l in enumerate(grown)},
                f'-> {recorded}')

    suite.section('labels that cannot tell two seeds apart')
    exc = raised(lambda: make_farmer(work, 'dupes', fns,
                                     seed_labels=[labels[0]] * n_seeds))
    suite.check('a repeated label is refused', isinstance(exc, ValueError),
                f'-> {exc!r}'[:70])
    suite.check('the message names the repeat',
                repr(labels[0]) in str(exc), f'-> {str(exc)[:70]}')

    suite.section('a label list of the wrong length')
    exc = raised(lambda: make_farmer(work, 'short', fns,
                                     seed_labels=labels[:1]))
    suite.check('too few labels are refused', isinstance(exc, ValueError),
                f'-> {exc!r}'[:70])
    suite.check('the message names seed_labels',
                'seed_labels' in str(exc), f'-> {str(exc)[:70]}')


def check_guard_is_optional(suite, work, labels=LABELS):
    suite.section('seed_labels stays optional')
    fns = seed_inputs(work)
    farm = 'unlabelled'
    farmer, log = captured(lambda: make_farmer(work, farm, fns))
    suite.check('a campaign without labels boots', farmer.seed_labels is None)
    suite.check('no record is written',
                not seed_map(work, farm).is_file())
    suite.check('and nothing is said about one',
                'seed_labels' not in log)

    farm = 'dropped'
    captured(lambda: make_farmer(work, farm, fns, seed_labels=labels))
    _, log = captured(lambda: make_farmer(work, farm, fns))
    suite.check('dropping the labels from a guarded campaign warns',
                'WARNING' in log and 'seed_labels' in log,
                f'-> {log.strip().splitlines()[-1][:60] if log.strip() else ""}')
    suite.check('the warning names the record it is not checking',
                str(seed_map(work, farm)) in log)


def check_guard_precedes_setup(suite, work, labels=LABELS):
    suite.section('the guard refuses before any seed directory is made')
    fns = seed_inputs(work)
    farm = 'ordering'
    captured(lambda: make_farmer(work, farm, fns, seed_labels=labels))
    fresh = work / 'ordering-fresh'
    template = make_template(work, 'ordering')
    template['traj_dir_top_level'] = str(fresh)
    util.write_seed_map(fresh / util.SEED_MAP_NAME,
                        dict(enumerate(reversed(labels))))
    exc = raised(lambda: fm.Farmer(
        n_seeds=len(labels), n_clones=N_CLONES, n_gens=N_GENS,
        config_template=template,
        seed_structure_fns=fns[0], system_fns=fns[1], top_fns=fns[2],
        scheduler='sbatch',
        scheduler_fstring=util.basic_scheduler_fstrings['slurm'],
        scheduler_kws=dict(gpu_line='', queue_name='gpu', exclude_nodes='',
                           cpus=CPUS, run_script_name='run.py'),
        scheduler_report_cmd='true', scheduler_assoc_rep_cmd='true',
        sep='_', dirname_pad=2, runner=gs.gmx_generation, dry_run=True,
        overwrite=True, jids_file=work / 'jids.txt', seed_labels=labels))
    suite.check('the re-indexed boot is refused', isinstance(exc, ValueError),
                f'-> {exc!r}'[:70])
    suite.check('no seed directory was made first',
                not list(fresh.glob('seed_[0-9]*')),
                f'-> {sorted(p.name for p in fresh.iterdir())}')




def main():
    suite = Suite('seed_identity')
    work = harness.workdir('seed_identity')
    check_seed_labels(suite, work)
    check_guard_is_optional(suite, work)
    check_guard_precedes_setup(suite, work)
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
