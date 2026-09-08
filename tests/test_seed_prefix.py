"""How much of a half-built campaign can start today.

`missing_seed_inputs` and `ready_seed_count` say how far down the seed lists
the inputs actually exist, and a Farmer asked for fewer seeds than the lists
hold runs that ready prefix. Every Farmer here is a dry run over placeholder
files: no GROMACS, no scheduler.
"""
import contextlib
import io
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
# A campaign whose inputs stop part way down the lists.
READY_SEEDS = 2


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


def seed_inputs(work, n_seeds=N_SEEDS, present=None):
    """Three parallel n_seeds-long lists, with only `present` of them on disk."""
    present = n_seeds if present is None else present
    work.mkdir(parents=True, exist_ok=True)
    structures, systems, tops = [], [], []
    for seed_index in range(n_seeds):
        names = [f'seed{seed_index}.gro', f'seed{seed_index}.mdp',
                 f'seed{seed_index}.top']
        for name, fns in zip(names, (structures, systems, tops)):
            path = work / name
            if seed_index < present:
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





def check_missing_inputs(suite, work, n_seeds=N_SEEDS, ready=READY_SEEDS):
    suite.section('which seeds are missing what')
    structures, systems, tops = seed_inputs(work / 'preflight',
                                            n_seeds=n_seeds, present=ready)
    missing = fm.missing_seed_inputs(structures, systems, tops)
    suite.check('only the unbuilt seeds are listed',
                sorted(missing) == list(range(ready, n_seeds)),
                f'-> {sorted(missing)}')
    suite.check('all three of an unbuilt seed\'s files are named',
                sorted(missing[ready]) == sorted(
                    [structures[ready], systems[ready], tops[ready]]),
                f'-> {missing[ready]}')
    suite.check('a fully built campaign reports nothing',
                fm.missing_seed_inputs(structures[:ready], systems[:ready],
                                       tops[:ready]) == {})
    suite.check('no lists at all is not an error',
                fm.missing_seed_inputs([], [], []) == {})

    suite.section('a list shorter than the others')
    missing = fm.missing_seed_inputs(structures[:ready], systems[:ready],
                                     tops[:ready - 1])
    suite.check('the seed with no entry is reported missing',
                list(missing) == [ready - 1], f'-> {list(missing)}')
    suite.check('the absent entry is named rather than indexed',
                missing[ready - 1] == [fm.MISSING_ENTRY],
                f'-> {missing[ready - 1]}')

    suite.section('how many seeds are ready to run')
    suite.check('the leading run of built seeds is counted',
                fm.ready_seed_count(structures, systems, tops) == ready,
                f'-> {fm.ready_seed_count(structures, systems, tops)}')
    suite.check('a fully built campaign is entirely ready',
                fm.ready_seed_count(structures[:ready], systems[:ready],
                                    tops[:ready]) == ready)
    gapped = [str(work / 'preflight' / 'nope.gro')] + structures[1:]
    suite.check('a gap at seed 0 leaves nothing ready',
                fm.ready_seed_count(gapped, systems, tops) == 0)
    late = list(structures[:ready])
    late[-1] = str(work / 'preflight' / 'nope.gro')
    suite.check('seeds past a gap are not counted',
                fm.ready_seed_count(late + structures[ready:], systems,
                                    tops) == ready - 1)


def check_prefix_boot(suite, work, n_seeds=N_SEEDS, ready=READY_SEEDS,
                      n_clones=N_CLONES):
    suite.section('fail-fast is still what a full campaign gets')
    fns = seed_inputs(work / 'partial', n_seeds=n_seeds, present=ready)
    exc = raised(lambda: make_farmer(work, 'partial-farm', fns,
                                     n_seeds=n_seeds))
    suite.check('a campaign asking for every seed raises on the missing one',
                isinstance(exc, FileNotFoundError), f'-> {exc!r}'[:70])
    suite.check('the raise names one of that seed\'s missing files',
                f'seed{ready}.' in str(exc), f'-> {exc}'[:70])

    suite.section('a ready prefix is opt-in, by asking for fewer seeds')
    exc = raised(lambda: make_farmer(work, 'prefix-farm', fns, n_seeds=ready))
    suite.check('asking for fewer seeds does not raise over the unbuilt ones',
                exc is None, f'-> {exc!r}'[:70])
    farmer, log = captured(lambda: make_farmer(work, 'prefix-farm', fns,
                                               n_seeds=ready))
    suite.check('the prefix boots over the seeds that exist',
                len(farmer.seed_state_fns) == ready,
                f'-> {len(farmer.seed_state_fns)}')
    suite.check('the unbuilt seeds are not resolved',
                len(farmer.system_fns) == ready == len(farmer.top_fns))
    suite.check('boot says it is running a prefix',
                'NOTE' in log and f'running the first {ready}' in log)
    built = sum(len(queue) for queue in farmer.priority_ordered_clones)
    suite.check('only the prefix\'s clones are built',
                built == ready * n_clones, f'-> {built}')

    suite.section('a seed inside the prefix that is still missing')
    holed = ([str(work / 'partial' / 'nope.gro')] + fns[0][1:], fns[1], fns[2])
    exc = raised(lambda: make_farmer(work, 'holed-farm', holed, n_seeds=ready))
    suite.check('a prefix does not excuse a gap inside it',
                isinstance(exc, FileNotFoundError), f'-> {exc!r}'[:70])

    suite.section('seed-indexed lists too short for n_seeds')
    exc = raised(lambda: make_farmer(work, 'short-farm',
                                     (fns[0][:1], fns[1], fns[2]),
                                     n_seeds=ready))
    suite.check('a short seed_structure_fns is named, not left to IndexError',
                isinstance(exc, ValueError) and 'seed_structure_fns' in str(exc),
                f'-> {exc!r}'[:70])
    exc = raised(lambda: make_farmer(work, 'short-ovr', fns, n_seeds=ready,
                                     seed_config_overrides=[{}]))
    suite.check('a short seed_config_overrides is still refused',
                isinstance(exc, ValueError)
                and 'seed_config_overrides' in str(exc),
                f'-> {exc!r}'[:70])
    farmer, _ = captured(lambda: make_farmer(
        work, 'long-ovr', fns, n_seeds=ready,
        seed_config_overrides=[{'mdrun_args': ['-nb', 'cpu']}] * n_seeds))
    suite.check('a longer seed_config_overrides is cut to the prefix',
                len(farmer.seed_config_overrides) == ready,
                f'-> {len(farmer.seed_config_overrides)}')


def main():
    suite = Suite('seed_prefix')
    work = harness.workdir('seed_prefix')
    check_missing_inputs(suite, work)
    check_prefix_boot(suite, work)
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
