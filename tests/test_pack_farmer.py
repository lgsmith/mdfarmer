"""Farmer builds ClonePacks, and the per-seed differences a pack can hold.

Needs no GROMACS: every launch is a dry run, so the referenced files only have
to exist.
"""
import io
import json
import sys
from contextlib import redirect_stdout
from pathlib import Path

import harness
from harness import Suite

from mdfarmer import utilities as util
from mdfarmer import gmx_pack as gp
from mdfarmer import gmx_simulate as gs
from mdfarmer import farmer as fm
from mdfarmer.seeder import ClonePack

CPUS = 16
N_SEEDS = 2
N_CLONES = 2
PACK_SIZE = 2
STEPS_PER_GEN = 1000
WRITE_INTERVAL = 100
MEMBER_CORES = [12, 4]


def make_farmer(work, n_seeds=N_SEEDS, n_clones=N_CLONES, cpus=CPUS,
                steps_per_gen=STEPS_PER_GEN, write_interval=WRITE_INTERVAL,
                **kwargs):
    template = gs.gmx_config_template(
        traj_dir_top_level=str(work / 'farm'), title='pk',
        structure_fn=str(work / 'a.gro'), mdp_fn=str(work / 'base.mdp'),
        dirname_pad=2, sep='_', traj_name='prod', traj_suffix='.xtc',
        restart_name='state.cpt', steps=steps_per_gen,
        steps_per_gen=steps_per_gen, write_interval=write_interval,
        temperature=300, gen_seed_base=1, mdrun_args=['-nb', 'gpu'],
        traj_list=str(work / 'tl.txt'))
    return fm.Farmer(
        n_seeds=n_seeds, n_clones=n_clones, n_gens=2, config_template=template,
        seed_structure_fns=[str(work / 'a.gro'), str(work / 'b.gro')][:n_seeds],
        system_fns=[str(work / 'base.mdp')] * n_seeds,
        top_fns=[str(work / 'topol.top')] * n_seeds,
        scheduler='sbatch',
        scheduler_fstring=util.basic_scheduler_fstrings['slurm'],
        scheduler_kws=dict(gpu_line='', queue_name='gpu', exclude_nodes='',
                           cpus=cpus, run_script_name='run.py'),
        scheduler_report_cmd='true', scheduler_assoc_rep_cmd='true',
        sep='_', dirname_pad=2, runner=gs.gmx_generation, dry_run=True,
        overwrite=True, jids_file=work / 'jids.txt', **kwargs)


def packs_of(farmer):
    return [pack for queue in farmer.priority_ordered_clones for pack in queue]


def member_keys(pack):
    return [(c.config['seed_index'], c.config['clone_index'])
            for c in pack.clones]


def main(cpus=CPUS, pack_size=PACK_SIZE, member_cores=MEMBER_CORES,
         n_seeds=N_SEEDS, n_clones=N_CLONES, write_interval=WRITE_INTERVAL):
    suite = Suite('pack_farmer')
    work = harness.workdir('pack_farmer')
    for name in ('a.gro', 'b.gro', 'topol.top', 'base.mdp'):
        (work / name).write_text('placeholder\n')

    suite.section('the default grouping cuts the priority order into runs')
    farmer = make_farmer(work, pack_size=pack_size)
    packs = packs_of(farmer)
    suite.check('the Farmer built packs, not clones',
                all(isinstance(p, ClonePack) for p in packs))
    suite.check(f'{n_seeds * n_clones} clones became '
                f'{n_seeds * n_clones // pack_size} packs',
                len(packs) == n_seeds * n_clones // pack_size,
                f'-> {len(packs)}')
    suite.check('every clone is packed exactly once',
                sorted(k for p in packs for k in member_keys(p))
                == [(s, c) for s in range(n_seeds) for c in range(n_clones)])
    suite.check('a pack answers the Clone interface launch uses',
                all(hasattr(p, name) for p in packs for name in
                    ('check_start_gen', 'current_gen', 'get_tag', '__hash__')))

    suite.section('a custom grouping and an uneven core split')

    def pair_seeds(clones):
        by_key = {(c.config['seed_index'], c.config['clone_index']): c
                  for c in clones}
        return [[by_key[(0, i)], by_key[(1, i)]] for i in range(n_clones)]

    farmer = make_farmer(work, pack_size=pack_size, pack_grouping=pair_seeds,
                         pack_member_cores=member_cores)
    packs = packs_of(farmer)
    suite.check('each pack spans both seeds',
                all({c.config['seed_index'] for c in p.clones} == {0, 1}
                    for p in packs))
    suite.check('member_cores reaches the pack',
                all(p.member_cores == member_cores for p in packs))
    packs[0].check_start_gen(set(), overwrite=True)
    manifest = json.loads((packs[0].pack_dir / 'pack.json').read_text())
    suite.check('member_cores reaches the manifest',
                manifest['member_cores'] == member_cores)
    layout = gp.member_core_layout(cpus, pack_size,
                                   member_cores=manifest['member_cores'])
    suite.check('the layout is contiguous, non-overlapping and fits',
                layout == [(12, 0), (4, 12)], f'-> {layout}')
    script = (packs[0].pack_dir / 'sbatch.sh').read_text()
    suite.check('the pack submit script is the MPS template',
                'nvidia-cuda-mps-control' in script)
    suite.check('the submit script holds a pack lock',
                'flock -n 9' in script and 'exec 9>pack.lock' in script)
    name = packs[0]._job_name()
    suite.check('the pack job name parses for re-association',
                all(field.isdigit() for field in name.split('_')[-3:]),
                f'-> {name}')

    suite.check('each seed grompps generation 0 from its own structure',
                {s: Path(c['structure_fn']).name
                 for s, c in {c.config['seed_index']: c.config
                              for p in packs for c in p.clones}.items()}
                == {0: 'a.gro', 1: 'b.gro'})
    clone = packs[0].clones[0]
    clone.job_number = 4242
    adopted = ClonePack(
        packs[0].clones, work / 'pack-adopt', 'sbatch',
        util.basic_scheduler_fstrings_mps['slurm'],
        dict(gpu_line='', queue_name='gpu', exclude_nodes='', cpus=cpus,
             run_script_name='run.py'),
        run_script=gp.default_gmx_pack_run_script, cpus_per_task=cpus,
        dry_run=True)
    suite.check("a pack adopts a member's re-associated job number",
                adopted.job_number == 4242, f'-> {adopted.job_number}')

    suite.section('a grouping that loses or repeats a clone is refused')
    for label, grouping in (('drops a clone', lambda cs: [cs[:pack_size]]),
                            ('repeats a clone',
                             lambda cs: [cs[:pack_size], cs[:pack_size]])):
        try:
            make_farmer(work, pack_size=pack_size, pack_grouping=grouping)
            suite.check(f'refuses a grouping that {label}', False,
                        '-> no exception')
        except ValueError as exc:
            suite.check(f'refuses a grouping that {label}', True,
                        f'-> {str(exc)[:50]}')

    suite.section('packs of different conditions get different core budgets')

    def cpus_for(group):
        return 24 if group[0].config['seed_index'] == 0 else 8

    def cores_for(group):
        return [12, 12] if group[0].config['seed_index'] == 0 else [4, 4]

    def by_seed_pairs(clones):
        by_key = {(c.config['seed_index'], c.config['clone_index']): c
                  for c in clones}
        return [[by_key[(s, 0)], by_key[(s, 1)]] for s in range(n_seeds)]

    farmer = make_farmer(work, pack_size=pack_size, pack_grouping=by_seed_pairs,
                         pack_cpus_per_task=cpus_for,
                         pack_member_cores=cores_for)
    budgets = {p.clones[0].config['seed_index']:
               (p.cpus_per_task, p.member_cores) for p in packs_of(farmer)}
    print('   seed -> (cpus, member_cores):', budgets, flush=True)
    suite.check('each pack asks for its own core budget',
                budgets == {0: (24, [12, 12]), 1: (8, [4, 4])}, f'-> {budgets}')
    arm_b = [p for p in packs_of(farmer)
             if p.clones[0].config['seed_index'] == 1][0]
    arm_b.check_start_gen(set(), overwrite=True)
    script = (arm_b.pack_dir / 'sbatch.sh').read_text()
    suite.check('the submit script asks for that many, not the global figure',
                '--cpus-per-task=8' in script and '--cpus-per-task=24' not in script)

    suite.section('per-seed config overrides')
    farmer = make_farmer(
        work, pack_size=pack_size, pack_grouping=pair_seeds,
        seed_config_overrides=[
            dict(mdrun_args=['-nb', 'gpu', '-update', 'cpu']),
            dict(mdrun_args=['-nb', 'gpu', '-update', 'gpu'],
                 write_interval=write_interval * 2)])
    by_seed = {c.config['seed_index']: c.config
               for p in packs_of(farmer) for c in p.clones}
    suite.check('each seed keeps its own -update setting',
                by_seed[0]['mdrun_args'][-1] == 'cpu'
                and by_seed[1]['mdrun_args'][-1] == 'gpu')
    suite.check('a seed can override write_interval',
                by_seed[0]['write_interval'] == write_interval
                and by_seed[1]['write_interval'] == write_interval * 2)
    suite.check('overrides do not leak between seeds',
                by_seed[0]['mdrun_args'] is not by_seed[1]['mdrun_args'])
    try:
        make_farmer(work, pack_size=pack_size,
                    seed_config_overrides=[dict(steps=5), dict()])
        suite.check('an override of a clone-derived key is refused', False,
                    '-> no exception')
    except ValueError as exc:
        suite.check('an override of a clone-derived key is refused',
                    'steps' in str(exc), f'-> {str(exc)[:50]}')
    try:
        make_farmer(work, pack_size=pack_size, seed_config_overrides=[dict()])
        suite.check('a wrong-length override list is refused', False,
                    '-> no exception')
    except ValueError as exc:
        suite.check('a wrong-length override list is refused', True,
                    f'-> {str(exc)[:50]}')

    suite.section('unequal steps per generation')
    members = packs_of(farmer)[0].clones
    members[1].total_steps = members[0].total_steps * 2
    pack_kwargs = dict(
        run_script=gp.default_gmx_pack_run_script, cpus_per_task=cpus,
        dry_run=True)
    scheduler_kws = dict(gpu_line='', queue_name='gpu', exclude_nodes='',
                         cpus=cpus, run_script_name='run.py')
    # Noted, not refused: the cost is a straggler and a harvest that needs its
    # own Harvester, neither of which is a reason to stop the campaign.
    noted = io.StringIO()
    with redirect_stdout(noted):
        uneven = ClonePack(members, work / 'pack-uneven', 'sbatch',
                           util.basic_scheduler_fstrings_mps['slurm'],
                           scheduler_kws, **pack_kwargs)
    said = noted.getvalue()
    suite.check('unequal steps still build a pack', uneven is not None)
    suite.check('unequal steps are called out at boot',
                'NOTE' in said and 'steps per generation' in said,
                f'-> {said.strip()[:60]}')
    suite.check('and the note says which Harvester problem it causes',
                'Harvester' in said)
    pack = ClonePack(members, work / 'pack-wallclock', 'sbatch',
                     util.basic_scheduler_fstrings_mps['slurm'], scheduler_kws,
                     wallclock_matched=True, **pack_kwargs)
    suite.check('wallclock_matched=True permits them',
                len(pack.clones) == pack_size)
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
