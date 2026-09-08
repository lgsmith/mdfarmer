"""Three generation-counting bugs: an off-by-one that submits one job too
many, a resumed generation 0 that redraws velocities, and a ClonePack that
can never retire a permanently failed member. Plus one recovery-hardening
fix: a corrupt DCD must cascade out of _try_recover_gen, not raise through
it and cost the whole clone.

Everything here is a dry run against hand-built Clones and a minimal Farmer,
so it needs no GROMACS.
"""
import json
import sys

import harness
from harness import Suite

from mdfarmer import utilities as util
from mdfarmer import farmer as fm
from mdfarmer.seeder import Clone, ClonePack, _try_recover_gen

N_GENS = 3
PACK_N_GENS = 2
STEPS_PER_GEN = 1000
WRITE_INTERVAL = 100
PACK_RESTARTS = 2


def base_config(work, seed_index=0, clone_index=0,
                steps_per_gen=STEPS_PER_GEN, write_interval=WRITE_INTERVAL):
    return dict(
        traj_dir_top_level=str(work / 'farm'), seed_index=seed_index,
        clone_index=clone_index, gen_index=0, title='budget',
        structure_fn=str(work / 'seed.gro'), dirname_pad=2, sep='_',
        traj_name='prod', traj_suffix='.xtc', restart_name='state.cpt',
        steps=steps_per_gen, steps_per_gen=steps_per_gen,
        write_interval=write_interval, new_velocities=True, append=False,
        mdrun_args=[])


class ReapCounter:
    """Stands in for a Harvester, so no scheduler job is ever submitted."""

    def __init__(self):
        self.reaped = []

    def reap(self, gen_dir, dry_run=False):
        self.reaped.append(gen_dir)


class AlwaysDoneClone(Clone):
    """Every check-in reports this generation as already finished."""

    def gen_remaining_steps(self):
        return 0

    def was_preempted(self):
        return False


class NeverRunsClone(Clone):
    """Every check-in reports zero progress, burning the restart budget."""

    def gen_remaining_steps(self):
        return self.total_steps

    def was_preempted(self):
        return False


def preseed(clone):
    """Places a placeholder checkpoint in the clone's current gen dir, standing
    in for what a real launch would have written there, so an AlwaysDoneClone
    (which never actually launches its first generation) can still advance."""
    (clone.current_gen_dir / clone.config['restart_name']).write_text('ckpt\n')


def make_clone(work, cls, clone_index, last_gen_index, harvester=None,
              restarts_per_gen=3, scheduler_kws=None):
    scheduler_kws = scheduler_kws or dict(
        gpu_line='', queue_name='gpu', exclude_nodes='', run_script_name='run.py')
    return cls(base_config(work, clone_index=clone_index), 'sbatch',
              util.basic_scheduler_fstrings['slurm'], scheduler_kws,
              seed_fn=str(work / 'seed.gro'), sep='_', dirname_pad=2,
              steps_per_gen=STEPS_PER_GEN, dry_run=True, harvester=harvester,
              restarts_per_gen=restarts_per_gen, last_gen_index=last_gen_index)


def gen_dirs(work, clone_index=0):
    clone_dir = util.dir_seeds_clones(work / 'farm', 0, clone_index, 2,
                                      sep='_', mkdir=False)
    return sorted(clone_dir.iterdir()) if clone_dir.is_dir() else []


def probe_recover(gen_path, *, state_step, header_info, truncate,
                  restart_name='state.cpt', traj_name='prod',
                  traj_suffix='.dcd', write_interval=100, total_steps=1000):
    """Calls _try_recover_gen with util's DCD/state.xml readers replaced by
    fakes, so the branch under test doesn't need a real OpenMM state.xml or
    DCD file. Returns (result, exception raised while recovering)."""
    saved = (util.is_state_xml_usable, util.state_xml_step_count,
            util.dcd_header_info, util.truncate_dcd_to_nframes)
    util.is_state_xml_usable = lambda p: True
    util.state_xml_step_count = lambda p: state_step
    util.dcd_header_info = header_info
    util.truncate_dcd_to_nframes = truncate
    try:
        result = _try_recover_gen(
            gen_path, append_mode=True, restart_name=restart_name,
            traj_name=traj_name, traj_suffix=traj_suffix,
            write_interval=write_interval, total_steps=total_steps,
            top_fn='unused.top')
        return result, None
    except Exception as exc:
        return None, exc
    finally:
        (util.is_state_xml_usable, util.state_xml_step_count,
         util.dcd_header_info, util.truncate_dcd_to_nframes) = saved


def main(n_gens=N_GENS, pack_n_gens=PACK_N_GENS, steps_per_gen=STEPS_PER_GEN):
    suite = Suite('generation_budget')
    work = harness.workdir('generation_budget')
    (work / 'seed.gro').write_text('placeholder\n')

    suite.section('a clone submits exactly n_gens generations, not one more')
    reaper = ReapCounter()
    clone = make_clone(work, AlwaysDoneClone, 0, n_gens - 1, harvester=reaper)
    preseed(clone)

    farm = object.__new__(fm.Farmer)
    farm.n_gens = n_gens
    farm.finished_clones = set()
    farm.active_set = {clone}
    farm.failed_clone_set = set()
    farm.priority_ordered_clones = [[clone]]
    farm.current_jids = set()
    farm.overwrite = True
    farm.launch_interval = 0

    for _ in range(n_gens + 3):
        if clone in farm.finished_clones:
            break
        farm.launch(update_jids=False)

    suite.check('the clone reaches finished_clones',
                clone in farm.finished_clones)
    suite.check('it is dropped from active_set',
                clone not in farm.active_set)
    suite.check(f'exactly {n_gens} generations were harvested, not one more',
                len(reaper.reaped) == n_gens, f'-> {len(reaper.reaped)}')
    suite.check(f'exactly {n_gens} generation directories were created',
                len(gen_dirs(work)) == n_gens, f'-> {len(gen_dirs(work))}')
    suite.check('no gen dir was created for the phantom extra generation',
                not any(p.name.endswith(f'{n_gens:02d}') for p in gen_dirs(work)))

    suite.section('resuming a partial generation 0 does not redraw velocities')
    resume_config = base_config(work, clone_index=1)
    resume_config['new_velocities'] = True
    resume_config['append'] = False

    class PartialGenZero(Clone):
        """Reports generation 0 as half finished on its only check-in."""

        def gen_remaining_steps(self):
            return self.total_steps // 2

        def was_preempted(self):
            return False

    resuming = PartialGenZero(
        resume_config, 'sbatch', util.basic_scheduler_fstrings['slurm'],
        dict(gpu_line='', queue_name='gpu', exclude_nodes='',
            run_script_name='run.py'),
        seed_fn=str(work / 'seed.gro'), sep='_', dirname_pad=2,
        steps_per_gen=steps_per_gen, dry_run=True)
    resuming.check_start_gen(set())
    suite.check('the resume is recorded as an append',
                resuming.config['append'] is True)
    suite.check('velocities are not redrawn on the resume',
                resuming.config['new_velocities'] is False,
                f"-> {resuming.config['new_velocities']}")

    suite.section('a ClonePack retires a member that exhausts its restarts')
    scheduler_kws = dict(gpu_line='', queue_name='gpu', exclude_nodes='',
                        cpus=4, run_script_name='run.py')
    healthy = make_clone(work, AlwaysDoneClone, 2, pack_n_gens - 1,
                        scheduler_kws=dict(scheduler_kws))
    preseed(healthy)
    doomed = make_clone(work, NeverRunsClone, 3, pack_n_gens - 1,
                       restarts_per_gen=PACK_RESTARTS,
                       scheduler_kws=dict(scheduler_kws))

    pack = ClonePack([healthy, doomed], work / 'pack', 'sbatch',
                     util.basic_scheduler_fstrings_mps['slurm'], scheduler_kws,
                     run_script='', cpus_per_task=4, dry_run=True)

    outcomes = [pack.check_start_gen(set(), overwrite=True) for _ in range(4)]
    suite.check('the pack keeps launching while the doomed member retries',
                outcomes[0] is True, f'-> {outcomes}')
    suite.check('the doomed member is retired once its restart budget is spent',
                1 in pack.retired, f'-> {pack.retired}')
    suite.check('the healthy member is never retired',
                0 not in pack.retired)
    suite.check('the pack keeps reporting progress, not failure',
                all(outcomes), f'-> {outcomes}')
    suite.check("current_gen tracks the live member, not the retired one's "
               'frozen value',
                pack.current_gen == pack_n_gens, f'-> {pack.current_gen}')
    suite.check(f'the healthy member ran exactly {pack_n_gens} generations',
                len(gen_dirs(work, clone_index=2)) == pack_n_gens,
                f'-> {len(gen_dirs(work, clone_index=2))}')
    suite.check('the doomed member never advanced past generation 0',
                doomed.config['gen_index'] == 0)

    suite.section('a ClonePack fails once every member is retired')
    doomed_a = make_clone(work, NeverRunsClone, 4, pack_n_gens - 1,
                         restarts_per_gen=PACK_RESTARTS,
                         scheduler_kws=dict(scheduler_kws))
    doomed_b = make_clone(work, NeverRunsClone, 5, pack_n_gens - 1,
                         restarts_per_gen=PACK_RESTARTS,
                         scheduler_kws=dict(scheduler_kws))
    all_fail_pack = ClonePack(
        [doomed_a, doomed_b], work / 'pack-allfail', 'sbatch',
        util.basic_scheduler_fstrings_mps['slurm'], scheduler_kws,
        run_script='', cpus_per_task=4, dry_run=True)
    fail_outcomes = [all_fail_pack.check_start_gen(set(), overwrite=True)
                    for _ in range(PACK_RESTARTS + 1)]
    suite.check('the pack fails once both members exhaust their restarts',
                fail_outcomes[-1] is False, f'-> {fail_outcomes}')
    suite.check('both members are recorded as retired',
                all_fail_pack.retired == {0, 1}, f'-> {all_fail_pack.retired}')

    suite.section('_try_recover_gen cascades instead of raising on a corrupt DCD')
    position_dir = work / 'corrupt_position'
    position_dir.mkdir()
    (position_dir / 'config.json').write_text(json.dumps(base_config(work)))
    (position_dir / 'state.cpt').write_text('ckpt\n')
    (position_dir / 'prod.dcd').write_bytes(b'not a real dcd, just nonempty')

    def raises(path, n):
        raise ValueError('simulated corrupt DCD header')

    result, exc = probe_recover(
        position_dir, state_step=500,
        header_info=lambda p: {'nset': 7, 'nsavc': 100}, truncate=raises)
    suite.check('a corrupt position DCD cascades rather than raising',
                exc is None and result is None, f'-> exc={exc!r} result={result!r}')

    tandem_dir = work / 'corrupt_tandem'
    tandem_dir.mkdir()
    (tandem_dir / 'config.json').write_text(json.dumps(dict(
        base_config(work), velocity_traj_suffix='.dcd', velocity_name='vel')))
    (tandem_dir / 'state.cpt').write_text('ckpt\n')
    (tandem_dir / 'prod.dcd').write_bytes(b'placeholder positions')
    (tandem_dir / 'vel.dcd').write_bytes(b'placeholder velocities')

    def truncate_mixed(path, n):
        if path.name == 'vel.dcd':
            raise ValueError('simulated corrupt tandem DCD header')
        return n  # position truncation succeeds, so the tandem loop is reached

    result, exc = probe_recover(
        tandem_dir, state_step=500,
        header_info=lambda p: {'nset': 7, 'nsavc': 100}, truncate=truncate_mixed)
    suite.check('a corrupt tandem DCD cascades rather than raising',
                exc is None and result is None, f'-> exc={exc!r} result={result!r}')

    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
