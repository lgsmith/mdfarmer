"""A generation halted between two frames resumes back onto the frame grid.

mdrun -cpi restarts from the checkpoint's step, which for a walltime or preempt
stop is wherever the run happened to be -- almost never a multiple of
nstxout-compressed. If the resumed run wrote a frame there, the generation would
carry one frame more than its budget buys and one that does not sit on the
stride, which is what harvester.resolve_seam counts on. That error is raised
during a harvest, where seeder swallows it after marking the generation reaped,
so the harvest would be lost silently.
"""
import sys

import numpy as np
import mdtraj as md

import harness
from harness import Suite

import mdfarmer
from mdfarmer import gmx_simulate as gs

STEPS_PER_GEN = 1000
WRITE_INTERVAL = 100
# Deliberately not a multiple of WRITE_INTERVAL: the point of the test.
HALT_STEP = 250
SEED_INDEX = 0
CLONE_INDEX = 0
DIRNAME_PAD = 2
SEP = '-'
CPU_ARGS = ('-nb', 'cpu', '-pme', 'cpu', '-ntomp', '2')


def main(gmx_bin=harness.GMX_BIN, halt_step=HALT_STEP,
         steps_per_gen=STEPS_PER_GEN, write_interval=WRITE_INTERVAL):
    suite = Suite('offgrid_restart')
    work = harness.workdir('offgrid_restart')
    structure, topology, mdp = harness.build_water_system(work, gmx_bin=gmx_bin)
    farm = work / 'farm'

    common = dict(
        traj_dir_top_level=str(farm), top_fn=str(topology),
        seed_index=SEED_INDEX, clone_index=CLONE_INDEX, gen_index=0,
        title='offgrid', seed_fn=str(structure), structure_fn=str(structure),
        mdp_fn=str(mdp), dirname_pad=DIRNAME_PAD, sep=SEP, traj_name='prod',
        traj_suffix='.xtc', restart_name='state.cpt', steps=steps_per_gen,
        steps_per_gen=steps_per_gen, write_interval=write_interval,
        target_step=harness.target_step(0, steps_per_gen),
        temperature=300, gen_seed_base=4, gmx_bin=gmx_bin, grompp_maxwarn=3,
        new_velocities=True, append=False, mdrun_args=CPU_ARGS)

    suite.section(f'halt between two frames, at step {halt_step}')
    try:
        gs.gmx_generation(**dict(
            common, mdrun_args=CPU_ARGS + ('-nsteps', str(halt_step))))
    except gs.GenIncomplete:
        pass
    gen_dir = mdfarmer.dir_seeds_clones_gens(
        farm, SEED_INDEX, CLONE_INDEX, 0, DIRNAME_PAD, sep=SEP)
    stopped_at = gs.checkpoint_step(gen_dir / common['restart_name'],
                                    gmx_bin=gmx_bin)
    suite.check('the checkpoint really is off the frame grid',
                stopped_at == halt_step and stopped_at % write_interval,
                f'-> step {stopped_at}, interval {write_interval}')

    suite.section('the resumed generation lands back on the grid')
    traj = gs.gmx_generation(**common)
    with md.open(str(traj)) as handle:
        steps = np.asarray(handle.read()[2]).astype(int)
    expected = steps_per_gen // write_interval + 1
    off_grid = steps[steps % write_interval != 0]
    suite.check('no frame is written at the restart step',
                len(off_grid) == 0, f'-> {off_grid.tolist()}')
    suite.check('the generation holds exactly the frames its budget buys',
                len(steps) == expected, f'-> {len(steps)} vs {expected}')
    suite.check('the steps are the exact stride sequence',
                np.array_equal(steps,
                               np.arange(0, steps_per_gen + 1, write_interval)),
                f'-> {steps.tolist()}')
    suite.check('and rise strictly, so nothing is written twice',
                bool(np.all(np.diff(steps) > 0)))
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
