"""nsteps is the absolute cumulative step target, so the tpr has to count from
zero.

grompp keeps whatever init-step a base .mdp sets, and mdrun then stops at
init-step + nsteps. The completion test compares the checkpoint's step counter
against (gen_index + 1) * steps_per_gen, so an inherited init-step makes
generation 0 look finished the moment it starts, and the generation after it
finds a seed already past its own target and has nothing to merge.
"""
import sys
from pathlib import Path

import mdtraj as md

import harness
from harness import Suite

import mdfarmer
from mdfarmer import gmx_simulate as gs

STEPS_PER_GEN = 400
WRITE_INTERVAL = 100
INHERITED_INIT_STEP = 500       # what the base .mdp carries in
SEED_INDEX = 0
CLONE_INDEX = 0
DIRNAME_PAD = 2
SEP = '-'
CPU_ARGS = ('-nb', 'cpu', '-pme', 'cpu', '-ntomp', '2')


def main(gmx_bin=harness.GMX_BIN, init_step=INHERITED_INIT_STEP):
    suite = Suite('init_step')
    work = harness.workdir('init_step')
    structure, topology, mdp = harness.build_water_system(work, gmx_bin=gmx_bin)
    mdp.write_text(mdp.read_text() + f'init-step = {init_step}\n')
    farm = work / 'farm'

    common = dict(
        traj_dir_top_level=str(farm), top_fn=str(topology),
        seed_index=SEED_INDEX, clone_index=CLONE_INDEX, title='init-step',
        structure_fn=str(structure), mdp_fn=str(mdp), dirname_pad=DIRNAME_PAD,
        sep=SEP, traj_name='prod', traj_suffix='.xtc', restart_name='state.cpt',
        steps=STEPS_PER_GEN, steps_per_gen=STEPS_PER_GEN,
        write_interval=WRITE_INTERVAL, temperature=300, gen_seed_base=3,
        gmx_bin=gmx_bin, grompp_maxwarn=3, mdrun_args=CPU_ARGS, append=False)

    def gen_dir(gen_index):
        return mdfarmer.dir_seeds_clones_gens(
            farm, SEED_INDEX, CLONE_INDEX, gen_index, DIRNAME_PAD, sep=SEP)

    suite.section('generation 0 runs its whole budget, not none of it')
    traj = gs.gmx_generation(gen_index=0, seed_fn=str(structure),
                             new_velocities=True, **common)
    reached = gs.checkpoint_step(gen_dir(0) / common['restart_name'],
                                 gmx_bin=gmx_bin)
    suite.check('the checkpoint stops at the generation budget',
                reached == STEPS_PER_GEN, f'-> {reached} vs {STEPS_PER_GEN}')
    status = gs.read_gen_status(gen_dir(0))
    suite.check('the status it reports agrees',
                status['reached_step'] == STEPS_PER_GEN and status['complete'],
                f'-> {status}')
    frames = md.load(str(traj), top=str(structure))
    expected = STEPS_PER_GEN // WRITE_INTERVAL + 1
    suite.check('the trajectory holds the frames that budget buys',
                len(frames) == expected, f'-> {len(frames)} vs {expected}')

    suite.section('and generation 1 has somewhere to start from')
    gen_one = gen_dir(1)
    gen_one.mkdir(parents=True, exist_ok=True)
    seed = gen_one / common['restart_name']
    seed.write_bytes((gen_dir(0) / common['restart_name']).read_bytes())
    traj = gs.gmx_generation(gen_index=1, seed_fn=str(seed),
                             new_velocities=False, **common)
    reached = gs.checkpoint_step(gen_one / common['restart_name'],
                                 gmx_bin=gmx_bin)
    suite.check('generation 1 reaches its own cumulative target',
                reached == 2 * STEPS_PER_GEN,
                f'-> {reached} vs {2 * STEPS_PER_GEN}')
    suite.check('and it produced a trajectory to merge',
                Path(traj).is_file())
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
