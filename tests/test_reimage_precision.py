"""Reimaging must not quantise a trajectory written at a finer precision.

LOOS's XTCWriter defaults to precision 1000, the same grid GROMACS defaults to.
A run that asked for a finer compressed-x-precision therefore came back through
the LOOS reimaging path on the coarser grid, losing the digits it paid for on
every frame of the campaign, with nothing said.
"""
import sys

import harness
from harness import Suite

from mdfarmer import gmx_simulate as gs, reimage

STEPS = 200
WRITE_INTERVAL = 100
FINE_PRECISION = 10000          # 10x finer than the GROMACS default
CPU_ARGS = ('-nb', 'cpu', '-pme', 'cpu', '-ntomp', '2')


def main(gmx_bin=harness.GMX_BIN, fine_precision=FINE_PRECISION):
    suite = Suite('reimage_precision')
    work = harness.workdir('reimage_precision')
    include_dir = str(harness.require_gmx(gmx_bin))
    structure, topology, mdp = harness.build_water_system(work, gmx_bin=gmx_bin)
    mdp.write_text(f'{mdp.read_text()}'
                   f'compressed-x-precision = {fine_precision}\n')

    suite.section(f'a run written at compressed-x-precision {fine_precision}')
    traj = gs.gmx_generation(
        traj_dir_top_level=str(work / 'farm'), top_fn=str(topology),
        seed_index=0, clone_index=0, gen_index=0, title='fine',
        seed_fn=str(structure), structure_fn=str(structure), mdp_fn=str(mdp),
        dirname_pad=2, sep='-', traj_name='prod', traj_suffix='.xtc',
        restart_name='state.cpt', steps=STEPS, steps_per_gen=STEPS,
        target_step=harness.target_step(0, STEPS),
        write_interval=WRITE_INTERVAL, temperature=300, gen_seed_base=1,
        gmx_bin=gmx_bin, grompp_maxwarn=3, new_velocities=True, append=False,
        mdrun_args=CPU_ARGS)
    suite.check('the source carries the finer precision',
                gs.xtc_precision(traj) == fine_precision,
                f'-> {gs.xtc_precision(traj)}')

    suite.section('reimaging keeps it')
    out, n_written = reimage.reimage_with_loos(
        traj, str(structure), str(work / 'whole.xtc'), top_fn=str(topology),
        include_dir=include_dir, verify=False)
    suite.check('every frame is written', n_written == STEPS // WRITE_INTERVAL + 1,
                f'-> {n_written}')
    suite.check('the reimaged trajectory holds the same precision',
                gs.xtc_precision(out) == fine_precision,
                f'-> {gs.xtc_precision(out)}')
    suite.check('which is not the writer default',
                fine_precision != reimage.XTC_PRECISION,
                f'-> default {reimage.XTC_PRECISION:g}')
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
