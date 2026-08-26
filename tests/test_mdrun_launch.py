"""An OMP_NUM_THREADS inherited from the submitting shell must not kill mdrun.

mdrun refuses to start when OMP_NUM_THREADS and -ntomp disagree, and sbatch
exports the submitting environment by default. gmx_pack gives every replica its
own -ntomp, so one export in a site profile would fatal every replica of every
packed job at startup. The launcher therefore sets OMP_NUM_THREADS to whatever
-ntomp asks for, and leaves it alone when no -ntomp is passed, since there it
is what sets mdrun's thread count.
"""
import os
import sys

import harness
from harness import Suite

import mdfarmer
from mdfarmer import gmx_simulate as gs

STEPS_PER_GEN = 200
WRITE_INTERVAL = 100
HOSTILE_OMP = '2'               # what the submitting shell exported
REPLICA_CORES = '1'             # what this replica's -ntomp asks for
SEED_INDEX = 0
CLONE_INDEX = 0
DIRNAME_PAD = 2
SEP = '-'


def main(gmx_bin=harness.GMX_BIN, hostile_omp=HOSTILE_OMP,
         replica_cores=REPLICA_CORES):
    suite = Suite('mdrun_launch')

    suite.section('what the launcher puts in the environment')
    quiet = {'PATH': '/usr/bin'}
    env = gs._mdrun_env(['gmx', 'mdrun', '-ntomp', replica_cores], environ=quiet)
    suite.check('-ntomp sets OMP_NUM_THREADS to match',
                env['OMP_NUM_THREADS'] == replica_cores,
                f'-> {env.get("OMP_NUM_THREADS")}')
    loud = dict(quiet, OMP_NUM_THREADS=hostile_omp)
    env = gs._mdrun_env(['gmx', 'mdrun', '-ntomp', replica_cores], environ=loud)
    suite.check('and overrides one the shell had already exported',
                env['OMP_NUM_THREADS'] == replica_cores,
                f'-> {env.get("OMP_NUM_THREADS")}')
    env = gs._mdrun_env(['gmx', 'mdrun', '-nb', 'cpu'], environ=loud)
    suite.check('without -ntomp the inherited value is left alone',
                env['OMP_NUM_THREADS'] == hostile_omp,
                f'-> {env.get("OMP_NUM_THREADS")}')
    suite.check('the rest of the environment comes through',
                env['PATH'] == quiet['PATH'])

    suite.section('a real generation survives a hostile OMP_NUM_THREADS')
    work = harness.workdir('mdrun_launch')
    structure, topology, mdp = harness.build_water_system(work, gmx_bin=gmx_bin)
    farm = work / 'farm'
    os.environ['OMP_NUM_THREADS'] = hostile_omp
    try:
        traj = gs.gmx_generation(
            traj_dir_top_level=str(farm), top_fn=str(topology),
            seed_index=SEED_INDEX, clone_index=CLONE_INDEX, gen_index=0,
            title='omp', seed_fn=str(structure), structure_fn=str(structure),
            mdp_fn=str(mdp), dirname_pad=DIRNAME_PAD, sep=SEP,
            traj_name='prod', traj_suffix='.xtc', restart_name='state.cpt',
            steps=STEPS_PER_GEN, steps_per_gen=STEPS_PER_GEN,
            write_interval=WRITE_INTERVAL, temperature=300, gen_seed_base=5,
            gmx_bin=gmx_bin, grompp_maxwarn=3, new_velocities=True, append=False,
            mdrun_args=('-nb', 'cpu', '-pme', 'cpu', '-ntomp', replica_cores))
        started = True
    except RuntimeError as exc:
        traj, started = None, False
        print(f'   mdrun refused to start: {exc}', flush=True)
    finally:
        os.environ.pop('OMP_NUM_THREADS', None)
    suite.check('mdrun starts and the generation finishes',
                started and traj is not None)
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
