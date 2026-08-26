"""A preempted generation still has to say how far it got.

The orchestrator decides whether to spend one of a generation's restarts by
asking whether it made progress since last time. gmx_gen_progress reads
gen_status.json, and with no status file it reports that nothing has run. So a
generation that is preempted an hour into a long launch, over and over, gets
charged a restart every time and is eventually failed outright with a
perfectly good checkpoint sitting next to it.
"""
import sys

import harness
from harness import Suite

import mdfarmer
from mdfarmer import gmx_simulate as gs

STEPS_PER_GEN = 400
WRITE_INTERVAL = 100
PARTIAL_STEPS = 200             # steps the first launch gets through
SEED_INDEX = 0
CLONE_INDEX = 0
DIRNAME_PAD = 2
SEP = '-'
CPU_ARGS = ('-nb', 'cpu', '-pme', 'cpu', '-ntomp', '2')


def preempt_immediately(cmd, cwd, handle_preempt, **kwargs):
    """Stand in for _run_mdrun, as if the sentinel appeared the moment it
    started -- the checkpoint already in the directory is all this launch has."""
    raise gs.Preempted('preempt sentinel (test)')


def main(gmx_bin=harness.GMX_BIN):
    suite = Suite('preempt_status')
    work = harness.workdir('preempt_status')
    structure, topology, mdp = harness.build_water_system(work, gmx_bin=gmx_bin)
    farm = work / 'farm'

    common = dict(
        traj_dir_top_level=str(farm), top_fn=str(topology),
        seed_index=SEED_INDEX, clone_index=CLONE_INDEX, gen_index=0,
        title='preempt', seed_fn=str(structure), structure_fn=str(structure),
        mdp_fn=str(mdp), dirname_pad=DIRNAME_PAD, sep=SEP, traj_name='prod',
        traj_suffix='.xtc', restart_name='state.cpt', steps=STEPS_PER_GEN,
        steps_per_gen=STEPS_PER_GEN, write_interval=WRITE_INTERVAL,
        temperature=300, gen_seed_base=11, gmx_bin=gmx_bin, grompp_maxwarn=3,
        new_velocities=True, append=False)
    gen_dir = mdfarmer.dir_seeds_clones_gens(
        farm, SEED_INDEX, CLONE_INDEX, 0, DIRNAME_PAD, sep=SEP)

    suite.section('a launch gets partway and is stopped')
    try:
        gs.gmx_generation(**dict(
            common, mdrun_args=CPU_ARGS + ('-nsteps', str(PARTIAL_STEPS))))
    except gs.GenIncomplete:
        pass
    reached = gs.checkpoint_step(gen_dir / common['restart_name'],
                                 gmx_bin=gmx_bin)
    suite.check('its checkpoint holds the steps it ran',
                reached == PARTIAL_STEPS, f'-> {reached}')
    (gen_dir / gs.GEN_STATUS_NAME).unlink()
    suite.check('and the status file is out of the way for the next launch',
                gs.read_gen_status(gen_dir) is None)

    suite.section('the next launch is preempted before it runs anything')
    real_run_mdrun = gs._run_mdrun
    gs._run_mdrun = preempt_immediately
    try:
        gs.gmx_generation(**dict(common, mdrun_args=CPU_ARGS))
        preempted = False
    except gs.Preempted:
        preempted = True
    finally:
        gs._run_mdrun = real_run_mdrun
    suite.check('the preempt reaches the orchestrator', preempted)

    status = gs.read_gen_status(gen_dir)
    suite.check('and it left a status file behind', status is not None,
                f'-> {status}')
    suite.check('reporting the progress the earlier launch made',
                status is not None and status['reached_step'] == PARTIAL_STEPS,
                f'-> {status}')
    suite.check('and not claiming the generation is finished',
                status is not None and not status['complete'])

    progress = gs.gmx_gen_progress(gen_dir, total_steps=STEPS_PER_GEN,
                                   gen_index=0)
    suite.check('so the orchestrator sees steps still to run, not a fresh start',
                progress == STEPS_PER_GEN - PARTIAL_STEPS,
                f'-> {progress} vs {STEPS_PER_GEN - PARTIAL_STEPS}')
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
