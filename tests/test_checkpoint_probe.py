"""A checkpoint gmx cannot read must not look like a file that is not a
checkpoint at all.

The two cases need opposite handling. At generation 0 the seed really is a
.gro, and the runner has to go on without -cpi. A generation's own state.cpt
that is truncated or corrupt is a different thing: treating it as "not a
checkpoint" falls back to the seed, which silently rewinds the whole
generation and then reports success. is_checkpoint therefore asks only whether
the file carries the GROMACS checkpoint magic number, and leaves reading it to
checkpoint_part_step, which raises.
"""
import shutil
import sys
from pathlib import Path

import harness
from harness import Suite

import mdfarmer
from mdfarmer import gmx_simulate as gs

STEPS_PER_GEN = 400             # one generation's full step budget
WRITE_INTERVAL = 100
PARTIAL_STEPS = 200             # steps a deliberately short launch runs
TRUNCATE_BYTES = 64             # enough to keep the magic, far too few to read
SEED_INDEX = 0
CLONE_INDEX = 0
DIRNAME_PAD = 2
SEP = '-'
CPU_ARGS = ('-nb', 'cpu', '-pme', 'cpu', '-ntomp', '2')


def run_gen(config, **overrides):
    """gmx_generation with these overrides; None if it left the generation
    incomplete, else the trajectory it returned."""
    try:
        return gs.gmx_generation(**dict(config, **overrides))
    except gs.GenIncomplete:
        return None


def truncate(path, n_bytes=TRUNCATE_BYTES):
    """Keep the first n_bytes of a file, so its magic number survives."""
    head = Path(path).read_bytes()[:n_bytes]
    Path(path).write_bytes(head)
    return path


def raises_reading(call):
    """(did it raise, what it said) for a call expected to refuse to guess."""
    try:
        call()
    except (RuntimeError, ValueError) as exc:
        return True, f'{type(exc).__name__}: {str(exc).splitlines()[0][:70]}'
    return False, 'no exception'


def main(gmx_bin=harness.GMX_BIN):
    suite = Suite('checkpoint_probe')
    work = harness.workdir('checkpoint_probe')
    structure, topology, mdp = harness.build_water_system(work, gmx_bin=gmx_bin)
    farm = work / 'farm'

    common = dict(
        traj_dir_top_level=str(farm), top_fn=str(topology),
        seed_index=SEED_INDEX, clone_index=CLONE_INDEX, title='cpt-probe',
        structure_fn=str(structure), mdp_fn=str(mdp), dirname_pad=DIRNAME_PAD,
        sep=SEP, traj_name='prod', traj_suffix='.xtc', restart_name='state.cpt',
        steps=STEPS_PER_GEN, steps_per_gen=STEPS_PER_GEN,
        write_interval=WRITE_INTERVAL, temperature=300, gen_seed_base=7,
        gmx_bin=gmx_bin, grompp_maxwarn=3, mdrun_args=CPU_ARGS, append=False)

    def gen_dir(gen_index):
        return mdfarmer.dir_seeds_clones_gens(
            farm, SEED_INDEX, CLONE_INDEX, gen_index, DIRNAME_PAD, sep=SEP)

    suite.section('what the magic number can and cannot tell apart')
    traj = run_gen(common, gen_index=0, seed_fn=str(structure),
                   new_velocities=True,
                   mdrun_args=CPU_ARGS + ('-nsteps', str(PARTIAL_STEPS)))
    own_cpt = gen_dir(0) / common['restart_name']
    suite.check('a short launch leaves an incomplete generation and a checkpoint',
                traj is None and own_cpt.is_file())
    suite.check('a real checkpoint is recognised', gs.is_checkpoint(own_cpt))
    suite.check('the generation-0 seed .gro is not a checkpoint',
                not gs.is_checkpoint(structure))
    suite.check('a file that is not there is not a checkpoint',
                not gs.is_checkpoint(gen_dir(0) / 'absent.cpt'))
    empty = work / 'empty.cpt'
    empty.write_bytes(b'')
    suite.check('an empty file is not a checkpoint', not gs.is_checkpoint(empty))

    stump = truncate(shutil.copy(own_cpt, work / 'stump.cpt'))
    suite.check('a truncated checkpoint is still a checkpoint',
                gs.is_checkpoint(stump))
    raised, detail = raises_reading(
        lambda: gs.checkpoint_part_step(stump, gmx_bin=gmx_bin))
    suite.check('and reading it raises rather than reporting a step',
                raised, f'-> {detail}')

    suite.section('generation 0 refuses to restart over a damaged checkpoint')
    truncate(own_cpt)
    raised, detail = raises_reading(
        lambda: run_gen(common, gen_index=0, seed_fn=str(structure),
                        new_velocities=True))
    suite.check('an unreadable state.cpt stops the run instead of starting over',
                raised, f'-> {detail}')

    suite.section('a later generation refuses to rewind to its seed')
    clean = work / 'farm2'
    common = dict(common, traj_dir_top_level=str(clean))

    def gen_dir(gen_index):
        return mdfarmer.dir_seeds_clones_gens(
            clean, SEED_INDEX, CLONE_INDEX, gen_index, DIRNAME_PAD, sep=SEP)

    traj = run_gen(common, gen_index=0, seed_fn=str(structure),
                   new_velocities=True)
    suite.check('generation 0 finishes', traj is not None)

    gen_one = gen_dir(1)
    gen_one.mkdir(parents=True, exist_ok=True)
    seed = shutil.copy(gen_dir(0) / common['restart_name'],
                       gen_one / common['restart_name'])
    traj = run_gen(common, gen_index=1, seed_fn=str(seed), new_velocities=False,
                   mdrun_args=CPU_ARGS + ('-nsteps', str(PARTIAL_STEPS)))
    own_cpt = gen_one / common['restart_name']
    suite.check('generation 1 runs partway and keeps its seed separately',
                traj is None and own_cpt.is_file()
                and (gen_one / gs.SEED_CPT_NAME).is_file())

    truncate(own_cpt)
    raised, detail = raises_reading(
        lambda: run_gen(common, gen_index=1, seed_fn=str(seed),
                        new_velocities=False))
    suite.check('a damaged state.cpt is refused, not swapped for seed.cpt',
                raised, f'-> {detail}')
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
