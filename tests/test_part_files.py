"""A relaunch from a checkpoint earlier than the parts already on disk (a
rewind) must not leave the abandoned branch's higher-numbered parts where
concat_parts' glob can still find them.

The scenario: a generation writes part0001, then a branch runs on from there
to part0002 and part0003. Something then rewinds the checkpoint back to the
point after part0001 (the same failure mode as a cascade that restores an
older checkpoint over state.cpt) and the generation is relaunched. The
relaunch's own part0002 overwrites the old one, but nothing about mdrun
touches the abandoned part0003 -- it is still sitting there, still matching
the glob, and the merge prefers its frames wherever the two branches' steps
overlap.
"""
import shutil
import sys
from pathlib import Path

import numpy as np
import mdtraj as md

import harness
from harness import Suite

import mdfarmer
from mdfarmer import gmx_simulate as gs

# Generation 0's full step budget. Bigger than the abandoned branch reaches
# (400 + 1000 + 1000), so that branch never finishes -- a finished generation
# merges its parts and deletes them, and then there is no stale part left to
# outrank anything. Lingering parts are exactly the unfinished case.
STEPS_PER_GEN = 3400
WRITE_INTERVAL = 100            # -> 0.2 ps frame spacing at dt=0.002
DT_PS = 0.002                   # must match harness.water_mdp's dt
REWIND_STEP = 400               # step of the checkpoint the rewind restores
BRANCH_STEP = 1000              # steps each abandoned-branch launch adds
SEED_INDEX = 0
CLONE_INDEX = 0
GEN_INDEX = 0
DIRNAME_PAD = 2
SEP = '-'

# Different thread counts give the two branches genuinely different rounding,
# which for a chaotic system the length of this test's run decorrelates into
# two unrelated trajectories -- exactly what makes "did the wrong branch's
# frames survive" a checkable question rather than a coin flip.
BRANCH_A_ARGS = ('-nb', 'cpu', '-pme', 'cpu', '-ntomp', '2')
BRANCH_B_ARGS = ('-nb', 'cpu', '-pme', 'cpu', '-ntomp', '1')

# The merge round-trips through compressed .xtc, which costs a bit of precision.
COORD_ATOL_NM = 2e-3
# Far above that round-trip noise, far below the ~0.3-1 nm the two branches
# reach once chaos has had BRANCH_STEP*2 steps to amplify a rounding
# difference into two unrelated configurations.
DIVERGED_NM = 0.02


def pbc_max_displacement(a, b, box):
    """Largest per-atom distance between two frames, the short way round the
    box, so atoms GROMACS wrapped to the far side don't dominate the answer.
    """
    delta = a - b
    delta -= box * np.round(delta / box)
    return float(np.linalg.norm(delta, axis=-1).max())


def run_gen(config, mdrun_args):
    """gmx_generation with these mdrun_args; None if it left the generation
    incomplete, else the trajectory it returned."""
    try:
        return gs.gmx_generation(**dict(config, mdrun_args=mdrun_args))
    except gs.GenIncomplete:
        return None


def main(gmx_bin=harness.GMX_BIN):
    suite = Suite('part_files')
    work = harness.workdir('part_files')
    structure, topology, mdp = harness.build_water_system(work, gmx_bin=gmx_bin)

    common = dict(
        traj_dir_top_level=str(work / 'farm'), top_fn=str(topology),
        seed_index=SEED_INDEX, clone_index=CLONE_INDEX, gen_index=GEN_INDEX,
        title='part-rewind', seed_fn=str(structure), structure_fn=str(structure),
        mdp_fn=str(mdp), dirname_pad=DIRNAME_PAD, sep=SEP, traj_name='prod',
        traj_suffix='.xtc', restart_name='state.cpt', steps=STEPS_PER_GEN,
        steps_per_gen=STEPS_PER_GEN, write_interval=WRITE_INTERVAL,
        target_step=harness.target_step(GEN_INDEX, STEPS_PER_GEN),
        temperature=300, gen_seed_base=99, gmx_bin=gmx_bin, grompp_maxwarn=3,
        new_velocities=True, append=False)

    gen_dir = mdfarmer.dir_seeds_clones_gens(
        Path(work / 'farm'), SEED_INDEX, CLONE_INDEX, GEN_INDEX, DIRNAME_PAD,
        sep=SEP)
    own_cpt = gen_dir / common['restart_name']

    suite.section('building the rewind: part0001, then an abandoned branch')
    traj = run_gen(common, BRANCH_A_ARGS + ('-nsteps', str(REWIND_STEP)))
    suite.check('the first launch stops at the rewind point, incomplete',
               traj is None and gs.checkpoint_step(own_cpt, gmx_bin=gmx_bin)
               == REWIND_STEP)

    early_cpt = work / 'early.cpt'
    shutil.copy(own_cpt, early_cpt)

    traj = run_gen(common, BRANCH_A_ARGS + ('-nsteps', str(BRANCH_STEP)))
    suite.check('the abandoned branch writes a second part, still incomplete',
               traj is None and gs.checkpoint_step(own_cpt, gmx_bin=gmx_bin)
               == REWIND_STEP + BRANCH_STEP)

    traj = run_gen(common, BRANCH_A_ARGS + ('-nsteps', str(BRANCH_STEP)))
    suite.check('the abandoned branch writes a third part, still incomplete',
               traj is None and gs.checkpoint_step(own_cpt, gmx_bin=gmx_bin)
               == REWIND_STEP + 2 * BRANCH_STEP)

    parts_before_rewind = gs.part_files(gen_dir)
    suite.check('three parts have accumulated before the rewind',
               len(parts_before_rewind) == 3,
               f'-> {[p.name for p in parts_before_rewind]}')
    # Keep the abandoned branch's own record of what it wrote, independent of
    # however the fix ends up moving the on-disk files aside.
    abandoned_reference = work / 'abandoned-branch-reference.xtc'
    shutil.copy(parts_before_rewind[-1], abandoned_reference)

    suite.section('rewinding: restore the early checkpoint and relaunch')
    shutil.copy(early_cpt, own_cpt)
    traj = run_gen(common, BRANCH_B_ARGS)
    suite.check('the relaunch resumes from the rewound checkpoint and finishes',
               traj is not None)

    # The generation finished, so its own parts were merged and deleted: what a
    # finished generation holds is one trajectory, not the launches it took.
    parts_after_rewind = gs.part_files(gen_dir)
    suite.check('a finished generation keeps no parts on the merge glob',
               not parts_after_rewind,
               f'-> {[p.name for p in parts_after_rewind]}')
    suite.check('and holds exactly one trajectory', Path(traj).is_file())
    abandoned_on_disk = sorted(gen_dir.glob(
        f'{gs.ABANDONED_PART_PREFIX}*{common["traj_suffix"]}'))
    # Moved aside before the relaunch, so never on the glob the merge deleted.
    suite.check('the abandoned parts are kept on disk, just moved aside',
               len(abandoned_on_disk) == 2,
               f'-> {[p.name for p in abandoned_on_disk]}')

    suite.section('the merged trajectory is exactly contiguous')
    merged = md.load(str(traj), top=str(structure))
    times = merged.time
    n_expected = STEPS_PER_GEN // WRITE_INTERVAL + 1
    suite.check('frame count matches the full step budget',
               len(times) == n_expected, f'-> {len(times)} vs {n_expected}')
    gaps = np.diff(times)
    suite.check('frame times are strictly increasing', bool(np.all(gaps > 0)),
               f'-> min gap {gaps.min():.4f} ps')
    expected_gap = WRITE_INTERVAL * DT_PS
    suite.check('frame spacing is uniform',
               bool(np.allclose(gaps, expected_gap, atol=1e-4)),
               f'-> {gaps.min():.4f}-{gaps.max():.4f} vs {expected_gap} ps')

    suite.section('no frame in the merge came from the abandoned branch')
    abandoned = md.load(str(abandoned_reference), top=str(structure))
    checked = 0
    for a_idx in range(len(abandoned.time)):
        t = abandoned.time[a_idx]
        m_idx = int(np.argmin(np.abs(times - t)))
        if abs(times[m_idx] - t) > 1e-4:
            continue                        # not one of the overlapping frames
        checked += 1
        gap = pbc_max_displacement(merged.xyz[m_idx], abandoned.xyz[a_idx],
                                   merged.unitcell_lengths[m_idx])
        suite.check(f'merged frame at t={t:.2f} ps differs from the '
                   f'abandoned branch', gap > DIVERGED_NM, f'-> {gap:.4f} nm')
    suite.check('the abandoned branch\'s time range was actually checked',
               checked == len(abandoned.time), f'-> {checked}/{len(abandoned.time)}')

    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
