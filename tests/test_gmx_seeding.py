"""Generation 0 starts from this generation's own seed, not a shared structure.

Seeds are allowed to differ in topology, and an adaptive scheme may reseed a
later generation from a configuration it picked, so grompp -c has to follow
seed_fn rather than the campaign-wide structure_fn.
"""
import shutil
import sys

import numpy as np

import harness
from harness import Suite

from mdfarmer import gmx_simulate as gs

STEPS_PER_GEN = 200
WRITE_INTERVAL = 100
SHIFT_NM = 0.37          # how far the alternative seed is displaced


def shifted_copy(source, dest, shift_nm=SHIFT_NM):
    """A .gro with every atom moved, so its coordinates are recognisable."""
    lines = source.read_text().splitlines()
    header, count, body, box = lines[0], lines[1], lines[2:-1], lines[-1]
    moved = []
    for line in body:
        x, y, z = (float(line[20:28]), float(line[28:36]), float(line[36:44]))
        moved.append(f'{line[:20]}{x + shift_nm:8.3f}{y:8.3f}{z:8.3f}')
    dest.write_text('\n'.join([header, count] + moved + [box]) + '\n')
    return dest


def first_frame(traj, structure):
    import mdtraj as md
    frame = md.load(str(traj), top=str(structure))
    return frame.xyz[0], frame.unitcell_lengths[0]


def image_gap(coords, reference, box):
    """Largest per-atom distance, measured the short way round the box.

    Plain subtraction would be dominated by atoms GROMACS wrapped to the other
    side, which says nothing about which file the run started from.
    """
    delta = coords - reference
    delta -= box * np.round(delta / box)
    return float(np.abs(delta).max())


def main(steps_per_gen=STEPS_PER_GEN, write_interval=WRITE_INTERVAL,
         shift_nm=SHIFT_NM, gmx_bin=harness.GMX_BIN):
    suite = Suite('gmx_seeding')
    work = harness.workdir('gmx_seeding')
    structure, topology, mdp = harness.build_water_system(work, gmx_bin=gmx_bin)
    other = shifted_copy(structure, work / 'other.gro', shift_nm=shift_nm)

    common = dict(
        traj_dir_top_level=str(work / 'farm'), top_fn=str(topology),
        seed_index=0, clone_index=0, title='seeding',
        structure_fn=str(structure), mdp_fn=str(mdp),
        dirname_pad=2, sep='-', traj_name='prod', traj_suffix='.xtc',
        restart_name='state.cpt', steps=steps_per_gen,
        steps_per_gen=steps_per_gen, write_interval=write_interval,
        temperature=300, gen_seed_base=42, gmx_bin=gmx_bin, grompp_maxwarn=3,
        mdrun_args=('-nb', 'cpu', '-pme', 'cpu', '-ntomp', '2'), append=False)

    suite.section('generation 0 follows seed_fn, not structure_fn')
    traj = gs.gmx_generation(gen_index=0, seed_fn=str(other),
                             new_velocities=True, **common)
    started_at, box = first_frame(traj, structure)
    from_seed, _ = first_frame(other, other)
    from_structure, _ = first_frame(structure, structure)
    gap_seed = image_gap(started_at, from_seed, box)
    gap_structure = image_gap(started_at, from_structure, box)
    print(f'   first frame is {gap_seed:.4f} nm from seed_fn, '
          f'{gap_structure:.4f} nm from structure_fn', flush=True)
    suite.check('the run starts from the coordinates in seed_fn',
                gap_seed < shift_nm / 2 < gap_structure,
                f'-> {gap_seed:.4f} vs {gap_structure:.4f} nm')

    suite.section('structure_fn is still the fallback')
    traj = gs.gmx_generation(
        gen_index=0, seed_fn=str(work / 'farm' / 'nothing.cpt'),
        new_velocities=True,
        **dict(common, clone_index=1, traj_dir_top_level=str(work / 'farm2')))
    started_at, box = first_frame(traj, structure)
    gap_seed = image_gap(started_at, from_seed, box)
    gap_structure = image_gap(started_at, from_structure, box)
    print(f'   first frame is {gap_seed:.4f} nm from seed_fn, '
          f'{gap_structure:.4f} nm from structure_fn', flush=True)
    suite.check('a checkpoint seed_fn falls back to structure_fn',
                gap_structure < shift_nm / 2 < gap_seed,
                f'-> {gap_structure:.4f} vs {gap_seed:.4f} nm')

    suite.section('neither one available')
    try:
        gs.gmx_generation(
            gen_index=0, seed_fn='', new_velocities=True,
            **dict(common, clone_index=2, structure_fn=None,
                   traj_dir_top_level=str(work / 'farm3')))
        suite.check('a missing structure is refused', False, '-> no exception')
    except ValueError as exc:
        suite.check('a missing structure is refused', True, f'-> {str(exc)[:60]}')
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
