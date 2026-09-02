"""A rhombic dodecahedron is reimaged per molecule, not per atom.

LOOS cannot hold a triclinic cell at all -- its periodic box is three numbers --
so the mdtraj backend is the only one that can answer here, and everything it
claims is checked against first principles rather than against another tool.

The molecules of a short run come out of mdrun already whole, so they are broken
deliberately first: every atom is wrapped into the cell on its own, in fractional
coordinates, which is exactly what imaging a triclinic box atom by atom does.
The headline is the pair -- the per-molecule path leaves no overlong bond where
the per-atom path leaves thousands -- and then that the molecules it hands back
are, bond vector for bond vector, the ones mdrun wrote.
"""
import shutil
import subprocess as sp
import sys

import numpy as np
import mdtraj as md
from mdtraj.formats import XTCTrajectoryFile

import harness
from harness import Suite

from mdfarmer import reimage

# Edge of the dodecahedron, nm. Big enough for solvate to fit several hundred
# waters in, small enough that mdrun on two threads stays quick.
DODEC_BOX_NM = 3.0
BOX_TYPE = 'dodecahedron'

# Energy minimisation before the run: without it mdrun dies on "water molecules
# can not be settled" for a freshly solvated dodecahedron.
EM_STEPS = 500
EM_TOL = 200

STEPS = 400
WRITE_INTERVAL = 100

# The fixture runs on CPU, so it needs neither a GPU nor MPI ranks.
MDRUN_ARGS = ('-nb', 'cpu', '-pme', 'cpu', '-ntomp', '2')

# How far a bond vector may move, nm. An .xtc stores 1e-3 nm, and a bond vector
# is a difference of two coordinates from two separately written files.
BOND_VECTOR_ATOL_NM = 5e-3

# Slack, nm, on the half-edge a wrapped molecule has to sit inside. A centroid
# is an average of coordinates an .xtc stores to 1e-3 nm, and the wrap leaves
# molecules right against the bound.
WRAP_ATOL_NM = 2e-3


def gmx(cmd, cwd, gmx_bin=harness.GMX_BIN):
    """Run a gmx subcommand that builds the fixture, or Skip if it fails."""
    result = sp.run([str(gmx_bin), *[str(c) for c in cmd]], cwd=str(cwd),
                    text=True, capture_output=True)
    if result.returncode != 0:
        raise harness.Skip(
            f'gmx {cmd[0]} failed, so there is no dodecahedron to reimage:'
            f'\n{result.stderr[-2000:]}')
    return result


def em_mdp(steps=EM_STEPS, tol=EM_TOL, cutoff_nm=harness.WATER_CUTOFF_NM,
           nstlist=harness.WATER_NSTLIST):
    return '\n'.join((
        'integrator    = steep',
        f'nsteps        = {steps}',
        f'emtol         = {tol}',
        'cutoff-scheme = Verlet',
        f'nstlist       = {nstlist}',
        'coulombtype   = PME',
        f'rcoulomb      = {cutoff_nm}',
        'vdw-type      = cut-off',
        f'rvdw          = {cutoff_nm}',
        '',
    ))


def build_dodec_system(work, gmx_bin=harness.GMX_BIN, box_nm=DODEC_BOX_NM,
                       box_type=BOX_TYPE, steps=STEPS, mdrun_args=MDRUN_ARGS):
    """Solvate a dodecahedron and run it, returning (traj, structure, topology)."""
    top_dir = harness.require_gmx(gmx_bin=gmx_bin)
    shutil.copy(top_dir / harness.WATER_BOX_NAME, work / 'spc216.gro')
    (work / 'topol.top').write_text(harness.water_top())

    gmx(['editconf', '-f', 'spc216.gro', '-o', 'empty.gro', '-bt', box_type,
         '-box', box_nm, '-c'], work, gmx_bin=gmx_bin)
    gmx(['solvate', '-cp', 'empty.gro', '-cs', 'spc216.gro', '-o', 'dodec.gro',
         '-p', 'topol.top'], work, gmx_bin=gmx_bin)

    (work / 'em.mdp').write_text(em_mdp())
    gmx(['grompp', '-f', 'em.mdp', '-c', 'dodec.gro', '-p', 'topol.top',
         '-o', 'em.tpr', '-maxwarn', '3'], work, gmx_bin=gmx_bin)
    gmx(['mdrun', '-s', 'em.tpr', '-deffnm', 'em', '-nb', 'cpu', '-ntomp', '2'],
        work, gmx_bin=gmx_bin)

    (work / 'base.mdp').write_text(harness.water_mdp())
    gmx(['grompp', '-f', 'base.mdp', '-c', 'em.gro', '-p', 'topol.top',
         '-o', 'prod.tpr', '-maxwarn', '3'], work, gmx_bin=gmx_bin)
    gmx(['mdrun', '-s', 'prod.tpr', '-deffnm', 'prod', '-nsteps', steps,
         *mdrun_args], work, gmx_bin=gmx_bin)
    return work / 'prod.xtc', work / 'dodec.gro', work / 'topol.top'


def write_xtc(out_fn, frames, step):
    """Write frames keeping the source's MD step numbers, which save_xtc drops."""
    with XTCTrajectoryFile(str(out_fn), 'w') as handle:
        handle.write(frames.xyz, time=frames.time, step=step,
                     box=frames.unitcell_vectors)
    return out_fn


def wrap_each_atom(frames):
    """Wrap every atom into the cell on its own, breaking every molecule.

    Done in fractional coordinates, so it is a real triclinic wrap rather than
    a rectangular one that would leave a dodecahedron's molecules alone.
    """
    for i, box in enumerate(frames.unitcell_vectors):
        fractional = frames.xyz[i] @ np.linalg.inv(box)
        frames.xyz[i] = (fractional - np.floor(fractional)) @ box
    return frames


def image_per_atom(frames):
    """Image frames with every atom its own molecule, the way mdtraj would.

    On a topology with no bonds find_molecules() returns one molecule per atom
    and nothing is made whole; this is that, and it is the failure the
    per-molecule path exists to prevent.
    """
    atoms = [[atom] for atom in frames.topology.atoms]
    return frames.image_molecules(inplace=True, anchor_molecules=atoms[:1],
                                  other_molecules=atoms[1:], make_whole=False)


def read_clock(path):
    """(time, step) for a whole .xtc."""
    with md.open(str(path)) as handle:
        _, time, step, _ = handle.read()
    return np.asarray(time), np.asarray(step)


def bond_vectors(frames, pairs):
    """(n_frames, n_bonds, 3) end-to-end vector of every bond, unimaged."""
    return frames.xyz[:, pairs[:, 0], :] - frames.xyz[:, pairs[:, 1], :]


def outside_half_box(frames, ranges, anchor_index, atol=WRAP_ATOL_NM):
    """Count molecule-frames whose centroid is not wrapped in beside the anchor.

    Wrapping per molecule means every centroid ends up within half a box edge
    of the anchor's along each axis; a molecule left where it was is a whole
    box out.
    """
    missed = 0
    for i, box in enumerate(frames.unitcell_vectors):
        centroids = np.array([frames.xyz[i, start:stop].mean(axis=0)
                              for start, stop in ranges])
        offset = np.abs(centroids - centroids[anchor_index])
        missed += int((offset > np.diag(box) / 2 + atol).any(axis=1).sum())
    return missed


def chosen_backend(traj_fn, **kwargs):
    """Which backend reimage_trajectory's 'auto' picks, without running one."""
    seen = []

    def spy(name):
        return lambda traj, *a, **k: (seen.append(name), (traj, 0))[1]

    saved = reimage.reimage_with_loos, reimage.reimage_with_mdtraj
    reimage.reimage_with_loos = spy(reimage.BACKEND_LOOS)
    reimage.reimage_with_mdtraj = spy(reimage.BACKEND_MDTRAJ)
    try:
        reimage.reimage_trajectory(traj_fn, **kwargs)
    finally:
        reimage.reimage_with_loos, reimage.reimage_with_mdtraj = saved
    return seen[0]


def main(gmx_bin=harness.GMX_BIN, write_interval=WRITE_INTERVAL, steps=STEPS,
         bond_vector_atol_nm=BOND_VECTOR_ATOL_NM):
    suite = Suite('reimage_triclinic')
    work = harness.workdir('reimage_triclinic')
    gmx_top = harness.require_gmx(gmx_bin=gmx_bin)

    traj, structure, topology = build_dodec_system(work, gmx_bin=gmx_bin)
    ranges = reimage.molecule_ranges(str(topology), include_dir=str(gmx_top))
    pairs = reimage.bond_pairs(str(topology), include_dir=str(gmx_top))
    mdtop = reimage.mdtraj_topology(str(topology), include_dir=str(gmx_top))
    source = md.load(str(traj), top=mdtop)
    _, source_step = read_clock(traj)

    suite.section('the cell really is triclinic')
    box = reimage.box_vectors(traj_fn=traj)
    off_diagonal = np.abs(box - np.diag(np.diag(box))).max()
    print(f'   box\n{np.array2string(box, precision=4)}', flush=True)
    print(f'   {len(ranges)} molecules, {len(pairs)} bonds, '
          f'{source.n_frames} frames', flush=True)
    suite.check('the box is not rectangular',
                not reimage.is_orthorhombic(box),
                f'-> largest off-diagonal {off_diagonal:.4f} nm')
    try:
        reimage.reimage_with_loos(traj, str(structure), work / 'loos.xtc',
                                  ranges=ranges, verify=False)
        suite.check('LOOS refuses the cell rather than keeping its diagonal',
                    False, '-> no exception')
    except reimage.BoxTypeError as exc:
        suite.check('LOOS refuses the cell rather than keeping its diagonal',
                    reimage.BACKEND_MDTRAJ in str(exc),
                    f'-> {str(exc).splitlines()[0][:50]}')

    suite.section('molecules broken on purpose, then put back together')
    broken = write_xtc(work / 'broken.xtc',
                       wrap_each_atom(md.load(str(traj), top=mdtop)),
                       source_step)
    n_broken, _ = reimage.check_bond_lengths(broken, pairs=pairs)
    print(f'   wrapping every atom on its own: {n_broken} bonds over '
          f'{reimage.MAX_BOND} nm', flush=True)
    suite.check('the fixture really is broken to start with', n_broken > 0,
                f'-> {n_broken}')

    whole, n_written = reimage.reimage_with_mdtraj(
        broken, str(topology), work / 'broken-whole.xtc', ranges=ranges,
        pairs=pairs, include_dir=str(gmx_top), verify=True)
    n_whole, _ = reimage.check_bond_lengths(whole, pairs=pairs)
    suite.check('every frame is written', n_written == source.n_frames,
                f'-> {n_written}')
    suite.check('no bond is longer than a bond can be', n_whole == 0,
                f'-> {n_whole}')

    suite.section('the per-molecule guarantee, against per-atom imaging')
    per_atom = write_xtc(
        work / 'per-atom.xtc',
        image_per_atom(md.load(str(broken), top=mdtop)), source_step)
    n_per_atom, _ = reimage.check_bond_lengths(per_atom, pairs=pairs)
    print(f'   per-atom imaging: {n_per_atom} long bonds; per-molecule: '
          f'{n_whole}', flush=True)
    suite.check('imaging atom by atom leaves molecules broken', n_per_atom > 0,
                f'-> {n_per_atom}')
    suite.check('the check that passes per molecule fails per atom',
                n_per_atom > n_whole == 0)

    suite.section('the molecules handed back are the ones mdrun wrote')
    imaged = md.load(str(whole), top=mdtop)
    moved = np.abs(bond_vectors(imaged, pairs)
                   - bond_vectors(source, pairs)).max()
    print(f'   worst bond vector moved {moved:.2e} nm', flush=True)
    suite.check('every bond vector survives, not just its length',
                moved < bond_vector_atol_nm, f'-> {moved:.2e} nm')
    suite.check('the box of every frame is carried through unchanged',
                np.allclose(imaged.unitcell_vectors, source.unitcell_vectors,
                            atol=1e-6))

    suite.section('the molecules are wrapped, not merely made whole')
    anchor = reimage.largest_molecule(ranges)
    missed = outside_half_box(imaged, ranges, anchor)
    before = outside_half_box(md.load(str(broken), top=mdtop), ranges, anchor)
    print(f'   molecules more than half a box from the anchor: {missed} after '
          f'imaging, {before} before', flush=True)
    suite.check('every molecule is wrapped in beside the anchor',
                missed == 0, f'-> {missed}')
    suite.check('which the source trajectory was not', before > 0,
                f'-> {before}')

    suite.section('the step and time axis survives')
    source_time, _ = read_clock(broken)
    whole_time, whole_step = read_clock(whole)
    wanted = np.arange(0, steps + 1, write_interval)
    suite.check("the MD steps are exactly the run's write interval",
                np.array_equal(whole_step, wanted), f'-> {whole_step}')
    suite.check('frame times survive reimaging',
                np.allclose(whole_time, source_time, atol=1e-5))

    suite.section('a grouping that contradicts the bonds is refused')
    per_atom_ranges = [(i, i + 1) for i in range(mdtop.n_atoms)]
    try:
        reimage.check_bonds_within_molecules(pairs, per_atom_ranges)
        suite.check('one molecule per atom is refused', False,
                    '-> no exception')
    except ValueError as exc:
        suite.check('one molecule per atom is refused', True,
                    f'-> {str(exc)[:60]}')
    suite.check("the topology's own blocks are accepted",
                reimage.check_bonds_within_molecules(pairs, ranges) is None)

    suite.section('auto picks the backend from the box')
    rectangular = md.load(str(traj), top=mdtop)
    rectangular.unitcell_vectors = np.array(
        [np.diag(np.diag(b)) for b in rectangular.unitcell_vectors])
    write_xtc(work / 'rectangular.xtc', rectangular, source_step)
    picked = chosen_backend(traj, structure_fn=str(structure),
                            top_fn=str(topology), include_dir=str(gmx_top))
    suite.check('a triclinic box goes to mdtraj',
                picked == reimage.BACKEND_MDTRAJ, f'-> {picked}')
    picked = chosen_backend(work / 'rectangular.xtc',
                            structure_fn=str(structure),
                            top_fn=str(topology), include_dir=str(gmx_top))
    suite.check('the same system in a rectangular box goes to LOOS',
                picked == reimage.BACKEND_LOOS, f'-> {picked}')

    suite.section('the anchor is the largest molecule')
    suite.check('largest_molecule picks the biggest block',
                reimage.largest_molecule([(0, 3), (3, 10), (10, 12)]) == 1)
    suite.check('a system of one molecule type anchors on the first',
                reimage.largest_molecule(ranges) == 0)
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
