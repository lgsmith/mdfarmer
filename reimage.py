"""Reimage trajectories so molecules are not broken by periodic boundaries.

GROMACS writes whatever coordinates the integrator holds. Depending on rank
count, update-group size and molecule size, molecules in the raw ``.xtc`` may be
split across a periodic boundary (or may simply have wandered outside the box).
Either way, any per-molecule geometry computed on the raw trajectory -- radius of
gyration, RMSD, contacts, a picture -- is wrong. This module makes molecules
whole again.

Two backends, because the right tool depends on the unit cell:

  * ``'loos'`` -- LOOS, for **orthorhombic** ("box") cells only. LOOS's periodic
    box is a ``GCoord``: three numbers (``src/PeriodicBox.hpp``). It has no
    triclinic representation at all, and its readers silently keep only the
    diagonal of a triclinic box -- ``xtc.cpp`` builds the box from elements
    0/4/8 of the 3x3 and discards the rest, and ``gro.cpp`` parses only the
    first three fields of the nine-field box line. Fed a rhombic dodecahedron,
    LOOS reports a rectangular cell, raises nothing, and every minimum-image
    calculation downstream is quietly wrong. So this backend *refuses* to run on
    a non-orthorhombic cell rather than producing plausible garbage.

  * ``'trjconv'`` -- shells out to ``gmx trjconv -pbc mol -ur compact``, which
    handles triclinic cells correctly. This is the right backend for
    dodecahedral / truncated-octahedral / general triclinic boxes.

``reimage_trajectory`` picks between them from the actual box read off the
trajectory, so the common case needs no decision from the caller.

Molecule membership never comes from LOOS. A GROMACS ``.gro`` carries no bonds,
and a bondless LOOS model makes ``splitByMolecule()`` return **one group holding
the entire system** (``AtomicGroup.cpp``, "If no connectivity, just return the
entire group") -- reimaging then degenerates into a single system-wide
translation that looks like it worked. Bond connectivity alone is not enough
either: a TIP4P-ice virtual site is bonded to nothing, so connected components
over bonds strand every ``MW`` in its own "molecule" and fling virtual sites
across the box. Instead, molecule blocks are read from the GROMACS topology via
``openmm.app.GromacsTopFile``, whose chains reproduce the ``[ molecules ]``
section exactly (verified on a 4-site-water + ion system: sizes {304: 1, 4:
3240, 1: 1}, tiling the atom index range with no gaps).

Nothing here ever writes over its input. The orchestrator counts frames in the
raw trajectory to decide whether a generation is finished, so clobbering it --
or replacing it with a symlink to a stripped copy -- would corrupt the restart
bookkeeping of a run that is still in flight.
"""

import shutil
import subprocess as sp
from pathlib import Path

import numpy as np


# A box is treated as triclinic when any off-diagonal element of the 3x3 box
# matrix exceeds this, relative to the largest diagonal element. GROMACS writes
# exact zeros for a rectangular cell, so this only has to clear float noise.
TRICLINIC_RTOL = 1e-6

# Backend selectors for `reimage_trajectory`.
BACKEND_AUTO = 'auto'
BACKEND_LOOS = 'loos'
BACKEND_TRJCONV = 'trjconv'

# `gmx trjconv` flags. '-pbc mol' puts each molecule's centre of mass in the box
# and makes the molecule whole (it needs a .tpr for molecule definitions);
# '-ur compact' is what makes a dodecahedral / octahedral cell come out as the
# compact unit cell rather than the triclinic parallelepiped.
TRJCONV_PBC = 'mol'
TRJCONV_UR = 'compact'

# Default GROMACS binary. Sites that build an MPI-only GROMACS have `gmx_mpi`.
GMX_BIN = 'gmx'

# Appended to the input stem to name the output, e.g. prod.xtc -> prod-whole.xtc.
OUTPUT_TAG = '-whole'

# Longest plausible chemical bond, nm. 0.25 nm == the 2.5 Angstrom default of
# LOOS's long-bond-finder, which is the tool this check reimplements.
MAX_BOND = 0.25

# Frames read per chunk when scanning a trajectory, so a multi-GB .xtc never
# has to be resident.
SCAN_CHUNK = 200

# Index group fed to trjconv on stdin. 0 is 'System' in the default group set.
TRJCONV_OUTPUT_GROUP = '0'

# LOOS reads and writes Angstroms; GROMACS .gro/.xtc are nm.
ANGSTROM_PER_NM = 10.0


class BoxTypeError(ValueError):
    """Raised when a backend is asked to handle a cell it cannot represent."""


def _read_gro_box_vectors(structure_fn):
    """Parse the final line of a .gro file into a 3x3 box matrix (nm).

    The GROMACS .gro box line is ``v1x v2y v3z [v1y v1z v2x v2z v3x v3y]`` --
    three fields for a rectangular cell, nine for a triclinic one.
    """
    last = None
    with open(structure_fn) as f:
        for line in f:
            if line.strip():
                last = line
    if last is None:
        raise ValueError(f'{structure_fn} is empty; no box line to read')
    fields = [float(x) for x in last.split()]
    if len(fields) not in (3, 9):
        raise ValueError(
            f'{structure_fn}: box line has {len(fields)} fields, expected 3 '
            f'(rectangular) or 9 (triclinic): {last.strip()!r}')
    box = np.zeros((3, 3))
    box[0, 0], box[1, 1], box[2, 2] = fields[:3]
    if len(fields) == 9:
        box[0, 1], box[0, 2] = fields[3], fields[4]
        box[1, 0], box[1, 2] = fields[5], fields[6]
        box[2, 0], box[2, 1] = fields[7], fields[8]
    return box


def box_vectors(traj_fn=None, structure_fn=None,
                angstrom_per_nm=ANGSTROM_PER_NM):
    """Return the 3x3 box matrix (nm) for a trajectory, falling back to a
    structure file.

    The trajectory is authoritative -- under a barostat the cell drifts away
    from whatever the starting structure recorded -- so the first frame is read
    when a trajectory is given. mdtraj is used because it reads the full 3x3;
    LOOS deliberately is not, since it would hand back the diagonal of a
    triclinic cell without complaint and defeat the whole point of this check.
    """
    if traj_fn is not None:
        vectors = _box_from_first_frame(
            Path(traj_fn), angstrom_per_nm=angstrom_per_nm)
        if vectors is not None:
            return vectors
    if structure_fn is None:
        raise ValueError(
            f'could not read a periodic box from {traj_fn} and no '
            'structure_fn was given to fall back on')
    return _read_gro_box_vectors(structure_fn)


def _box_from_first_frame(traj_p, angstrom_per_nm=ANGSTROM_PER_NM):
    """3x3 box (nm) from a trajectory's first frame, or None if it carries none.

    The raw mdtraj readers return different tuples per format and in different
    units, so dispatch explicitly rather than guessing from array shapes: XTC
    gives ``(xyz, time, step, box)`` with the box already a 3x3 in nm, while DCD
    gives ``(xyz, cell_lengths, cell_angles)`` with lengths in Angstroms.
    """
    import mdtraj
    from mdtraj.utils import lengths_and_angles_to_box_vectors
    suffix = traj_p.suffix.lower()
    with mdtraj.open(str(traj_p)) as fh:
        frame = fh.read(1)
    if suffix == '.xtc':
        box = np.asarray(frame[3])[0]
        return None if np.allclose(box, 0.0) else box
    if suffix == '.dcd':
        lengths = np.asarray(frame[1])
        angles = np.asarray(frame[2])
        if lengths.size == 0 or np.allclose(lengths, 0.0):
            return None
        lengths = lengths[0] / angstrom_per_nm
        angles = angles[0]
        return np.array(lengths_and_angles_to_box_vectors(*lengths, *angles))
    raise ValueError(
        f'{traj_p.suffix} is not a trajectory format mdfarmer reads boxes '
        'from; pass structure_fn instead')


def is_orthorhombic(box, triclinic_rtol=TRICLINIC_RTOL):
    """True when the 3x3 box matrix is rectangular to within `triclinic_rtol`."""
    box = np.asarray(box, dtype=float)
    scale = np.abs(np.diag(box)).max()
    if scale == 0.0:
        raise ValueError(f'degenerate box with zero diagonal: {box!r}')
    off_diagonal = box - np.diag(np.diag(box))
    return bool(np.abs(off_diagonal).max() <= triclinic_rtol * scale)


def molecule_ranges(top_fn, include_dir=None):
    """Return ``[(start, stop), ...]`` atom index ranges, one per molecule.

    Read from the GROMACS topology through ``openmm.app.GromacsTopFile``: its
    chains are one-per-molecule-instance and reproduce the ``[ molecules ]``
    section, including molecules whose atoms carry no bonds at all (bare ions,
    4-site-water virtual sites). Ranges are returned sorted and are verified to
    tile ``[0, n_atoms)`` without gaps or overlaps, which is the invariant the
    LOOS backend relies on.
    """
    from openmm import app
    top_p = Path(top_fn).resolve()
    if top_p.suffix != '.top':
        raise ValueError(
            f'molecule_ranges needs a GROMACS .top, got {top_p.name}. The '
            'topology is the only place molecule blocks are recorded; a .gro '
            'carries no connectivity.')
    kwargs = {} if include_dir is None else {'includeDir': str(include_dir)}
    # A .top's `#include "./ff/..."` lines resolve relative to the process cwd,
    # so parse from the topology's own directory.
    topology = app.GromacsTopFile(str(top_p), **kwargs).topology

    ranges = []
    for chain in topology.chains():
        indices = [a.index for a in chain.atoms()]
        if not indices:
            continue
        start, stop = min(indices), max(indices) + 1
        if stop - start != len(indices):
            raise ValueError(
                f'{top_p.name}: molecule spanning atoms [{start}, {stop}) is '
                f'not a contiguous index range ({len(indices)} atoms); '
                'mdfarmer cannot reimage this topology.')
        ranges.append((start, stop))
    ranges.sort()

    n_atoms = topology.getNumAtoms()
    cursor = 0
    for start, stop in ranges:
        if start != cursor:
            raise ValueError(
                f'{top_p.name}: molecule blocks do not tile the atom range -- '
                f'expected next molecule to start at {cursor}, got {start}.')
        cursor = stop
    if cursor != n_atoms:
        raise ValueError(
            f'{top_p.name}: molecule blocks cover {cursor} atoms but the '
            f'topology has {n_atoms}.')
    return ranges


def bond_pairs(top_fn, include_dir=None):
    """``(n_bonds, 2)`` array of bonded atom index pairs from a GROMACS topology.

    Connectivity comes from the topology rather than from a distance cutoff on
    the first frame: inferring bonds geometrically is exactly backwards when the
    thing you are trying to detect is a frame whose geometry is wrong.
    """
    from openmm import app
    top_p = Path(top_fn).resolve()
    kwargs = {} if include_dir is None else {'includeDir': str(include_dir)}
    topology = app.GromacsTopFile(str(top_p), **kwargs).topology
    pairs = np.array([[a.index, b.index] for a, b in topology.bonds()],
                     dtype=int)
    if not len(pairs):
        raise ValueError(
            f'{top_p.name} yielded no bonds; a bond-length check against it '
            'would pass vacuously.')
    return pairs


def check_bond_lengths(traj_fn, top_fn=None, pairs=None, max_bond=MAX_BOND,
                       stride=1, scan_chunk=SCAN_CHUNK, stop_early=True):
    """Find bonded pairs stretched further than a chemical bond can reach.

    This is the physical invariant that says whether a trajectory is imaged
    correctly, and it is deliberately independent of *which* tool did the
    imaging -- a reimaging pass agreeing with some other reimaging pass proves
    only that both made the same choice. A molecule split across a periodic
    boundary shows up here as a bond roughly one box-length long.

    Distances are computed **without** the minimum-image convention, on purpose:
    a PBC-aware distance is short for a split molecule by construction, which is
    precisely the defect being hunted. (Same reasoning as LOOS's
    ``long-bond-finder``, whose 2.5 Angstrom default `max_bond` this matches.)

    Returns ``(n_violations, violations)`` where `violations` is a list of
    ``(frame_index, atom_i, atom_j, length_nm)``. With `stop_early` the scan
    returns on the first bad frame, which is much faster for a pass/fail gate.
    """
    import mdtraj
    traj_p = Path(traj_fn)
    if pairs is None:
        if top_fn is None:
            raise ValueError('check_bond_lengths needs either pairs or top_fn')
        pairs = bond_pairs(top_fn)
    left, right = pairs[:, 0], pairs[:, 1]

    violations = []
    frame_index = 0
    with mdtraj.open(str(traj_p)) as fh:
        while True:
            chunk = fh.read(scan_chunk)
            xyz = np.asarray(chunk[0])
            if xyz.size == 0:
                break
            if xyz.shape[1] <= max(left.max(), right.max()):
                raise ValueError(
                    f'{traj_p} has {xyz.shape[1]} atoms but the topology '
                    f'describes at least {max(left.max(), right.max()) + 1}')
            sl = slice(None, None, stride)
            block = xyz[sl]
            lengths = np.linalg.norm(block[:, left, :] - block[:, right, :],
                                     axis=-1)
            bad = np.argwhere(lengths > max_bond)
            for row, bond in bad:
                violations.append((frame_index + int(row) * stride,
                                   int(left[bond]), int(right[bond]),
                                   float(lengths[row, bond])))
            if violations and stop_early:
                return len(violations), violations
            frame_index += xyz.shape[0]
    return len(violations), violations


def check_anchor_distances(traj_fn, ranges, structure_fn=None,
                           scan_chunk=SCAN_CHUNK):
    """How far the furthest atom sits from its own molecule's first atom.

    This measures exactly the assumption LOOS's ``mergeImage()`` makes. It
    minimum-images every atom of a molecule against that molecule's *first*
    atom, so it is correct only while no atom is more than half a box edge from
    that anchor. Past that, an atom that is legitimately far from the anchor
    gets wrapped to the wrong image and the molecule is quietly mangled -- with
    no error, and often with bond lengths that still look fine because the
    *bonded* neighbours moved together.

    Note this is not a correctness check on the trajectory: a molecule may
    legitimately be larger than half the box. It is a check on whether the LOOS
    backend is inside its safe regime for that trajectory. When it is not, use
    the trjconv backend, which walks the bond graph instead of anchoring.

    It also catches atoms carrying no bonds at all -- a TIP4P-ice ``MW`` virtual
    site stranded a box-length from its own water passes every bond check there
    is, because it has no bonds to be long, but its anchor distance is enormous.
    """
    import mdtraj
    traj_p = Path(traj_fn)
    starts = np.array([start for start, _ in ranges])
    box = box_vectors(traj_fn=traj_p, structure_fn=structure_fn)
    # reimageByAtom() wraps each Cartesian component independently, so the
    # limit is per axis against that axis' own edge -- not the vector norm
    # against the shortest edge, which would flag molecules that are actually
    # fine in a box with unequal edges.
    edges = np.abs(np.diag(box))
    limits = edges / 2.0

    # Repeating the anchor for every atom turns "displacement from my molecule's
    # first atom" into one vectorised subtraction over the whole frame.
    counts = np.array([stop - start for start, stop in ranges])
    anchor_of = np.repeat(starts, counts)

    worst, worst_frame, worst_atom, worst_axis = 0.0, -1, -1, -1
    frame_index = 0
    with mdtraj.open(str(traj_p)) as fh:
        while True:
            xyz = np.asarray(fh.read(scan_chunk)[0])
            if xyz.size == 0:
                break
            # fraction of each axis' half-edge used up, so axes compare directly
            offsets = np.abs(xyz - xyz[:, anchor_of, :])
            usage = offsets / limits
            flat = int(np.argmax(usage))
            here = float(usage.flat[flat])
            if here > worst:
                worst = here
                n_atoms, n_axes = usage.shape[1], usage.shape[2]
                worst_frame = frame_index + flat // (n_atoms * n_axes)
                worst_atom = (flat // n_axes) % n_atoms
                worst_axis = flat % n_axes
            frame_index += xyz.shape[0]
    return {'worst_half_box_fraction': worst,
            'max_anchor_offset': worst * float(limits[worst_axis])
            if worst_axis >= 0 else 0.0,
            'limit': float(limits[worst_axis]) if worst_axis >= 0 else 0.0,
            'frame': worst_frame, 'atom': worst_atom, 'axis': worst_axis,
            'loos_safe': worst < 1.0}


def reimage_with_loos(traj_fn, structure_fn, out_fn, top_fn=None,
                      ranges=None, center_selection=None,
                      skip_first_frame=False, verify=True,
                      max_bond=MAX_BOND,
                      triclinic_rtol=TRICLINIC_RTOL,
                      output_tag=OUTPUT_TAG):
    """Make molecules whole and wrap them into the primary image, with LOOS.

    Orthorhombic cells only -- see the module docstring for why a triclinic cell
    raises instead of being silently mishandled.

    `ranges` are ``(start, stop)`` atom index ranges, one per molecule; if not
    given they are read from `top_fn`. They are injected into the LOOS model as
    bonds, which is what makes ``splitByMolecule()`` return real molecules whose
    subgroups share -- and therefore track -- the parent's per-frame periodic
    box. Groups assembled by hand instead do *not* inherit the box: they report
    ``isPeriodic() == False`` with a sentinel 99999 cell, and reimaging against
    that is silently wrong.

    Per frame: ``mergeImage()`` on each molecule (unbreak it), then ``reimage()``
    (wrap its centroid into the cell). That order matters -- ``reimage()`` alone
    wraps a broken molecule by its meaningless centroid and leaves it broken.
    """
    import loos
    from loos import pyloos

    traj_p, out_p = Path(traj_fn), Path(out_fn)
    if out_p.resolve() == traj_p.resolve():
        raise ValueError(
            f'refusing to reimage {traj_p} onto itself; the orchestrator counts '
            f'frames in the raw trajectory to decide whether a generation is '
            f'complete. Write to a new file (e.g. {traj_p.stem}{output_tag}'
            f'{traj_p.suffix}).')

    box = box_vectors(traj_fn=traj_p, structure_fn=structure_fn)
    if not is_orthorhombic(box, triclinic_rtol=triclinic_rtol):
        raise BoxTypeError(
            f'{traj_p} has a non-orthorhombic box:\n{np.array2string(box, precision=4)}\n'
            'LOOS represents a periodic box as three numbers and would silently '
            'keep only the diagonal, making every minimum-image result wrong. '
            f'Use the {BACKEND_TRJCONV!r} backend (gmx trjconv -pbc mol -ur '
            'compact) for this cell.')

    if ranges is None:
        if top_fn is None:
            raise ValueError('reimage_with_loos needs either ranges or top_fn')
        ranges = molecule_ranges(top_fn)

    model = loos.createSystem(str(structure_fn))
    if len(model) != ranges[-1][1]:
        raise ValueError(
            f'{structure_fn} has {len(model)} atoms but the topology describes '
            f'{ranges[-1][1]}; they are not the same system.')

    # Star-bond each molecule so splitByMolecule() recovers exactly these
    # groups. A star is enough: nothing here walks the bond graph, it only needs
    # the connected components to match the molecule blocks.
    for start, stop in ranges:
        first = model[start]
        for i in range(start + 1, stop):
            other = model[i]
            first.addBond(other)
            other.addBond(first)

    molecules = model.splitByMolecule()
    if len(molecules) != len(ranges):
        raise ValueError(
            f'splitByMolecule() produced {len(molecules)} groups but the '
            f'topology has {len(ranges)} molecules; refusing to reimage '
            'against a grouping that does not match the topology.')

    center = None
    if center_selection:
        center = loos.selectAtoms(model, center_selection)
        if not len(center):
            raise ValueError(
                f'center_selection {center_selection!r} matched no atoms')

    writer = _loos_writer(out_p)
    traj = pyloos.Trajectory(str(traj_p), model)
    n_written = 0
    for index, _ in enumerate(traj):
        if skip_first_frame and index == 0:
            continue
        for molecule in molecules:
            molecule.mergeImage()
        if center is not None:
            # Anchor on a single atom before taking any centroid: the centroid
            # of a selection that is itself split across the boundary points
            # somewhere meaningless.
            model.translate(-center[0].coords())
            for molecule in molecules:
                molecule.reimage()
            model.translate(-center.centroid())
        for molecule in molecules:
            molecule.reimage()
        writer.writeFrame(model)
        n_written += 1
    del writer

    if verify:
        _verify_reimaged(out_p, top_fn=top_fn, ranges=ranges,
                         structure_fn=structure_fn, max_bond=max_bond)
    return out_p, n_written


def _verify_reimaged(out_p, top_fn=None, ranges=None, structure_fn=None,
                     max_bond=MAX_BOND):
    """Check a freshly-reimaged trajectory against the physical invariant.

    Reimaging by atom has many quiet failure modes, so the backend does not get
    to assume it succeeded. A bond stretched past `max_bond` means the output is
    wrong, and it is far better to hear that here than to find it in an
    analysis three weeks later.
    """
    if top_fn is not None:
        n_bad, violations = check_bond_lengths(
            out_p, top_fn=top_fn, max_bond=max_bond, stop_early=True)
        if n_bad:
            frame, i, j, length = violations[0]
            hint = ''
            if ranges is not None:
                margin = check_anchor_distances(
                    out_p, ranges, structure_fn=structure_fn)
                if not margin['loos_safe']:
                    hint = (
                        f' The furthest atom sits {margin["max_anchor_offset"]:.3f} nm '
                        f'along axis {margin["axis"]} from its molecule\'s anchor '
                        f'atom, past the {margin["limit"]:.3f} nm half-edge limit '
                        f'that LOOS\'s mergeImage() assumes, so this system is '
                        f'outside the LOOS backend\'s safe regime -- use the '
                        f'{BACKEND_TRJCONV!r} backend.')
            raise RuntimeError(
                f'{out_p} still has an overlong bond after reimaging: atoms '
                f'{i}-{j} are {length:.3f} nm apart in frame {frame} (limit '
                f'{max_bond} nm).{hint}')
    elif ranges is not None:
        margin = check_anchor_distances(out_p, ranges,
                                        structure_fn=structure_fn)
        if not margin['loos_safe']:
            print(f'[reimage] WARNING: cannot verify bond lengths without a '
                  f'.top, and the furthest atom is '
                  f'{margin["max_anchor_offset"]:.3f} nm from its anchor along '
                  f'axis {margin["axis"]} (half-edge limit '
                  f'{margin["limit"]:.3f} nm) -- LOOS\'s mergeImage() may have '
                  f'mis-wrapped atoms. Pass top_fn to get a real check, or use '
                  f'the trjconv backend.', flush=True)


def _loos_writer(out_p):
    suffix = out_p.suffix.lower()
    if suffix == '.xtc':
        return loos_xtc_writer(out_p)
    if suffix == '.dcd':
        import loos
        return loos.DCDWriter(str(out_p))
    raise ValueError(
        f'{out_p.suffix} is not a format the LOOS backend writes; use .xtc or '
        '.dcd')


def loos_xtc_writer(out_p):
    import loos
    return loos.XTCWriter(str(out_p))


def reimage_with_trjconv(traj_fn, tpr_fn, out_fn, gmx_bin=GMX_BIN,
                         pbc=TRJCONV_PBC, ur=TRJCONV_UR,
                         output_group=TRJCONV_OUTPUT_GROUP,
                         center_group=None, index_fn=None,
                         skip_first_frame=False, output_tag=OUTPUT_TAG):
    """Make molecules whole with ``gmx trjconv``, for any cell including triclinic.

    ``-pbc mol -ur compact`` is the combination GROMACS's own help recommends
    for a triclinic cell; ``-pbc mol`` needs the ``.tpr`` because that is where
    molecule definitions live (a ``.gro`` would silently give ``-pbc atom``
    behaviour for molecule purposes). ``-fit`` is deliberately not offered here:
    GROMACS documents that it cannot be combined with ``-pbc`` in one pass, and
    doing both at once silently produces the wrong answer -- fit afterwards, in
    a second call.

    trjconv reads its group selection from stdin, so the groups are piped rather
    than typed: `center_group` first when centring, then `output_group`.
    """
    traj_p, out_p = Path(traj_fn), Path(out_fn)
    if out_p.resolve() == traj_p.resolve():
        raise ValueError(
            f'refusing to reimage {traj_p} onto itself; write to a new file '
            f'(e.g. {traj_p.stem}{output_tag}{traj_p.suffix}).')
    if not Path(tpr_fn).is_file():
        raise FileNotFoundError(
            f'{tpr_fn}: gmx trjconv -pbc {pbc} needs the run .tpr for molecule '
            'definitions')

    cmd = [gmx_bin, 'trjconv', '-s', str(tpr_fn), '-f', str(traj_p),
           '-o', str(out_p), '-pbc', pbc, '-ur', ur]
    if index_fn:
        cmd += ['-n', str(index_fn)]
    if center_group is not None:
        cmd += ['-center']
    if skip_first_frame:
        # -b is a time, and trjconv's own frame 0 is the duplicate of the
        # previous generation's last frame; ask for everything strictly after it.
        first_time = _first_frame_time(traj_p)
        cmd += ['-b', repr(first_time + _frame_spacing(traj_p))]

    groups = [] if center_group is None else [str(center_group)]
    groups.append(str(output_group))
    stdin = '\n'.join(groups) + '\n'

    print('[reimage]', ' '.join(cmd), f'<<< {groups}', flush=True)
    result = sp.run(cmd, input=stdin, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            f'gmx trjconv exited {result.returncode}\n'
            f'--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}')
    if not out_p.is_file() or out_p.stat().st_size == 0:
        raise RuntimeError(
            f'gmx trjconv reported success but {out_p} is missing or empty\n'
            f'{result.stderr}')
    return out_p


def _frame_times(traj_p, limit=2):
    """First `limit` frame times (ps) of an .xtc."""
    import mdtraj
    if traj_p.suffix.lower() != '.xtc':
        raise ValueError(
            f'{traj_p.suffix}: frame times for skip_first_frame are only read '
            'from .xtc')
    with mdtraj.open(str(traj_p)) as fh:
        frame = fh.read(limit)
    return np.asarray(frame[1], dtype=float)


def _first_frame_time(traj_p):
    times = _frame_times(traj_p, limit=1)
    if times.size < 1:
        raise ValueError(f'{traj_p} has no frames')
    return float(times[0])


def _frame_spacing(traj_p):
    times = _frame_times(traj_p, limit=2)
    if times.size < 2:
        raise ValueError(
            f'{traj_p} has fewer than two frames; cannot infer frame spacing '
            'to skip the duplicate boundary frame')
    return float(times[1] - times[0])


def reimage_trajectory(traj_fn, out_fn=None, structure_fn=None, top_fn=None,
                       tpr_fn=None, backend=BACKEND_AUTO,
                       center_selection=None, center_group=None,
                       index_fn=None, skip_first_frame=False,
                       gmx_bin=GMX_BIN, triclinic_rtol=TRICLINIC_RTOL,
                       output_tag=OUTPUT_TAG, pbc=TRJCONV_PBC, ur=TRJCONV_UR):
    """Reimage `traj_fn`, choosing the backend from the trajectory's own box.

    With ``backend='auto'`` an orthorhombic cell goes to LOOS and anything else
    to ``gmx trjconv``, which is the split the two tools' capabilities actually
    dictate. Pass ``backend='trjconv'`` to force GROMACS even for a rectangular
    box (e.g. when no ``.top`` is available, or to match a colleague's pipeline).

    Returns the output path. Never modifies `traj_fn`.
    """
    traj_p = Path(traj_fn)
    if out_fn is None:
        out_fn = traj_p.with_name(f'{traj_p.stem}{output_tag}{traj_p.suffix}')

    if backend == BACKEND_AUTO:
        box = box_vectors(traj_fn=traj_p, structure_fn=structure_fn)
        orthorhombic = is_orthorhombic(box, triclinic_rtol=triclinic_rtol)
        backend = BACKEND_LOOS if orthorhombic else BACKEND_TRJCONV
        print(f'[reimage] box is '
              f'{"orthorhombic" if orthorhombic else "triclinic"}; '
              f'using the {backend} backend', flush=True)
        if backend == BACKEND_LOOS and structure_fn is None:
            print('[reimage] no structure_fn given for the LOOS backend; '
                  'falling back to trjconv', flush=True)
            backend = BACKEND_TRJCONV

    if backend == BACKEND_LOOS:
        out_p, _ = reimage_with_loos(
            traj_p, structure_fn, out_fn, top_fn=top_fn,
            center_selection=center_selection,
            skip_first_frame=skip_first_frame,
            triclinic_rtol=triclinic_rtol, output_tag=output_tag)
        return out_p
    if backend == BACKEND_TRJCONV:
        return reimage_with_trjconv(
            traj_p, tpr_fn, out_fn, gmx_bin=gmx_bin, pbc=pbc, ur=ur,
            center_group=center_group, index_fn=index_fn,
            skip_first_frame=skip_first_frame, output_tag=output_tag)
    raise ValueError(
        f'unknown backend {backend!r}; choose from '
        f'{[BACKEND_AUTO, BACKEND_LOOS, BACKEND_TRJCONV]}')


def reimage_gen_dir(gen_dir, config=None, backend=BACKEND_AUTO,
                    center_selection=None, center_group=None,
                    skip_first_frame=None, gmx_bin=GMX_BIN,
                    output_tag=OUTPUT_TAG, triclinic_rtol=TRICLINIC_RTOL):
    """Reimage one mdfarmer generation directory in place-adjacent fashion.

    Reads ``config.json`` from `gen_dir` for the trajectory name, the topology
    and the starting structure, and writes ``<traj_name><output_tag><suffix>``
    beside the raw trajectory. The raw trajectory is left untouched so a
    generation that is still being resumed keeps its frame accounting intact.

    `skip_first_frame` defaults to "drop it for every generation after the
    first": GROMACS writes an output frame at step 0 of every run, so generation
    N>0 opens with an exact duplicate of generation N-1's final frame.
    """
    import json
    gen_p = Path(gen_dir)
    if config is None:
        config = json.loads((gen_p / 'config.json').read_text())

    traj_p = (gen_p / config['traj_name']).with_suffix(config['traj_suffix'])
    if not traj_p.is_file():
        raise FileNotFoundError(f'{traj_p}: nothing to reimage')
    if skip_first_frame is None:
        skip_first_frame = config.get('gen_index', 0) > 0

    tpr_p = gen_p / 'prod.tpr'
    return reimage_trajectory(
        traj_p,
        structure_fn=config.get('structure_fn'),
        top_fn=config.get('top_fn'),
        tpr_fn=tpr_p if tpr_p.is_file() else None,
        backend=backend, center_selection=center_selection,
        center_group=center_group, skip_first_frame=skip_first_frame,
        gmx_bin=gmx_bin, triclinic_rtol=triclinic_rtol, output_tag=output_tag)


def have_gmx(gmx_bin=GMX_BIN):
    """True when `gmx_bin` is runnable, so callers can pick a backend up front."""
    return shutil.which(gmx_bin) is not None
