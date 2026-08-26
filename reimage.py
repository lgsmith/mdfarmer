"""Put molecules back together after periodic boundaries have split them.

GROMACS writes whatever coordinates the integrator is holding, so a molecule can
come out cut in half across a box face. Anything you then measure per molecule
is wrong: radius of gyration, RMSD, contacts, a picture.

Two backends. LOOS handles rectangular boxes and is refused anything else, since
its periodic box is three numbers and it keeps only the diagonal of a triclinic
one without complaining. Anything else goes to gmx trjconv, which does triclinic
correctly. reimage_trajectory reads the box and picks for you.

Which atoms make up a molecule always comes from the GROMACS topology, never
from LOOS. A .gro carries no bonds, and LOOS treats a model with no bonds as one
single molecule, so reimaging would quietly become one system-wide shift.

Nothing here overwrites its input: the orchestrator counts frames in the raw
trajectory to decide whether a generation has finished.
"""

import subprocess as sp
from pathlib import Path

import numpy as np

from . import utilities as util


# Largest off-diagonal box element, relative to the diagonal, still counted as
# rectangular. GROMACS writes exact zeros, so this only clears float noise.
TRICLINIC_RTOL = 1e-6

# Backend selectors for reimage_trajectory.
BACKEND_AUTO = 'auto'
BACKEND_LOOS = 'loos'
BACKEND_TRJCONV = 'trjconv'

# trjconv flags. '-pbc mol' makes each molecule whole and puts its centre of
# mass in the box; '-ur compact' gives a dodecahedron its compact shape.
TRJCONV_PBC = 'mol'
TRJCONV_UR = 'compact'

# Default GROMACS binary. Sites that build an MPI-only GROMACS have gmx_mpi.
GMX_BIN = 'gmx'

# What the runner calls a generation's tpr, when its config does not say.
TPR_NAME = 'prod.tpr'

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
    """Parse the last line of a .gro file into a 3x3 box matrix, in nm."""
    # The box line is v1x v2y v3z, plus six more fields when triclinic.
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
    """The 3x3 box in nm, from the trajectory if there is one, else the structure.

    The trajectory wins because a barostat moves the cell away from whatever the
    starting structure recorded. Read with mdtraj, which gives the full 3x3.
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
    """3x3 box in nm from a trajectory's first frame, or None if it has none."""
    # Each format is handled by name, not by guessing at array shapes: XTC gives
    # a 3x3 in nm, DCD gives lengths and angles in Angstroms.
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


def length_scale(traj_fn, angstrom_per_nm=ANGSTROM_PER_NM):
    """What to multiply this format's raw coordinates by to get nanometres.

    mdtraj.open() hands back whatever the file holds: nm for an .xtc, Angstroms
    for a .dcd. Everything here compares against nanometre thresholds.
    """
    return 1.0 if Path(traj_fn).suffix.lower() == '.xtc' else 1 / angstrom_per_nm


def is_orthorhombic(box, triclinic_rtol=TRICLINIC_RTOL):
    """True when the 3x3 box matrix is rectangular to within triclinic_rtol."""
    box = np.asarray(box, dtype=float)
    scale = np.abs(np.diag(box)).max()
    if scale == 0.0:
        raise ValueError(f'degenerate box with zero diagonal: {box!r}')
    off_diagonal = box - np.diag(np.diag(box))
    return bool(np.abs(off_diagonal).max() <= triclinic_rtol * scale)


def gromacs_topology(top_fn, include_dir=None):
    """The openmm Topology for a GROMACS .top.

    Its chains match the [ molecules ] section, including molecules whose atoms
    carry no bonds, which is what makes them usable as molecule blocks.
    """
    from openmm import app
    top_p = Path(top_fn).resolve()
    if top_p.suffix != '.top':
        raise ValueError(
            f'expected a GROMACS .top, got {top_p.name}. The topology is the '
            'only place molecule blocks are recorded; a .gro carries no '
            'connectivity.')
    # A .top's #include lines resolve relative to the process cwd, so a force
    # field that lives elsewhere needs include_dir.
    kwargs = {} if include_dir is None else {'includeDir': str(include_dir)}
    return app.GromacsTopFile(str(top_p), **kwargs).topology


def molecule_ranges(top_fn, include_dir=None):
    """[(start, stop), ...] atom index ranges, one per molecule in the topology.

    Read through openmm.app.GromacsTopFile, whose chains match the [ molecules ]
    section even for molecules with no bonds, such as bare ions. Checked to
    cover every atom exactly once, which the LOOS backend relies on.
    """
    top_p = Path(top_fn).resolve()
    topology = gromacs_topology(top_p, include_dir=include_dir)

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
                f'{top_p.name}: molecule blocks do not tile the atom range, '
                f'expected next molecule to start at {cursor}, got {start}.')
        cursor = stop
    if cursor != n_atoms:
        raise ValueError(
            f'{top_p.name}: molecule blocks cover {cursor} atoms but the '
            f'topology has {n_atoms}.')
    return ranges


def bond_pairs(top_fn, include_dir=None):
    """(n_bonds, 2) array of bonded atom index pairs, read from the topology."""
    # Not guessed from distances in the first frame: the thing being looked for
    # is a frame whose geometry is wrong, so its distances cannot be trusted.
    top_p = Path(top_fn).resolve()
    topology = gromacs_topology(top_p, include_dir=include_dir)
    pairs = np.array([[a.index, b.index] for a, b in topology.bonds()],
                     dtype=int)
    if not len(pairs):
        raise ValueError(
            f'{top_p.name} yielded no bonds; a bond-length check against it '
            'would pass vacuously.')
    return pairs


def check_bond_lengths(traj_fn, top_fn=None, pairs=None, max_bond=MAX_BOND,
                       scan_chunk=SCAN_CHUNK, stop_early=True):
    """Find bonds longer than a chemical bond can be, one per split molecule.

    This is what actually says whether a trajectory is imaged correctly, and it
    does not care which tool did the imaging. A molecule cut across a boundary
    shows up as a bond about one box length long.

    Returns (n_violations, [(frame, atom_i, atom_j, length_nm), ...]).
    stop_early returns on the first bad frame, for a quick pass/fail.
    """
    # Distances ignore the minimum image convention on purpose: it would make a
    # split molecule's bonds short, hiding the very thing being looked for.
    import mdtraj
    traj_p = Path(traj_fn)
    if pairs is None:
        if top_fn is None:
            raise ValueError('check_bond_lengths needs either pairs or top_fn')
        pairs = bond_pairs(top_fn)
    left, right = pairs[:, 0], pairs[:, 1]

    violations = []
    frame_index = 0
    to_nm = length_scale(traj_p)
    with mdtraj.open(str(traj_p)) as fh:
        while True:
            chunk = fh.read(scan_chunk)
            xyz = np.asarray(chunk[0]) * to_nm
            if xyz.size == 0:
                break
            if xyz.shape[1] <= max(left.max(), right.max()):
                raise ValueError(
                    f'{traj_p} has {xyz.shape[1]} atoms but the topology '
                    f'describes at least {max(left.max(), right.max()) + 1}')
            lengths = np.linalg.norm(xyz[:, left, :] - xyz[:, right, :],
                                     axis=-1)
            for row, bond in np.argwhere(lengths > max_bond):
                violations.append((frame_index + int(row),
                                   int(left[bond]), int(right[bond]),
                                   float(lengths[row, bond])))
            if violations and stop_early:
                return len(violations), violations
            frame_index += xyz.shape[0]
    return len(violations), violations


def check_anchor_distances(traj_fn, ranges, structure_fn=None,
                           scan_chunk=SCAN_CHUNK):
    """How far the furthest atom sits from its own molecule's first atom.

    LOOS's mergeImage() measures every atom of a molecule against that
    molecule's first atom, so it is right only while no atom is more than half a
    box edge away. This says whether that holds, which is a question about the
    LOOS backend rather than about the trajectory: a molecule is allowed to be
    bigger than half the box, but then trjconv has to do the reimaging.

    It also catches an unbonded atom stranded far from its molecule, which every
    bond-length check misses because it has no bonds to be long.
    """
    import mdtraj
    traj_p = Path(traj_fn)
    starts = np.array([start for start, _ in ranges])
    box = box_vectors(traj_fn=traj_p, structure_fn=structure_fn)
    # Checked per axis, since LOOS wraps each component on its own. A vector
    # norm against the shortest edge would flag molecules that are really fine.
    edges = np.abs(np.diag(box))
    limits = edges / 2.0

    # Repeating the anchor for every atom turns "displacement from my molecule's
    # first atom" into one vectorised subtraction over the whole frame.
    counts = np.array([stop - start for start, stop in ranges])
    anchor_of = np.repeat(starts, counts)

    worst, worst_frame, worst_atom, worst_axis = 0.0, -1, -1, -1
    frame_index = 0
    to_nm = length_scale(traj_p)
    with mdtraj.open(str(traj_p)) as fh:
        while True:
            xyz = np.asarray(fh.read(scan_chunk)[0]) * to_nm
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
                      ranges=None, include_dir=None, center_selection=None,
                      skip_first_frame=False, verify=True,
                      max_bond=MAX_BOND,
                      triclinic_rtol=TRICLINIC_RTOL,
                      output_tag=OUTPUT_TAG):
    """Make molecules whole and wrap them back into the box, with LOOS.

    Rectangular boxes only. ranges are the atom index ranges of each molecule,
    read from top_fn if not given, and are added to the LOOS model as bonds so
    that splitByMolecule() returns real molecules whose groups share the
    parent's box. A group built by hand instead reports no box at all.

    Per frame: mergeImage() on each molecule (unbreak it), then reimage()
    (wrap its centroid into the cell). That order matters: reimage() alone
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
        ranges = molecule_ranges(top_fn, include_dir=include_dir)

    model = loos.createSystem(str(structure_fn))
    if len(model) != ranges[-1][1]:
        raise ValueError(
            f'{structure_fn} has {len(model)} atoms but the topology describes '
            f'{ranges[-1][1]}; they are not the same system.')

    # Bond every atom of a molecule to its first atom, so splitByMolecule()
    # finds exactly these groups. Nothing here walks the bonds themselves.
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
    timing = util.frame_timing(traj_p)
    traj = pyloos.Trajectory(str(traj_p), model)
    n_written = 0
    for index, _ in enumerate(traj):
        if skip_first_frame and index == 0:
            continue
        for molecule in molecules:
            molecule.mergeImage()
        if center is not None:
            # Move to one atom first: the centroid of a selection that is
            # itself split across the boundary points nowhere useful.
            model.translate(-center[0].coords())
            for molecule in molecules:
                molecule.reimage()
            model.translate(-center.centroid())
        for molecule in molecules:
            molecule.reimage()
        if timing is None:
            writer.writeFrame(model)
        else:
            step0, steps_per_frame, time0, time_per_frame = timing
            writer.writeFrame(model, step0 + index * steps_per_frame,
                              time0 + index * time_per_frame)
        n_written += 1
    del writer

    if verify:
        _verify_reimaged(out_p, top_fn=top_fn, include_dir=include_dir,
                         ranges=ranges,
                         structure_fn=structure_fn, max_bond=max_bond)
    return out_p, n_written


def _verify_reimaged(out_p, top_fn=None, include_dir=None, ranges=None,
                     structure_fn=None, max_bond=MAX_BOND):
    """Check a just-reimaged trajectory, since LOOS can fail quietly.

    Two independent checks. A bond longer than max_bond says a molecule is still
    split. The anchor margin says whether mergeImage() was even entitled to an
    answer, and is the only one that catches an atom with no bonds stranded from
    its own molecule, which no bond length can be long enough to reveal.
    """
    margin = None
    if ranges is not None:
        margin = check_anchor_distances(out_p, ranges,
                                        structure_fn=structure_fn)
    unsafe = '' if margin is None or margin['loos_safe'] else (
        f' The furthest atom sits {margin["max_anchor_offset"]:.3f} nm along '
        f'axis {margin["axis"]} from its molecule\'s anchor atom, past the '
        f'{margin["limit"]:.3f} nm half-edge limit mergeImage() assumes, so '
        f'this system is outside the LOOS backend\'s safe regime. Use the '
        f'{BACKEND_TRJCONV!r} backend.')

    if top_fn is not None:
        n_bad, violations = check_bond_lengths(
            out_p, pairs=bond_pairs(top_fn, include_dir=include_dir),
            max_bond=max_bond, stop_early=True)
        if n_bad:
            frame, i, j, length = violations[0]
            raise RuntimeError(
                f'{out_p} still has an overlong bond after reimaging: atoms '
                f'{i}-{j} are {length:.3f} nm apart in frame {frame} (limit '
                f'{max_bond} nm).{unsafe}')
    if unsafe:
        raise RuntimeError(
            f'{out_p} was reimaged outside the LOOS backend\'s safe '
            f'regime.{unsafe}')
    if top_fn is None and margin is None:
        print(f'[reimage] {out_p} was written unchecked: pass top_fn or ranges '
              'to have it verified.', flush=True)


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
    """Make molecules whole with gmx trjconv, for any box including triclinic.

    -pbc mol needs the .tpr, which is where molecule definitions live. There is
    deliberately no -fit option: GROMACS cannot combine it with -pbc in one pass
    and gives a wrong answer if asked to, so fit in a second call.
    """
    traj_p, out_p = Path(traj_fn), Path(out_fn)
    if out_p.resolve() == traj_p.resolve():
        raise ValueError(
            f'refusing to reimage {traj_p} onto itself; write to a new file '
            f'(e.g. {traj_p.stem}{output_tag}{traj_p.suffix}).')
    if tpr_fn is None or not Path(tpr_fn).is_file():
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
    """First limit frame times (ps) of an .xtc."""
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
                       tpr_fn=None, backend=BACKEND_AUTO, include_dir=None,
                       center_selection=None, center_group=None,
                       index_fn=None, skip_first_frame=False,
                       gmx_bin=GMX_BIN, triclinic_rtol=TRICLINIC_RTOL,
                       output_tag=OUTPUT_TAG, pbc=TRJCONV_PBC, ur=TRJCONV_UR):
    """Reimage a trajectory, choosing the backend from its own box.

    'auto' sends a rectangular box to LOOS and anything else to gmx trjconv.
    Force 'trjconv' when there is no .top, or to match someone else's pipeline.
    Returns the output path, and never touches the input.
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
            include_dir=include_dir, center_selection=center_selection,
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
                    include_dir=None, center_selection=None, center_group=None,
                    skip_first_frame=None, gmx_bin=GMX_BIN,
                    output_tag=OUTPUT_TAG, triclinic_rtol=TRICLINIC_RTOL,
                    tpr_name=TPR_NAME):
    """Reimage one generation directory, writing beside the raw trajectory.

    Reads config.json for the trajectory name, topology and structure. The raw
    trajectory is left alone so a generation still being resumed keeps its frame
    count. skip_first_frame defaults to dropping it for every generation after
    the first, since GROMACS repeats the previous generation's last frame.
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

    tpr_p = gen_p / config.get('tpr_name', tpr_name)
    return reimage_trajectory(
        traj_p,
        structure_fn=config.get('structure_fn'),
        top_fn=config.get('top_fn'),
        tpr_fn=tpr_p if tpr_p.is_file() else None,
        backend=backend, include_dir=include_dir,
        center_selection=center_selection,
        center_group=center_group, skip_first_frame=skip_first_frame,
        gmx_bin=gmx_bin, triclinic_rtol=triclinic_rtol, output_tag=output_tag)


