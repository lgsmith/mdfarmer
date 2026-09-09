"""Put molecules back together after periodic boundaries have split them.

GROMACS writes whatever coordinates the integrator is holding, so a molecule can
come out cut in half across a box face. Anything you then measure per molecule
is wrong: radius of gyration, RMSD, contacts, a picture.

Two backends. LOOS handles rectangular boxes and is refused anything else, since
its periodic box is three numbers and it keeps only the diagonal of a triclinic
one without complaining. Anything else goes to mdtraj, which carries the full
3x3 cell. reimage_trajectory reads the box and picks for you.

Which atoms make up a molecule always comes from the GROMACS topology, never
from the structure file. A .gro carries no bonds; LOOS treats a model with no
bonds as one single molecule, so reimaging would quietly become one system-wide
shift, and mdtraj's find_molecules() makes every atom its own molecule, so
reimaging would become a per-atom wrap that breaks molecules irrecoverably.

Nothing here overwrites its input: the orchestrator counts frames in the raw
trajectory to decide whether a generation has finished.
"""

import os
from pathlib import Path

import numpy as np

from . import utilities as util


# Largest off-diagonal box element, relative to the diagonal, still counted as
# rectangular. GROMACS writes exact zeros, so this only clears float noise.
TRICLINIC_RTOL = 1e-6

# Backend selectors for reimage_trajectory.
BACKEND_AUTO = 'auto'
BACKEND_LOOS = 'loos'
BACKEND_MDTRAJ = 'mdtraj'

# Appended to the input stem to name the output, e.g. prod.xtc -> prod-whole.xtc.
OUTPUT_TAG = '-whole'

# Longest plausible chemical bond, nm. 0.25 nm == the 2.5 Angstrom default of
# LOOS's long-bond-finder, which is the tool this check reimplements.
MAX_BOND = 0.25

# Frames read per chunk when scanning a trajectory, so a multi-GB .xtc never
# has to be resident.
SCAN_CHUNK = 200

# LOOS reads and writes Angstroms; GROMACS .gro/.xtc are nm.
ANGSTROM_PER_NM = 10.0

# Coordinate precision an .xtc is written at, in reciprocal nm. GROMACS and
# LOOS both default to this; a run that asked for finer says so in its frames.
XTC_PRECISION = 1000.0

# Placeholders for LOOS's own clock, which every frame here overrides.
LOOS_WRITER_DT = 1.0
LOOS_WRITER_STEPS_PER_FRAME = 1


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
    try:
        return app.GromacsTopFile(str(top_p), **kwargs).topology
    except Exception as exc:
        if include_dir is not None:
            raise
        # openmm guesses an include directory from $GMXDATA, then $GMXBIN.
        # $GMXBIN holding a binary NAME rather than the directory that binary
        # lives in, which is how a module-installed GROMACS is often set up,
        # sends that guess somewhere that does not exist.
        raise ValueError(
            f'{top_p.name}: {type(exc).__name__}: {exc}. No include_dir was '
            f'given, so openmm guessed one from GMXDATA='
            f'{os.environ.get("GMXDATA")!r} and GMXBIN='
            f'{os.environ.get("GMXBIN")!r}. GMXBIN must be the directory the '
            f'gmx binary lives in, not its name. Pass include_dir naming the '
            f"force field's share/gromacs/top instead.") from exc


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
    bigger than half the box, but then the mdtraj backend has to do the
    reimaging, since it walks the bonds instead.

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
                local_frame, worst_atom, worst_axis = (
                    int(i) for i in np.unravel_index(flat, usage.shape))
                worst_frame = frame_index + local_frame
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
            f'Use the {BACKEND_MDTRAJ!r} backend for this cell.')

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

    # A run that asked for a finer compressed-x-precision must not be
    # quantised back to the default on its way through LOOS.
    writer = _loos_writer(out_p, precision=source_precision(traj_p))
    timing = util.frame_timing(traj_p)
    axis = WrittenAxis()
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
        axis.took(index)
        n_written += 1
    del writer
    stamp_written_axis(out_p, timing, axis)

    if verify:
        _verify_reimaged(out_p, top_fn=top_fn, include_dir=include_dir,
                         ranges=ranges,
                         structure_fn=structure_fn, max_bond=max_bond)
    return out_p, n_written


def _verify_reimaged(out_p, top_fn=None, include_dir=None, ranges=None,
                     structure_fn=None, max_bond=MAX_BOND):
    """Check a just-reimaged trajectory, since imaging can fail quietly.

    A bond longer than max_bond says a molecule is still split, whichever
    backend wrote the file. Given ranges it also measures the anchor margin,
    which says whether LOOS's mergeImage() was entitled to an answer and is the
    only check that catches an unbonded atom stranded from its own molecule.
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
        f'{BACKEND_MDTRAJ!r} backend.')

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


def source_precision(traj_p, default=XTC_PRECISION):
    """The .xtc precision a trajectory was written at, or the default."""
    from . import gmx_simulate
    found = gmx_simulate.xtc_precision(traj_p)
    return default if found is None else found


class _DcdWriter:
    """loos.DCDWriter, taking the same three arguments an XTCWriter takes.

    A DCD keeps one timing rule in its header rather than a stamp per frame, so
    the step and time handed in here are recorded for stamp_written_axis to
    write at close, not passed to LOOS -- whose writeFrame takes a group and
    nothing else, and whose header hardcodes istart and nsavc to 1.
    """

    def __init__(self, out_p):
        import loos
        self.inner = loos.DCDWriter(str(out_p))

    def writeFrame(self, group, step=None, time=None):
        self.inner.writeFrame(group)


class WrittenAxis:
    """Where the frames a writer kept sit on the axis of the source it read.

    A stream that keeps every Nth frame has its own spacing, not its source's,
    and a DCD can stated only one uniform rule -- so what was kept is tracked
    as it is written rather than assumed from the policy that chose it.
    """

    def __init__(self):
        self.first = self.stride = self.last = None
        self.count = 0

    def took(self, index):
        if self.count == 0:
            self.first = index
        elif self.count == 1:
            self.stride = index - self.first
        self.last = index
        self.count += 1

    def uniform(self):
        """Whether one linear rule reaches every frame that was kept."""
        if self.count < 2:
            return True
        return self.last == self.first + (self.count - 1) * self.stride


def stamp_written_axis(out_p, timing, axis):
    """Carry a source's axis onto a DCD output, in that output's own spacing.

    Does nothing for the formats that stamp their frames as they go, or when
    the source had no axis to carry. True when the header was written.
    """
    out_p = Path(out_p)
    if timing is None or out_p.suffix.lower() != '.dcd' or axis.count < 2:
        return False
    if not axis.uniform():
        raise ValueError(
            f'{out_p} keeps frames {axis.first} to {axis.last} of its source '
            f'unevenly, and a DCD states one rule for the whole file. Refusing '
            'to stamp an axis that skips frames it claims to cover.')
    step0, steps_per_frame, time0, time_per_frame = timing
    return util.stamp_dcd_timing(
        out_p,
        step0 + axis.first * steps_per_frame, axis.stride * steps_per_frame,
        time0 + axis.first * time_per_frame, axis.stride * time_per_frame)


def _loos_writer(out_p, precision=XTC_PRECISION):
    suffix = out_p.suffix.lower()
    if suffix == '.xtc':
        return loos_xtc_writer(out_p, precision=precision)
    if suffix == '.dcd':
        return _DcdWriter(out_p)
    raise ValueError(
        f'{out_p.suffix} is not a format the LOOS backend writes; use .xtc or '
        '.dcd')


def loos_xtc_writer(out_p, precision=XTC_PRECISION,
                    steps_per_frame=LOOS_WRITER_STEPS_PER_FRAME,
                    dt=LOOS_WRITER_DT):
    """An .xtc writer holding this precision. Frames carry their own step and
    time, so dt and steps_per_frame only feed the unused 1-argument form."""
    import loos
    return loos.XTCWriter(str(out_p), dt, steps_per_frame, float(precision))


class _MdtrajWriter:
    """Writer that appends, which Trajectory.save() cannot do."""

    def __init__(self, out_p, angstrom_per_nm=ANGSTROM_PER_NM):
        import mdtraj as md
        self.suffix = out_p.suffix.lower()
        self.angstrom_per_nm = angstrom_per_nm
        if self.suffix == '.xtc':
            self.fh = md.formats.XTCTrajectoryFile(str(out_p), 'w')
        elif self.suffix == '.dcd':
            self.fh = md.formats.DCDTrajectoryFile(str(out_p), 'w')
        else:
            raise ValueError(
                f'{out_p.suffix} is not a format the mdtraj backend writes; '
                'use .xtc or .dcd')

    # The harvest writes through this class too, so a format either backend
    # cannot write fails the same way whichever one the box shape picked.

    def write(self, traj, step=None):
        if self.suffix == '.xtc':
            # Without an explicit step mdtraj writes the frame index, which
            # silently replaces the MD step counter with a small integer.
            self.fh.write(traj.xyz, time=traj.time, step=step,
                          box=traj.unitcell_vectors)
        else:
            self.fh.write(traj.xyz * self.angstrom_per_nm,
                          cell_lengths=(None if traj.unitcell_lengths is None
                                        else traj.unitcell_lengths
                                        * self.angstrom_per_nm),
                          cell_angles=traj.unitcell_angles)

    def close(self):
        self.fh.close()


def mdtraj_topology(top_fn, include_dir=None):
    """The mdtraj Topology for a GROMACS .top, carrying its molecules and bonds."""
    import mdtraj as md
    return md.Topology.from_openmm(
        gromacs_topology(top_fn, include_dir=include_dir))


def molecule_atoms(topology, ranges):
    """One list of mdtraj Atom objects per molecule block in ranges."""
    atoms = list(topology.atoms)
    if len(atoms) != ranges[-1][1]:
        raise ValueError(
            f'the topology has {len(atoms)} atoms but the molecule ranges '
            f'cover {ranges[-1][1]}; they do not describe the same system.')
    return [atoms[start:stop] for start, stop in ranges]


def largest_molecule(ranges):
    """Index of the molecule block with the most atoms, the solute if there is one."""
    return int(np.argmax([stop - start for start, stop in ranges]))


def check_bonds_within_molecules(pairs, ranges):
    """Raise unless every bond joins two atoms of the same molecule block.

    The analogue of the LOOS backend's splitByMolecule() count check: a bond
    crossing a block means the blocks being imaged are not the molecules, and
    wrapping them per block would tear a real molecule apart.
    """
    counts = [stop - start for start, stop in ranges]
    molecule_of = np.repeat(np.arange(len(ranges)), counts)
    if pairs.max() >= len(molecule_of):
        raise ValueError(
            f'a bond names atom {pairs.max()} but the molecule ranges cover '
            f'only {len(molecule_of)} atoms.')
    crossing = molecule_of[pairs[:, 0]] != molecule_of[pairs[:, 1]]
    if crossing.any():
        i, j = pairs[int(np.argmax(crossing))]
        raise ValueError(
            f'bond {i}-{j} joins two different molecule blocks '
            f'({int(molecule_of[i])} and {int(molecule_of[j])}) of the '
            f'{len(ranges)} the topology declares; refusing to reimage against '
            'a grouping that does not match the bonds.')


def _step_numbers(timing, index):
    """MD step of each frame at index, or None for a format that keeps its own."""
    return None if timing is None else timing[0] + index * timing[1]


def reimage_with_mdtraj(traj_fn, top_fn, out_fn, ranges=None, pairs=None,
                        include_dir=None, anchor_index=None,
                        skip_first_frame=False, verify=True,
                        max_bond=MAX_BOND, output_tag=OUTPUT_TAG,
                        scan_chunk=SCAN_CHUNK):
    """Make molecules whole and wrap them per molecule, with mdtraj.

    Any cell, triclinic included, which is the reason this backend exists.
    Molecules and bonds are read from top_fn and handed to image_molecules
    explicitly: left to itself mdtraj calls find_molecules(), which strands
    every unbonded atom in a molecule of its own and images the system atom by
    atom, and no later pass can put those molecules back together.

    anchor_index names the molecule centred in the box, the largest by default.
    One anchor is what keeps the wrap tight and is far cheaper than anchoring
    on every molecule.

    Returns (out_p, n_written).
    """
    import mdtraj as md
    traj_p, out_p = Path(traj_fn), Path(out_fn)
    if out_p.resolve() == traj_p.resolve():
        raise ValueError(
            f'refusing to reimage {traj_p} onto itself; the orchestrator counts '
            f'frames in the raw trajectory to decide whether a generation is '
            f'complete. Write to a new file (e.g. {traj_p.stem}{output_tag}'
            f'{traj_p.suffix}).')
    if top_fn is None:
        raise ValueError(
            'reimage_with_mdtraj needs top_fn: the GROMACS topology is the '
            'only place molecule blocks and bonds are recorded, and mdtraj '
            'would otherwise image this trajectory atom by atom.')

    topology = mdtraj_topology(top_fn, include_dir=include_dir)
    if ranges is None:
        ranges = molecule_ranges(top_fn, include_dir=include_dir)
    if pairs is None:
        pairs = bond_pairs(top_fn, include_dir=include_dir)
    check_bonds_within_molecules(pairs, ranges)

    molecules = molecule_atoms(topology, ranges)
    if anchor_index is None:
        anchor_index = largest_molecule(ranges)
    anchor = molecules[anchor_index]
    others = molecules[:anchor_index] + molecules[anchor_index + 1:]

    with md.open(str(traj_p)) as fh:
        timing = util.frame_timing(traj_p, n_frames=len(fh))

    writer = _MdtrajWriter(out_p)
    n_read = n_written = 0
    try:
        for chunk in md.iterload(str(traj_p), top=topology, chunk=scan_chunk):
            index = np.arange(n_read, n_read + chunk.n_frames)
            keep = (index > 0) if skip_first_frame else np.ones(len(index), bool)
            chunk.image_molecules(inplace=True, anchor_molecules=[anchor],
                                  other_molecules=others, make_whole=True)
            if keep.any():
                writer.write(chunk[keep],
                             step=_step_numbers(timing, index[keep]))
                n_written += int(keep.sum())
            n_read += chunk.n_frames
    finally:
        writer.close()

    if verify:
        # No ranges: the anchor margin measures the half-edge assumption in
        # LOOS's mergeImage(), which walking the bonds does not make.
        _verify_reimaged(out_p, top_fn=top_fn, include_dir=include_dir,
                         max_bond=max_bond)
    return out_p, n_written


def reimage_trajectory(traj_fn, out_fn=None, structure_fn=None, top_fn=None,
                       backend=BACKEND_AUTO, include_dir=None,
                       center_selection=None, anchor_index=None,
                       skip_first_frame=False, triclinic_rtol=TRICLINIC_RTOL,
                       output_tag=OUTPUT_TAG):
    """Reimage a trajectory, choosing the backend from its own box.

    'auto' sends a rectangular box to LOOS and anything else to mdtraj, which
    is the only backend here that can represent a triclinic cell. Force
    'mdtraj' when there is no structure file, or when a molecule is larger than
    half a box edge. Returns the output path, and never touches the input.
    """
    traj_p = Path(traj_fn)
    if out_fn is None:
        out_fn = traj_p.with_name(f'{traj_p.stem}{output_tag}{traj_p.suffix}')

    if backend == BACKEND_AUTO:
        box = box_vectors(traj_fn=traj_p, structure_fn=structure_fn)
        orthorhombic = is_orthorhombic(box, triclinic_rtol=triclinic_rtol)
        backend = BACKEND_LOOS if orthorhombic else BACKEND_MDTRAJ
        print(f'[reimage] box is '
              f'{"orthorhombic" if orthorhombic else "triclinic"}; '
              f'using the {backend} backend', flush=True)
        if backend == BACKEND_LOOS and structure_fn is None:
            print('[reimage] no structure_fn given for the LOOS backend; '
                  'falling back to mdtraj', flush=True)
            backend = BACKEND_MDTRAJ

    if backend == BACKEND_LOOS:
        out_p, _ = reimage_with_loos(
            traj_p, structure_fn, out_fn, top_fn=top_fn,
            include_dir=include_dir, center_selection=center_selection,
            skip_first_frame=skip_first_frame,
            triclinic_rtol=triclinic_rtol, output_tag=output_tag)
        return out_p
    if backend == BACKEND_MDTRAJ:
        out_p, _ = reimage_with_mdtraj(
            traj_p, top_fn, out_fn, include_dir=include_dir,
            anchor_index=anchor_index, skip_first_frame=skip_first_frame,
            output_tag=output_tag)
        return out_p
    raise ValueError(
        f'unknown backend {backend!r}; choose from '
        f'{[BACKEND_AUTO, BACKEND_LOOS, BACKEND_MDTRAJ]}')


def reimage_gen_dir(gen_dir, config=None, backend=BACKEND_AUTO,
                    include_dir=None, center_selection=None, anchor_index=None,
                    skip_first_frame=None, output_tag=OUTPUT_TAG,
                    triclinic_rtol=TRICLINIC_RTOL):
    """Reimage one generation directory, writing beside the raw trajectory.

    Reads config.json for the trajectory name, topology, structure and, when
    the caller names none, the include_dir its #include lines need. The raw
    trajectory is left alone so a generation still being resumed keeps its frame
    count. skip_first_frame defaults to dropping it for every generation after
    the first, since GROMACS repeats the previous generation's last frame.
    """
    import json
    gen_p = Path(gen_dir)
    if config is None:
        config = json.loads((gen_p / 'config.json').read_text())
    if include_dir is None:
        include_dir = config.get('include_dir')

    traj_p = (gen_p / config['traj_name']).with_suffix(config['traj_suffix'])
    if not traj_p.is_file():
        raise FileNotFoundError(f'{traj_p}: nothing to reimage')
    if skip_first_frame is None:
        skip_first_frame = config.get('gen_index', 0) > 0

    return reimage_trajectory(
        traj_p,
        structure_fn=config.get('structure_fn'),
        top_fn=config.get('top_fn'),
        backend=backend, include_dir=include_dir,
        center_selection=center_selection, anchor_index=anchor_index,
        skip_first_frame=skip_first_frame,
        triclinic_rtol=triclinic_rtol, output_tag=output_tag)


