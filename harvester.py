"""Harvesting: reduce a finished generation to the streams we keep, then --
and only then -- delete the original.

A harvest produces two outputs from one generation's raw trajectory:

  * the **dry** stream, every frame, solute only;
  * the **downsampled** stream, every Nth frame, solvent kept.

and, if both check out, replaces the original with a symlink to the dry stream.
That last step is irreversible, so everything here is built around making it
provably safe:

**Counts, not file sizes.** The original guard was ``st_size > 0`` on both
outputs. A harvest killed partway through the write loop -- preemption,
walltime, OOM -- leaves both outputs nonzero and truncated, and the guard passes.
Frames are counted in the source first, the expected output counts are computed
from a single frame plan, and the written files are re-counted off disk before
anything is unlinked.

**Idempotence.** After a successful harvest the original name is a symlink to
the dry stream. Re-running the old code opened that symlink as input *and* the
dry file as output -- same file, read and written at once, dry copy destroyed,
original already gone. A requeueing scheduler makes re-runs routine, so a
``.harvested`` sentinel carrying the counts is written last and short-circuits
the whole thing; a symlinked input with no sentinel is recognised as a harvest
that died in that window and is repaired rather than repeated.

**Streaming.** Generations here run to ~1 microsecond. ``md.load`` is
whole-trajectory: 240k frames of 10k atoms is 29 GB of coordinates before
imaging or slicing. The LOOS backend streams frame by frame; the mdtraj backend
uses ``md.iterload``. Neither holds the trajectory.

**Backend by box shape, decided by looking.** LOOS cannot represent a triclinic
cell -- and does not raise, it silently keeps the diagonal (see ``reimage`` for
the gory details) -- so the choice cannot be a try/except around the LOOS call;
the exception would never fire on the case it exists to catch. The box is read
off frame 0 and the off-diagonals decide: rectangular goes to LOOS, anything
else to mdtraj.

**One frame plan.** Which frames land in each stream is computed once, in
``frame_plan``, from the generation's *global* index -- so the downsample phase
is continuous across generations instead of resetting at every seam, and the
frame each engine writes at the restart step is dropped exactly once. Both
backends consume that plan, so they cannot disagree.
"""

import json
import subprocess as sp
from pathlib import Path

from . import utilities as util
from . import reimage


# Written last, after both outputs are verified and the original is gone. Its
# presence is what makes a harvest idempotent under a requeueing scheduler.
SENTINEL_NAME = '.harvested'

# Output name prefixes, joined to the trajectory name with the config's `sep`.
DRY_PREFIX = 'dry'
DOWNSAMPLE_PREFIX = 'downsample'

# Solute-only topology written beside the dry stream, so everything downstream
# -- including get_traj_len's LOOS fallback, which would otherwise build the
# model from the wet topology and choke on the atom count -- can read it.
DRY_TOPOLOGY_NAME = 'dry-top.pdb'

# Backend selectors for `harvest_generation`.
BACKEND_AUTO = 'auto'
BACKEND_LOOS = 'loos'
BACKEND_MDTRAJ = 'mdtraj'

# Frames per md.iterload chunk. Only the mdtraj backend uses this; the LOOS
# backend is already frame-at-a-time.
ITERLOAD_CHUNK = 100

# What to do with the frame an engine writes at the step it restarted from,
# which duplicates the previous generation's last frame. 'auto' decides from the
# frame count: a generation holding steps_per_gen/write_interval + 1 frames has
# one, a generation holding exactly steps_per_gen/write_interval does not.
# GROMACS writes it; the OpenMM reporters do not.
SEAM_AUTO = 'auto'
SEAM_DROP = 'drop'
SEAM_KEEP = 'keep'

# Selection-language dialects understood for `harvester_subset`.
SYNTAX_LOOS = 'loos'
SYNTAX_MDTRAJ = 'mdtraj'

# LOOS reads and writes Angstroms; GROMACS .gro/.xtc are nm.
ANGSTROM_PER_NM = reimage.ANGSTROM_PER_NM


class HarvestError(RuntimeError):
    """Raised when a harvest cannot be completed safely.

    Always raised *before* the original trajectory is removed, so a generation
    that fails to harvest still holds everything it started with.
    """


class Harvester:
    __slots__ = ('run_config', 'main_path', 'run_config_name',
                 'harvester_template', 'scheduler', 'scriptname')

    def __init__(self,  harvester_template: str,
                 scheduler: str, run_config=None, scriptname='harvest.sh',
                 run_config_name='hconfig.json'):
        # Expects a dict to be fed to harvester_template.format(**run_config)
        self.run_config = run_config
        self.run_config_name = run_config_name
        # Template shell script, optionally with str.format slots for run_dict.
        self.harvester_template = harvester_template
        self.scheduler = scheduler
        self.scriptname = scriptname

    def prep_and_write_inputs(self, current_dir: Path):
        if self.run_config:
            harvest_script = self.harvester_template.format(**self.run_config)
            harvest_config_p = current_dir/self.run_config_name
            harvest_config_p.write_text(json.dumps(self.run_config))
        else:
            harvest_script = self.harvester_template
        harvest_script_p = current_dir/self.scriptname
        harvest_script_p.write_text(harvest_script)
        return harvest_script_p

    def reap(self, current_dir, dry_run=False):
        harvest_script_p = self.prep_and_write_inputs(current_dir)
        if not dry_run:
            try:
                with harvest_script_p.open() as f:
                    scheduler_output = sp.check_output(
                        self.scheduler, stdin=f, cwd=current_dir, text=True
                    )
                    print('harvester scheduler return:', scheduler_output)
            except sp.CalledProcessError as err:
                print(f'{self.scheduler} call threw error', err.stdout, err.stderr)
                raise
        else:
            scheduler_output = None
        return scheduler_output


# ---------------------------------------------------------------------------
# Config-time checks
# ---------------------------------------------------------------------------

def check_commensurability(steps_per_gen, write_interval, downsample_frq):
    """Both spacing conditions, checked where a violation is still free to fix.

    ``steps_per_gen % write_interval == 0`` puts a generation's last written
    frame exactly on the checkpoint the next generation restarts from. Miss it
    and every seam silently drops the sliver of trajectory between the last
    frame and the restart state -- invisible in the timestamps, real in the
    sampling.

    ``frames_per_gen % downsample_frq == 0`` keeps the downsampled stream evenly
    spaced across the concatenation. Generations are harvested independently, so
    an incommensurate frame count breaks the spacing at every boundary.

    Returns frames_per_gen, which is the number of *new* frames a generation
    contributes -- one less than the file holds when the engine writes a frame
    at its restart step.
    """
    if write_interval <= 0:
        raise ValueError(f'write_interval must be positive, got {write_interval}')
    if steps_per_gen % write_interval:
        raise ValueError(
            f'steps_per_gen={steps_per_gen} is not a whole number of write '
            f'intervals ({write_interval}); the last frame of each generation '
            f'would not land on the checkpoint the next one restarts from, so '
            f'every seam would quietly drop '
            f'{steps_per_gen % write_interval} steps of trajectory.')
    frames_per_gen = steps_per_gen // write_interval
    if downsample_frq and frames_per_gen % downsample_frq:
        raise ValueError(
            f'{frames_per_gen} frames per generation is not a whole number of '
            f'downsample periods ({downsample_frq}); the downsampled stream '
            f'would change spacing at every generation boundary.')
    return frames_per_gen


# ---------------------------------------------------------------------------
# The frame plan -- the single place that decides what goes where
# ---------------------------------------------------------------------------

def resolve_seam(n_orig, frames_per_gen, gen_index, seam=SEAM_AUTO):
    """True when this generation's frame 0 duplicates the previous one's last.

    Generation N+1 restarts from generation N's checkpoint; GROMACS writes a
    frame at that step, so the two files share a time point. Concatenating
    without dropping one puts a repeated frame at every seam, which never fails
    loudly -- it just biases lag times and kinetics, which is the entire point of
    the dataset.
    """
    if seam == SEAM_KEEP:
        return False
    if seam == SEAM_DROP:
        return gen_index > 0
    if seam != SEAM_AUTO:
        raise ValueError(
            f'seam={seam!r} is not one of {SEAM_AUTO!r}, {SEAM_DROP!r}, '
            f'{SEAM_KEEP!r}')
    if n_orig == frames_per_gen + 1:
        return gen_index > 0
    if n_orig == frames_per_gen:
        return False
    raise HarvestError(
        f'generation {gen_index} holds {n_orig} frames, but a generation of '
        f'{frames_per_gen} write intervals should hold either {frames_per_gen} '
        f'(no frame written at the restart step) or {frames_per_gen + 1} (one '
        f'written). Refusing to guess which frames are duplicates -- pass '
        f'seam={SEAM_DROP!r} or {SEAM_KEEP!r} explicitly if this count is '
        f'expected.')


def keeps_frame(local_index, first_global_index, downsample_frq, skip_first):
    """``(write_dry, write_downsample)`` for one frame. The rule, in one place.

    `first_global_index` is where this generation's frames sit in the
    concatenated stream, which is what keeps the downsample phase continuous
    across generations rather than resetting at every seam.
    """
    if skip_first and local_index == 0:
        return False, False
    return True, ((first_global_index + local_index) % downsample_frq == 0)


def frame_plan(n_orig, first_global_index, downsample_frq, skip_first):
    """Yield ``(local_index, write_dry, write_downsample)`` for every frame."""
    for local in range(n_orig):
        dry, down = keeps_frame(local, first_global_index, downsample_frq,
                                skip_first)
        yield local, dry, down


def expected_counts(n_orig, first_global_index, downsample_frq, skip_first):
    """(n_dry, n_downsample) the plan will produce. Computed, never observed."""
    n_dry = n_down = 0
    for _, dry, down in frame_plan(n_orig, first_global_index, downsample_frq,
                                   skip_first):
        n_dry += dry
        n_down += down
    return n_dry, n_down


# ---------------------------------------------------------------------------
# Backend selection and subset resolution
# ---------------------------------------------------------------------------

def source_frame_timing(traj_fn, n_orig):
    """``(step0, steps_per_frame, time0, time_per_frame)`` read off the source.

    Neither writer preserves this on its own, and the harvest deletes the
    original, so getting it wrong destroys the time axis permanently rather
    than merely inconveniently:

      * LOOS's XTCWriter numbers frames from its own counters -- ``dt_ = 1.0``,
        ``step_ = 0``, ``steps_per_frame_ = 1`` (``src/xtcwriter.hpp``) -- so
        every harvested frame came out 1 ps apart and stamped with its frame
        index instead of its MD step, whatever the trajectory actually held.
      * mdtraj carries ``time`` but fills ``step`` with the frame index unless
        told otherwise.

    Both are fixed by reading the source's frame 0 and frame 1 and writing every
    frame with an explicit step and time. The extrapolation is checked against
    the *last* frame before it is used, because it is only valid if the source
    is evenly spaced -- which a generation always is, but which is exactly the
    assumption you want to hear about rather than trust.

    Returns None for a format that carries no per-frame timing (DCD keeps it in
    the header, and mdtraj does not surface it on read).
    """
    import numpy as np
    import mdtraj as md

    traj_p = Path(traj_fn)
    if traj_p.suffix.lower() != '.xtc':
        return None
    with md.open(str(traj_p)) as fh:
        _, time, step, _ = fh.read(min(2, n_orig))
        if n_orig < 2:
            return int(step[0]), 0, float(time[0]), 0.0
        step0, time0 = int(step[0]), float(time[0])
        steps_per_frame = int(step[1]) - step0
        time_per_frame = float(time[1]) - time0
        fh.seek(n_orig - 1)
        _, last_time, last_step, _ = fh.read(1)
    predicted_step = step0 + (n_orig - 1) * steps_per_frame
    predicted_time = time0 + (n_orig - 1) * time_per_frame
    if int(last_step[0]) != predicted_step or not np.isclose(
            float(last_time[0]), predicted_time, rtol=1e-5,
            atol=1e-5 * max(abs(predicted_time), 1.0)):
        raise HarvestError(
            f'{traj_p} is not evenly spaced: frames 0 and 1 are '
            f'{steps_per_frame} steps / {time_per_frame} ps apart, which puts '
            f'frame {n_orig - 1} at step {predicted_step} / {predicted_time} ps, '
            f'but it is at step {int(last_step[0])} / {float(last_time[0])} ps. '
            'Refusing to restamp frames from an assumption the trajectory '
            'contradicts.')
    return step0, steps_per_frame, time0, time_per_frame


def select_backend(traj_fn, structure_fn=None, backend=BACKEND_AUTO,
                   triclinic_rtol=reimage.TRICLINIC_RTOL):
    """Pick a backend by looking at the box, not by catching an exception.

    LOOS does not raise on a triclinic cell; it keeps the diagonal and carries
    on. A try/except around the LOOS call would therefore never fire on the one
    case it was written for, so the off-diagonals are read first and the choice
    is made before either engine is touched.
    """
    if backend != BACKEND_AUTO:
        if backend not in (BACKEND_LOOS, BACKEND_MDTRAJ):
            raise ValueError(
                f'backend={backend!r} is not one of {BACKEND_AUTO!r}, '
                f'{BACKEND_LOOS!r}, {BACKEND_MDTRAJ!r}')
        return backend
    box = reimage.box_vectors(traj_fn=traj_fn, structure_fn=structure_fn)
    if reimage.is_orthorhombic(box, triclinic_rtol=triclinic_rtol):
        return BACKEND_LOOS
    print(f'[harvest] {traj_fn} has a non-orthorhombic box; LOOS would keep '
          f'only its diagonal, so harvesting with mdtraj instead.', flush=True)
    return BACKEND_MDTRAJ


def subset_indices(structure_fn, selection, syntax=SYNTAX_LOOS):
    """0-based atom indices matched by `selection`, in file order."""
    if syntax == SYNTAX_LOOS:
        import loos
        model = loos.createSystem(str(structure_fn))
        group = loos.selectAtoms(model, selection)
        return [atom.index() for atom in group]
    if syntax == SYNTAX_MDTRAJ:
        import mdtraj as md
        return [int(i) for i in md.load(str(structure_fn)).top.select(selection)]
    raise ValueError(
        f'syntax={syntax!r} is not one of {SYNTAX_LOOS!r}, {SYNTAX_MDTRAJ!r}')


def indices_to_loos_selection(indices):
    """A LOOS selection string matching exactly `indices`.

    ``index`` is a numeric selector in the LOOS grammar (``src/grammar.yy``,
    ``pushAtomIndex``) reading ``Atom::index()`` -- the 0-based position in the
    model, which is what both engines agree on. Contiguous runs are collapsed so
    a solute subset stays a handful of clauses rather than thousands.

    Going through a selection string rather than assembling a group atom by atom
    is not stylistic: ``AtomicGroup::select`` shares the parent's
    SharedPeriodicBox, and a hand-built group does not. A dry trajectory written
    from a group with its own box gets a frozen (or absent) cell, which is
    exactly the sort of thing found two years later.
    """
    runs, start, previous = [], None, None
    for index in sorted(indices):
        if start is None:
            start = previous = index
        elif index == previous + 1:
            previous = index
        else:
            runs.append((start, previous))
            start = previous = index
    if start is not None:
        runs.append((start, previous))
    return ' || '.join(
        f'(index == {lo})' if lo == hi
        else f'(index >= {lo} && index <= {hi})' for lo, hi in runs)


def resolve_subset(structure_fn, selection, syntax=SYNTAX_LOOS):
    """Resolve a subset once, into a form each backend can use.

    The two engines speak different selection languages, and the backend is
    chosen from the box rather than by the user -- so a config carrying only a
    LOOS string must still work when the box turns out to be triclinic, and vice
    versa. Resolving once and handing each backend its own view of the *same*
    atoms is what keeps the choice of backend from changing what gets written.

    Returns ``None`` for "keep everything", else a dict with ``indices`` (for
    mdtraj) and ``loos_selection`` (for LOOS).
    """
    if not selection:
        return None
    indices = subset_indices(structure_fn, selection, syntax=syntax)
    if not indices:
        raise HarvestError(
            f'harvester_subset {selection!r} ({syntax} syntax) matched no atoms '
            f'in {structure_fn}; refusing to write an empty dry trajectory.')
    return dict(
        indices=indices,
        loos_selection=(selection if syntax == SYNTAX_LOOS
                        else indices_to_loos_selection(indices)))


# ---------------------------------------------------------------------------
# Backends
# ---------------------------------------------------------------------------

def _harvest_loos(traj_fn, structure_fn, subset_spec, dry_out, down_out,
                  first_global_index, downsample_frq, skip_first,
                  timing=None, dry_topology_name=DRY_TOPOLOGY_NAME):
    """Stream the trajectory once, writing both outputs. Orthorhombic only.

    The subset comes from ``AtomicGroup::select``, which shares -- rather than
    copies -- the parent's SharedPeriodicBox (``src/AtomicGroup.cpp:349``). That
    is what makes the dry stream carry a live per-frame box instead of a frozen
    one, and it is why the subset is taken with a selection string rather than
    assembled atom by atom.
    """
    import loos
    from loos import pyloos

    model = loos.createSystem(str(structure_fn))
    if subset_spec is None:
        subset = model
    else:
        subset = loos.selectAtoms(model, subset_spec['loos_selection'])
        if len(subset) != len(subset_spec['indices']):
            raise HarvestError(
                f'the LOOS selection matched {len(subset)} atoms but the subset '
                f'resolves to {len(subset_spec["indices"])}; the two backends '
                'would not write the same atoms.')

    dry_writer = reimage._loos_writer(Path(dry_out))
    down_writer = reimage._loos_writer(Path(down_out))
    traj = pyloos.Trajectory(str(traj_fn), model)
    n_orig = n_dry = n_down = 0
    for local, _ in enumerate(traj):
        n_orig += 1
        dry, down = keeps_frame(local, first_global_index, downsample_frq,
                                skip_first)
        if not (dry or down):
            continue
        if timing is None:
            # No per-frame timing to carry (DCD); let the writer count.
            if dry:
                dry_writer.writeFrame(subset)
            if down:
                down_writer.writeFrame(model)
        else:
            step0, steps_per_frame, time0, time_per_frame = timing
            step = step0 + local * steps_per_frame
            time = time0 + local * time_per_frame
            if dry:
                dry_writer.writeFrame(subset, step, time)
            if down:
                down_writer.writeFrame(model, step, time)
        n_dry += dry
        n_down += down
    del dry_writer, down_writer

    # Solute-only topology for anything that reads the dry stream later.
    subset.pruneBonds()
    Path(dry_out).parent.joinpath(dry_topology_name).write_text(
        str(loos.PDB.fromAtomicGroup(subset)))
    return n_orig, n_dry, n_down


def _harvest_mdtraj(traj_fn, structure_fn, subset_spec, dry_out, down_out,
                    first_global_index, downsample_frq, skip_first,
                    timing=None, dry_topology_name=DRY_TOPOLOGY_NAME,
                    iterload_chunk=ITERLOAD_CHUNK,
                    angstrom_per_nm=ANGSTROM_PER_NM):
    """Same plan, chunked. Handles the cells LOOS cannot represent.

    ``md.iterload`` is what keeps this off the heap, but it introduces its own
    trap: the downsample must be driven by a running global frame index, not by
    slicing each chunk ``[::N]``. Per-chunk slicing resets the phase at every
    chunk boundary -- the same bug as a per-generation phase reset, just finer
    grained and harder to see.
    """
    import numpy as np
    import mdtraj as md

    model = md.load(str(structure_fn))
    indices = None if subset_spec is None else subset_spec['indices']
    dry_writer = _MdtrajWriter(Path(dry_out), angstrom_per_nm=angstrom_per_nm)
    down_writer = _MdtrajWriter(Path(down_out), angstrom_per_nm=angstrom_per_nm)
    n_orig = n_dry = n_down = 0
    try:
        for chunk in md.iterload(str(traj_fn), top=model.top,
                                 chunk=iterload_chunk):
            local = np.arange(n_orig, n_orig + chunk.n_frames)
            global_index = first_global_index + local
            keep_dry = np.ones(chunk.n_frames, dtype=bool)
            keep_down = (global_index % downsample_frq) == 0
            if skip_first:
                keep_dry &= local > 0
                keep_down &= local > 0
            step = None
            if timing is not None:
                step0, steps_per_frame, time0, time_per_frame = timing
                step = step0 + local * steps_per_frame
                predicted = time0 + local * time_per_frame
                if not np.allclose(predicted, chunk.time, rtol=1e-5,
                                   atol=1e-5 * max(abs(time0), 1.0)):
                    raise HarvestError(
                        f'{traj_fn}: frame times drift from the even spacing '
                        f'read off frames 0 and 1 (worst difference '
                        f'{np.abs(predicted - chunk.time).max():g} ps in the '
                        f'chunk starting at frame {n_orig}).')
            if keep_dry.any():
                sub = chunk[keep_dry]
                dry_writer.write(sub if indices is None
                                 else sub.atom_slice(indices),
                                 step=None if step is None else step[keep_dry])
                n_dry += int(keep_dry.sum())
            if keep_down.any():
                down_writer.write(
                    chunk[keep_down],
                    step=None if step is None else step[keep_down])
                n_down += int(keep_down.sum())
            n_orig += chunk.n_frames
    finally:
        dry_writer.close()
        down_writer.close()

    dry_model = model if indices is None else model.atom_slice(indices)
    dry_model.save_pdb(str(Path(dry_out).parent / dry_topology_name))
    return n_orig, n_dry, n_down


class _MdtrajWriter:
    """Append-as-you-go writer, because Trajectory.save() cannot append.

    XTC carries nm and a 3x3 box; DCD carries Angstroms and lengths/angles.
    Getting that wrong writes a trajectory scaled by ten, which looks fine until
    someone measures something.
    """

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
                f'{out_p.suffix} is not a format the mdtraj harvest backend '
                'writes; use .xtc or .dcd')

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


HARVEST_BACKENDS = {BACKEND_LOOS: _harvest_loos, BACKEND_MDTRAJ: _harvest_mdtraj}


# ---------------------------------------------------------------------------
# The harvest itself
# ---------------------------------------------------------------------------

def harvest_generation(config_fn, harvester_config_fn,
                       backend=BACKEND_AUTO, seam=SEAM_AUTO,
                       sentinel_name=SENTINEL_NAME,
                       dry_prefix=DRY_PREFIX,
                       downsample_prefix=DOWNSAMPLE_PREFIX,
                       dry_topology_name=DRY_TOPOLOGY_NAME,
                       iterload_chunk=ITERLOAD_CHUNK,
                       triclinic_rtol=reimage.TRICLINIC_RTOL,
                       backends=None):
    """Harvest one generation directory. Idempotent, and safe to interrupt.

    `harvester_config_fn` keys:

    ``harvester_subset``      selection string for the dry stream. Omitted or
                              empty keeps every atom.
    ``harvester_subset_syntax``  ``'loos'`` (default) or ``'mdtraj'``
    ``harvester_structure``   structure file to build the model from. Defaults
                              to ``config['top_fn']``, which is right for OpenMM
                              runs and wrong for GROMACS ones -- a ``.top`` is a
                              force-field topology and neither engine can build
                              a model from it, so GROMACS runs must point this
                              at a ``.gro`` / ``.pdb``.
    ``downsample_frq``        keep every Nth frame in the solvated stream
    ``harvester_unlink``      delete the original once verified (default True)
    """
    backends = HARVEST_BACKENDS if backends is None else backends
    config = json.loads(Path(config_fn).read_text())
    hconfig = json.loads(Path(harvester_config_fn).read_text())
    gen_dir = Path(config_fn).resolve().parent

    sep = config['sep']
    traj_fn = f"{config['traj_name']}{config['traj_suffix']}"
    traj_p = gen_dir / traj_fn
    dry_p = gen_dir / f'{dry_prefix}{sep}{traj_fn}'
    down_p = gen_dir / f'{downsample_prefix}{sep}{traj_fn}'
    sentinel_p = gen_dir / sentinel_name
    dry_top_p = gen_dir / dry_topology_name

    if sentinel_p.is_file():
        record = json.loads(sentinel_p.read_text())
        print(f'[harvest] {sentinel_p} exists; generation already harvested '
              f'({record.get("n_dry")} dry frames). Nothing to do.', flush=True)
        return dict(record, status='already-harvested')

    structure_fn = hconfig.get('harvester_structure') or config['top_fn']
    downsample_frq = hconfig['downsample_frq']
    gen_index = config['gen_index']
    write_interval = config['write_interval']
    steps_per_gen = _steps_per_gen(config, hconfig)
    frames_per_gen = check_commensurability(
        steps_per_gen, write_interval, downsample_frq)
    first_global_index = gen_index * frames_per_gen

    # A symlinked input means a previous harvest unlinked the original but died
    # before writing the sentinel. Re-running the write loop from here would
    # read the dry stream and write it at the same time; repair instead.
    if traj_p.is_symlink():
        return _repair_from_symlink(
            traj_p, dry_p, down_p, sentinel_p, dry_top_p, structure_fn,
            frames_per_gen=frames_per_gen, gen_index=gen_index,
            first_global_index=first_global_index,
            downsample_frq=downsample_frq, seam=seam)

    n_orig = util.get_traj_len(traj_p, structure_fn)
    if not n_orig:
        raise HarvestError(
            f'{traj_p} holds no readable frames; refusing to harvest a '
            'generation that has not run.')

    skip_first = resolve_seam(n_orig, frames_per_gen, gen_index, seam=seam)
    n_dry_expected, n_down_expected = expected_counts(
        n_orig, first_global_index, downsample_frq, skip_first)

    chosen = select_backend(traj_p, structure_fn=structure_fn, backend=backend,
                            triclinic_rtol=triclinic_rtol)
    subset_spec = resolve_subset(
        structure_fn, hconfig.get('harvester_subset'),
        syntax=hconfig.get('harvester_subset_syntax', SYNTAX_LOOS))

    print(f'[harvest] {gen_dir}: {n_orig} frames, backend={chosen}, '
          f'global frames {first_global_index}..'
          f'{first_global_index + n_orig - 1}, '
          f'{"dropping" if skip_first else "keeping"} the restart-step frame; '
          f'expecting {n_dry_expected} dry and {n_down_expected} downsampled.',
          flush=True)

    kwargs = dict(dry_topology_name=dry_topology_name,
                  timing=source_frame_timing(traj_p, n_orig))
    if chosen == BACKEND_MDTRAJ:
        kwargs['iterload_chunk'] = iterload_chunk
    n_seen, n_dry, n_down = backends[chosen](
        traj_p, structure_fn, subset_spec, dry_p, down_p,
        first_global_index, downsample_frq, skip_first, **kwargs)

    _verify_counts(traj_p, dry_p, down_p, dry_top_p, structure_fn,
                   n_orig=n_orig, n_seen=n_seen, n_dry=n_dry, n_down=n_down,
                   n_dry_expected=n_dry_expected,
                   n_down_expected=n_down_expected)

    record = dict(
        status='harvested', backend=chosen, n_orig=n_orig, n_dry=n_dry,
        n_down=n_down, first_global_index=first_global_index,
        frames_per_gen=frames_per_gen, gen_index=gen_index,
        downsample_frq=downsample_frq, skip_first=skip_first,
        write_interval=write_interval, steps_per_gen=steps_per_gen,
        structure_fn=str(structure_fn), dry=dry_p.name, downsample=down_p.name,
        original=traj_p.name,
        subset=hconfig.get('harvester_subset'),
        n_subset_atoms=(None if subset_spec is None
                        else len(subset_spec['indices'])))

    if hconfig.get('harvester_unlink', True):
        traj_p.unlink()
        # Leave a symlink so frame counting on the original name still works.
        traj_p.symlink_to(dry_p.name)
        record['unlinked'] = True
    else:
        record['unlinked'] = False

    # Written last: everything above is redoable, this is what says not to.
    sentinel_p.write_text(json.dumps(record, indent=2))
    return record


def _steps_per_gen(config, hconfig):
    """The full generation length, which is not always ``config['steps']``.

    On a resumed generation ``config['steps']`` has been narrowed to the steps
    still owed, so reading it would shorten frames_per_gen and put every
    subsequent generation's global frame index -- and therefore the downsample
    phase -- in the wrong place. `Clone.from_disk` records the untouched value
    as ``steps_per_gen``; hand-built configs may predate it.
    """
    for source, key in ((hconfig, 'steps_per_gen'), (config, 'steps_per_gen')):
        if source.get(key) is not None:
            return source[key]
    print(f'[harvest] WARNING: neither the harvester config nor the run config '
          f'records steps_per_gen, falling back to steps={config["steps"]}. '
          f'That is the full generation length only if this generation ran in '
          f'one launch; if it was resumed, the downsample phase for every later '
          f'generation will be wrong. Add steps_per_gen to the harvester '
          f'config.', flush=True)
    return config['steps']


def _verify_counts(traj_p, dry_p, down_p, dry_top_p, structure_fn, n_orig,
                   n_seen, n_dry, n_down, n_dry_expected, n_down_expected):
    """Check written frames against the plan, on disk, before anything is lost.

    Both the in-loop tallies and the files themselves are checked. The tallies
    catch a plan/backend disagreement; re-reading the files catches a writer
    that returned normally having flushed less than it was handed, which is what
    a job killed inside the write loop looks like.
    """
    problems = []
    if n_seen != n_orig:
        problems.append(
            f'read {n_seen} frames from {traj_p.name} but counted {n_orig} '
            'before starting; the trajectory changed under us')
    if n_dry != n_dry_expected:
        problems.append(f'wrote {n_dry} dry frames, planned {n_dry_expected}')
    if n_down != n_down_expected:
        problems.append(
            f'wrote {n_down} downsampled frames, planned {n_down_expected}')
    # Each stream gets the topology that actually describes it: the dry stream
    # is the subset, the downsampled stream is still the whole system. Only the
    # LOOS fallback consults them, but handing it the wrong one turns a healthy
    # harvest into a spurious refusal.
    dry_on_disk = util.get_traj_len(
        dry_p, dry_top_p if dry_top_p.is_file() else None)
    down_on_disk = util.get_traj_len(down_p, structure_fn)
    if dry_on_disk != n_dry_expected:
        problems.append(
            f'{dry_p.name} holds {dry_on_disk} frames on disk, planned '
            f'{n_dry_expected}')
    if down_on_disk != n_down_expected:
        problems.append(
            f'{down_p.name} holds {down_on_disk} frames on disk, planned '
            f'{n_down_expected}')
    if problems:
        raise HarvestError(
            f'harvest of {traj_p.parent} did not produce what was planned, so '
            f'{traj_p.name} is being left alone:\n  ' + '\n  '.join(problems))


def _repair_from_symlink(traj_p, dry_p, down_p, sentinel_p, dry_top_p,
                         structure_fn, frames_per_gen, gen_index,
                         first_global_index, downsample_frq, seam=SEAM_AUTO):
    """Finish a harvest that unlinked the original and then died.

    The window is small but a requeueing scheduler will find it. The outputs
    are checked against the same plan the original harvest would have used; if
    they match, the sentinel that was never written gets written now.
    """
    if not (dry_p.is_file() and down_p.is_file()):
        raise HarvestError(
            f'{traj_p} is a symlink -- a previous harvest removed the original '
            f'-- but {dry_p.name} and {down_p.name} are not both present. The '
            'raw trajectory for this generation is gone and cannot be rebuilt.')
    dry_on_disk = util.get_traj_len(
        dry_p, dry_top_p if dry_top_p.is_file() else None)
    down_on_disk = util.get_traj_len(down_p, structure_fn)
    for candidate in (frames_per_gen + 1, frames_per_gen):
        skip_first = resolve_seam(candidate, frames_per_gen, gen_index,
                                  seam=seam)
        n_dry, n_down = expected_counts(candidate, first_global_index,
                                        downsample_frq, skip_first)
        if dry_on_disk == n_dry and down_on_disk == n_down:
            record = dict(
                status='repaired', backend=None, n_orig=candidate, n_dry=n_dry,
                n_down=n_down, first_global_index=first_global_index,
                frames_per_gen=frames_per_gen, gen_index=gen_index,
                downsample_frq=downsample_frq, skip_first=skip_first,
                dry=dry_p.name, downsample=down_p.name, original=traj_p.name,
                unlinked=True)
            sentinel_p.write_text(json.dumps(record, indent=2))
            print(f'[harvest] {traj_p.parent}: a previous harvest completed but '
                  f'never wrote {sentinel_p.name}; counts match the plan '
                  f'({n_dry} dry, {n_down} downsampled), sentinel written.',
                  flush=True)
            return record
    raise HarvestError(
        f'{traj_p} is a symlink, so the original is gone, but {dry_p.name} '
        f'({dry_on_disk} frames) and {down_p.name} ({down_on_disk} frames) do '
        f'not match any expected count for a {frames_per_gen}-interval '
        f'generation. This generation was harvested incompletely and the raw '
        f'trajectory cannot be rebuilt.')


# ---------------------------------------------------------------------------
# Auditing a campaign after the fact
# ---------------------------------------------------------------------------

def unharvested_gen_dirs(top_level, sentinel_name=SENTINEL_NAME,
                         config_name='config.json'):
    """Generation directories that ran but carry no harvest sentinel.

    Harvest failures are deliberately swallowed by the tender -- losing a
    harvest must not stop a campaign -- and the harvest job's own exit status is
    never observed by anything. This is how you find out.
    """
    top = Path(top_level)
    stale = []
    for config_p in sorted(top.glob(f'*/*/*/{config_name}')):
        gen_dir = config_p.parent
        if not (gen_dir / sentinel_name).is_file():
            stale.append(gen_dir)
    return stale


def _frame_times(traj_p, scan_chunk=reimage.SCAN_CHUNK):
    """Every frame's time, without ever holding the coordinates."""
    import numpy as np
    import mdtraj as md
    times = []
    with md.open(str(traj_p)) as fh:
        while True:
            frame = fh.read(scan_chunk)
            time = np.asarray(frame[1])
            if not time.size:
                break
            times.append(time)
            if time.size < scan_chunk:
                break
    return np.concatenate(times) if times else np.array([])


def verify_dry_chain(gen_dirs, structure_fn=None, sentinel_name=SENTINEL_NAME,
                     dry_prefix=DRY_PREFIX, scan_chunk=reimage.SCAN_CHUNK):
    """Check the harvested stream of a clone is contiguous and unduplicated.

    The seam convention never fails loudly, so it gets asserted: across N
    harvested generations the dry stream must hold ``N * frames_per_gen + 1``
    frames, and no two consecutive frames may carry the same time.
    """
    import numpy as np
    import mdtraj as md

    times, records = [], []
    for gen_dir in gen_dirs:
        gen_p = Path(gen_dir)
        sentinel_p = gen_p / sentinel_name
        if not sentinel_p.is_file():
            raise HarvestError(f'{gen_p} has not been harvested')
        record = json.loads(sentinel_p.read_text())
        records.append(record)
        times.append(_frame_times(gen_p / record['dry'], scan_chunk=scan_chunk))
    order = sorted(range(len(records)), key=lambda i: records[i]['gen_index'])
    records = [records[i] for i in order]
    times = [times[i] for i in order]
    gen_indices = [r['gen_index'] for r in records]
    if gen_indices != list(range(gen_indices[0], gen_indices[0] + len(records))):
        raise HarvestError(
            f'generations {gen_indices} are not a contiguous run; concatenating '
            'them would splice across a gap.')
    frames_per_gen = set(r['frames_per_gen'] for r in records)
    if len(frames_per_gen) != 1:
        raise HarvestError(
            f'generations disagree about frames_per_gen ({sorted(frames_per_gen)}); '
            'the concatenated stream would not be evenly spaced.')
    frames_per_gen = frames_per_gen.pop()

    time = np.concatenate(times)
    # Each generation contributes frames_per_gen new frames. The engine's frame
    # at the very start of gen 0 is the only extra one -- for every later
    # generation it duplicates a seam and was dropped at harvest.
    writes_step_zero = records[0]['n_orig'] == frames_per_gen + 1
    expected = len(records) * frames_per_gen + (1 if writes_step_zero else 0)
    spacing = np.diff(time)
    return dict(
        n_frames=len(time), expected=expected, contiguous=len(time) == expected,
        strictly_increasing=bool((spacing > 0).all()) if len(spacing) else True,
        n_duplicate_times=int((spacing == 0).sum()),
        uniform_spacing=bool(np.allclose(spacing, spacing[0])) if len(spacing)
        else True,
        first_time=float(time[0]), last_time=float(time[-1]))
