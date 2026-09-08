"""Reduce a finished generation to the two trajectories we keep.

A harvest writes a dry trajectory (every frame, solute only) and a downsampled
one (every Nth frame, solvent kept), then replaces the original with a symlink
to the dry one. Deleting the original cannot be undone, so the frames written
are counted and compared against a plan first, and a .harvested sentinel is
written at the end so a re-run does nothing instead of reading and writing the
same file at once.

Both backends stream, since a generation can be tens of GB. LOOS is used when
the box is rectangular and mdtraj when it is not.

A generation that finished while the tender was down is never reaped, and
recovery advances the clone past it, so nothing comes back for it. Those are
found by unharvested_gen_dirs and sorted by classify_gen_dir, which reads only:
a generation is called safe to harvest late when the engine's own record says
it reached its step target, that record's target is the one the chain implies,
and the trajectory holds exactly the frames the configs predict. Everything
else is reported for a human, because an unharvested generation is evidence
something went wrong and the harvest deletes the only copy of the original.
"""

import json
import subprocess as sp
from pathlib import Path

from . import utilities as util
from . import reimage
from . import gmx_simulate


# Written last, after both outputs are verified and the original is gone. Its
# presence means the generation is harvested; a re-run stops on it.
SENTINEL_NAME = '.harvested'

# Output name prefixes, joined to the trajectory name with the config's sep.
DRY_PREFIX = 'dry'
DOWNSAMPLE_PREFIX = 'downsample'

# Solute-only topology written beside the dry stream, so anything reading it
# later, get_traj_len's LOOS fallback included, has a matching atom count.
DRY_TOPOLOGY_NAME = 'dry-top.pdb'

# What every generation directory calls its run record.
CONFIG_NAME = util.CONFIG_NAME

# What Harvester.reap writes into a generation directory at submission time.
HARVESTER_CONFIG_NAME = 'hconfig.json'

# The GROMACS runner's own progress record, read as a completion witness.
GEN_STATUS_NAME = gmx_simulate.GEN_STATUS_NAME

# Backend selectors for harvest_generation.
BACKEND_AUTO = 'auto'
BACKEND_LOOS = 'loos'
BACKEND_MDTRAJ = 'mdtraj'

# Frames per md.iterload chunk. Only the mdtraj backend uses this; the LOOS
# backend is already frame-at-a-time.
ITERLOAD_CHUNK = 100

# What to do with the frame that repeats the previous generation's last one.
# GROMACS writes it, the OpenMM reporters do not, and 'auto' tells which by
# whether the file holds one frame more than the generation is long.
SEAM_AUTO = 'auto'
SEAM_DROP = 'drop'
SEAM_KEEP = 'keep'

# Selection-language dialects understood for harvester_subset.
SYNTAX_LOOS = 'loos'
SYNTAX_MDTRAJ = 'mdtraj'

# What classify_gen_dir can conclude about one generation directory.
CATEGORY_HARVESTED = 'harvested'
CATEGORY_COMPLETE = 'complete'
CATEGORY_REPAIRABLE = 'repairable'
CATEGORY_PARTLY_HARVESTED = 'partly-harvested'
CATEGORY_UNFINISHED = 'unfinished'
CATEGORY_INCONSISTENT = 'inconsistent'
CATEGORY_UNPROVEN = 'unproven'
CATEGORY_UNREADABLE = 'unreadable'

# The only two a harvest may act on: one is provably finished, the other has
# already lost its original and needs nothing but its sentinel.
SAFE_CATEGORIES = (CATEGORY_COMPLETE, CATEGORY_REPAIRABLE)

# Loudest first, so a report ends on what a human may safely act on.
CATEGORY_ORDER = (CATEGORY_UNREADABLE, CATEGORY_INCONSISTENT,
                  CATEGORY_PARTLY_HARVESTED, CATEGORY_UNPROVEN,
                  CATEGORY_UNFINISHED, CATEGORY_REPAIRABLE,
                  CATEGORY_COMPLETE, CATEGORY_HARVESTED)

# A recovery harvest keeps the original by default; see harvest_recovered.
RECOVERY_UNLINK = False


class HarvestError(RuntimeError):
    """Raised before the original is removed, so a failed harvest loses nothing."""


class Harvester:
    __slots__ = ('run_config', 'main_path', 'run_config_name',
                 'harvester_template', 'scheduler', 'scriptname')

    def __init__(self,  harvester_template: str,
                 scheduler: str, run_config=None, scriptname='harvest.sh',
                 run_config_name=HARVESTER_CONFIG_NAME):
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
            util.write_json_atomic(current_dir/self.run_config_name,
                                   self.run_config)
        else:
            harvest_script = self.harvester_template
        harvest_script_p = current_dir/self.scriptname
        harvest_script_p.write_text(harvest_script)
        return harvest_script_p

    def reap(self, current_dir, dry_run=False, sentinel_name=SENTINEL_NAME):
        current_dir = Path(current_dir)
        if (current_dir / sentinel_name).is_file():
            # A second harvest job would read and write the same files as the
            # first, at the same time.
            print(f'{current_dir} already carries {sentinel_name}; not '
                  f'submitting another harvest.')
            return None
        harvest_script_p = self.prep_and_write_inputs(current_dir)
        if dry_run:
            return None
        with harvest_script_p.open() as f:
            result = sp.run(self.scheduler, stdin=f, cwd=current_dir,
                            text=True, capture_output=True)
        if result.returncode != 0:
            print(f'{self.scheduler} for {current_dir} exited '
                  f'{result.returncode}')
            print('  stdout:', result.stdout)
            print('  stderr:', result.stderr)
            raise sp.CalledProcessError(result.returncode, self.scheduler,
                                        result.stdout, result.stderr)
        print('harvester scheduler return:', result.stdout)
        return result.stdout


def frames_before(config, downsample_frq, config_name=CONFIG_NAME):
    """Frames the generations before this one wrote, from their own configs."""
    return chain_offsets(config, downsample_frq, config_name=config_name)[0]


def chain_offsets(config, downsample_frq, config_name=CONFIG_NAME):
    """(frames, steps) the generations before this one contributed.

    Counted, not multiplied by this generation's length: a seed may be given a
    different generation length between boots, and a wallclock-matched scheme
    ends generations wherever the clock ran out. The step total is where this
    generation starts, so its target step is that plus its own steps_per_gen.
    """
    frames = steps = 0
    for earlier in util.earlier_gen_configs(
            config['traj_dir_top_level'], config['seed_index'],
            config['clone_index'], config['gen_index'], config['dirname_pad'],
            sep=config['sep'], config_name=config_name):
        frames += check_commensurability(earlier['steps_per_gen'],
                                         earlier['write_interval'],
                                         downsample_frq)
        steps += earlier['steps_per_gen']
    return frames, steps


def check_commensurability(steps_per_gen, write_interval, downsample_frq):
    """Return frames per generation, checking both spacings divide evenly.

    This is the count of NEW frames, one less than the file holds when the
    engine also writes a frame at the step it restarted from.
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


def resolve_seam(n_orig, frames_per_gen, gen_index, seam=SEAM_AUTO):
    """True when this generation's first frame repeats the previous one's last."""
    # GROMACS writes a frame at the step it restarts from; OpenMM does not.
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
        f'written). Refusing to guess which frames are duplicates; pass '
        f'seam={SEAM_DROP!r} or {SEAM_KEEP!r} explicitly if this count is '
        f'expected.')


def keeps_frame(local_index, first_global_index, downsample_frq, skip_first):
    """(write_dry, write_downsample) for one frame.

    The LOOS backend calls this; the mdtraj one works the same rule out in
    numpy over a whole chunk. _verify_counts is what proves the two agreed, by
    comparing what each wrote against what expected_counts predicts.
    """
    # Counting from the frame's place in the whole trajectory, not in this
    # generation, is what keeps the downsample phase running across seams.
    if skip_first and local_index == 0:
        return False, False
    return True, ((first_global_index + local_index) % downsample_frq == 0)


def frame_plan(n_orig, first_global_index, downsample_frq, skip_first):
    """Yield (local_index, write_dry, write_downsample) for every frame.

    Answers "which frames would this harvest keep?" without harvesting, which
    is the cheapest way to check a downsample phase before spending a run on it.
    """
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


def select_backend(traj_fn, structure_fn=None, backend=BACKEND_AUTO,
                   triclinic_rtol=reimage.TRICLINIC_RTOL):
    """Return 'loos' for a rectangular box, 'mdtraj' for anything else."""
    # Looked up rather than caught: LOOS keeps only the diagonal of a triclinic
    # box and raises nothing, so a try/except would never fire.
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
    """0-based atom indices matched by selection, in file order."""
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
    """A LOOS selection string matching exactly these atom indices."""
    # Runs are collapsed so a solute is a few clauses, not thousands. It has
    # to be a selection: LOOS shares the box with a selected group, but gives a
    # hand-built one its own, which then freezes in the dry trajectory.
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
    """Resolve a subset into indices for mdtraj and a selection for LOOS.

    None means keep every atom. Both forms come from one resolution so that
    which backend runs cannot change which atoms get written.
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


def _harvest_loos(traj_fn, structure_fn, subset_spec, dry_out, down_out,
                  first_global_index, downsample_frq, skip_first,
                  timing=None, dry_topology_name=DRY_TOPOLOGY_NAME):
    """Stream the trajectory once, writing both outputs. Rectangular boxes only."""
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
                    iterload_chunk=ITERLOAD_CHUNK):
    """Same, in chunks, for the boxes LOOS cannot represent."""
    # The downsample follows the running frame index, not a [::N] slice of each
    # chunk, which would restart the phase at every chunk boundary.
    import numpy as np
    import mdtraj as md

    model = md.load(str(structure_fn))
    indices = None if subset_spec is None else subset_spec['indices']
    dry_writer = _MdtrajWriter(Path(dry_out))
    down_writer = _MdtrajWriter(Path(down_out))
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
    """Writer that appends, which Trajectory.save() cannot do."""

    def __init__(self, out_p):
        import mdtraj as md
        self.suffix = out_p.suffix.lower()
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
            self.fh.write(traj.xyz * reimage.ANGSTROM_PER_NM,
                          cell_lengths=(None if traj.unitcell_lengths is None
                                        else traj.unitcell_lengths
                                        * reimage.ANGSTROM_PER_NM),
                          cell_angles=traj.unitcell_angles)

    def close(self):
        self.fh.close()


HARVEST_BACKENDS = {BACKEND_LOOS: _harvest_loos, BACKEND_MDTRAJ: _harvest_mdtraj}


def harvest_generation(config_fn, harvester_config_fn,
                       backend=BACKEND_AUTO, seam=SEAM_AUTO,
                       sentinel_name=SENTINEL_NAME,
                       dry_prefix=DRY_PREFIX,
                       downsample_prefix=DOWNSAMPLE_PREFIX,
                       dry_topology_name=DRY_TOPOLOGY_NAME,
                       iterload_chunk=ITERLOAD_CHUNK,
                       triclinic_rtol=reimage.TRICLINIC_RTOL,
                       default_syntax=SYNTAX_LOOS,
                       backends=None,
                       unlink=None):
    """Harvest one generation directory. Safe to re-run and safe to interrupt.

    unlink overrides the harvester config's harvester_unlink; None, the
    default, leaves the decision to that file.

    Keys read from the harvester config file:
      harvester_subset         which atoms the dry trajectory keeps (all, if unset)
      harvester_subset_syntax  'loos' (default) or 'mdtraj'
      harvester_structure      structure file to build the model from. GROMACS
                               runs MUST set this: top_fn is a force field
                               topology, and no reader builds a model from one.
      downsample_frq           keep every Nth frame in the solvated trajectory
      steps_per_gen            full generation length, if not in the run config
      harvester_unlink         delete the original once verified (default True)
    """
    backends = HARVEST_BACKENDS if backends is None else backends
    config = json.loads(Path(config_fn).read_text())
    hconfig = json.loads(Path(harvester_config_fn).read_text())
    gen_dir = Path(config_fn).resolve().parent

    traj_p, dry_p, down_p, dry_top_p, sentinel_p = harvest_paths(
        gen_dir, config, dry_prefix=dry_prefix,
        downsample_prefix=downsample_prefix,
        dry_topology_name=dry_topology_name, sentinel_name=sentinel_name)

    if sentinel_p.is_file():
        record = json.loads(sentinel_p.read_text())
        print(f'[harvest] {sentinel_p} exists; generation already harvested '
              f'({record.get("n_dry")} dry frames). Nothing to do.', flush=True)
        return dict(record, status='already-harvested')

    structure_fn = hconfig.get('harvester_structure') or config['top_fn']
    downsample_frq = hconfig['downsample_frq']
    if downsample_frq < 1:
        raise HarvestError(
            f'downsample_frq is {downsample_frq}; keeping every Nth frame '
            'needs N of at least 1.')
    gen_index = config['gen_index']
    write_interval = config['write_interval']
    steps_per_gen = _steps_per_gen(config, hconfig)
    frames_per_gen = check_commensurability(
        steps_per_gen, write_interval, downsample_frq)
    first_global_index = frames_before(config, downsample_frq)

    # A symlink with no sentinel is a harvest that deleted the original and
    # then died. Writing again here would read and write one file at once.
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
        syntax=hconfig.get('harvester_subset_syntax', default_syntax))

    print(f'[harvest] {gen_dir}: {n_orig} frames, backend={chosen}, '
          f'global frames {first_global_index}..'
          f'{first_global_index + n_orig - 1}, '
          f'{"dropping" if skip_first else "keeping"} the restart-step frame; '
          f'expecting {n_dry_expected} dry and {n_down_expected} downsampled.',
          flush=True)

    kwargs = dict(dry_topology_name=dry_topology_name,
                  timing=util.frame_timing(traj_p, n_orig))
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

    if unlink is None:
        unlink = hconfig.get('harvester_unlink', True)
    if unlink:
        traj_p.unlink()
        # Leave a symlink so frame counting on the original name still works.
        traj_p.symlink_to(dry_p.name)
        record['unlinked'] = True
    else:
        record['unlinked'] = False

    # Written last and atomically: everything above is redoable, and a torn
    # read of this one would say the generation still needs harvesting.
    util.write_json_atomic(sentinel_p, record)
    return record


def harvest_paths(gen_dir, config, dry_prefix=DRY_PREFIX,
                  downsample_prefix=DOWNSAMPLE_PREFIX,
                  dry_topology_name=DRY_TOPOLOGY_NAME,
                  sentinel_name=SENTINEL_NAME):
    """(original, dry, downsampled, dry topology, sentinel) for a generation.

    Named in one place so the classifier looks for exactly the files
    harvest_generation would write.
    """
    gen_p = Path(gen_dir)
    traj_fn = f"{config['traj_name']}{config['traj_suffix']}"
    return (gen_p / traj_fn,
            gen_p / f"{dry_prefix}{config['sep']}{traj_fn}",
            gen_p / f"{downsample_prefix}{config['sep']}{traj_fn}",
            gen_p / dry_topology_name,
            gen_p / sentinel_name)


def _steps_per_gen(config, hconfig):
    """The full generation length, which is not always config['steps'].

    On a resumed generation config['steps'] is the steps still owed, which
    would shorten frames_per_gen and misplace every later generation's global
    frame index, and with it the downsample phase. Clone records the untouched
    value as steps_per_gen.
    """
    from_run, from_harvester = config.get('steps_per_gen'), hconfig.get('steps_per_gen')
    if None not in (from_run, from_harvester) and from_run != from_harvester:
        raise HarvestError(
            f'the run config says steps_per_gen={from_run} and the harvester '
            f'config says {from_harvester}. They decide where this generation '
            'sits in the whole trajectory, so guessing between them would put '
            'the downsample phase in the wrong place.')
    if from_run is not None:
        return from_run          # written by Clone, and validated there
    if from_harvester is not None:
        return from_harvester
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

    The in-loop tallies catch a plan/backend disagreement; re-reading the files
    catches a writer that returned normally having flushed less than it was
    handed, which is what a job killed inside the write loop looks like.
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
    # Each stream gets the topology that matches it: the dry one is the
    # subset, the downsampled one is still the whole system.
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


def repair_candidates(n_dry_on_disk, n_down_on_disk, frames_per_gen, gen_index,
                      first_global_index, downsample_frq, seam=SEAM_AUTO):
    """Original frame counts whose plan produces the outputs already on disk.

    Empty means the two outputs are not what any harvest of this generation
    would have written, so the deleted original cannot be accounted for.
    """
    matches = []
    for candidate in (frames_per_gen + 1, frames_per_gen):
        skip_first = resolve_seam(candidate, frames_per_gen, gen_index,
                                  seam=seam)
        n_dry, n_down = expected_counts(candidate, first_global_index,
                                        downsample_frq, skip_first)
        if (n_dry_on_disk, n_down_on_disk) == (n_dry, n_down):
            matches.append((candidate, skip_first, n_dry, n_down))
    return matches


def _repair_from_symlink(traj_p, dry_p, down_p, sentinel_p, dry_top_p,
                         structure_fn, frames_per_gen, gen_index,
                         first_global_index, downsample_frq, seam=SEAM_AUTO):
    """Finish a harvest that unlinked the original and then died.

    The outputs are checked against the same plan the harvest would have used;
    if they match, the missing sentinel is written now.
    """
    if not (dry_p.is_file() and down_p.is_file()):
        raise HarvestError(
            f'{traj_p} is a symlink, so a previous harvest removed the original, '
            f'-- but {dry_p.name} and {down_p.name} are not both present. The '
            'raw trajectory for this generation is gone and cannot be rebuilt.')
    dry_on_disk = util.get_traj_len(
        dry_p, dry_top_p if dry_top_p.is_file() else None)
    down_on_disk = util.get_traj_len(down_p, structure_fn)
    matches = repair_candidates(dry_on_disk, down_on_disk, frames_per_gen,
                                gen_index, first_global_index, downsample_frq,
                                seam=seam)
    if matches:
        # After generation 0 both candidates produce the same two counts, so
        # how many frames the deleted original held cannot be recovered from
        # them. The outputs are right either way; ambiguous marks the fields
        # of this record that are a guess.
        candidate, skip_first, n_dry, n_down = matches[0]
        ambiguous = len(matches) > 1
        record = dict(
            status='repaired', backend=None, n_orig=candidate, n_dry=n_dry,
            n_down=n_down, first_global_index=first_global_index,
            frames_per_gen=frames_per_gen, gen_index=gen_index,
            downsample_frq=downsample_frq, skip_first=skip_first,
            dry=dry_p.name, downsample=down_p.name, original=traj_p.name,
            unlinked=True, ambiguous=ambiguous)
        util.write_json_atomic(sentinel_p, record)
        print(f'[harvest] {traj_p.parent}: a previous harvest completed but '
              f'never wrote {sentinel_p.name}; counts match the plan '
              f'({n_dry} dry, {n_down} downsampled), sentinel written.'
              + (f' Both {frames_per_gen} and {frames_per_gen + 1} original '
                 f'frames fit those counts, so n_orig={candidate} and '
                 f'skip_first={skip_first} are recorded as a guess '
                 f'(ambiguous).' if ambiguous else ''), flush=True)
        return record
    raise HarvestError(
        f'{traj_p} is a symlink, so the original is gone, but {dry_p.name} '
        f'({dry_on_disk} frames) and {down_p.name} ({down_on_disk} frames) do '
        f'not match any expected count for a {frames_per_gen}-interval '
        f'generation. This generation was harvested incompletely and the raw '
        f'trajectory cannot be rebuilt.')


def unharvested_gen_dirs(top_level, sentinel_name=SENTINEL_NAME,
                         config_name='config.json', skip_newest=False):
    """Generation directories that ran but carry no harvest sentinel.

    The tender swallows harvest failures so a lost harvest cannot stop a
    campaign, and nothing observes the harvest job's exit status. This does.

    On a campaign that is still running, each clone's newest generation is
    normally the one in flight and has no sentinel yet; skip_newest leaves
    those out. On a finished campaign leave it False, since the newest
    generation is exactly the one whose lost harvest you want to hear about.
    """
    top = Path(top_level)
    by_clone = {}
    for config_p in sorted(top.glob(f'*/*/*/{config_name}')):
        gen_dir = config_p.parent
        by_clone.setdefault(gen_dir.parent, []).append(gen_dir)
    stale = []
    for clone_dir, gen_dirs in sorted(by_clone.items()):
        gen_dirs.sort(key=lambda p: _gen_sort_key(p))
        if skip_newest:
            gen_dirs = gen_dirs[:-1]
        stale += [d for d in gen_dirs if not (d / sentinel_name).is_file()]
    return stale


def _gen_sort_key(gen_dir):
    """Order generation directories by their number, not by their name."""
    tail = gen_dir.name.rsplit('-', 1)[-1].rsplit('_', 1)[-1]
    return (0, int(tail)) if tail.isdigit() else (1, 0)


def _frame_times(traj_p, scan_chunk=reimage.SCAN_CHUNK):
    """Every frame's time, without ever holding the coordinates.

    None for a DCD, whose second field is cell lengths and which carries no
    per-frame time at all.
    """
    import numpy as np
    import mdtraj as md
    if Path(traj_p).suffix.lower() != '.xtc':
        return None
    times = []
    with md.open(str(traj_p)) as fh:
        while True:
            time = np.asarray(fh.read(scan_chunk)[1])
            if not time.size:
                break
            times.append(time)
            if time.size < scan_chunk:
                break
    return np.concatenate(times) if times else np.array([])


def verify_dry_chain(gen_dirs, sentinel_name=SENTINEL_NAME,
                     scan_chunk=reimage.SCAN_CHUNK):
    """Check the harvested stream of a clone is contiguous and unduplicated.

    Across N harvested generations the dry stream holds
    N * frames_per_gen + 1 frames, and no two consecutive frames carry the
    same time. Neither fails loudly on its own, hence the check.
    """
    import numpy as np

    if not gen_dirs:
        raise HarvestError('verify_dry_chain was given no generations')
    times, records = [], []
    for gen_dir in gen_dirs:
        gen_p = Path(gen_dir)
        sentinel_p = gen_p / sentinel_name
        if not sentinel_p.is_file():
            raise HarvestError(f'{gen_p} has not been harvested')
        record = json.loads(sentinel_p.read_text())
        records.append(record)
        frame_times = _frame_times(gen_p / record['dry'],
                                   scan_chunk=scan_chunk)
        if frame_times is None:
            raise HarvestError(
                f'{gen_p / record["dry"]} carries no per-frame time, so the '
                'spacing of the chain cannot be checked.')
        times.append(frame_times)
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
    # Every generation contributes frames_per_gen new frames. Only generation
    # 0's frame at step 0 is extra; the later ones were dropped as seams, so a
    # chain that starts partway through a campaign has no extra frame at all.
    writes_step_zero = (records[0]['gen_index'] == 0
                        and records[0]['n_orig'] == frames_per_gen + 1)
    expected = len(records) * frames_per_gen + (1 if writes_step_zero else 0)
    spacing = np.diff(time)
    return dict(
        n_frames=len(time), expected=expected, contiguous=len(time) == expected,
        strictly_increasing=bool((spacing > 0).all()) if len(spacing) else True,
        n_duplicate_times=int((spacing == 0).sum()),
        uniform_spacing=bool(np.allclose(spacing, spacing[0])) if len(spacing)
        else True,
        first_time=float(time[0]), last_time=float(time[-1]))


def _read_json(path):
    """The JSON object at path, or None if it is missing or will not parse."""
    p = Path(path)
    if not p.is_file():
        return None
    try:
        return json.loads(p.read_text())
    except (json.JSONDecodeError, UnicodeDecodeError):
        print(f'[harvest] {p} will not parse as JSON; ignoring.', flush=True)
        return None


def _verdict(row, category, reason):
    """Stamp a category and the evidence for it onto a classification row."""
    row['category'] = category
    row['reason'] = reason
    return row


def _resolve_hconfig(gen_p, hconfig, harvester_config_name):
    """(harvester config, where it came from) for one generation.

    A generation the tender never reaped has no harvester config of its own,
    since Harvester.reap writes that at submission time. A dict or path given
    by the caller wins; otherwise the nearest sibling generation's is borrowed,
    and _steps_per_gen refuses it if it belongs to a different generation
    length.
    """
    if isinstance(hconfig, dict):
        return hconfig, 'caller'
    if hconfig is not None:
        return _read_json(hconfig), str(hconfig)
    own = _read_json(gen_p / harvester_config_name)
    if own is not None:
        return own, harvester_config_name
    return _sibling_hconfig(gen_p, harvester_config_name)


def _sibling_hconfig(gen_p, harvester_config_name):
    """The harvester config of the nearest other generation in this clone."""
    here = _gen_sort_key(gen_p)[1]
    candidates = sorted(
        (p for p in gen_p.parent.glob(f'*/{harvester_config_name}')
         if p.parent != gen_p),
        key=lambda p: abs(_gen_sort_key(p.parent)[1] - here))
    for candidate in candidates:
        found = _read_json(candidate)
        if found is not None:
            return found, str(candidate)
    return None, None


def _gromacs_witness(gen_p, restart_p, target_step, gen_status_name):
    """gen_status.json, plus the checkpoint the next generation is seeded from."""
    status = gmx_simulate.read_gen_status(gen_p,
                                          gen_status_name=gen_status_name)
    if status is None:
        return gen_status_name, None, f'{gen_status_name} will not parse'
    reached, recorded = status.get('reached_step'), status.get('target_step')
    if not status.get('complete'):
        return gen_status_name, False, (
            f'{gen_status_name} says complete=false, at step {reached} of '
            f'{recorded}')
    if recorded is not None and recorded != target_step:
        return gen_status_name, None, (
            f'{gen_status_name} targets step {recorded}, but the generations '
            f'before this one plus its own length end at {target_step}')
    if reached is not None and reached < target_step:
        return gen_status_name, None, (
            f'{gen_status_name} says complete but reached only {reached} of '
            f'{target_step}')
    if not gmx_simulate.is_checkpoint(restart_p):
        return gen_status_name, None, (
            f'{gen_status_name} says complete but {restart_p.name} is not a '
            'GROMACS checkpoint, so nothing was seeded from this generation')
    return gen_status_name, True, (
        f'{gen_status_name}: complete at step {reached} of {target_step}')


def _openmm_witness(restart_p, target_step):
    """state.xml's stepCount, which accumulates across generations."""
    if not util.is_state_xml_usable(restart_p):
        return restart_p.name, None, f'{restart_p.name} will not deserialize'
    try:
        state_step = util.state_xml_step_count(restart_p)
    except ValueError as exc:
        return restart_p.name, None, str(exc)
    detail = f'{restart_p.name}: stepCount {state_step} of {target_step}'
    if state_step == target_step:
        return restart_p.name, True, detail
    if state_step < target_step:
        return restart_p.name, False, detail
    return restart_p.name, None, detail + ', past the end of this generation'


def completion_witness(gen_dir, config, target_step,
                       gen_status_name=GEN_STATUS_NAME):
    """(witness, finished, detail): does the engine's own record say it ended?

    witness names the record that was read and is None when the generation left
    none, which is not the same as saying it did not finish. finished is None
    when the record contradicts itself or the rest of the chain.

    Both engines are asked the same question in the terms each records it:
    GROMACS writes gen_status.json only after its last mdrun, and OpenMM leaves
    a state.xml whose stepCount is absolute across generations.
    """
    gen_p = Path(gen_dir)
    restart_p = gen_p / config['restart_name']
    if (gen_p / gen_status_name).is_file():
        return _gromacs_witness(gen_p, restart_p, target_step, gen_status_name)
    if gmx_simulate.is_checkpoint(restart_p):
        return None, False, (
            f'a GROMACS checkpoint but no {gen_status_name}, so nothing '
            'records whether the run reached its target')
    if restart_p.is_file():
        return _openmm_witness(restart_p, target_step)
    return None, False, f'no {gen_status_name} and no {restart_p.name}'


def classify_gen_dir(gen_dir, hconfig=None, seam=SEAM_AUTO,
                     sentinel_name=SENTINEL_NAME, config_name=CONFIG_NAME,
                     harvester_config_name=HARVESTER_CONFIG_NAME,
                     gen_status_name=GEN_STATUS_NAME,
                     dry_prefix=DRY_PREFIX,
                     downsample_prefix=DOWNSAMPLE_PREFIX,
                     dry_topology_name=DRY_TOPOLOGY_NAME):
    """Sort one generation directory into a category, with its evidence.

    Reads; writes nothing. The row carries gen_dir, gen_index, category and
    reason, plus whichever of steps_per_gen, frames_per_gen, target_step,
    first_global_index, n_orig, skip_first, n_dry, n_down and witness were
    established before the verdict.

    Only CATEGORY_COMPLETE says a late harvest is safe, and it needs all of:
    both configs parse; the generation length is recorded rather than guessed;
    every earlier generation's config is present, so the chain places this one;
    the engine's own record says it reached the step the chain implies; the
    checkpoint that record depends on is readable; and the trajectory holds
    exactly the frames that length predicts, seam resolved. Anything short of
    that is a human's decision, not this function's.
    """
    gen_p = Path(gen_dir)
    row = dict(gen_dir=str(gen_p), gen_index=None, category=None, reason='')
    if (gen_p / sentinel_name).is_file():
        return _verdict(row, CATEGORY_HARVESTED, f'{sentinel_name} is present')
    config = _read_json(gen_p / config_name)
    if config is None:
        return _verdict(row, CATEGORY_UNREADABLE, f'no readable {config_name}')
    row['gen_index'] = config.get('gen_index')
    hconfig, row['hconfig_source'] = _resolve_hconfig(
        gen_p, hconfig, harvester_config_name)
    if hconfig is None:
        return _verdict(row, CATEGORY_UNREADABLE,
                        f'no {harvester_config_name} here or in a sibling '
                        'generation, so the harvest plan is unknown')
    try:
        return _classify_against_plan(
            gen_p, config, hconfig, row, seam=seam,
            sentinel_name=sentinel_name, gen_status_name=gen_status_name,
            dry_prefix=dry_prefix, downsample_prefix=downsample_prefix,
            dry_topology_name=dry_topology_name)
    except KeyError as exc:
        return _verdict(row, CATEGORY_UNREADABLE,
                        f'a config in this chain has no {exc} entry')
    except FileNotFoundError as exc:
        return _verdict(row, CATEGORY_UNREADABLE, str(exc))
    except (HarvestError, ValueError) as exc:
        return _verdict(row, CATEGORY_INCONSISTENT, str(exc))


def _classify_against_plan(gen_p, config, hconfig, row, seam, sentinel_name,
                           gen_status_name, dry_prefix, downsample_prefix,
                           dry_topology_name):
    """The part of classify_gen_dir that needs both configs to have parsed."""
    traj_p, dry_p, down_p, dry_top_p, _ = harvest_paths(
        gen_p, config, dry_prefix=dry_prefix,
        downsample_prefix=downsample_prefix,
        dry_topology_name=dry_topology_name, sentinel_name=sentinel_name)
    structure_fn = hconfig.get('harvester_structure') or config['top_fn']
    downsample_frq = hconfig['downsample_frq']
    if downsample_frq < 1:
        raise HarvestError(
            f'downsample_frq is {downsample_frq}; keeping every Nth frame '
            'needs N of at least 1.')
    if (config.get('steps_per_gen') is None
            and hconfig.get('steps_per_gen') is None):
        return _verdict(row, CATEGORY_UNPROVEN,
                        'neither config records steps_per_gen, so this '
                        "generation's length would have to be guessed")
    steps_per_gen = _steps_per_gen(config, hconfig)
    frames_per_gen = check_commensurability(
        steps_per_gen, config['write_interval'], downsample_frq)
    first_global_index, chain_steps = chain_offsets(config, downsample_frq)
    row.update(steps_per_gen=steps_per_gen, frames_per_gen=frames_per_gen,
               first_global_index=first_global_index,
               target_step=chain_steps + steps_per_gen)

    if traj_p.is_symlink():
        return _classify_symlinked(
            row, traj_p, dry_p, down_p, dry_top_p, structure_fn,
            frames_per_gen=frames_per_gen, gen_index=config['gen_index'],
            first_global_index=first_global_index,
            downsample_frq=downsample_frq, seam=seam,
            sentinel_name=sentinel_name)
    leftovers = [p.name for p in (dry_p, down_p, dry_top_p) if p.exists()]
    if leftovers:
        return _verdict(row, CATEGORY_PARTLY_HARVESTED,
                        f'{", ".join(leftovers)} present with no '
                        f'{sentinel_name}, but {traj_p.name} is still the '
                        'original: a harvest died before it swapped them')
    if not traj_p.is_file():
        return _verdict(row, CATEGORY_UNFINISHED, f'no {traj_p.name}')
    row['n_orig'] = util.get_traj_len(traj_p, structure_fn)
    if not row['n_orig']:
        return _verdict(row, CATEGORY_UNREADABLE,
                        f'{traj_p.name} holds no readable frames')
    # Raises unless the count is one of the two a generation of this length can
    # hold, which is what proves the trajectory against the config.
    row['skip_first'] = resolve_seam(row['n_orig'], frames_per_gen,
                                     config['gen_index'], seam=seam)

    witness, finished, detail = completion_witness(
        gen_p, config, row['target_step'], gen_status_name=gen_status_name)
    row['witness'] = witness
    if witness is None:
        return _verdict(row, CATEGORY_UNPROVEN, detail)
    if finished is None:
        return _verdict(row, CATEGORY_INCONSISTENT, detail)
    if not finished:
        return _verdict(row, CATEGORY_UNFINISHED, detail)
    return _verdict(row, CATEGORY_COMPLETE,
                    f'{detail}; {traj_p.name} holds {row["n_orig"]} frames for '
                    f'a {frames_per_gen}-frame generation')


def _classify_symlinked(row, traj_p, dry_p, down_p, dry_top_p, structure_fn,
                        frames_per_gen, gen_index, first_global_index,
                        downsample_frq, seam, sentinel_name):
    """A generation whose original is already a symlink: a harvest died late."""
    if not (dry_p.is_file() and down_p.is_file()):
        return _verdict(row, CATEGORY_PARTLY_HARVESTED,
                        f'{traj_p.name} is a symlink, so the original is gone, '
                        f'but {dry_p.name} and {down_p.name} are not both here')
    row['n_dry'] = util.get_traj_len(
        dry_p, dry_top_p if dry_top_p.is_file() else None)
    row['n_down'] = util.get_traj_len(down_p, structure_fn)
    counts = (f'{dry_p.name} holds {row["n_dry"]} frames and {down_p.name} '
              f'{row["n_down"]}')
    if repair_candidates(row['n_dry'], row['n_down'], frames_per_gen, gen_index,
                         first_global_index, downsample_frq, seam=seam):
        return _verdict(row, CATEGORY_REPAIRABLE,
                        f'the original is already gone and {counts}, which '
                        f'matches the plan; only {sentinel_name} is missing')
    return _verdict(row, CATEGORY_PARTLY_HARVESTED,
                    f'{traj_p.name} is a symlink, so the original is gone, but '
                    f'{counts}, which matches no plan for a '
                    f'{frames_per_gen}-frame generation')


def classify_campaign(top_level, hconfig=None, seam=SEAM_AUTO,
                      skip_newest=False, sentinel_name=SENTINEL_NAME,
                      config_name=CONFIG_NAME,
                      harvester_config_name=HARVESTER_CONFIG_NAME,
                      gen_status_name=GEN_STATUS_NAME):
    """Classify every unharvested generation in a campaign. Reads only.

    skip_newest leaves out each clone's newest generation, which on a running
    campaign is the one in flight rather than one whose harvest was lost.
    """
    return [classify_gen_dir(gen_dir, hconfig=hconfig, seam=seam,
                             sentinel_name=sentinel_name,
                             config_name=config_name,
                             harvester_config_name=harvester_config_name,
                             gen_status_name=gen_status_name)
            for gen_dir in unharvested_gen_dirs(
                top_level, sentinel_name=sentinel_name,
                config_name=config_name, skip_newest=skip_newest)]


def format_report(rows, category_order=CATEGORY_ORDER,
                  safe_categories=SAFE_CATEGORIES,
                  harvester_config_name=HARVESTER_CONFIG_NAME):
    """The per-generation table a human reads before deciding anything."""
    counts = ', '.join(f'{c} {sum(r["category"] == c for r in rows)}'
                       for c in category_order
                       if any(r['category'] == c for r in rows))
    lines = [f'unharvested generations: {counts or "none"}']
    for category in category_order:
        in_category = [r for r in rows if r['category'] == category]
        if not in_category:
            continue
        lines.append(f'\n{category} ({len(in_category)})')
        for row in in_category:
            lines.append(f'  {row["gen_dir"]}')
            lines.append(f'      {row["reason"]}')
            borrowed = row.get('hconfig_source')
            if borrowed not in (None, harvester_config_name):
                lines.append(f'      judged against the harvester config at '
                             f'{borrowed}')
    n_safe = sum(r['category'] in safe_categories for r in rows)
    lines.append(f'\n{n_safe} of {len(rows)} are safe to harvest without a '
                 f'human looking first ({", ".join(safe_categories)}); the '
                 'rest are left alone.')
    return '\n'.join(lines)


def harvest_recovered(rows, hconfig=None, unlink=RECOVERY_UNLINK,
                      seam=SEAM_AUTO, safe_categories=SAFE_CATEGORIES,
                      config_name=CONFIG_NAME,
                      harvester_config_name=HARVESTER_CONFIG_NAME):
    """Harvest the classified generations that were proved safe, and no others.

    unlink defaults to False here, unlike a harvest the tender submits: these
    generations went unwitnessed, so the original is kept and the disk is
    reclaimed by hand once the dry copies have been looked at. Pass unlink=True
    to accept _verify_counts alone, as the tender's own harvest does.

    The harvest runs in this process, on the same code path the harvest job
    runs; a campaign with many generations to recover is better handed to the
    scheduler one directory at a time.
    """
    results = []
    for row in rows:
        if row['category'] not in safe_categories:
            continue
        gen_p = Path(row['gen_dir'])
        hconfig_p = gen_p / harvester_config_name
        if not hconfig_p.is_file():
            # harvest_generation reads the plan from this directory, and a
            # generation the tender never reaped has none of its own.
            borrowed, source = _resolve_hconfig(gen_p, hconfig,
                                                harvester_config_name)
            print(f'[harvest] {gen_p}: no {harvester_config_name} of its own; '
                  f'writing the one from {source}.', flush=True)
            util.write_json_atomic(hconfig_p, borrowed)
        results.append(harvest_generation(gen_p / config_name, hconfig_p,
                                          seam=seam, unlink=unlink))
    return results


def _main(argv=None):
    """Report what a campaign never harvested, and only then act on it."""
    import argparse
    ap = argparse.ArgumentParser(
        prog='python -m mdfarmer harvest',
        description='Report the generations a campaign never harvested, and '
                    'optionally harvest the ones that are provably safe.')
    ap.add_argument('top_level', help="the campaign's traj_dir_top_level")
    ap.add_argument('--skip-newest', action='store_true',
                    help="leave out each clone's newest generation, which on a "
                         'running campaign is the one in flight')
    ap.add_argument('--hconfig', default=None,
                    help='harvester config to judge every generation against; '
                         "by default each borrows a sibling's")
    ap.add_argument('--harvest', action='store_true',
                    help=f'harvest the {"/".join(SAFE_CATEGORIES)} generations. '
                         'Without this nothing is modified.')
    ap.add_argument('--unlink', action='store_true', default=RECOVERY_UNLINK,
                    help='let the harvest replace each original with a symlink '
                         'to its dry copy; off on this path')
    args = ap.parse_args(argv)
    rows = classify_campaign(args.top_level, hconfig=args.hconfig,
                             skip_newest=args.skip_newest)
    print(format_report(rows), flush=True)
    if args.harvest:
        harvest_recovered(rows, hconfig=args.hconfig, unlink=args.unlink)
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(_main())
