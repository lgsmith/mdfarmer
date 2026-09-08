"""Run one GROMACS generation, in place of simulate.omm_generation.

It takes the same arguments the OpenMM runner does, so Farmer and Clone drive
GROMACS without knowing which engine they have. Nothing here imports OpenMM; it
shells out to the gmx binary.

Generations chain with gmx convert-tpr, to extend the step budget, plus
mdrun -cpi to resume from the last generation's checkpoint. This is the only
exact continuation: grompp -t reads coordinates, velocities and box size and
nothing else, which silently zeroes the Nose-Hoover and Parrinello-Rahman
integrals at every generation boundary. The chain also keeps step and time
counting from the start of the run rather than restarting at zero, which is what
lets the generations be put back together in order.

Each launch writes its own prod.partNNNN.xtc, because mdrun -cpi will not append
into a directory that does not already hold the files its checkpoint names.
concat_parts merges them when the generation finishes, reading with mdtraj so
every frame keeps its own step, time and box, and writing with mdtraj unless the
run's compressed-x-precision is finer than the 1000 that writer is fixed at, in
which case LOOS writes instead. Where two parts cover the same steps the later
one's frames survive, so a part left behind by a rewound relaunch is moved aside
before mdrun runs rather than left to outrank the branch that actually
continued.

Whether a generation is finished is read from the checkpoint's step counter, not
by counting frames. GROMACS writes a frame at step 0 too, so frame counting is
off by one write interval, and a generation killed inside its last interval
would look finished when it is not.

Per-seed Farmer inputs map onto GROMACS as:
    seed_structure_fns -> the .gro structure, also the generation-0 seed
    top_fns            -> the .top topology
    system_fns         -> the .mdp run parameters, see mdp_fn
"""

import json
import os
import re
import shutil
import signal
import struct
import subprocess as sp
from pathlib import Path

import numpy as np

from . import utilities as util
from .reimage import ANGSTROM_PER_NM, TRICLINIC_RTOL, BoxTypeError, is_orthorhombic


# Preempt sentinel the SIGTERM trap touches; same name as simulate's.
PREEMPT_SENTINEL_NAME = 'PREEMPT_SIGTERM'

# Written after each mdrun so the orchestrator can judge progress without gmx.
GEN_STATUS_NAME = 'gen_status.json'

# Incoming seed checkpoint, moved off the name mdrun writes its own one to.
SEED_CPT_NAME = 'seed.cpt'

# Per-generation tpr; gen N>0 is convert-tpr'd from gen N-1's file of this name.
TPR_NAME = 'prod.tpr'

# mdrun -deffnm stem; also the prefix of the .partNNNN outputs.
DEFFNM = 'prod'

# Prefix that hides a stale part from part_files' glob without deleting it.
ABANDONED_PART_PREFIX = 'abandoned-'

# First bytes of any GROMACS checkpoint: tells one from a .gro without gmx.
CHECKPOINT_MAGIC = b'\x00\x02\x9f\x29'

# mdrun checkpoint period. The GROMACS default of 15 is what a hard kill costs.
CHECKPOINT_MINUTES = 5

# gen-seed = base + GEN_SEED_STRIDE * seed + clone; must exceed n_clones.
GEN_SEED_STRIDE = 1000

# Structure formats grompp -c reads; a continuing gen is seeded with a .cpt.
GROMPP_STRUCTURE_SUFFIXES = ('.gro', '.g96', '.pdb', '.brk', '.ent')

# gmx_pack passes these at run time; a config template must not also carry them.
RUNTIME_ONLY_KEYS = ('fleet', 'fleet_key', 'grompp_lock')

# Filled in per generation by a Clone; gmx_config_template leaves placeholders.
CLONE_FILLED_KEYS = ('seed_index', 'clone_index', 'gen_index', 'seed_fn',
                     'top_fn')

# Seconds between preempt-sentinel polls while mdrun runs.
PREEMPT_POLL_SECONDS = 5

# Default binary. Sites with an MPI-only build have gmx_mpi instead.
GMX_BIN = 'gmx'

# Trajectory suffixes the part merge carries over without losing anything; a
# .trr would arrive with velocities and forces mdtraj does not read back out.
MERGEABLE_SUFFIXES = ('.xtc',)

# Frames read per chunk while merging, so a multi-GB part is never resident.
CONCAT_CHUNK = 200

# Coordinate agreement the merge demands of its own output, nm. Far under the
# 5e-4 nm a silent drop to mdtraj's fixed 1e-3 xtc precision would cost.
CONCAT_ATOL_NM = 1e-6

# The only precision mdtraj's .xtc writer emits; it takes no parameter.
MDTRAJ_XTC_PRECISION = 1000.0

# Byte offset of the precision in an .xtc frame: magic, atom count, step and
# time, then nine box floats, then the count that opens the coordinate block.
XTC_PRECISION_OFFSET = 56

# Atoms at or under which GROMACS stores coordinates raw, with no precision.
XTC_UNCOMPRESSED_ATOMS = util.XTC_UNCOMPRESSED_ATOMS

# Nominal clock LOOS's XTCWriter stamps frames with when it is not told one.
# Every frame here is written with its own step and time, so it is never used.
LOOS_WRITER_DT = 1.0
LOOS_WRITER_STEPS_PER_FRAME = 1


class Preempted(Exception):
    """Raised when a preempt sentinel is seen mid-mdrun; the gen is incomplete
    and will resume via -cpi on the next launch."""
    pass


class GenIncomplete(Exception):
    """Raised when mdrun returned cleanly but the generation has not reached its
    step target, for instance because mdrun stopped itself at -maxh. It
    resumes on the next launch; this is a normal event, not a failure."""
    pass


# Default run.py the Farmer writes into each gen dir; the batch script runs it.
default_gmx_run_script = """
from mdfarmer.gmx_simulate import gmx_basic_sim_block_json as runner
runner('config.json')
"""


def _norm_mdp_key(k):
    """Canonical .mdp key: GROMACS is case-insensitive and reads '_' as '-'."""
    return k.strip().lower().replace('_', '-')


# Spellings GROMACS removed: grompp fatals on them, so map to the modern name.
LEGACY_MDP_KEYS = {
    'nstxtcout': 'nstxout-compressed',
    'xtc-precision': 'compressed-x-precision',
    'xtc-grps': 'compressed-x-grps',
    'unconstrained-start': 'continuation',
}


def write_gen_mdp(base_mdp, out_mdp, *, nsteps, nstxout_compressed,
                  gen_vel, continuation, gen_seed=None, gen_temp=None,
                  ld_seed=None, legacy_mdp_keys=LEGACY_MDP_KEYS):
    """Copy base_mdp to out_mdp, changing only the per-generation control keys.

    Everything else is inherited verbatim. Only generation 0 needs this; later
    generations inherit their parameters from the previous tpr through
    convert-tpr, which is what keeps them exact.
    """
    overrides = {
        'nsteps': str(int(nsteps)),
        # nsteps is the absolute cumulative target, so this run counts from 0.
        'init-step': '0',
        'nstxout-compressed': str(int(nstxout_compressed)),
        'gen-vel': 'yes' if gen_vel else 'no',
        'continuation': 'yes' if continuation else 'no',
    }
    # Set even without gen-vel: ld-seed drives the thermostat for the whole run.
    if ld_seed is not None:
        overrides['ld-seed'] = str(int(ld_seed))
    if gen_vel:
        if gen_seed is not None:
            overrides['gen-seed'] = str(int(gen_seed))
        if gen_temp is not None:
            overrides['gen-temp'] = str(gen_temp)
    targets = set(overrides)
    seen = set()
    out_lines = []
    for line in Path(base_mdp).read_text().splitlines():
        code = line.split(';', 1)[0]
        if '=' in code:
            key = _norm_mdp_key(code.split('=', 1)[0])
            key = legacy_mdp_keys.get(key, key)
            if key in targets:
                if key == 'init-step' and code.split('=', 1)[1].strip() != '0':
                    print(f'[gmx] {base_mdp} sets {line.strip()}; this runner '
                          'needs nsteps to be the absolute step target, so '
                          'init-step is pinned to 0.', flush=True)
                if key in seen:
                    # A later duplicate would override ours, so drop it.
                    continue
                out_lines.append(f'{key} = {overrides[key]}')
                seen.add(key)
                continue
        out_lines.append(line)
    for key in targets - seen:
        out_lines.append(f'{key} = {overrides[key]}')
    Path(out_mdp).write_text('\n'.join(out_lines) + '\n')


def _run(cmd, cwd):
    print('[gmx]', ' '.join(map(str, cmd)), flush=True)
    result = sp.run([str(c) for c in cmd], cwd=str(cwd), text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f'command exited {result.returncode}: {" ".join(map(str, cmd))}')


def _run_capture(cmd, cwd=None):
    result = sp.run([str(c) for c in cmd], cwd=None if cwd is None else str(cwd),
                    text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            f'command exited {result.returncode}: {" ".join(map(str, cmd))}\n'
            f'{result.stderr[-2000:]}')
    return result.stdout


_STEP_RE = re.compile(r'^\s*step\s*=\s*(\d+)', re.MULTILINE)
_PART_RE = re.compile(r'^\s*simulation part\s*#\s*=\s*(\d+)', re.MULTILINE)


def checkpoint_part_step(cpt_fn, gmx_bin=GMX_BIN):
    """(simulation part #, step) recorded in a GROMACS checkpoint.

    The part number is the part mdrun was writing when the checkpoint was
    saved; a relaunch with -cpi on this checkpoint writes part number + 1.
    """
    out = _run_capture([gmx_bin, 'dump', '-cp', str(cpt_fn)])
    step_match = _STEP_RE.search(out)
    if step_match is None:
        raise ValueError(f'no step counter in checkpoint {cpt_fn}')
    part_match = _PART_RE.search(out)
    if part_match is None:
        raise ValueError(f'no simulation part counter in checkpoint {cpt_fn}')
    return int(part_match.group(1)), int(step_match.group(1))


def checkpoint_step(cpt_fn, gmx_bin=GMX_BIN):
    """The step counter in a GROMACS checkpoint: how far this generation got.

    Frame counts cannot answer that. GROMACS writes a frame at step 0, and the
    checkpoint is allowed to lag the last frame written.
    """
    return checkpoint_part_step(cpt_fn, gmx_bin=gmx_bin)[1]


def is_checkpoint(cpt_fn, magic=CHECKPOINT_MAGIC):
    """True when this file carries the GROMACS checkpoint magic number.

    Only asks whether the file is a checkpoint at all, which at generation 0 it
    is not: the seed is a .gro, and handing that to -cpi is fatal. A checkpoint
    gmx cannot read is a different problem, and checkpoint_part_step raises on
    it rather than letting it look like a generation that never started.
    """
    path = Path(cpt_fn)
    if not path.is_file():
        return False
    with path.open('rb') as handle:
        return handle.read(len(magic)) == magic


def part_files(gen_dir, deffnm=DEFFNM, traj_suffix='.xtc'):
    """The prod.partNNNN.<suffix> files a generation has accumulated, in order."""
    gen_p = Path(gen_dir)
    parts = sorted(gen_p.glob(f'{deffnm}.part[0-9][0-9][0-9][0-9]{traj_suffix}'))
    return parts


def _part_number(part_fn, deffnm=DEFFNM):
    """The NNNN in a deffnm.partNNNN.<suffix> filename part_files returned."""
    stem = Path(part_fn).name
    return int(stem[len(deffnm) + len('.part'):len(deffnm) + len('.partNNNN')])


def _move_aside_stale_parts(gen_dir, resume_part, deffnm=DEFFNM,
                            traj_suffix='.xtc', prefix=ABANDONED_PART_PREFIX):
    """Move aside any part numbered past resume_part.

    The launch about to happen writes part resume_part + 1, so a higher part
    already on disk was written by a branch this checkpoint has rewound past.
    Left in place it would still match part_files' glob, and concat_parts has
    no way to tell it apart from the branch that actually continued.
    """
    for part in part_files(gen_dir, deffnm=deffnm, traj_suffix=traj_suffix):
        if _part_number(part, deffnm=deffnm) > resume_part:
            part.rename(part.with_name(prefix + part.name))


def _first_frame_step(part_p):
    """The MD step stamped on a part's first frame."""
    import mdtraj
    with mdtraj.open(str(part_p)) as handle:
        step = handle.read(1)[2]
    if not len(step):
        raise ValueError(f'{part_p} holds no frames to merge')
    return int(step[0])


def part_cutoffs(parts):
    """The step each part hands over to the next at; None for the last one.

    Where two parts cover the same steps the later one wins, so a part's frames
    are kept only while their step is below the step the next part begins at.
    """
    starts = [_first_frame_step(p) for p in parts]
    for index in range(1, len(starts)):
        if starts[index] < starts[index - 1]:
            raise ValueError(
                f'{parts[index].name} starts at step {starts[index]}, before '
                f'{parts[index - 1].name} starts at {starts[index - 1]}. Parts '
                'are merged in the order mdrun wrote them, so one that begins '
                'earlier than the part before it means a branch this run '
                'rewound past is still on the merge glob.')
    return starts[1:] + [None]


def _kept_frames(parts, chunk=CONCAT_CHUNK):
    """(xyz, time, step, box) chunks holding the frames a merge of parts keeps.

    Reads one part at a time, so a multi-GB generation never has to be resident.
    """
    import mdtraj
    for part_p, cutoff in zip(parts, part_cutoffs(parts)):
        with mdtraj.open(str(part_p)) as handle:
            while True:
                xyz, time, step, box = handle.read(chunk)
                if not len(step):
                    break
                n_kept = len(step) if cutoff is None else int(
                    np.count_nonzero(np.asarray(step) < cutoff))
                if n_kept:
                    yield (xyz[:n_kept], time[:n_kept], step[:n_kept],
                           box[:n_kept])
                if n_kept < len(step):
                    break        # steps only rise, so the rest is overlap too


def xtc_precision(traj_fn, offset=XTC_PRECISION_OFFSET,
                  uncompressed_atoms=XTC_UNCOMPRESSED_ATOMS):
    """The compressed-x-precision stamped on an .xtc's first frame, or None.

    GROMACS writes the precision into every frame, at the head of the
    compressed coordinate block. None means there is no precision to read: a
    system of nine atoms or fewer is stored raw, and so loses nothing whichever
    writer the merge picks.
    """
    with open(traj_fn, 'rb') as handle:
        head = handle.read(offset + 4)
    if len(head) < offset + 4:
        return None
    n_atoms = struct.unpack('>i', head[offset - 4:offset])[0]
    if n_atoms <= uncompressed_atoms:
        return None
    return float(struct.unpack('>f', head[offset:offset + 4])[0])


# Precisions already remarked on, so a packed job says it once, not per gen.
_warned_precisions = set()


def warn_fine_precision(precision, mdtraj_precision=MDTRAJ_XTC_PRECISION,
                        warned=_warned_precisions):
    """Print, once per precision per process, that a finer grid is being paid for."""
    if precision in warned:
        return
    warned.add(precision)
    print(f'[gmx] compressed-x-precision is {precision:g}, '
          f'{precision / mdtraj_precision:g}x finer than the '
          f'{mdtraj_precision:g} GROMACS defaults to. Are you sure about this? '
          f'A {1 / precision:g} nm grid is far below the accuracy the force '
          f'field itself has, so the extra digits store integrator noise, not '
          f'signal, and every frame of the campaign pays for them. It also '
          f'costs the merge its default writer: mdtraj cannot hold this '
          f'precision, so parts are merged with LOOS, which needs a structure '
          f'file and refuses a triclinic box. Set compressed-x-precision = '
          f'{mdtraj_precision:g} unless you truly need the finer grid.',
          flush=True)


def _refuse_triclinic(box, frame, precision, mdtraj_precision=MDTRAJ_XTC_PRECISION,
                      triclinic_rtol=TRICLINIC_RTOL):
    """Raise unless this cell is one the fine-precision writer can represent."""
    if is_orthorhombic(box, triclinic_rtol=triclinic_rtol):
        return
    raise BoxTypeError(
        f'frame {frame} has a non-orthorhombic box:\n'
        f'{np.array2string(np.asarray(box), precision=4)}\n'
        f'and the parts were written at compressed-x-precision {precision:g}. '
        f'No writer holds both: mdtraj writes an .xtc only at '
        f'{mdtraj_precision:g}, and LOOS, which does take a precision, stores '
        f'a periodic box as three numbers and would silently keep just the '
        f'diagonal of this one. Set compressed-x-precision = '
        f'{mdtraj_precision:g} for a triclinic cell.')


def _loos_model(structure_fn, n_atoms):
    """The LOOS AtomicGroup each merged frame is stamped onto before writing."""
    import loos
    if structure_fn is None:
        raise ValueError(
            'merging these parts needs the LOOS writer, which builds every '
            'frame out of a model, so concat_parts needs structure_fn. Pass '
            'the .gro the run was built from, or set compressed-x-precision to '
            f'{MDTRAJ_XTC_PRECISION:g} so mdtraj can do the merge.')
    model = loos.createSystem(str(structure_fn))
    if len(model) != n_atoms:
        raise ValueError(
            f'{structure_fn} has {len(model)} atoms but the parts hold '
            f'{n_atoms}; LOOS writes a frame out of the model, so the two have '
            'to be the same system. A compressed-x-grps subset cannot be '
            'merged at this precision.')
    return model


def _write_merged_loos(parts, out_p, structure_fn, precision, chunk=CONCAT_CHUNK,
                       angstrom_per_nm=ANGSTROM_PER_NM, writer_dt=LOOS_WRITER_DT,
                       steps_per_frame=LOOS_WRITER_STEPS_PER_FRAME):
    """Stream the kept frames into out_p through LOOS; returns how many.

    LOOS's XTCWriter takes a precision, which mdtraj's does not, so this is the
    path for parts written finer than 1e-3 nm. mdtraj still does the reading:
    LOOS's XTC reader is not wrapped for Python and cannot report the per-frame
    step and time the generation chain is spliced on.

    LOOS works in Angstroms and divides by ten again on the way out, so
    coordinates and box go in scaled up -- in float64, because scaling the
    float32 arrays mdtraj hands back would return a box a rounding step away
    from the one the source part carries.
    """
    import loos
    model, n_written = None, 0
    writer = loos.XTCWriter(str(out_p), writer_dt, steps_per_frame,
                            float(precision))
    for xyz, time, step, box in _kept_frames(parts, chunk=chunk):
        if model is None:
            model = _loos_model(structure_fn, xyz.shape[1])
        for index in range(len(step)):
            _refuse_triclinic(box[index], n_written, precision)
            model.setCoords(np.asarray(xyz[index], dtype=float) * angstrom_per_nm)
            model.periodicBox(loos.GCoord(
                *(np.diag(box[index]).astype(float) * angstrom_per_nm)))
            writer.writeFrame(model, int(step[index]), float(time[index]))
            n_written += 1
    return n_written


def _write_merged(parts, out_p, chunk=CONCAT_CHUNK):
    """Stream the kept frames of these parts into out_p; returns how many."""
    from mdtraj.formats import XTCTrajectoryFile
    n_written = 0
    with XTCTrajectoryFile(str(out_p), 'w') as handle:
        for xyz, time, step, box in _kept_frames(parts, chunk=chunk):
            handle.write(xyz, time=time, step=step, box=box)
            n_written += len(step)
    return n_written


def _check_clock(merged_p, what, got, expected, first_frame):
    """Raise unless this axis of a merged chunk came back exactly as written."""
    got, expected = np.asarray(got), np.asarray(expected)
    if np.array_equal(got, expected):
        return
    row = int(np.argwhere(got != expected)[0][0])
    raise RuntimeError(
        f'{merged_p}: the {what} of frame {first_frame + row} reads '
        f'{got[row]!r}, not the {expected[row]!r} its source part carries. '
        'The generation chain is spliced on that clock, so a merge that does '
        'not preserve it exactly is refused.')


def verify_merge(parts, merged_p, chunk=CONCAT_CHUNK, atol_nm=CONCAT_ATOL_NM):
    """Read a merge back and check every frame against the part it came from.

    Coordinates that come back moved mean the rewrite requantised them, which
    is what a writer coarser than the parts' own compressed-x-precision does.
    Returns the frames checked.
    """
    import mdtraj
    n_checked = 0
    with mdtraj.open(str(merged_p)) as handle:
        for xyz, time, step, box in _kept_frames(parts, chunk=chunk):
            got_xyz, got_time, got_step, got_box = handle.read(len(step))
            if len(got_step) != len(step):
                raise RuntimeError(
                    f'{merged_p} stops after {n_checked + len(got_step)} '
                    f'frames, but the parts it was merged from still have '
                    f'{len(step) - len(got_step)} more to contribute.')
            _check_clock(merged_p, 'MD step', got_step, step, n_checked)
            _check_clock(merged_p, 'frame time', got_time, time, n_checked)
            _check_clock(merged_p, 'box', got_box, box, n_checked)
            moved = float(np.abs(np.asarray(got_xyz)
                                 - np.asarray(xyz)).max())
            if moved > atol_nm:
                raise RuntimeError(
                    f'{merged_p}: coordinates moved by up to {moved:.2e} nm '
                    f'between the source parts and the merge, from frame '
                    f'{n_checked} on, so the rewrite requantised them or kept '
                    f'the wrong frames. The parts carry '
                    f'compressed-x-precision {xtc_precision(parts[0])}; the '
                    f'merge has to be written by a writer that holds it.')
            n_checked += len(step)
        if len(handle.read(1)[2]):
            raise RuntimeError(
                f'{merged_p} carries frames past the {n_checked} its source '
                'parts account for.')
    return n_checked


def concat_parts(gen_dir, out_fn, deffnm=DEFFNM, traj_suffix='.xtc',
                 structure_fn=None, verify=True, chunk=CONCAT_CHUNK,
                 atol_nm=CONCAT_ATOL_NM, mergeable_suffixes=MERGEABLE_SUFFIXES,
                 mdtraj_precision=MDTRAJ_XTC_PRECISION):
    """Merge a generation's parts into the one trajectory the orchestrator wants.

    Every frame keeps its own step, time and box: the generation chain is
    spliced on that clock. Where two parts cover the same steps the later part's
    frames are the ones that survive, which is only correct because callers move
    any part left behind by a rewound relaunch aside before it reaches this glob.

    mdtraj writes the merge, at its fixed precision of 1000. Parts written finer
    than that go to LOOS instead, which takes a precision but builds each frame
    out of structure_fn and cannot represent a triclinic box.

    verify reads the merge back and compares it against the parts it came from.
    """
    if traj_suffix.lower() not in mergeable_suffixes:
        raise ValueError(
            f'{traj_suffix} parts cannot be merged without losing something: '
            f'mdfarmer reads them with mdtraj, which does not carry over the '
            f'velocities and forces a .trr holds. Run with traj_suffix in '
            f'{list(mergeable_suffixes)}.')
    parts = part_files(gen_dir, deffnm=deffnm, traj_suffix=traj_suffix)
    # A launch killed before its first write leaves a part with nothing in it.
    parts = [p for p in parts if p.stat().st_size]
    if not parts:
        raise FileNotFoundError(
            f'no {deffnm}.partNNNN{traj_suffix} files in {gen_dir} to merge')
    out_p = Path(out_fn)
    # Temp name, suffix kept so the writer reads the format from it.
    tmp_p = out_p.with_name(f'{out_p.stem}.concat-tmp{out_p.suffix}')
    precision = xtc_precision(parts[0])
    fine = precision is not None and precision > mdtraj_precision
    if fine:
        warn_fine_precision(precision, mdtraj_precision=mdtraj_precision)
    try:
        if len(parts) == 1:
            # Copy, not rename, so re-running is idempotent and the part survives.
            shutil.copy(parts[0], tmp_p)
        else:
            if fine:
                _write_merged_loos(parts, tmp_p, structure_fn, precision,
                                   chunk=chunk)
            else:
                _write_merged(parts, tmp_p, chunk=chunk)
            if verify:
                verify_merge(parts, tmp_p, chunk=chunk, atol_nm=atol_nm)
    except BaseException:
        tmp_p.unlink(missing_ok=True)   # a refused merge leaves nothing behind
        raise
    tmp_p.replace(out_p)
    return out_p


def write_gen_status(gen_dir, *, target_step, reached_step, complete,
                     gen_status_name=GEN_STATUS_NAME, traj_fn=None):
    """Record how far this generation got, for the orchestrator to read."""
    status = {'target_step': int(target_step),
              'reached_step': int(reached_step),
              'complete': bool(complete)}
    if traj_fn is not None:
        status['traj'] = str(traj_fn)
    path = Path(gen_dir) / gen_status_name
    util.write_json_atomic(path, status)
    return status


def read_gen_status(gen_dir, gen_status_name=GEN_STATUS_NAME):
    """Status dict for a generation, or None if it has never reported."""
    path = Path(gen_dir) / gen_status_name
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError:
        print(f'read_gen_status: malformed {path}; ignoring.')
        return None


class MdrunFleet:
    """One stop signal shared by every mdrun in a packed job.

    A preempt or walltime warning has to reach all of them and let all of them
    checkpoint, or a replica loses work a lone job would have kept. One watcher
    passes SIGTERM on to every registered process, so they do not race.
    """

    def __init__(self, sentinel_path, poll_seconds=PREEMPT_POLL_SECONDS):
        import threading
        self.sentinel = Path(sentinel_path)
        self.poll_seconds = poll_seconds
        self._procs = {}
        self._lock = threading.Lock()
        self._stopping = threading.Event()

    def clear_sentinel(self):
        if self.sentinel.exists():
            self.sentinel.unlink()

    def register(self, key, proc):
        with self._lock:
            self._procs[key] = proc
            # A replica starting after the signal fired still has to be told.
            if self._stopping.is_set():
                proc.send_signal(signal.SIGTERM)

    def unregister(self, key):
        with self._lock:
            self._procs.pop(key, None)

    @property
    def stopping(self):
        return self._stopping.is_set()

    def poll_and_signal(self):
        """True once the sentinel has been seen and everyone has been told."""
        if self._stopping.is_set():
            return True
        if not self.sentinel.is_file():
            return False
        self._stopping.set()
        with self._lock:
            targets = list(self._procs.items())
        print(f'[gmx] preempt sentinel seen; SIGTERM -> {len(targets)} mdrun(s) '
              '(each writes a final checkpoint and stops).', flush=True)
        for key, proc in targets:
            try:
                proc.send_signal(signal.SIGTERM)
            except ProcessLookupError:
                pass                      # already exited on its own
        return True


def _mdrun_env(cmd, environ=None):
    """Environment for an mdrun, with OMP_NUM_THREADS agreeing with -ntomp.

    mdrun refuses to start when the two disagree, and a batch job inherits
    whatever the submitting shell had exported.
    """
    env = dict(os.environ if environ is None else environ)
    cmd = [str(c) for c in cmd]
    if '-ntomp' in cmd and cmd.index('-ntomp') + 1 < len(cmd):
        env['OMP_NUM_THREADS'] = cmd[cmd.index('-ntomp') + 1]
    return env


def _wait_in_fleet(proc, fleet, fleet_key):
    """Wait for a packed mdrun. The fleet owner watches and signals; a member
    only registers itself and reports whether a stop was called."""
    fleet.register(fleet_key, proc)
    try:
        rc = proc.wait()
    finally:
        fleet.unregister(fleet_key)
    if fleet.stopping:
        # mdrun exits 0 after a clean SIGTERM stop, so rc cannot tell us.
        raise Preempted(f'preempt sentinel at {fleet.sentinel}')
    return rc


def _wait_alone(proc, sentinel, poll_seconds=PREEMPT_POLL_SECONDS):
    """Wait for a solo mdrun, watching for a sentinel. None means do not watch."""
    if sentinel is None:
        return proc.wait()
    while True:
        try:
            return proc.wait(timeout=poll_seconds)
        except sp.TimeoutExpired:
            if not sentinel.is_file():
                continue
            print('[gmx] preempt sentinel seen; SIGTERM -> mdrun '
                  '(it will write a final checkpoint and stop).', flush=True)
            proc.send_signal(signal.SIGTERM)
            proc.wait()          # mdrun stops at next NS step and checkpoints
            raise Preempted(f'preempt sentinel at {sentinel}')


def _run_mdrun(cmd, cwd, handle_preempt, poll_seconds=PREEMPT_POLL_SECONDS,
               fleet=None, fleet_key=None,
               sentinel_name=PREEMPT_SENTINEL_NAME):
    """Run mdrun, passing on a preempt signal so it checkpoints before it dies."""
    cwd = Path(cwd)
    sentinel = None
    if fleet is None and handle_preempt:
        # A solo gen clears its own stale sentinel; a pack's owner clears the
        # pack's once, before any member starts.
        sentinel = cwd / sentinel_name
        if sentinel.exists():
            sentinel.unlink()
    print('[gmx mdrun]', ' '.join(map(str, cmd)), flush=True)
    proc = sp.Popen([str(c) for c in cmd], cwd=str(cwd), text=True,
                    env=_mdrun_env(cmd))
    if fleet is not None:
        rc = _wait_in_fleet(proc, fleet, fleet_key)
    else:
        rc = _wait_alone(proc, sentinel, poll_seconds=poll_seconds)
    if rc != 0:
        raise RuntimeError(f'gmx mdrun exited {rc}')


def gmx_generation(traj_dir_top_level: str,
                   top_fn: str,
                   seed_index: int,
                   clone_index: int,
                   gen_index: int,
                   title: str,
                   # gen 0: the starting .gro; gen N: the seed state.cpt.
                   seed_fn: str,
                   # Constant starting structure (.gro): grompp -c at gen 0, and
                   # the model a fine-precision part merge writes frames out of.
                   structure_fn: str = None,
                   # base .mdp; only generation 0 uses it.
                   mdp_fn: str = None,
                   # Farmer sets config['system_fn'] per seed; repurposed as the .mdp.
                   system_fn: str = None,
                   append: bool = False,  # unread; mdrun always -noappends
                   dirname_pad: int = 2,
                   sep: str = '-',
                   traj_name: str = 'prod',
                   traj_suffix: str = '.xtc',
                   restart_name: str = 'state.cpt',
                   # A resume shrinks this to the remainder still owed.
                   steps: int = 500000,
                   # Full generation length. Unread here, but recorded in
                   # config.json, which is what the harvest counts frames from.
                   steps_per_gen: int = None,
                   # The absolute step this generation ends at, counted by
                   # Clone from what the earlier generations recorded.
                   target_step: int = None,
                   # xtc stride; steps must be a whole number of these.
                   write_interval: int = 50000,
                   temperature=None,            # gen-temp for gen-0 velocities (K)
                   new_velocities: bool = False,  # True only on gen 0
                   gen_seed_base: int = 1,
                   # Overrides GEN_SEED_STRIDE; Farmer checks it > n_clones.
                   gen_seed_stride: int = GEN_SEED_STRIDE,
                   # ld-seed when set; None keeps whatever the mdp holds.
                   ld_seed: int = None,
                   maxh: float = 23.5,           # mdrun -maxh backstop
                   checkpoint_minutes: float = CHECKPOINT_MINUTES,
                   gmx_bin: str = GMX_BIN,
                   ndx_fn: str = None,
                   grompp_maxwarn: int = 2,
                   # -update cpu is MANDATORY with TIP4P-ice virtual sites.
                   mdrun_args=('-nb', 'gpu', '-bonded', 'gpu', '-pme', 'gpu',
                               '-update', 'cpu', '-pin', 'on', '-nstlist', '200'),
                   handle_preempt: bool = False,
                   deffnm: str = DEFFNM,
                   tpr_name: str = TPR_NAME,
                   seed_cpt_name: str = SEED_CPT_NAME,
                   gen_status_name: str = GEN_STATUS_NAME,
                   # Shared stop signal when generations are packed in one job.
                   fleet=None,
                   fleet_key=None,
                   # Held while the tpr is built, so K grompps do not contend.
                   grompp_lock=None,
                   **_unused):
    """Run generation gen_index to completion and return its merged trajectory.

    Raises Preempted if the scheduler asked the job to stop, and GenIncomplete
    if mdrun returned cleanly short of the step target. Both leave a checkpoint
    and a gen_status.json behind, so the next launch resumes rather than
    restarts; neither is a failure.
    """
    if target_step is None:
        raise ValueError(
            'gmx_generation needs target_step, the absolute step this '
            'generation ends at. Clone.target_step counts it from the chain; '
            'a direct caller has to work it out and pass it.')
    target_step = int(target_step)

    print('starting', title, seed_index, clone_index, gen_index, flush=True)
    gen_dir = util.dir_seeds_clones_gens(
        Path(traj_dir_top_level), seed_index, clone_index, gen_index,
        dirname_pad, sep=sep).resolve()
    traj = (gen_dir / traj_name).with_suffix(traj_suffix)
    tpr = gen_dir / tpr_name
    own_cpt = gen_dir / restart_name
    seed_cpt = gen_dir / seed_cpt_name

    # Absolute: -cpi resumes at the checkpoint's step and runs to the tpr's.
    # -cpi on another run's checkpoint, or on a .gro, is fatal: move seed aside.
    if not seed_cpt.exists() and own_cpt.exists() and not tpr.is_file():
        # No tpr yet => mdrun has not run here => restart_name is the seed.
        own_cpt.replace(seed_cpt)

    # ------------------------------- build the tpr --------------------------
    if not tpr.is_file():
        _build_gen_tpr(
            tpr=tpr, gen_dir=gen_dir, gen_index=gen_index,
            new_velocities=new_velocities, target_step=target_step,
            write_interval=write_interval, mdp_fn=mdp_fn, system_fn=system_fn,
            structure_fn=structure_fn, seed_fn=seed_fn, top_fn=top_fn,
            ndx_fn=ndx_fn,
            temperature=temperature, gen_seed_base=gen_seed_base,
            gen_seed_stride=gen_seed_stride, ld_seed=ld_seed,
            clone_index=clone_index, seed_index=seed_index,
            traj_dir_top_level=traj_dir_top_level, dirname_pad=dirname_pad,
            sep=sep, tpr_name=tpr_name, gmx_bin=gmx_bin,
            grompp_maxwarn=grompp_maxwarn, grompp_lock=grompp_lock)

    # -cpi takes our own checkpoint if mdrun has run here, else the seed.
    resume_from = None
    if is_checkpoint(own_cpt):
        resume_from = own_cpt
    elif gen_index > 0 and is_checkpoint(seed_cpt):
        resume_from = seed_cpt
    elif gen_index > 0:
        raise FileNotFoundError(
            f'generation {gen_index} has no usable checkpoint to continue from '
            f'(looked at {own_cpt} and {seed_cpt}). Its predecessor did not '
            'leave a readable state.cpt.')

    # This launch writes part resume_part + 1; hide higher, rewound-past parts.
    if resume_from is not None:
        resume_part, already = checkpoint_part_step(resume_from, gmx_bin=gmx_bin)
    else:
        resume_part, already = 0, None
    _move_aside_stale_parts(gen_dir, resume_part, deffnm=deffnm,
                            traj_suffix=traj_suffix)

    # mdrun aborts on a checkpoint at or past its nsteps, so finalise instead.
    if already is not None and already >= target_step:
        print(f'[gmx] generation {gen_index} is already at step {already} '
              f'of {target_step}; finalising without running mdrun.',
              flush=True)
        if resume_from != own_cpt:
            shutil.copy(resume_from, own_cpt)
        reached = already
    else:
        reached = None

    if reached is None:
        # -noappend: mdrun cannot append into a dir lacking its cpt's own files.
        mdrun = [gmx_bin, 'mdrun', '-s', tpr, '-deffnm', deffnm,
                 '-cpo', own_cpt, '-maxh', maxh, '-cpt', checkpoint_minutes,
                 '-noappend', *mdrun_args]
        if resume_from is not None:
            mdrun += ['-cpi', str(resume_from)]
        try:
            _run_mdrun(mdrun, gen_dir, handle_preempt,
                       fleet=fleet, fleet_key=fleet_key)
        except Preempted:
            # Record what it reached, or it is charged a restart for real work.
            if is_checkpoint(own_cpt):
                write_gen_status(
                    gen_dir, target_step=target_step, complete=False,
                    reached_step=checkpoint_step(own_cpt, gmx_bin=gmx_bin),
                    gen_status_name=gen_status_name)
            raise
        reached = checkpoint_step(own_cpt, gmx_bin=gmx_bin)

    # ------------------------------- assess ---------------------------------
    complete = reached >= target_step
    if not complete:
        write_gen_status(gen_dir, target_step=target_step, reached_step=reached,
                         complete=False, gen_status_name=gen_status_name)
        raise GenIncomplete(
            f'generation {gen_index} stopped at step {reached} of {target_step} '
            f'(mdrun -maxh, or the scheduler stopped it). It will resume from '
            f'{own_cpt} on the next launch.')

    concat_parts(gen_dir, traj, deffnm=deffnm, traj_suffix=traj_suffix,
                 structure_fn=structure_fn)
    write_gen_status(gen_dir, target_step=target_step, reached_step=reached,
                     complete=True, gen_status_name=gen_status_name,
                     traj_fn=traj)
    print('Done!', flush=True)
    return traj.resolve()


def _build_gen_tpr(*, tpr, gen_dir, gen_index, new_velocities, target_step,
                   write_interval, mdp_fn, system_fn, structure_fn, seed_fn,
                   top_fn,
                   ndx_fn, temperature, gen_seed_base, gen_seed_stride,
                   ld_seed, clone_index, seed_index,
                   traj_dir_top_level, dirname_pad, sep, tpr_name, gmx_bin,
                   grompp_maxwarn, grompp_lock=None,
                   grompp_structure_suffixes=GROMPP_STRUCTURE_SUFFIXES):
    """Build this generation's tpr, holding grompp_lock if one was given.

    Generation 0 is grompp'd from the .mdp with fresh velocities. Later ones are
    convert-tpr'd from the one before, which extends the step budget and carries
    the parameters over rather than rebuilding them.
    """
    import contextlib
    guard = grompp_lock if grompp_lock is not None else contextlib.nullcontext()
    with guard:
        if tpr.is_file():
            return tpr                    # another replica may have won the race
        if gen_index == 0 or new_velocities:
            mdp_fn = mdp_fn or system_fn
            if mdp_fn is None:
                raise ValueError(
                    'gmx_generation needs an .mdp via mdp_fn (or system_fn) to '
                    'build generation 0.')
            # Prefer seed_fn, so seeds may differ and be adaptively reseeded.
            start_fn = seed_fn if (
                seed_fn and Path(seed_fn).suffix.lower()
                in grompp_structure_suffixes) else structure_fn
            if start_fn is None:
                raise ValueError(
                    'gmx_generation needs a structure for grompp -c, as either '
                    'seed_fn or structure_fn.')
            gen_mdp = gen_dir / 'gen.mdp'
            write_gen_mdp(str(Path(mdp_fn).resolve()), str(gen_mdp),
                          nsteps=target_step,
                          nstxout_compressed=write_interval,
                          gen_vel=True, continuation=False,
                          gen_seed=(gen_seed_base
                                    + gen_seed_stride * seed_index
                                    + clone_index),
                          ld_seed=ld_seed,
                          gen_temp=temperature)
            grompp = [gmx_bin, 'grompp', '-f', gen_mdp,
                      '-c', str(Path(start_fn).resolve()),
                      '-p', str(Path(top_fn).resolve()),
                      '-o', tpr, '-po', gen_dir / 'mdout.mdp',
                      '-maxwarn', grompp_maxwarn]
            if ndx_fn:
                grompp += ['-n', str(Path(ndx_fn).resolve())]
            # Run from the topology's dir so relative #includes still resolve.
            _run(grompp, str(Path(top_fn).resolve().parent))
        else:
            prev_tpr = _previous_gen_tpr(
                Path(traj_dir_top_level), seed_index, clone_index, gen_index,
                dirname_pad, sep, tpr_name)
            _run([gmx_bin, 'convert-tpr', '-s', str(prev_tpr),
                  '-nsteps', str(target_step), '-o', str(tpr)], gen_dir)
    return tpr


def _previous_gen_tpr(top_level, seed_index, clone_index, gen_index,
                      dirname_pad, sep, tpr_name=TPR_NAME):
    prev_dir = util.dir_seeds_clones_gens(
        Path(top_level), seed_index, clone_index, gen_index - 1, dirname_pad,
        sep=sep, mkdir=False)
    prev_tpr = prev_dir / tpr_name
    if not prev_tpr.is_file():
        raise FileNotFoundError(
            f'{prev_tpr} is missing, so generation {gen_index} cannot extend '
            'its predecessor. An exact continuation needs the previous '
            "generation's tpr.")
    return prev_tpr


def gmx_config_template(clone_filled_keys=CLONE_FILLED_KEYS,
                        runtime_only_keys=RUNTIME_ONLY_KEYS, **overrides):
    """A Farmer config_template recording the whole gmx_generation call.

    Fills in placeholders for the arguments a Clone supplies per generation and
    leaves out the ones gmx_pack passes at run time, so a driver need not know
    either list. Everything else comes from overrides or the defaults.
    """
    placeholders = {key: (0 if key.endswith('_index') else '')
                    for key in clone_filled_keys}
    template = util.merge_args_defaults_dict(
        gmx_generation, **{**placeholders, **overrides})
    for key in runtime_only_keys:
        template.pop(key, None)
    return template


def gmx_basic_sim_block_json(config):
    """What the run.py in each generation directory calls.

    traj_list gains a line only when the generation finished, so it holds one
    per generation rather than one per launch. An unfinished generation exits 0:
    its work is checkpointed and will be resumed, so the job did not fail.
    """
    conf = json.loads(Path(config).read_text())
    traj_list = Path(conf.pop('traj_list'))
    try:
        new_traj = gmx_generation(**conf)
    except (Preempted, GenIncomplete) as exc:
        print(f'generation not finished: {exc}', flush=True)
        return
    with traj_list.open('a') as tl:
        tl.write(str(new_traj) + '\n')


def gmx_gen_progress(gen_path, *, total_steps, gen_index=None,
                     gen_status_name=GEN_STATUS_NAME, **_unused):
    """Steps a generation still owes, read from the status the runner wrote.

    Returns total_steps when it has never reported, which the orchestrator reads
    as nothing having run yet.
    """
    # Not counted from frames, which for GROMACS is out by one write interval.
    gen_p = Path(gen_path)
    status = read_gen_status(gen_p, gen_status_name=gen_status_name)
    if status is None:
        return total_steps
    if status.get('complete'):
        return 0
    target = status.get('target_step')
    reached = status.get('reached_step', 0)
    if target is None:
        return total_steps
    return max(0, int(target) - int(reached))


def gmx_try_recover_gen(gen_path: Path, *,
                        append_mode: bool,
                        restart_name: str,
                        traj_name: str,
                        traj_suffix: str,
                        write_interval: int,
                        total_steps: int,
                        top_fn: str,
                        gen_status_name: str = GEN_STATUS_NAME,
                        seed_cpt_name: str = SEED_CPT_NAME,
                        tpr_name: str = TPR_NAME):
    """The GROMACS version of seeder._try_recover_gen.

    Sorts a generation directory into done, partial or never started, and
    returns (gen_index, seed_fn, steps_to_run, append), or None to fall back to
    an older generation. No trajectory has to be trimmed: mdrun -cpi resumes
    from the checkpoint and concat_parts drops the overlap. The one thing it must
    never do is call a generation complete before its checkpoint reaches the
    step target, which would rewind the trajectory at the boundary.
    """
    config_p = gen_path / 'config.json'
    if not config_p.is_file():
        return None
    try:
        prev = json.loads(config_p.read_text())
    except json.JSONDecodeError:
        print(f'gmx_try_recover_gen: malformed config at {config_p}; skipping.')
        return None
    gen_index = prev['gen_index']

    status = read_gen_status(gen_path, gen_status_name=gen_status_name)
    own_cpt = gen_path / restart_name
    have_own_cpt = own_cpt.is_file() and own_cpt.stat().st_size > 0
    tpr = gen_path / tpr_name

    if status is None:
        # Never reported: nothing ran, or it died before its first checkpoint.
        seed_cpt = gen_path / seed_cpt_name
        if have_own_cpt and tpr.is_file():
            return gen_index, str(own_cpt.resolve()), total_steps, True
        if seed_cpt.is_file() and seed_cpt.stat().st_size > 0:
            return gen_index, str(seed_cpt.resolve()), total_steps, False
        return None

    if status.get('complete'):
        # Advance, seeding the next generation from this one's final checkpoint.
        if not have_own_cpt:
            print(f'gmx_try_recover_gen: {gen_path} reports complete but '
                  f'{own_cpt} is missing; cascading.')
            return None
        return gen_index + 1, str(own_cpt.resolve()), total_steps, False

    # Partial. Resume this generation from its own checkpoint.
    if not have_own_cpt:
        print(f'gmx_try_recover_gen: {gen_path} is partial but has no usable '
              f'{own_cpt}; cascading.')
        return None
    target = status.get('target_step', total_steps)
    reached = status.get('reached_step', 0)
    return gen_index, str(own_cpt.resolve()), max(0, int(target) - int(reached)), True
