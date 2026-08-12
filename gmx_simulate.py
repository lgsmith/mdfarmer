"""GROMACS per-generation runner for mdfarmer.

Mirrors the contract of ``simulate.omm_generation`` / ``omm_basic_sim_block_json``
so the engine-agnostic ``Farmer``/``Clone`` orchestration (priority queue,
active-clone throttling, disk-driven resume, preempt budget, bad-node registry)
drives GROMACS exactly as it drives OpenMM.

HOW GENERATIONS CHAIN
---------------------
Generations are chained the way GROMACS documents for an *exact* continuation:
``gmx convert-tpr`` extends the run's step budget and ``gmx mdrun -cpi`` resumes
from the previous generation's checkpoint. The obvious-looking alternative --
``grompp -t prev/state.cpt`` -- is **not** an exact continuation and was measured
to lose real state:

    checkpoint written by gen N     tpr built by `grompp -t` from it
    ---------------------------     --------------------------------
    nosehoover-xi   = -9.30e-02     nosehoover_xi: not available
    nosehoover-vxi  = -2.93e-01     (absent)
    pres_prev       = <nonzero>     pres_prev = 0
                                    boxv      = 0

grompp itself only claims "Reading Coordinates, Velocities and Box size from old
trajectory". With ``tcoupl = nose-hoover`` and ``pcoupl = Parrinello-Rahman`` --
deterministic, memory-bearing couplings -- that zeroes the thermostat and
barostat integrals at every generation boundary. ``-cpi`` carries them.

The chain also keeps the step and time counters **globally continuous**
(gen N ends at t=10 ps, gen N+1's frames start at t=10 ps, not t=0), which is
what makes ``gmx trjcat`` able to assemble the dataset: trjcat orders and
de-duplicates by time, so per-generation files that all restart at t=0 would be
silently mangled by it.

WHY THERE ARE ``partNNNN`` FILES
-------------------------------
``mdrun -cpi`` refuses to append into a directory that does not already contain
the exact output files named in the checkpoint, so a cross-directory resume must
use ``-noappend``, and every mdrun invocation writes its own
``prod.partNNNN.xtc``. A generation therefore accumulates one part per launch
(first run, plus one per walltime/preempt continuation). When the generation
finishes, the parts are merged with ``gmx trjcat`` into the single
``<traj_name><traj_suffix>`` the orchestrator expects. Because the parts carry
globally-correct times, trjcat removes the overlap a resumed part re-covers;
this was verified to turn parts of 11 and 1 frames into a merged 11 frames with
monotonic, uniformly-spaced times.

COMPLETION IS DECIDED BY THE CHECKPOINT, NOT BY COUNTING FRAMES
--------------------------------------------------------------
GROMACS writes an output frame at step 0 of every run, so a complete generation
holds ``nsteps/nstxout-compressed + 1`` frames, not ``nsteps/nstxout-compressed``
(verified: 11, not 10). Deciding completion by ``total_steps - frames *
write_interval`` is therefore off by one interval: a generation killed anywhere
inside its final write interval scores exactly 0 and is misread as finished,
after which the next generation seeds from a checkpoint that does not match the
end of the saved trajectory -- a silent rewind at the boundary. Instead this
runner reads the step counter out of the checkpoint and records it, with the
generation's step target, in ``GEN_STATUS_NAME``. The orchestrator reads that
file; it never has to count frames or run ``gmx`` itself.

Per-seed Farmer inputs map to GROMACS as:
    seed_structure_fns -> the .gro structure (also the gen-0 seed)
    top_fns            -> the .top topology
    system_fns         -> the .mdp run parameters   (repurposed; see ``mdp_fn``)

This module does not import OpenMM; it only shells out to the ``gmx`` binary.
"""

import json
import re
import shutil
import signal
import subprocess as sp
from pathlib import Path

from . import utilities as util


# Sentinel the batch script's SIGTERM trap touches on preempt (matches
# simulate.PREEMPT_SENTINEL_NAME so the same scheduler fstrings work).
PREEMPT_SENTINEL_NAME = 'PREEMPT_SIGTERM'

# Written by the runner after every mdrun so the orchestrator can decide whether
# a generation is finished without counting frames or invoking gmx.
GEN_STATUS_NAME = 'gen_status.json'

# The incoming seed checkpoint is copied in by Clone under `restart_name`; the
# runner immediately moves it aside to this name so `restart_name` refers only
# to a checkpoint *this* generation's mdrun wrote. Sharing one name for both is
# what made `-cpi` fatal on a generation that died before its first checkpoint.
SEED_CPT_NAME = 'seed.cpt'

# Per-generation tpr. Generation N>0 is built from generation N-1's by
# convert-tpr, so this name is also how a generation finds its predecessor.
TPR_NAME = 'prod.tpr'

# mdrun -deffnm stem; also the prefix of the .partNNNN outputs.
DEFFNM = 'prod'

# How often mdrun writes a checkpoint, in minutes. GROMACS defaults to 15, which
# is how much work a hard kill can cost; a shorter period costs almost nothing.
CHECKPOINT_MINUTES = 5

# Seconds between preempt-sentinel polls while mdrun runs.
PREEMPT_POLL_SECONDS = 5

# Default binary. Sites with an MPI-only build have `gmx_mpi` instead.
GMX_BIN = 'gmx'


class Preempted(Exception):
    """Raised when a preempt sentinel is seen mid-mdrun; the gen is incomplete
    and will resume via ``-cpi`` on the next launch."""
    pass


class GenIncomplete(Exception):
    """Raised when mdrun returned cleanly but the generation has not reached its
    step target -- e.g. mdrun stopped itself at ``-maxh``. The generation
    resumes on the next launch; this is a normal event, not a failure."""
    pass


# Default run.py body the Farmer writes into each gen dir. The batch script runs
# `python run.py`; this dispatches to the GROMACS block runner.
default_gmx_run_script = """
from mdfarmer.gmx_simulate import gmx_basic_sim_block_json as runner
runner('config.json')
"""


# .mdp keys this runner controls; everything else is inherited verbatim from the
# base .mdp. GROMACS treats '-'/'_' as equivalent and is case-insensitive.
def _norm_mdp_key(k):
    return k.strip().lower().replace('_', '-')


# Legacy spellings GROMACS removed. Silently leaving one of these in a .mdp
# means the intended setting is ignored and grompp fatals on the unknown key, so
# they are normalised to the modern name rather than left to fail at run time.
LEGACY_MDP_KEYS = {
    'nstxtcout': 'nstxout-compressed',
    'xtc-precision': 'compressed-x-precision',
    'xtc-grps': 'compressed-x-grps',
    'unconstrained-start': 'continuation',
}


def write_gen_mdp(base_mdp, out_mdp, *, nsteps, nstxout_compressed,
                  gen_vel, continuation, gen_seed=None, gen_temp=None,
                  legacy_mdp_keys=LEGACY_MDP_KEYS):
    """Copy base_mdp to out_mdp, overriding only the per-gen control keys.

    Only generation 0 needs this: later generations inherit their parameters
    from the previous generation's tpr through ``convert-tpr``, which is what
    makes the continuation exact.
    """
    overrides = {
        'nsteps': str(int(nsteps)),
        'nstxout-compressed': str(int(nstxout_compressed)),
        'gen-vel': 'yes' if gen_vel else 'no',
        'continuation': 'yes' if continuation else 'no',
    }
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
                if key in seen:
                    # A duplicate assignment later in the file would override
                    # ours; drop it rather than emit a second, conflicting line.
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


def checkpoint_step(cpt_fn, gmx_bin=GMX_BIN):
    """Step counter recorded in a GROMACS checkpoint.

    This is the authoritative measure of how far a generation actually got.
    Frame counts are not: GROMACS writes a frame at step 0, and the checkpoint
    can legitimately lag the last written frame, which is exactly the gap that
    lets a short generation masquerade as a complete one.
    """
    out = _run_capture([gmx_bin, 'dump', '-cp', str(cpt_fn)])
    match = _STEP_RE.search(out)
    if match is None:
        raise ValueError(f'no step counter in checkpoint {cpt_fn}')
    return int(match.group(1))


def is_checkpoint(cpt_fn, gmx_bin=GMX_BIN):
    """True when `cpt_fn` is a checkpoint gmx can actually read.

    Guards the case where the file at ``restart_name`` is not a checkpoint at
    all -- at generation 0 the seed is the starting ``.gro``, and handing that
    to ``-cpi`` fatals with 'Start of file magic number mismatch'.
    """
    path = Path(cpt_fn)
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        checkpoint_step(path, gmx_bin=gmx_bin)
    except (RuntimeError, ValueError):
        return False
    return True


def part_files(gen_dir, deffnm=DEFFNM, traj_suffix='.xtc'):
    """The ``prod.partNNNN.<suffix>`` files a generation has accumulated, in order."""
    gen_p = Path(gen_dir)
    parts = sorted(gen_p.glob(f'{deffnm}.part[0-9][0-9][0-9][0-9]{traj_suffix}'))
    return parts


def concat_parts(gen_dir, out_fn, deffnm=DEFFNM, traj_suffix='.xtc',
                 gmx_bin=GMX_BIN):
    """Merge a generation's ``-noappend`` parts into one trajectory.

    ``gmx trjcat`` orders by time and drops the overlap that a resumed part
    re-covers, which is correct here precisely because the continuation keeps
    times globally monotonic.
    """
    parts = part_files(gen_dir, deffnm=deffnm, traj_suffix=traj_suffix)
    if not parts:
        raise FileNotFoundError(
            f'no {deffnm}.partNNNN{traj_suffix} files in {gen_dir} to merge')
    out_p = Path(out_fn)
    if len(parts) == 1:
        # Nothing to merge; copy rather than rename so a re-run of this step is
        # idempotent and the part stays as the provenance record.
        shutil.copy(parts[0], out_p)
        return out_p
    # The temp file must keep the trajectory SUFFIX: gmx picks the output format
    # from the extension, so a name like 'prod.xtc.tmp' is rejected outright.
    # Writing to a temp and renaming keeps a reader from ever seeing a partly
    # merged trajectory at the real path.
    tmp_p = out_p.with_name(f'{out_p.stem}.trjcat-tmp{out_p.suffix}')
    _run([gmx_bin, 'trjcat', '-f', *[str(p) for p in parts], '-o', str(tmp_p)],
         gen_dir)
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
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(status, indent=2))
    tmp.replace(path)            # atomic, so a reader never sees a torn file
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
    """Shared stop-signal for a group of mdruns running in one job.

    When K replicas are packed onto one GPU under MPS, the preempt/walltime
    handshake has to reach ALL of them and wait for ALL to checkpoint -- this is
    the code path that protects trajectory contiguity, so a replica that misses
    the signal loses work that a single-replica job would have kept. One watcher
    polls the sentinel and fans SIGTERM out to every registered process, rather
    than each replica polling independently and racing.
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
            # A replica that starts after the signal already fired still has to
            # be told, or it would run on alone until walltime kills it hard.
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


def _run_mdrun(cmd, cwd, handle_preempt, poll_seconds=PREEMPT_POLL_SECONDS,
               fleet=None, fleet_key=None):
    """Run mdrun; if handle_preempt, watch for the preempt sentinel and forward
    SIGTERM to mdrun (which writes a checkpoint) before raising Preempted.

    With a `fleet`, the sentinel is watched by whoever owns the fleet and this
    call only registers its process and reports whether the stop was signalled.
    """
    cwd = Path(cwd)
    sentinel = cwd / PREEMPT_SENTINEL_NAME
    if fleet is None and handle_preempt and sentinel.exists():
        # Clear a stale sentinel from a previous preempted attempt in this gen
        # dir. With a fleet the owner does this once, before any member starts.
        sentinel.unlink()
    print('[gmx mdrun]', ' '.join(map(str, cmd)), flush=True)
    proc = sp.Popen([str(c) for c in cmd], cwd=str(cwd), text=True)
    if fleet is not None:
        fleet.register(fleet_key, proc)
    try:
        while True:
            try:
                watching = handle_preempt or fleet is not None
                rc = proc.wait(timeout=poll_seconds if watching else None)
            except sp.TimeoutExpired:
                if fleet is not None:
                    # The fleet owner signals; just keep waiting for our exit.
                    continue
                if handle_preempt and sentinel.is_file():
                    print('[gmx] preempt sentinel seen; SIGTERM -> mdrun '
                          '(it will write a final checkpoint and stop).',
                          flush=True)
                    proc.send_signal(signal.SIGTERM)
                    proc.wait()  # mdrun stops at next NS step and checkpoints
                    raise Preempted(f'preempt sentinel at {sentinel}')
                continue
            break
    finally:
        if fleet is not None:
            fleet.unregister(fleet_key)
    if fleet is not None and fleet.stopping:
        # mdrun exits 0 after a clean SIGTERM stop, so the exit code alone
        # cannot distinguish "preempted" from "finished".
        raise Preempted(f'preempt sentinel at {fleet.sentinel}')
    if rc != 0:
        raise RuntimeError(f'gmx mdrun exited {rc}')


def gmx_generation(traj_dir_top_level: str,
                   top_fn: str,
                   seed_index: int,
                   clone_index: int,
                   gen_index: int,
                   title: str,
                   # gen 0: path to the starting .gro; gen N: path to the seed
                   # state.cpt (copied into this gen dir by Clone).
                   seed_fn: str,
                   # constant starting structure (.gro) for grompp -c at gen 0.
                   structure_fn: str = None,
                   # base .mdp; only generation 0 uses it.
                   mdp_fn: str = None,
                   # Farmer sets config['system_fn'] per seed; repurposed as the .mdp.
                   system_fn: str = None,
                   append: bool = False,
                   dirname_pad: int = 2,
                   sep: str = '-',
                   traj_name: str = 'prod',
                   traj_suffix: str = '.xtc',
                   restart_name: str = 'state.cpt',
                   # Full per-generation step count. On a resume the orchestrator
                   # narrows config['steps'], so the untouched full length is
                   # carried separately -- the tpr's nsteps must be the absolute
                   # cumulative target, never a remainder.
                   steps: int = 500000,
                   steps_per_gen: int = None,
                   # xtc stride; steps must be a whole number of these or the
                   # last frame of a generation does not land on its final step.
                   write_interval: int = 50000,
                   temperature=None,            # gen-temp for gen-0 velocities (K)
                   new_velocities: bool = False,  # True only on gen 0
                   gen_seed_base: int = 1,       # gen-seed = base + clone_index
                   maxh: float = 23.5,           # mdrun -maxh backstop
                   checkpoint_minutes: float = CHECKPOINT_MINUTES,
                   gmx_bin: str = GMX_BIN,
                   ndx_fn: str = None,
                   grompp_maxwarn: int = 2,
                   # mdrun hardware flags. -update cpu is MANDATORY with TIP4P-ice
                   # virtual sites.
                   mdrun_args=('-nb', 'gpu', '-bonded', 'gpu', '-pme', 'gpu',
                               '-update', 'cpu', '-pin', 'on', '-nstlist', '200'),
                   handle_preempt: bool = False,
                   deffnm: str = DEFFNM,
                   tpr_name: str = TPR_NAME,
                   seed_cpt_name: str = SEED_CPT_NAME,
                   gen_status_name: str = GEN_STATUS_NAME,
                   # Shared stop-signal when several generations run in one
                   # job (MPS packing); None for a solo generation.
                   fleet=None,
                   fleet_key=None,
                   # Held while the tpr is built. grompp is cheap but K of them
                   # at once just contend for cores at job startup.
                   grompp_lock=None,
                   **_unused):
    steps_per_gen = int(steps_per_gen if steps_per_gen is not None else steps)
    if steps_per_gen % write_interval:
        raise ValueError(
            f'steps_per_gen={steps_per_gen} is not a whole number of '
            f'write_interval={write_interval} steps. The last frame of a '
            'generation would not land on its final step, so the frame spacing '
            'across the generation boundary would be irregular.')

    print('starting', title, seed_index, clone_index, gen_index, flush=True)
    gen_dir = util.dir_seeds_clones_gens(
        Path(traj_dir_top_level), seed_index, clone_index, gen_index,
        dirname_pad, sep=sep).resolve()
    traj = (gen_dir / traj_name).with_suffix(traj_suffix)
    tpr = gen_dir / tpr_name
    own_cpt = gen_dir / restart_name
    seed_cpt = gen_dir / seed_cpt_name

    # The cumulative step the tpr must target. Absolute, because -cpi resumes at
    # the checkpoint's absolute step and runs until the tpr's nsteps.
    target_step = (gen_index + 1) * steps_per_gen

    # --- separate the incoming seed from this generation's own checkpoint ----
    # Clone copies the previous generation's checkpoint in under `restart_name`,
    # which is also mdrun's -cpo target. Left alone, `-cpi restart_name` would
    # hand mdrun a checkpoint from another run (or, at gen 0, the .gro) and
    # fatal. Moving it aside means anything later found at `restart_name` can
    # only be a checkpoint this generation's mdrun wrote.
    if not seed_cpt.exists() and own_cpt.exists() and not tpr.is_file():
        # No tpr yet => mdrun has not run here => restart_name is the seed.
        own_cpt.replace(seed_cpt)

    # ------------------------------- build the tpr --------------------------
    if not tpr.is_file():
        _build_gen_tpr(
            tpr=tpr, gen_dir=gen_dir, gen_index=gen_index,
            new_velocities=new_velocities, target_step=target_step,
            write_interval=write_interval, mdp_fn=mdp_fn, system_fn=system_fn,
            structure_fn=structure_fn, top_fn=top_fn, ndx_fn=ndx_fn,
            temperature=temperature, gen_seed_base=gen_seed_base,
            clone_index=clone_index, seed_index=seed_index,
            traj_dir_top_level=traj_dir_top_level, dirname_pad=dirname_pad,
            sep=sep, tpr_name=tpr_name, gmx_bin=gmx_bin,
            grompp_maxwarn=grompp_maxwarn, grompp_lock=grompp_lock)

    # ------------------------------- run ------------------------------------
    # -cpi must name a checkpoint this run can legitimately continue from: our
    # own if mdrun has already run here, otherwise the seed. -noappend because
    # mdrun refuses to append into a directory lacking the checkpoint's own
    # output files, which is always the case for a fresh generation directory.
    resume_from = None
    if is_checkpoint(own_cpt, gmx_bin=gmx_bin):
        resume_from = own_cpt
    elif gen_index > 0 and is_checkpoint(seed_cpt, gmx_bin=gmx_bin):
        resume_from = seed_cpt
    elif gen_index > 0:
        raise FileNotFoundError(
            f'generation {gen_index} has no usable checkpoint to continue from '
            f'(looked at {own_cpt} and {seed_cpt}). Its predecessor did not '
            'leave a readable state.cpt.')

    # Already there? Don't re-run. mdrun handed a checkpoint at or past the
    # tpr's nsteps aborts ("the checkpoint file has already reached step N"),
    # so a redundant launch -- which happens whenever the orchestrator acts on a
    # stale or missing status file, e.g. after the job was SIGKILLed between
    # mdrun finishing and the status being written -- would otherwise turn a
    # finished generation into a hard failure. Finalising instead makes the
    # relaunch self-healing.
    if resume_from is not None:
        already = checkpoint_step(resume_from, gmx_bin=gmx_bin)
        if already >= target_step:
            print(f'[gmx] generation {gen_index} is already at step {already} '
                  f'of {target_step}; finalising without running mdrun.',
                  flush=True)
            if resume_from != own_cpt:
                shutil.copy(resume_from, own_cpt)
            reached = already
        else:
            reached = None
    else:
        reached = None

    if reached is None:
        mdrun = [gmx_bin, 'mdrun', '-s', tpr, '-deffnm', deffnm,
                 '-cpo', own_cpt, '-maxh', maxh, '-cpt', checkpoint_minutes,
                 '-noappend', *mdrun_args]
        if resume_from is not None:
            mdrun += ['-cpi', str(resume_from)]
        _run_mdrun(mdrun, gen_dir, handle_preempt,
                   fleet=fleet, fleet_key=fleet_key)
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
                 gmx_bin=gmx_bin)
    write_gen_status(gen_dir, target_step=target_step, reached_step=reached,
                     complete=True, gen_status_name=gen_status_name,
                     traj_fn=traj)
    print('Done!', flush=True)
    return traj.resolve()


def _build_gen_tpr(*, tpr, gen_dir, gen_index, new_velocities, target_step,
                   write_interval, mdp_fn, system_fn, structure_fn, top_fn,
                   ndx_fn, temperature, gen_seed_base, clone_index, seed_index,
                   traj_dir_top_level, dirname_pad, sep, tpr_name, gmx_bin,
                   grompp_maxwarn, grompp_lock=None):
    """Build this generation's tpr, holding `grompp_lock` if one was supplied.

    Generation 0 is grompp'd from the .mdp with fresh per-clone velocities.
    Every later generation is ``convert-tpr``'d from its predecessor, extending
    the cumulative step budget -- that inheritance is what makes the
    continuation exact, since the parameters (and, via ``-cpi``, the coupling
    state) are carried rather than rebuilt.
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
            if structure_fn is None:
                raise ValueError(
                    'gmx_generation needs structure_fn (the .gro for grompp -c).')
            gen_mdp = gen_dir / 'gen.mdp'
            write_gen_mdp(str(Path(mdp_fn).resolve()), str(gen_mdp),
                          nsteps=target_step,
                          nstxout_compressed=write_interval,
                          gen_vel=True, continuation=False,
                          gen_seed=gen_seed_base + clone_index,
                          gen_temp=temperature)
            grompp = [gmx_bin, 'grompp', '-f', gen_mdp,
                      '-c', str(Path(structure_fn).resolve()),
                      '-p', str(Path(top_fn).resolve()),
                      '-o', tpr, '-po', gen_dir / 'mdout.mdp',
                      '-maxwarn', grompp_maxwarn]
            if ndx_fn:
                grompp += ['-n', str(Path(ndx_fn).resolve())]
            # grompp runs from the topology's directory so a .top with relative
            # force-field includes resolves regardless of the gen-dir cwd.
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


def gmx_basic_sim_block_json(config):
    """Entry point invoked by the per-gen run.py.

    Appends to traj_list only when the generation actually finished, so the list
    holds one line per generation rather than one per mdrun invocation. An
    incomplete generation (preempt, ``-maxh``) exits 0: the work is checkpointed
    and the orchestrator resumes it, so the job is not a failure.
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
    """Steps still owed by a generation, from its recorded status.

    The orchestrator's engine-agnostic default infers this from frame counts,
    which is wrong for GROMACS by one write interval (the step-0 frame) and
    cannot see that a checkpoint lags the last written frame. Reading the step
    the runner recorded avoids both.

    Returns `total_steps` when the generation has never reported, which is what
    the orchestrator reads as "nothing ran here yet".
    """
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
    """GROMACS analog of ``seeder._try_recover_gen``.

    Classifies a generation directory as done / partial / never-started and
    returns ``(gen_index, seed_fn, steps_to_run, append)`` or None to cascade to
    an older generation.

    There is no trajectory surgery to do: ``mdrun -cpi`` resumes from the
    checkpoint and ``trjcat`` drops the overlap, so a partial generation only
    has to be pointed back at its own checkpoint. What this must get right is
    never calling a generation complete when its checkpoint has not reached the
    step target -- that is the silent-rewind case.
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
        # Never reported. Either nothing ran, or it died before finishing its
        # first mdrun. Both are relaunchable provided we can still find a
        # checkpoint to continue from -- but only one this generation wrote, or
        # the seed that was staged for it.
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
