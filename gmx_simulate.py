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
lets gmx trjcat put the generations back together in order.

Each launch writes its own prod.partNNNN.xtc, because mdrun -cpi will not append
into a directory that does not already hold the files its checkpoint names. The
parts are merged with trjcat when the generation finishes; where two parts cover
the same time, trjcat keeps the later file's frames, so a part left behind by a
rewound relaunch is moved aside before mdrun runs rather than left for trjcat to
prefer over the branch that actually continued.

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

# Clone copies the incoming checkpoint in under restart_name, and the runner
# moves it here at once, so restart_name only ever names a checkpoint this
# generation's own mdrun wrote.
SEED_CPT_NAME = 'seed.cpt'

# Per-generation tpr. Generation N>0 is built from generation N-1's by
# convert-tpr, so this name is also how a generation finds its predecessor.
TPR_NAME = 'prod.tpr'

# mdrun -deffnm stem; also the prefix of the .partNNNN outputs.
DEFFNM = 'prod'

# Prefix a stale part is renamed to, so it stops matching part_files' glob but
# stays on disk as a record of the abandoned branch.
ABANDONED_PART_PREFIX = 'abandoned-'

# First four bytes of every GROMACS checkpoint, big-endian 171817. Reading them
# tells a checkpoint from a .gro without asking gmx.
CHECKPOINT_MAGIC = b'\x00\x02\x9f\x29'

# How often mdrun writes a checkpoint, in minutes. GROMACS defaults to 15, which
# is how much work a hard kill can cost; a shorter period costs almost nothing.
CHECKPOINT_MINUTES = 5

# Spacing between seeds in the gen-seed sequence. gen-seed is
# base + GEN_SEED_STRIDE * seed_index + clone_index, so this must exceed
# n_clones or two seeds draw the same initial velocities.
GEN_SEED_STRIDE = 1000

# Structure formats grompp -c will read. A generation that continues another is
# seeded with a checkpoint instead, which is why the suffix has to be checked.
GROMPP_STRUCTURE_SUFFIXES = ('.gro', '.g96', '.pdb', '.brk', '.ent')

# Parameters gmx_pack injects into every gmx_generation call at runtime. A
# config template records them too, so the file's copies must be dropped before
# it is splatted, or the call gets two values for the same keyword.
RUNTIME_ONLY_KEYS = ('fleet', 'fleet_key', 'grompp_lock')

# Arguments a Clone supplies per generation. A driver builds its template
# before it knows any of them, so gmx_config_template leaves placeholders.
CLONE_FILLED_KEYS = ('seed_index', 'clone_index', 'gen_index', 'seed_fn',
                     'top_fn')

# Seconds between preempt-sentinel polls while mdrun runs.
PREEMPT_POLL_SECONDS = 5

# Default binary. Sites with an MPI-only build have gmx_mpi instead.
GMX_BIN = 'gmx'


class Preempted(Exception):
    """Raised when a preempt sentinel is seen mid-mdrun; the gen is incomplete
    and will resume via -cpi on the next launch."""
    pass


class GenIncomplete(Exception):
    """Raised when mdrun returned cleanly but the generation has not reached its
    step target, for instance because mdrun stopped itself at -maxh. It
    resumes on the next launch; this is a normal event, not a failure."""
    pass


# Default run.py body the Farmer writes into each gen dir. The batch script runs
# python run.py, which dispatches to the GROMACS block runner.
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
                  ld_seed=None, legacy_mdp_keys=LEGACY_MDP_KEYS):
    """Copy base_mdp to out_mdp, changing only the per-generation control keys."""
    # Only generation 0 needs it. Later ones inherit their parameters from the
    # previous tpr through convert-tpr, which is what keeps them exact.
    overrides = {
        'nsteps': str(int(nsteps)),
        'nstxout-compressed': str(int(nstxout_compressed)),
        'gen-vel': 'yes' if gen_vel else 'no',
        'continuation': 'yes' if continuation else 'no',
    }
    # ld-seed is set whether or not velocities are generated: it drives a
    # stochastic thermostat for the whole run, not just the start.
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


def concat_parts(gen_dir, out_fn, deffnm=DEFFNM, traj_suffix='.xtc',
                 gmx_bin=GMX_BIN):
    """Merge a generation's parts into the one trajectory the orchestrator wants.

    trjcat sorts by time and, where two parts cover the same time, keeps the
    later file's frames. That is only correct because callers move any part
    left behind by a rewound relaunch aside before it ever reaches this glob.
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
    # Written to a temp name and renamed, so nothing ever reads a half-merged
    # trajectory. The temp name keeps the suffix, since gmx reads the format
    # from the extension.
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
    """Run mdrun, passing on a preempt signal so it checkpoints before it dies.

    With a fleet, the fleet owner watches the sentinel and this only registers
    its process and reports whether a stop was signalled.
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
                   # What this generation starts from: a checkpoint when it
                   # continues another, a structure when it starts fresh. Only
                   # append is unread, since mdrun always -noappends.
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
                   # Full generation length. config['steps'] shrinks to the
                   # remainder on a resume, but the tpr's nsteps has to be the
                   # total from the start of the run.
                   steps: int = 500000,
                   steps_per_gen: int = None,
                   # xtc stride; steps must be a whole number of these or the
                   # last frame of a generation does not land on its final step.
                   write_interval: int = 50000,
                   temperature=None,            # gen-temp for gen-0 velocities (K)
                   new_velocities: bool = False,  # True only on gen 0
                   gen_seed_base: int = 1,
                   # gen-seed = base + stride * seed + clone, so two seeds
                   # cannot draw the same velocities. Farmer checks the stride
                   # is bigger than n_clones.
                   gen_seed_stride: int = GEN_SEED_STRIDE,
                   # Written to ld-seed when set, making a stochastic thermostat
                   # reproducible. None keeps whatever the mdp holds.
                   ld_seed: int = None,
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

    # Move the incoming seed aside. restart_name is also where mdrun writes
    # its own checkpoint, and -cpi given a checkpoint from another run, or the
    # .gro at generation 0, is fatal.
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
    # -noappend because mdrun will not append into a directory that does not
    # already hold the output files its checkpoint names.
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

    # The launch about to happen writes part resume_part + 1 (part 1 when
    # nothing is resumed), so any higher part already here is a branch this
    # checkpoint has rewound past. Move it aside before concat_parts can see
    # it, whether or not mdrun actually runs below.
    if resume_from is not None:
        resume_part, already = checkpoint_part_step(resume_from, gmx_bin=gmx_bin)
    else:
        resume_part, already = 0, None
    _move_aside_stale_parts(gen_dir, resume_part, deffnm=deffnm,
                            traj_suffix=traj_suffix)

    # Already finished? Finalise instead of re-running. mdrun given a
    # checkpoint at or past its nsteps aborts, so a relaunch after a lost status
    # file would otherwise turn a finished generation into a failure.
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
            # grompp -c takes this generation's own starting structure, which
            # is seed_fn whenever the generation starts fresh. Using it rather
            # than structure_fn is what lets seeds differ in topology, and lets
            # an adaptive scheme reseed from a configuration it picked.
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
    from the checkpoint and trjcat drops the overlap. The one thing it must
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
        # Never reported: nothing ran, or it died before its first
        # checkpoint. Relaunchable either way, if a checkpoint can be found.
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
