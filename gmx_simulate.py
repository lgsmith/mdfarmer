"""GROMACS per-generation runner for mdfarmer.

Mirrors the contract of ``simulate.omm_generation`` / ``omm_basic_sim_block_json``
so the engine-agnostic ``Farmer``/``Clone`` orchestration (priority queue,
active-clone throttling, disk-driven resume, preempt budget, bad-node registry)
drives GROMACS exactly as it drives OpenMM. The differences from the OpenMM
runner are concentrated here:

  * State is carried between generations with GROMACS files, not an OpenMM State
    XML. The per-gen "checkpoint" is ``state.cpt`` (set ``restart_name='state.cpt'``
    in the config template); gen N seeds from gen N-1's final ``state.cpt`` via
    ``grompp -t`` (positions + velocities + box), so velocities flow across the
    gen boundary and the concatenated per-gen trajectories are one physically
    contiguous run. (Step/time counters reset per gen — analyze by frame order.)
  * Generation 0 of each clone starts from the SAME structure (.gro) with FRESH
    Maxwell-Boltzmann velocities and a DISTINCT ``gen-seed`` per clone
    (``gen_seed_base + clone_index``), so replicas diverge immediately. Later
    gens continue (``gen-vel=no, continuation=yes``).
  * Recovery leans on ``gmx mdrun -cpi`` (which truncates the trajectory back to
    the checkpoint on resume), so the OpenMM DCD-surgery in ``_try_recover_gen``
    is unnecessary — ``gmx_try_recover_gen`` only classifies a gen as
    done / partial / never-started from the xtc frame count + checkpoint.
  * Preemption: ``gmx mdrun`` checkpoints natively on SIGTERM. With
    ``handle_preempt=True`` (paired with ``utilities.basic_scheduler_fstrings_preempt``,
    whose bash trap touches ``PREEMPT_SIGTERM``), this runner polls for that
    sentinel and forwards SIGTERM to mdrun, which writes a final checkpoint and
    exits; the gen resumes via ``-cpi`` on the next launch.

Per-seed Farmer inputs map to GROMACS as:
    seed_structure_fns -> the .gro structure (also the gen-0 seed)
    top_fns            -> the .top topology
    system_fns         -> the .mdp run parameters   (repurposed; see ``mdp_fn``)
and the constant structure for ``grompp -c`` on later gens is taken from
``structure_fn`` in the config template (one structure per Farmer driver, i.e.
n_seeds == 1; for multi-seed GROMACS drivers, pass a per-gen structure instead).

This module does not import OpenMM; it only shells out to the ``gmx`` binary.
"""

import json
import signal
import subprocess as sp
from pathlib import Path

from . import utilities as util


# Sentinel the batch script's SIGTERM trap touches on preempt (matches
# simulate.PREEMPT_SENTINEL_NAME so the same scheduler fstrings work).
PREEMPT_SENTINEL_NAME = 'PREEMPT_SIGTERM'


class Preempted(Exception):
    """Raised when a preempt sentinel is seen mid-mdrun; the gen is incomplete
    and will resume via append/-cpi on the next launch."""
    pass


# Default run.py body the Farmer writes into each gen dir. The batch script runs
# `python run.py`; this dispatches to the GROMACS block runner.
default_gmx_run_script = """
from mdfarmer.gmx_simulate import gmx_basic_sim_block_json as runner
runner('config.json')
"""


# .mdp keys this runner controls per generation; everything else is inherited
# verbatim from the base .mdp (rcoulomb, rvdw, tcoupl, ref-t, constraints, ...).
# GROMACS treats '-'/'_' as equivalent and is case-insensitive in mdp keys.
def _norm_mdp_key(k):
    return k.strip().lower().replace('_', '-')


def write_gen_mdp(base_mdp, out_mdp, *, nsteps, nstxout_compressed,
                  gen_vel, continuation, gen_seed=None, gen_temp=None):
    """Copy base_mdp to out_mdp, overriding only the per-gen control keys."""
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
            if key in targets:
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


def _run_mdrun(cmd, cwd, handle_preempt, poll_seconds=5):
    """Run mdrun; if handle_preempt, watch for the preempt sentinel and forward
    SIGTERM to mdrun (which writes a checkpoint) before raising Preempted."""
    cwd = Path(cwd)
    sentinel = cwd / PREEMPT_SENTINEL_NAME
    # Clear a stale sentinel from a previous preempted attempt in this gen dir.
    if handle_preempt and sentinel.exists():
        sentinel.unlink()
    print('[gmx mdrun]', ' '.join(map(str, cmd)), flush=True)
    proc = sp.Popen([str(c) for c in cmd], cwd=str(cwd), text=True)
    while True:
        try:
            rc = proc.wait(timeout=poll_seconds if handle_preempt else None)
        except sp.TimeoutExpired:
            if handle_preempt and sentinel.is_file():
                print('[gmx] preempt sentinel seen; SIGTERM -> mdrun '
                      '(it will write a final checkpoint and stop).', flush=True)
                proc.send_signal(signal.SIGTERM)
                proc.wait()  # mdrun stops at next step and checkpoints
                raise Preempted(f'preempt sentinel at {sentinel}')
            continue
        break
    if rc != 0:
        raise RuntimeError(f'gmx mdrun exited {rc}')


def gmx_generation(traj_dir_top_level: str,
                   top_fn: str,
                   seed_index: int,
                   clone_index: int,
                   gen_index: int,
                   title: str,
                   # gen 0: path to the starting .gro; gen N: path to the seed
                   # state.cpt (copied into this gen dir by check_copy_set_restart_seed).
                   seed_fn: str,
                   # constant starting structure (.gro) used for grompp -c on every
                   # gen (the -t checkpoint overrides its coords/vels on gen N).
                   structure_fn: str,
                   # base .mdp; per-gen copies override nsteps/gen-vel/seed/continuation.
                   mdp_fn: str = None,
                   # Farmer sets config['system_fn'] per seed; we repurpose it as the .mdp.
                   system_fn: str = None,
                   append: bool = False,
                   dirname_pad: int = 2,
                   sep: str = '-',
                   traj_name: str = 'prod',
                   traj_suffix: str = '.xtc',
                   restart_name: str = 'state.cpt',
                   # full per-gen step count (Clone passes remaining on a resume,
                   # but on append we don't re-grompp, so only the fresh value matters).
                   steps: int = 500000,
                   # xtc/checkpoint stride; must match nstxout-compressed so
                   # calx_remaining_steps (frames * write_interval) is exact.
                   write_interval: int = 50000,
                   temperature=None,            # gen-temp for gen-0 velocities (K)
                   new_velocities: bool = False,  # True only on gen 0
                   gen_seed_base: int = 1,       # gen-seed = base + clone_index
                   maxh: float = 23.5,           # mdrun -maxh backstop before walltime
                   gmx_bin: str = 'gmx',
                   ndx_fn: str = None,
                   grompp_maxwarn: int = 2,
                   # mdrun hardware flags. -update cpu is MANDATORY with TIP4P-ice
                   # virtual sites; baked in here so it can't be misconfigured per-gen.
                   mdrun_args=('-nb', 'gpu', '-bonded', 'gpu', '-pme', 'gpu',
                               '-update', 'cpu', '-pin', 'on', '-nstlist', '200'),
                   handle_preempt: bool = False,
                   **_unused):
    mdp_fn = mdp_fn or system_fn
    if mdp_fn is None:
        raise ValueError('gmx_generation needs an .mdp via mdp_fn (or system_fn).')
    if structure_fn is None:
        raise ValueError('gmx_generation needs structure_fn (the .gro for grompp -c).')

    print('starting', title, seed_index, clone_index, gen_index, flush=True)
    gen_dir = util.dir_seeds_clones_gens(
        Path(traj_dir_top_level), seed_index, clone_index, gen_index,
        dirname_pad, sep=sep).resolve()
    # Absolute paths everywhere: grompp runs from the topology's directory (so a
    # .top with relative FF includes like `#include "./amber03w.ff/..."` resolves
    # regardless of the gen-dir cwd), while mdrun runs in the gen dir.
    top_fn = str(Path(top_fn).resolve())
    top_dir = str(Path(top_fn).parent)
    structure_fn = str(Path(structure_fn).resolve())
    mdp_fn = str(Path(mdp_fn).resolve())
    tpr = gen_dir / 'prod.tpr'
    traj = (gen_dir / traj_name).with_suffix(traj_suffix)
    cpt = gen_dir / restart_name

    # Build the .tpr unless we're resuming an existing one (append).
    if not (append and tpr.is_file()):
        gen_mdp = gen_dir / 'gen.mdp'
        if new_velocities:           # generation 0: fresh, per-clone velocities
            write_gen_mdp(mdp_fn, str(gen_mdp), nsteps=steps,
                          nstxout_compressed=write_interval,
                          gen_vel=True, continuation=False,
                          gen_seed=gen_seed_base + clone_index,
                          gen_temp=temperature)
            grompp = [gmx_bin, 'grompp', '-f', gen_mdp, '-c', structure_fn,
                      '-p', top_fn, '-o', tpr, '-po', gen_dir / 'mdout.mdp',
                      '-maxwarn', grompp_maxwarn]
        else:                        # generation N>0: continue from seed cpt
            write_gen_mdp(mdp_fn, str(gen_mdp), nsteps=steps,
                          nstxout_compressed=write_interval,
                          gen_vel=False, continuation=True)
            grompp = [gmx_bin, 'grompp', '-f', gen_mdp, '-c', structure_fn,
                      '-t', str(Path(seed_fn).resolve()), '-p', top_fn,
                      '-o', tpr, '-po', gen_dir / 'mdout.mdp',
                      '-maxwarn', grompp_maxwarn]
        if ndx_fn:
            grompp += ['-n', str(Path(ndx_fn).resolve())]
        _run(grompp, top_dir)

    mdrun = [gmx_bin, 'mdrun', '-s', tpr, '-deffnm', 'prod',
             '-x', traj, '-cpo', cpt, '-maxh', maxh, *mdrun_args]
    if append and cpt.is_file():
        # resume this gen's partial run; -cpi truncates the traj back to the cpt.
        mdrun += ['-cpi', cpt, '-append']
    _run_mdrun(mdrun, gen_dir, handle_preempt)

    print('Done!', flush=True)
    return traj.resolve()


def gmx_basic_sim_block_json(config):
    """Entry point invoked by the per-gen run.py. Mirrors
    simulate.omm_basic_sim_block_json: run the gen, then append the trajectory
    path to traj_list — unless preempted (incomplete), in which case exit 0 so
    Slurm records the job CANCELLED and the orchestrator resumes via append."""
    conf = json.loads(Path(config).read_text())
    traj_list = Path(conf.pop('traj_list'))
    try:
        new_traj = gmx_generation(**conf)
    except Preempted:
        return
    with traj_list.open('a') as tl:
        tl.write(str(new_traj) + '\n')


def gmx_try_recover_gen(gen_path: Path, *,
                        append_mode: bool,
                        restart_name: str,
                        traj_name: str,
                        traj_suffix: str,
                        write_interval: int,
                        total_steps: int,
                        top_fn: str):
    """GROMACS analog of seeder._try_recover_gen. Classifies a gen dir as
    done / partial / never-started from the xtc frame count + checkpoint, and
    returns (gen_index, seed_fn, steps_to_run, append) or None to cascade.

    No DCD/state surgery: ``gmx mdrun -cpi`` truncates the trajectory back to the
    checkpoint on resume, so we only need to decide what the next launch does.
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

    cpt = gen_path / restart_name
    have_cpt = cpt.is_file() and cpt.stat().st_size > 0
    traj = (gen_path / traj_name).with_suffix(traj_suffix)
    have_traj = traj.is_file() and traj.stat().st_size > 0

    if not have_traj:
        # Nothing ran yet. If a usable seed checkpoint is present, (re)launch
        # this gen from it; otherwise cascade to an older gen / the initial seed.
        if have_cpt:
            return gen_index, str(cpt.resolve()), total_steps, False
        return None

    if not have_cpt:
        # Trajectory but no checkpoint to resume from -> unrecoverable; cascade.
        print(f'gmx_try_recover_gen: {traj} present but no usable {cpt}; cascading.')
        return None

    remaining = util.calx_remaining_steps(
        str(traj), top_fn, total_steps, write_interval)
    seed_fn = str(cpt.resolve())
    if remaining > 0:
        if append_mode:
            return gen_index, seed_fn, remaining, True
        return None  # non-append: leave partial, cascade to older gen
    # Gen complete; advance, seeding the next gen from this gen's final cpt.
    return gen_index + 1, seed_fn, total_steps, False
