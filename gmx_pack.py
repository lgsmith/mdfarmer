"""Pack K GROMACS replicas onto one GPU in a single Slurm job, under MPS.

WHY THIS EXISTS
---------------
Running two replicas per card costs about 6% of the per-replica sampling rate
and halves the GPU-hour bill. The zero-code way to get that would be to let
Slurm co-schedule two independent jobs onto one GPU, but this cluster's Slurm
has ``GresTypes = gpu`` only -- no ``mps``, no ``shard`` -- so it cannot. The
packing therefore has to happen *inside* one job, which is what this module
does: one sbatch, K generations advancing together, K separate ``gmx mdrun``
processes sharing the card through the MPS daemon.

K separate processes, not ``mdrun -multidir``: GROMACS could run the ensemble
itself, but that makes the K trajectories a single failure domain, and one crash
taking out every replica is the wrong trade when contiguity is the point.

THE PART THAT SILENTLY EATS THE BENEFIT
---------------------------------------
CPU pinning. Two replicas in one job that both say a bare ``-pin on`` each pin
starting at core 0 and fight over the same cores. You keep the packing and lose
most of the retention, and it looks like node variance rather than a
misconfiguration. `replica_mdrun_args` therefore *strips* any inherited
``-ntomp``/``-pin*`` flags and derives this replica's own from its index within
the pack -- and refuses to run rather than overcommit the allocation.

FAILURE IS PER-REPLICA
----------------------
A packed job's outcome is K outcomes. One replica raising must not abort the
others, so every member runs to completion and the exceptions are collected;
only then is the per-member result reported, which lets the tender fail exactly
one clone. Preempt/walltime is the exception that must reach *everyone*: the
sentinel is watched once and SIGTERM is fanned out to all K mdruns (see
``gmx_simulate.MdrunFleet``), because that handshake is what protects
trajectory contiguity.
"""

import itertools
import json
import os
import subprocess as sp
import threading
import time
from pathlib import Path

from . import gmx_simulate as gmx


# Replicas per GPU. 2 is the measured sweet spot on Blackwell.
REPS_PER_CARD = 2

# mdrun pinning stride. 1 keeps a replica's threads on consecutive cores.
PIN_STRIDE = 1

# Manifest a packed job's run.py reads: which generation directories to advance.
PACK_MANIFEST_NAME = 'pack.json'

# Per-member outcome record written beside the manifest.
PACK_STATUS_NAME = 'pack_status.json'

# mdrun flags that must be derived per replica, never inherited from the shared
# template. Each maps to the number of VALUES that follow it.
PER_REPLICA_MDRUN_FLAGS = {'-ntomp': 1, '-pin': 1, '-pinoffset': 1,
                           '-pinstride': 1, '-ntmpi': 1, '-nt': 1}

# Thread-MPI ranks per replica. Without it a thread-MPI mdrun picks its own
# decomposition from every core it can SEE rather than the ones this replica was
# given -- several ranks of -ntomp threads each, and pinning offsets that then
# describe nothing. One GPU per member means one rank.
#
# Only a thread-MPI build accepts the flag; a real-MPI build takes its rank
# count from mpirun and makes it fatal. Such a build already gets one rank per
# replica, since the pack launches mdrun directly rather than under mpirun, so
# `gmx_supports_ntmpi` decides once per job whether to emit it at all.
#
# '-nt' is stripped for a related reason: it fixes TOTAL threads, so it cannot
# vary per replica and goes fatal the moment two members differ in size.
NTMPI = 1

# Line `gmx -version` prints for the MPI flavour, and the value a thread-MPI
# build reports there.
GMX_MPI_VERSION_KEY = 'MPI library:'
GMX_THREAD_MPI_VALUE = 'thread_mpi'

# Named by gmx_simulate, which owns the parameters.
RUNTIME_ONLY_KEYS = gmx.RUNTIME_ONLY_KEYS

# Seconds between preempt-sentinel polls while the pack runs.
POLL_SECONDS = 5

# Environment variables the MPS daemon is configured through.
MPS_PIPE_ENV = 'CUDA_MPS_PIPE_DIRECTORY'
MPS_LOG_ENV = 'CUDA_MPS_LOG_DIRECTORY'

MPS_CONTROL_BIN = 'nvidia-cuda-mps-control'


default_gmx_pack_run_script = """
from mdfarmer.gmx_pack import gmx_pack_sim_block_json as runner
runner('pack.json')
"""


def member_core_layout(cpus_per_task, n_replicas, member_cores=None):
    """``[(cores, offset), ...]``: one contiguous, non-overlapping core block
    per replica.

    `member_cores` sets each member's width independently, for a pack whose
    members have different core knees. Absent, the allocation is split evenly,
    which is what a pack of one condition wants.
    """
    cpus_per_task = int(cpus_per_task)
    if n_replicas < 1:
        raise ValueError(f'n_replicas must be >= 1, got {n_replicas}')
    if member_cores is None:
        cores = [cpus_per_task // n_replicas] * n_replicas
    else:
        cores = [int(c) for c in member_cores]
        if len(cores) != n_replicas:
            raise ValueError(
                f'member_cores has {len(cores)} entries for {n_replicas} '
                'replicas')
    if min(cores) < 1:
        raise ValueError(
            f'{cpus_per_task} cpus cannot be split as {cores} across '
            f'{n_replicas} replicas; ask for at least {n_replicas} '
            'cpus-per-task.')
    if sum(cores) > cpus_per_task:
        raise ValueError(
            f'{cores} sums to {sum(cores)} threads, over the {cpus_per_task} '
            'cpus this job holds')
    offsets = list(itertools.accumulate(cores, initial=0))[:-1]
    return list(zip(cores, offsets))


def replica_mdrun_args(base_args, replica_index, n_replicas, cpus_per_task,
                       pin_stride=PIN_STRIDE, ntmpi=NTMPI, member_cores=None,
                       per_replica_flags=PER_REPLICA_MDRUN_FLAGS):
    """This replica's mdrun flags: shared template minus pinning, plus its own.

    ``-ntomp <cores> -pin on -pinoffset <offset> -pinstride 1`` gives replica i
    a private, contiguous block of cores. Raises rather than overcommitting: K
    replicas each taking more cores than the allocation holds is exactly the
    silent slowdown this function exists to prevent.
    """
    if not 0 <= replica_index < n_replicas:
        raise ValueError(
            f'replica_index {replica_index} out of range for {n_replicas} replicas')
    ntomp, offset = member_core_layout(
        cpus_per_task, n_replicas, member_cores=member_cores)[replica_index]

    stripped = []
    args = list(base_args or ())
    i = 0
    while i < len(args):
        token = str(args[i])
        if token in per_replica_flags:
            i += 1 + per_replica_flags[token]
            continue
        stripped.append(token)
        i += 1
    # ntmpi=None means "this build cannot take the flag"; see NTMPI above.
    rank_args = [] if ntmpi is None else ['-ntmpi', str(ntmpi)]
    return stripped + rank_args + [
        '-ntomp', str(ntomp), '-pin', 'on',
        '-pinoffset', str(offset),
        '-pinstride', str(pin_stride)]


def gmx_supports_ntmpi(gmx_bin=gmx.GMX_BIN,
                       mpi_version_key=GMX_MPI_VERSION_KEY,
                       thread_mpi_value=GMX_THREAD_MPI_VALUE, timeout=60):
    """True when this GROMACS is a thread-MPI build, so -ntmpi is legal.

    `gmx -version` reports either ``MPI library: thread_mpi`` or ``MPI library:
    MPI (...)``. Only the first accepts -ntmpi; the second makes it fatal, so an
    unprobeable binary returns False.
    """
    try:
        result = sp.run([gmx_bin, '-version'], capture_output=True, text=True,
                        timeout=timeout)
    except (FileNotFoundError, sp.TimeoutExpired, OSError) as exc:
        print(f'[pack] could not probe {gmx_bin} for its MPI flavour ({exc}); '
              'not emitting -ntmpi.', flush=True)
        return False
    for line in (result.stdout + result.stderr).splitlines():
        if line.strip().startswith(mpi_version_key):
            return thread_mpi_value in line.lower()
    return False


def mps_is_running(mps_control_bin=MPS_CONTROL_BIN, timeout=10):
    """True when an MPS control daemon answers on this job's pipe directory.

    Worth checking rather than assuming: if the daemon failed to start the
    mdruns still run, just time-sliced at maybe 10-20% worse throughput. That
    has to be loud, or a week gets spent blaming the nodes.
    """
    try:
        result = sp.run([mps_control_bin], input='get_server_list\n',
                        text=True, capture_output=True, timeout=timeout)
    except (FileNotFoundError, sp.TimeoutExpired, OSError):
        return False
    return result.returncode == 0


def report_mps_state(mps_control_bin=MPS_CONTROL_BIN,
                     mps_pipe_env=MPS_PIPE_ENV, mps_log_env=MPS_LOG_ENV):
    """Log whether MPS is actually in effect. Returns the state dict."""
    pipe_dir = os.environ.get(mps_pipe_env)
    log_dir = os.environ.get(mps_log_env)
    running = mps_is_running(mps_control_bin=mps_control_bin)
    state = {'mps_running': running, 'pipe_dir': pipe_dir, 'log_dir': log_dir}
    if running:
        print(f'[pack] MPS daemon reachable (pipe dir {pipe_dir}).', flush=True)
    else:
        print('[pack] ' + '=' * 68, flush=True)
        print('[pack] WARNING: no MPS daemon is reachable. The replicas in this '
              'job will still run, but they will TIME-SLICE the GPU instead of '
              'sharing it, at roughly 10-20% worse throughput. This is a '
              'configuration failure, not node variance -- check that the '
              f'batch script started {mps_control_bin} and that '
              f'{mps_pipe_env} ({pipe_dir!r}) is writable.', flush=True)
        print('[pack] ' + '=' * 68, flush=True)
    return state


def write_pack_manifest(pack_dir, member_config_fns, *, cpus_per_task,
                        reps_per_card=REPS_PER_CARD, member_cores=None,
                        pack_manifest_name=PACK_MANIFEST_NAME):
    """Record which generation configs one packed job should advance.

    `member_cores` is one core count per member, in the same order; None splits
    the allocation evenly. Validated here so a bad split fails at submission
    rather than inside the job.
    """
    members = [str(Path(p).resolve()) for p in member_config_fns]
    if member_cores is not None:
        member_cores = [int(c) for c in member_cores]
        member_core_layout(cpus_per_task, len(members),
                           member_cores=member_cores)
    manifest = {
        'members': members,
        'cpus_per_task': int(cpus_per_task),
        'reps_per_card': int(reps_per_card),
        'member_cores': member_cores,
    }
    path = Path(pack_dir) / pack_manifest_name
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(manifest, indent=2))
    tmp.replace(path)
    return path


def _resolve_cpus_per_task(manifest):
    """Cores this job actually got, preferring what Slurm reports."""
    for var in ('SLURM_CPUS_PER_TASK', 'SLURM_CPUS_ON_NODE'):
        value = os.environ.get(var)
        if value and value.isdigit():
            return int(value)
    cpus = manifest.get('cpus_per_task')
    if cpus:
        return int(cpus)
    return os.cpu_count() or 1


def gmx_pack_sim_block_json(manifest_fn=PACK_MANIFEST_NAME,
                            poll_seconds=POLL_SECONDS,
                            reps_per_card=REPS_PER_CARD,
                            pin_stride=PIN_STRIDE,
                            pack_status_name=PACK_STATUS_NAME,
                            ntmpi=NTMPI,
                            runtime_only_keys=RUNTIME_ONLY_KEYS,
                            gmx_bin=None):
    """Entry point for a packed job's run.py: advance every member concurrently.

    Threads rather than processes: the work is all subprocess waiting, and one
    process keeps signal handling and the preempt sentinel in a single place.
    """
    manifest_p = Path(manifest_fn).resolve()
    manifest = json.loads(manifest_p.read_text())
    members = manifest['members']
    n_replicas = len(members)
    cpus_per_task = _resolve_cpus_per_task(manifest)
    pack_dir = manifest_p.parent

    # Probe the binary once, not once per replica: -ntmpi is fatal on a
    # real-MPI build and necessary on a thread-MPI one.
    if gmx_bin is None:
        gmx_bin = json.loads(Path(members[0]).read_text()).get(
            'gmx_bin', gmx.GMX_BIN)
    if ntmpi is not None and not gmx_supports_ntmpi(gmx_bin):
        print(f'[pack] {gmx_bin} is a real-MPI build; omitting -ntmpi (it gets '
              'one rank per replica from being launched without mpirun).',
              flush=True)
        ntmpi = None

    member_cores = manifest.get('member_cores')
    layout = member_core_layout(cpus_per_task, n_replicas,
                                member_cores=member_cores)
    print(f'[pack] {n_replicas} replicas, {cpus_per_task} cpus, '
          f'cores(offset) ' + ' '.join(f'{c}({o})' for c, o in layout),
          flush=True)
    mps_state = report_mps_state()

    fleet = gmx.MdrunFleet(pack_dir / gmx.PREEMPT_SENTINEL_NAME,
                           poll_seconds=poll_seconds)
    fleet.clear_sentinel()
    grompp_lock = threading.Lock()
    results = [None] * n_replicas

    def advance(index, config_fn):
        outcome = {'config': str(config_fn), 'replica': index}
        try:
            conf = json.loads(Path(config_fn).read_text())
            traj_list = Path(conf.pop('traj_list'))
            # Drop the file's copies of exactly the keys being injected,
            # derived from the injected dict so the two cannot drift apart.
            runtime_kwargs = dict(fleet=fleet, fleet_key=index,
                                  grompp_lock=grompp_lock)
            if set(runtime_kwargs) != set(runtime_only_keys):
                raise RuntimeError(
                    f'runtime kwargs {sorted(runtime_kwargs)} no longer match '
                    f'RUNTIME_ONLY_KEYS {sorted(runtime_only_keys)}; update the '
                    'constant so config.json copies keep getting dropped.')
            for key in runtime_kwargs:
                conf.pop(key, None)
            conf['mdrun_args'] = replica_mdrun_args(
                conf.get('mdrun_args'), index, n_replicas, cpus_per_task,
                pin_stride=pin_stride, ntmpi=ntmpi,
                member_cores=member_cores)
            traj = gmx.gmx_generation(**runtime_kwargs, **conf)
        except gmx.Preempted as exc:
            outcome.update(status='preempted', detail=str(exc))
        except gmx.GenIncomplete as exc:
            outcome.update(status='incomplete', detail=str(exc))
        except Exception as exc:
            # One replica's failure must not take the others down with it.
            outcome.update(status='failed',
                           detail=f'{type(exc).__name__}: {exc}')
            import traceback
            traceback.print_exc()
        else:
            outcome.update(status='complete', traj=str(traj))
            with traj_list.open('a') as tl:
                tl.write(str(traj) + '\n')
        results[index] = outcome

    threads = [threading.Thread(target=advance, args=(i, fn), daemon=False,
                                name=f'replica-{i}')
               for i, fn in enumerate(members)]
    for t in threads:
        t.start()

    # One watcher for the whole pack. Each replica waits only on its own mdrun.
    while any(t.is_alive() for t in threads):
        fleet.poll_and_signal()
        time.sleep(poll_seconds)
    for t in threads:
        t.join()

    status = {'members': results, 'cpus_per_task': cpus_per_task,
              'n_replicas': n_replicas, **mps_state}
    status_p = pack_dir / pack_status_name
    tmp = status_p.with_name(status_p.name + '.tmp')
    tmp.write_text(json.dumps(status, indent=2))
    tmp.replace(status_p)

    for outcome in results:
        print(f'[pack] replica {outcome["replica"]}: {outcome["status"]}'
              + (f' -- {outcome.get("detail", "")[:120]}'
                 if outcome['status'] != 'complete' else ''), flush=True)
    return status
