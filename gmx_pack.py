"""Run several GROMACS replicas on one GPU, in a single job, sharing it via MPS.

Two replicas per card costs about 6% of each one's sampling rate and halves the
GPU-hour bill. Slurm here offers only a plain gpu gres, with no mps or shard, so
it cannot put two jobs on one card and the packing has to happen inside one job:
one sbatch, K generations advancing together, K separate mdrun processes.

They are separate processes rather than mdrun -multidir so that one crash costs
one replica instead of all of them.

The thing that quietly wastes the benefit is CPU pinning. Two replicas that both
say a bare '-pin on' both start at core 0 and fight over the same cores, which
reads as a slow node rather than a mistake. replica_mdrun_args strips any
inherited thread and pinning flags and works out this replica's own from its
place in the pack, refusing to run if they would not fit.

A packed job has one outcome per replica. One replica failing does not stop the
others, so the tender can fail exactly one clone. Preemption is the exception
that has to reach everyone: the sentinel is watched once and SIGTERM is passed
on to every mdrun, since that handshake is what keeps trajectories contiguous.
"""

import itertools
import json
import os
import subprocess as sp
import threading
import time
from pathlib import Path

from . import gmx_simulate as gmx


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

# Ranks per replica. One GPU per member means one rank; without saying so a
# thread-MPI mdrun spreads over every core it can see, not the ones it was
# given, and the pinning offsets stop meaning anything. Only a thread-MPI build
# takes the flag at all, so gmx_supports_ntmpi asks before it is used.
# '-nt' is stripped for a related reason: it fixes total threads, so it cannot
# differ between replicas.
NTMPI = 1

# Line gmx -version prints for the MPI flavour, and the value a thread-MPI
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
    """[(cores, offset), ...], one block of neighbouring cores per replica.

    member_cores widens some members and narrows others, for a pack whose
    members stop scaling at different core counts. Left out, the split is even.
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
    """This replica's mdrun flags: the shared ones, plus its own core block.

    Raises rather than handing out more cores than the job holds, which is the
    silent slowdown this exists to prevent.
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
    """True when this GROMACS is a thread-MPI build, which is what takes -ntmpi.

    A real-MPI build makes the flag fatal, so a binary that cannot be asked
    returns False.
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
              'configuration failure, not node variance. Check that the '
              f'batch script started {mps_control_bin} and that '
              f'{mps_pipe_env} ({pipe_dir!r}) is writable.', flush=True)
        print('[pack] ' + '=' * 68, flush=True)
    return state


def write_pack_manifest(pack_dir, member_config_fns, *, cpus_per_task,
                        member_cores=None,
                        pack_manifest_name=PACK_MANIFEST_NAME):
    """Record which generation configs one packed job should advance.

    member_cores is one core count per member, in the same order; None splits
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

    # Checked once, out here: a mismatch is a mistake in this module, not one
    # replica's bad luck, and inside a member it would be reported as one.
    injected = {'fleet', 'fleet_key', 'grompp_lock'}
    if injected != set(runtime_only_keys):
        raise RuntimeError(
            f'this module injects {sorted(injected)} but RUNTIME_ONLY_KEYS '
            f'names {sorted(runtime_only_keys)}; they have to match or a '
            "config's own copies stop being dropped.")

    fleet = gmx.MdrunFleet(pack_dir / gmx.PREEMPT_SENTINEL_NAME,
                           poll_seconds=poll_seconds)
    fleet.clear_sentinel()
    grompp_lock = threading.Lock()
    results = [None] * n_replicas

    def advance(index, config_fn):
        # Published before anything can raise, and updated in place, so the
        # summary below always finds an outcome for every member.
        outcome = {'config': str(config_fn), 'replica': index,
                   'status': 'failed', 'detail': 'never reported'}
        results[index] = outcome
        try:
            conf = json.loads(Path(config_fn).read_text())
            traj_list = Path(conf.pop('traj_list'))
            # Drop the file's copies of exactly the keys being injected,
            # derived from the injected dict so the two cannot drift apart.
            runtime_kwargs = dict(fleet=fleet, fleet_key=index,
                                  grompp_lock=grompp_lock)
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
