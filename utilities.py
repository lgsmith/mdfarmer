"""Helpers shared by the orchestrator, the runners and the harvest.

Job scripts. basic_scheduler_fstrings and its variants are ready-made submit
scripts, keyed by scheduler. Anything you put in their place has to keep the
placeholders: {job_name}, which campaign_jobs below matches on, {queue_name},
{gpu_line}, and {exclude_nodes}, which expands to a directive line excluding the
nodes BadNodeRegistry has flagged and to nothing when none are. Keep the NODE:
and GPU: echoes as well, since that registry reads them out of the job's log.

Job names. A job is named title, seed, clone and gen joined by the config's sep,
so 'mycampaign-0-3-5'. The queue reports ask for the whole queue and
campaign_jobs picks ours out in Python, splitting the trailing indices off and
comparing the title that remains by string equality. Equality rather than a
prefix, because a second campaign titled 'mycampaign-long' names its jobs
'mycampaign-long-0-3-5', and a looser test would bind those ids to this
campaign's clone (0, 3, 5).

Bad nodes. A generation that aborts at 0 steps because the node is broken (a
stale CUDA driver, a GLIBC mismatch, no visible GPU) aborts the same way on
every retry, since the scheduler tends to hand the resubmission back to the same
node, and the farmer spends each clone's restart budget doing it.
BadNodeRegistry recognises those failures, remembers the host in a file that
survives a restart, and excludes it from later submissions.

Resuming. The last frame of a generation's DCD has to sit at exactly the step
its state.xml stopped at, or the next append lands at the wrong point in time
and gen-to-gen concatenation drifts. A run killed between the trajectory report
and the checkpoint report has one frame too many, so is_state_xml_usable,
state_xml_step_count, dcd_header_info and truncate_dcd_to_nframes are here to
check that and trim it.
"""

import inspect
import json
import os
import struct
from pathlib import Path
import subprocess as sp
import openmm as mm
from openmm import app


openmm_topology_readers = {
    '.top': app.GromacsTopFile,
    '.prmtop': app.AmberPrmtopFile,
    '.psf': app.CharmmPsfFile,
    '.pdb': app.PDBFile,
    '.cif': app.PDBxFile,
    '.pdbx': app.PDBxFile,
}


def read_openmm_top(top_fn):
    """The OpenMM Topology in a structure file, using the reader its suffix names.

    Only the suffix lookup is guarded: a reader raises a KeyError of its own for
    an atom type it was never given, and calling that an unsupported format
    would send the user after the wrong problem.
    """
    top_p = Path(top_fn)
    try:
        reader = openmm_topology_readers[top_p.suffix]
    except KeyError:
        raise ValueError(
            f'No topology reader for {top_p.suffix!r}. Choices are: '
            f'{", ".join(openmm_topology_readers)}') from None
    return reader(top_fn).topology


# Frame counting. mdtraj.open() measures a trajectory without reading its
# coordinates or needing a topology at all, so LOOS is only the fallback.
try:
    import mdtraj as _mdtraj
except ImportError:
    _mdtraj = None

try:
    import loos
    from loos import pyloos
except ImportError:
    loos = None
    pyloos = None


def _traj_len_mdtraj(traj_fn):
    with _mdtraj.open(str(traj_fn)) as fh:
        return len(fh)


def _traj_len_loos(traj_fn, top_fn):
    model = loos.createSystem(str(top_fn))
    return len(pyloos.Trajectory(str(traj_fn), model))


# Solute-only topology the harvester leaves beside a dry trajectory. Named here
# rather than imported from harvester, which imports this module.
DRY_TOPOLOGY_NAME = 'dry-top.pdb'


def get_traj_len(traj_fn, top_fn, dry_topology_name=DRY_TOPOLOGY_NAME):
    """Number of frames in a trajectory, or 0 if it is empty or unreadable.

    top_fn is only consulted by the LOOS fallback; the mdtraj path does not
    need it, which is what makes this work for GROMACS runs whose top_fn is a
    .top.

    The fallback tries the harvester's solute-only topology as well: after a
    harvest the trajectory name is a symlink to a stripped copy, and LOOS built
    from the wet top_fn would hit an atom-count mismatch, be swallowed by the
    broad except below, and report 0 frames, which reads as "never ran".

    With neither backend installed a trajectory that exists cannot be measured,
    so this raises rather than answering 0. Refusing here rather than at import
    keeps a scheduler-only install, which never counts a frame, working.
    """
    traj_p = Path(traj_fn)
    if not traj_p.is_file() or traj_p.stat().st_size == 0:
        return 0
    if _mdtraj is None and loos is None:
        # 0 for a trajectory that exists reads as a generation that never ran,
        # and the orchestrator deletes it.
        raise ImportError(
            f'Counting the frames in {traj_fn} needs mdtraj or LOOS, and '
            'neither is importable.')
    if _mdtraj is not None:
        try:
            return _traj_len_mdtraj(traj_p)
        except Exception as exc:
            print(f'mdtraj could not read {traj_fn}: '
                  f'{type(exc).__name__}: {exc}; trying LOOS.')
    if loos is not None:
        candidates = [top_fn] if top_fn is not None else []
        dry_top_p = traj_p.parent / dry_topology_name
        if dry_top_p.is_file():
            candidates.append(dry_top_p)
        for candidate in candidates:
            try:
                return _traj_len_loos(traj_p, candidate)
            except Exception as exc:
                # Broad: an unsupported topology is a RuntimeError, an
                # unreadable frame a LOOSError, and both have to be survivable.
                print(f'LOOS could not read {traj_fn} with topology '
                      f'{candidate}: {type(exc).__name__}: {exc}.')
        print(f'No usable topology for {traj_fn}; treating as empty.')
    return 0


def frame_timing(traj_fn, n_frames=None):
    """(step0, steps_per_frame, time0, time_per_frame), or None for a DCD.

    Neither LOOS nor mdtraj keeps a source's step and time on its own: LOOS
    numbers frames from zero at 1 ps apart, and mdtraj writes the frame index as
    the step. Anything that rewrites a trajectory has to read these and pass
    them back in.

    Given n_frames, the spacing read off frames 0 and 1 is checked against the
    last frame, since extrapolating from two frames is only right if the source
    is evenly spaced.
    """
    traj_p = Path(traj_fn)
    if traj_p.suffix.lower() != '.xtc':
        return None          # a DCD keeps its timing in the header
    with _mdtraj.open(str(traj_p)) as fh:
        available = len(fh)
        _, time, step, _ = fh.read(min(2, available))
        if available < 2:
            return int(step[0]), 0, float(time[0]), 0.0
        step0, time0 = int(step[0]), float(time[0])
        steps_per_frame = int(step[1]) - step0
        time_per_frame = float(time[1]) - time0
        if n_frames is None:
            return step0, steps_per_frame, time0, time_per_frame
        fh.seek(n_frames - 1)
        _, last_time, last_step, _ = fh.read(1)
    predicted_step = step0 + (n_frames - 1) * steps_per_frame
    predicted_time = time0 + (n_frames - 1) * time_per_frame
    if int(last_step[0]) != predicted_step or abs(
            float(last_time[0]) - predicted_time) > 1e-5 * max(
                abs(predicted_time), 1.0):
        raise ValueError(
            f'{traj_p} is not evenly spaced: frames 0 and 1 are '
            f'{steps_per_frame} steps apart, which puts frame {n_frames - 1} at '
            f'step {predicted_step}, but it is at {int(last_step[0])}. Refusing '
            'to restamp frames from an assumption the trajectory contradicts.')
    return step0, steps_per_frame, time0, time_per_frame


def strip_and_downsample(config_fn, harvester_config_fn):
    """The one old entry point, kept for the harvest.sh scripts already on disk.

    New scripts call harvester.harvest_generation, which this hands off to with
    no backend pinned, so the box shape picks one.

    hconfig keys:
     - harvester_subset: which atoms the solute trajectory keeps, in LOOS syntax
       unless harvester_subset_syntax says 'mdtraj'.
     - downsample_frq: keep every Nth frame in the solvated stream.
     - harvester_structure: structure file to build the model from. REQUIRED for
       GROMACS runs, whose top_fn is a force-field topology that neither LOOS nor
       mdtraj can build a model from. Defaults to config['top_fn'].
    """
    from . import harvester
    return harvester.harvest_generation(config_fn, harvester_config_fn)


# Ready-made job scripts, keyed by scheduler. The module docstring says what a
# replacement has to keep.
basic_scheduler_fstrings = {
    "lsf": inspect.cleandoc("""#!/bin/bash
                #BSUB -J {job_name}
                #BSUB -o lsf.out
                {gpu_line}
                #BSUB -q {queue_name}
                {exclude_nodes}

                echo "JOB_NAME: {job_name}"
                echo "NODE: $LSB_HOSTS"
                echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | paste -sd, -)"

                python {run_script_name}
                """),
    # Every #SBATCH has to stay above the echoes: Slurm stops reading
    # directives at the first real command.
    "slurm": inspect.cleandoc("""#!/bin/bash
                #SBATCH -J {job_name}
                #SBATCH -e slurm.out
                #SBATCH -o slurm.out
                {gpu_line}
                #SBATCH -p {queue_name}
                {exclude_nodes}

                echo "JOB_NAME: {job_name}"
                echo "SLURM_JOB_ID: $SLURM_JOB_ID"
                echo "NODE: $SLURMD_NODENAME"
                echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | paste -sd, -)"

                python {run_script_name}
                """)
}

# Preempt-aware variants for Farmer(handle_preempt=True): the trap touches the
# file SentinelReporter watches, python is backgrounded so bash can deliver the
# signal at all, and the sleep outlives the grace period so Slurm says CANCELLED.
basic_scheduler_fstrings_preempt = {
    "lsf": inspect.cleandoc("""#!/bin/bash
                #BSUB -J {job_name}
                #BSUB -o lsf.out
                {gpu_line}
                #BSUB -q {queue_name}
                {exclude_nodes}

                echo "JOB_NAME: {job_name}"
                echo "NODE: $LSB_HOSTS"
                echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | paste -sd, -)"

                preempt_handler() {{ touch PREEMPT_SIGTERM; sleep 70; }}
                trap preempt_handler SIGTERM

                python {run_script_name} &
                wait
                """),
    # --signal=B:TERM@120 warns the batch shell (B:) 120 seconds before the
    # allocation ends or a preemption lands, which is the time to checkpoint in.
    "slurm": inspect.cleandoc("""#!/bin/bash
                #SBATCH -J {job_name}
                #SBATCH -e slurm.out
                #SBATCH -o slurm.out
                {gpu_line}
                #SBATCH -p {queue_name}
                #SBATCH --signal=B:TERM@120
                {exclude_nodes}

                echo "JOB_NAME: {job_name}"
                echo "SLURM_JOB_ID: $SLURM_JOB_ID"
                echo "NODE: $SLURMD_NODENAME"
                echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | paste -sd, -)"

                preempt_handler() {{ touch PREEMPT_SIGTERM; sleep 70; }}
                trap preempt_handler SIGTERM

                python {run_script_name} &
                wait
                """)
}

basic_gpu_lines = {
    "lsf": "#BSUB -gpu 'num=1:j_exclusive=yes'",
    "slurm": "#SBATCH --gpus=1"
}

# One job, one GPU, several replicas sharing it through CUDA MPS. Goes with
# gmx_pack.gmx_pack_sim_block_json, which does the core pinning.
basic_scheduler_fstrings_mps = {
    "slurm": inspect.cleandoc("""#!/bin/bash
                #SBATCH -J {job_name}
                #SBATCH -e slurm.out
                #SBATCH -o slurm.out
                {gpu_line}
                #SBATCH -p {queue_name}
                #SBATCH --cpus-per-task={cpus}
                #SBATCH --signal=B:TERM@120
                {exclude_nodes}

                echo "JOB_NAME: {job_name}"
                echo "SLURM_JOB_ID: $SLURM_JOB_ID"
                echo "NODE: $SLURMD_NODENAME"
                echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | paste -sd, -)"

                # Refuse to be a second job in this pack directory. Two
                # would share one checkpoint and one set of part numbers, which
                # no checkpoint can undo. Held on fd 9 until the job ends.
                exec 9>pack.lock
                if ! flock -n 9; then
                    echo "PACK LOCK: another job already holds $(pwd)/pack.lock; exiting rather than putting a second mdrun on this checkpoint."
                    exit 0
                fi

                # Per-job MPS daemon. Keyed on the job id so two packed jobs on
                # one node never share or clobber each other's daemon.
                export CUDA_MPS_PIPE_DIRECTORY="/tmp/mps-$USER-$SLURM_JOB_ID/pipe"
                export CUDA_MPS_LOG_DIRECTORY="/tmp/mps-$USER-$SLURM_JOB_ID/log"
                mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
                if nvidia-cuda-mps-control -d; then
                    echo "MPS: daemon started ($CUDA_MPS_PIPE_DIRECTORY)"
                else
                    echo "MPS: WARNING daemon FAILED to start; replicas will time-slice the GPU at 10-20% worse throughput."
                fi

                stop_mps() {{
                    echo quit | nvidia-cuda-mps-control 2>/dev/null || true
                    rm -rf "/tmp/mps-$USER-$SLURM_JOB_ID"
                }}
                preempt_handler() {{ touch PREEMPT_SIGTERM; sleep 70; }}
                trap preempt_handler SIGTERM
                trap stop_mps EXIT

                python {run_script_name} &
                wait
                """)
}

# One 'jobid jobname' line per job in this user's whole queue. Which of them
# belong to a campaign is decided by campaign_jobs below, not by the shell.
basic_scheduler_reports = {
    # -o 'JOBID JOB_NAME', not the default, which pads out several columns.
    "lsf": "bjobs -o 'JOBID JOB_NAME' -noheader",
    # '%i %j' prints the id and the untruncated name; -O Name would truncate to
    # 8 characters and silently stop matching.
    "slurm": "squeue --me -h -o '%i %j'"
}

# The name Farmer's association argument takes. Both reports read the same
# lines now, and differ only in what the orchestrator takes out of them.
basic_scheduler_assoc_reports = dict(basic_scheduler_reports)

# A job name ends in this many integers: its seed, clone and gen index.
JOB_NAME_INDEX_COUNT = 3


def parse_scheduler_report(text):
    """The (job_id, job_name) pairs in a queue report, one per line.

    The id is the first whitespace-separated field and the name is the rest of
    the line, so a name holding a space survives. A line carrying no name is
    warned about, since a report without names can never match; a line whose id
    is not an integer is dropped quietly, being an array task or a header, and
    so never one of ours.
    """
    jobs = []
    for line in text.splitlines():
        fields = line.split(None, 1)
        if not fields:
            continue
        if len(fields) < 2:
            print(f'WARNING: scheduler report line {line!r} carries a job id '
                  'and no job name, so no clone can be bound to it.')
        elif fields[0].isascii() and fields[0].isdigit():
            jobs.append((int(fields[0]), fields[1].strip()))
    return jobs


def split_job_name(name, sep='-', index_count=JOB_NAME_INDEX_COUNT):
    """A job name split into (title, indices), or None.

    None when the name does not end in index_count integers. The title leads
    and may itself contain sep, so only the last index_count fields are read as
    indices and everything before them is rejoined as the title.
    """
    fields = name.split(sep)
    if len(fields) <= index_count:
        return None
    tail = fields[-index_count:]
    if not all(f.isascii() and f.isdigit() for f in tail):
        return None
    return sep.join(fields[:-index_count]), tuple(int(f) for f in tail)


def campaign_jobs(text, title, sep='-', index_count=JOB_NAME_INDEX_COUNT):
    """This campaign's jobs in a queue report, as (job_id, indices) pairs.

    Ordered as the report listed them. The title is compared by equality
    against the title each name was split into, never as a pattern and never as
    a bare prefix, so a campaign called 'sampling' does not claim the jobs of
    one called 'sampling-long'. A name that reads as ours but carries the wrong
    indices is warned about and left out: nothing can be bound to it, so a
    second job may land on top of it.
    """
    ours = []
    for jid, name in parse_scheduler_report(text):
        split = split_job_name(name, sep=sep, index_count=index_count)
        if split is None:
            if name == title or name.startswith(title + sep):
                print(f'WARNING: queued job {jid} is named {name!r}, which '
                      f'does not end in {index_count} {sep!r}-separated '
                      'indices. No clone will be bound to it, and one may '
                      'launch a second job on top of it.')
            continue
        if split[0] == title:
            ours.append((jid, split[1]))
    return ours


def slurm_was_preempted(jid):
    """Whether Slurm's accounting says this finished job was preempted.

    False whenever the question cannot be answered (no sacct, a timeout, a gap
    in the accounting), so a real failure still costs the clone a restart.
    """
    try:
        out = sp.check_output(
            ['sacct', '-j', str(jid), '-n', '-o', 'State', '-X'],
            text=True, timeout=30
        ).strip()
    except (sp.CalledProcessError, sp.TimeoutExpired, FileNotFoundError):
        return False
    for line in out.splitlines():
        if 'PREEMPTED' in line:
            return True
    return False


def lsf_was_preempted(jid):
    """Whether LSF recorded TERM_PREEMPT for this finished job. False on any
    error, as above."""
    try:
        out = sp.check_output(
            ['bjobs', '-d', '-o', 'exit_reason', '-noheader', str(jid)],
            text=True, timeout=30
        ).strip()
    except (sp.CalledProcessError, sp.TimeoutExpired, FileNotFoundError):
        return False
    return 'TERM_PREEMPT' in out


# Clone.check_start_gen uses these so a preemption is not charged against the
# generation's restart_attempts budget.
preemption_checkers = {
    'sbatch': slurm_was_preempted,
    'bsub': lsf_was_preempted,
}


# Log strings that mean the node is broken rather than the simulation.
default_bad_node_patterns = (
    'CUDA_ERROR_UNSUPPORTED_PTX_VERSION',
    'CUDA_ERROR_NO_DEVICE',
    'CUDA_ERROR_INVALID_DEVICE',
    'CUDA_ERROR_NOT_INITIALIZED',
    'CUDA driver version is insufficient',
    'No CUDA-capable device is detected',
    'Failed to initialize NVML',
    # GLIBC ABI mismatch; only the quoted half of the message is stable.
    "version `GLIBC_",
)

# Seconds a scheduler query may take before it counts as failed.
SCHEDULER_QUERY_TIMEOUT = 120

# Some schedulers say "nothing matched" with a non-zero exit rather than with
# empty output. LSF's bjobs does; an empty Slurm queue exits 0.
empty_query_messages = ('no unfinished job found', 'no matching job found',
                        'is not found')


def scheduler_query(command, timeout=SCHEDULER_QUERY_TIMEOUT,
                    empty_messages=empty_query_messages):
    """Ask the scheduler something. Returns (trusted, text).

    trusted is False when the query itself failed, which is not at all the same
    as the queue being empty: reading a squeue that died as "no jobs are
    running" relaunches every live clone on top of itself. The reports above
    are single commands, so their own exit status settles that, but pipefail
    stays for the pipelines a site may substitute, where the last stage would
    otherwise exit 0 over a dead first stage. A scheduler that reports an empty
    queue by exiting non-zero is recognised by what it says.
    """
    try:
        result = sp.run(f'set -o pipefail; {command}', shell=True,
                        executable='/bin/bash', text=True,
                        capture_output=True, timeout=timeout)
    except (sp.TimeoutExpired, OSError) as exc:
        print(f'WARNING: scheduler query {command!r} did not run: {exc}')
        return False, ''
    output = result.stdout.strip()
    if result.returncode == 0:
        return True, output
    said = (result.stderr or '').lower()
    if not output and any(m in said for m in empty_messages):
        return True, ''
    print(f'WARNING: scheduler query {command!r} exited {result.returncode}: '
          f'{(result.stderr or "").strip()[:200]}')
    return False, output


default_scheduler_log_names = {
    'sbatch': 'slurm.out',
    'slurm': 'slurm.out',
    'bsub': 'lsf.out',
    'lsf': 'lsf.out',
}

# Submit command -> the family the fstring tables are keyed by.
scheduler_families = {
    'sbatch': 'slurm',
    'slurm': 'slurm',
    'bsub': 'lsf',
    'lsf': 'lsf',
}


def _format_exclude_slurm(nodes):
    if not nodes:
        return ''
    return f'#SBATCH --exclude={",".join(sorted(nodes))}'


def _format_exclude_lsf(nodes):
    if not nodes:
        return ''
    selectors = ' && '.join(f"hname!='{n}'" for n in sorted(nodes))
    return f'#BSUB -R "select[{selectors}]"'


default_exclude_node_formatters = {
    'sbatch': _format_exclude_slurm,
    'slurm': _format_exclude_slurm,
    'bsub': _format_exclude_lsf,
    'lsf': _format_exclude_lsf,
}


class BadNodeRegistry:
    """Tracks nodes that produced node-local failures and excludes them
    from subsequent submissions.

    On boot: parses the persistence file (default bad_nodes.txt) and
    rebuilds the in-memory exclude set so a farmer restart doesn't
    re-learn the same bad nodes. Seeds scheduler_kws['exclude_nodes']
    with the corresponding directive line.

    Per detection: scan_and_record(gen_dir, clone_tag) reads the gen's
    scheduler log (slurm.out / lsf.out), looks for a bad_node_patterns
    hit, harvests the NODE: line (and GPU: if present), appends a breadcrumb
    row, and refreshes scheduler_kws['exclude_nodes'].

    Breadcrumb file is plain text, tab-separated. Comment-out (prefix
    with #) or delete rows to clear entries; the farmer rereads the file
    on boot.
    """

    BREADCRUMB_HEADER = (
        '# mdfarmer bad-nodes blocklist\n'
        '# This file is appended to whenever a clone aborts at 0 steps\n'
        "# on a node whose log matches a known fatal-on-this-node pattern.\n"
        '# The farmer excludes these nodes on subsequent submissions via\n'
        "# the '{exclude_nodes}' placeholder in scheduler_fstring.\n"
        '#\n'
        '# To clear an entry: comment it out, or delete the row, then\n'
        '# restart the farmer. Comment and blank lines are ignored.\n'
        '# are ignored on reload.\n'
        '#\n'
        "# If many nodes from one partition fail with the same pattern,\n"
        "# that's a hint about how to reconfigure: e.g. PTX-version errors\n"
        '# usually mean the partition has older driver/CUDA-toolkit nodes,\n'
        "# so narrowing your --constraint (Slurm) or queue is more durable\n"
        "# than relying on this exclude list to grow.\n"
        '#\n'
        '# Columns (tab-separated):\n'
        '#   timestamp\tnode\tgpu\tpattern\tlog_path\tclone_tag\n'
    )

    def __init__(self, persist_path, scheduler, scheduler_kws,
                 patterns=None, log_name=None, exclude_formatter=None):
        self.persist_path = Path(persist_path)
        self.scheduler = scheduler
        # Held by reference; mutating exclude_nodes here updates the dict
        # the Farmer hands to every Clone for str.format() at launch.
        self.scheduler_kws = scheduler_kws
        self.patterns = tuple(patterns) if patterns is not None \
            else default_bad_node_patterns
        self.log_name = log_name or default_scheduler_log_names.get(
            scheduler, 'slurm.out')
        self.exclude_formatter = (
            exclude_formatter
            or default_exclude_node_formatters.get(scheduler)
            or (lambda nodes: '')
        )
        self.bad_nodes = self._load_persisted()
        self._refresh_kws()

    def _load_persisted(self):
        nodes = set()
        if not self.persist_path.is_file():
            return nodes
        for raw in self.persist_path.read_text().splitlines():
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            # Older and hand-edited rows are tolerated by taking the first
            # field that looks like a hostname rather than column 2.
            for cand in parts:
                cand = cand.strip()
                if not cand:
                    continue
                # ISO timestamps start with a 4-digit year + '-'.
                if len(cand) >= 5 and cand[4] == '-' and cand[:4].isdigit():
                    continue
                nodes.add(cand)
                break
        return nodes

    def _refresh_kws(self):
        self.scheduler_kws['exclude_nodes'] = self.exclude_formatter(self.bad_nodes)

    def _ensure_header(self):
        if (not self.persist_path.is_file()
                or self.persist_path.stat().st_size == 0):
            self.persist_path.write_text(self.BREADCRUMB_HEADER)

    def _append_row(self, node, gpu, pattern, log_path, clone_tag):
        import datetime as _dt
        ts = _dt.datetime.now(_dt.timezone.utc).isoformat(timespec='seconds')
        row = f'{ts}\t{node}\t{gpu}\t{pattern}\t{log_path}\t{clone_tag}\n'
        with self.persist_path.open('a') as f:
            f.write(row)

    @staticmethod
    def _extract_value(text, key):
        for line in text.splitlines():
            if line.startswith(key):
                rest = line[len(key):].strip()
                return rest if rest else None
        return None

    def scan_and_record(self, gen_dir, clone_tag):
        """Blocklist the node behind a generation's failure, and name it.

        The generation's scheduler log is searched for a bad_node_patterns hit;
        on one, the NODE: line it echoed is recorded and excluded from here on.
        None when nothing matched, or when the script echoed no NODE: line.
        """
        log_p = Path(gen_dir) / self.log_name
        if not log_p.is_file():
            return None
        try:
            text = log_p.read_text(errors='replace')
        except OSError as exc:
            print(f'BadNodeRegistry: cannot read {log_p}: {exc}')
            return None
        matched = next((p for p in self.patterns if p in text), None)
        if matched is None:
            return None
        node = self._extract_value(text, 'NODE:')
        if node is None:
            print(f'BadNodeRegistry: pattern {matched!r} matched in {log_p} '
                  "but no 'NODE:' line found in scheduler log; cannot "
                  'exclude. Add an echo of NODE: $SLURMD_NODENAME (or LSF '
                  'equivalent) to your scheduler_fstring.')
            return None
        gpu = self._extract_value(text, 'GPU:') or 'unknown'
        self._ensure_header()
        self._append_row(node, gpu, matched, log_p.resolve(), clone_tag)
        if node not in self.bad_nodes:
            print(f'BadNodeRegistry: node {node!r} (gpu={gpu!r}) hit '
                  f'fatal-on-node pattern {matched!r}; adding to exclude '
                  f'list (logged to {self.persist_path.resolve()}).')
            self.bad_nodes.add(node)
            self._refresh_kws()
        else:
            print(f'BadNodeRegistry: node {node!r} (already excluded) hit '
                  f'pattern {matched!r} again; logged in {self.persist_path}.')
        return node


# Indent for the JSON files a runner and the orchestrator share.
JSON_INDENT = 2


def write_json_atomic(path, obj, indent=JSON_INDENT):
    """Write obj to path as JSON under a temp name, and return path.

    The rename is atomic and the temp name sits in the target's own directory,
    so a reader racing the write sees either the whole old file or the whole new
    one. Every file a running job and the orchestrator share is written this
    way: a generation's config and status, a pack manifest, the seed map.

    indent is a parameter rather than fixed because config.json has always been
    written at 4 and the rest at 2, and reshaping files already on disk would
    make a stored config differ from itself on the next boot.
    """
    path = Path(path)
    tmp_p = path.with_name(path.name + '.tmp')
    tmp_p.write_text(json.dumps(obj, indent=indent))
    tmp_p.replace(path)
    return path


# Where a campaign records what each seed index means.
SEED_MAP_NAME = 'seed_map.json'


def read_seed_map(seed_map_p: Path):
    """The recorded {seed index: label}, or {} when nothing is recorded yet."""
    seed_map_p = Path(seed_map_p)
    if not seed_map_p.is_file():
        return {}
    return {int(index): label
            for index, label in json.loads(seed_map_p.read_text()).items()}


def write_seed_map(seed_map_p: Path, seed_map: dict):
    """Write {seed index: label}, under a temp name so no boot reads a torn file."""
    seed_map_p = Path(seed_map_p)
    seed_map_p.parent.mkdir(parents=True, exist_ok=True)
    return write_json_atomic(
        seed_map_p,
        {str(index): label for index, label in sorted(seed_map.items())})


def changed_seed_labels(recorded: dict, seed_labels):
    """Recorded indices whose label has changed, as {index: (was, now)}."""
    return {index: (recorded[index], label)
            for index, label in enumerate(seed_labels)
            if index in recorded and recorded[index] != label}


def check_seed_map(seed_map_p: Path, seed_labels):
    """Bind each seed index to its label, refusing a boot that re-indexes one.

    A seed index is an on-disk identity: it names the seed directory, appears
    in job names, and is the key running jobs are re-associated by. Every label
    already recorded must still mean the same thing; indices past the end of
    the record are new seeds and are added to it. Returns the merged mapping.
    """
    seed_labels = list(seed_labels)
    duplicated = sorted({label for label in seed_labels
                         if seed_labels.count(label) > 1})
    if duplicated:
        raise ValueError(
            f'seed_labels repeats {duplicated}. A label is a seed\'s identity, '
            'so two seeds sharing one make a swap between them undetectable.')
    recorded = read_seed_map(seed_map_p)
    changed = changed_seed_labels(recorded, seed_labels)
    if changed:
        detail = '; '.join(f'seed {index} was {was!r} and is now {now!r}'
                           for index, (was, now) in sorted(changed.items()))
        raise ValueError(
            f'seed_labels disagrees with {seed_map_p}: {detail}. That index '
            'already names a directory of finished data and any job still '
            'queued for it, so booting would run this seed into another '
            "replica's trajectory. Put the seed lists back in their recorded "
            'order, or, if the re-index is deliberate, move the existing '
            f'data aside and delete {seed_map_p}.')
    merged = {**recorded, **dict(enumerate(seed_labels))}
    write_seed_map(seed_map_p, merged)
    return merged


def select_platform(platform_name=None, platform_properties=None):
    """(platform, properties) for a simulation to run on.

    A named platform is demanded and raises if it will not load. Otherwise the
    fastest platform that is not Reference is taken, and having only Reference
    raises too, since MD on it is a hang as far as the scheduler can tell.

    Properties the chosen platform does not expose are dropped, so a request
    like {'Precision': 'mixed'} can be left in a config that also runs on CPU.
    """
    if platform_name is not None:
        platform = mm.Platform.getPlatformByName(platform_name)
    else:
        candidates = []
        for i in range(mm.Platform.getNumPlatforms()):
            p = mm.Platform.getPlatform(i)
            if p.getName() != 'Reference':
                candidates.append(p)
        if not candidates:
            raise RuntimeError(
                'No non-Reference OpenMM platform available; refusing to launch '
                'a simulation that would silently hang on Reference.'
            )
        platform = max(candidates, key=lambda p: p.getSpeed())
    if platform_properties is None:
        filtered_properties = None
    else:
        supported = set(platform.getPropertyNames())
        filtered_properties = {}
        ignored = []
        for k, v in platform_properties.items():
            if k in supported:
                filtered_properties[k] = v
            else:
                ignored.append(k)
        if ignored:
            print(
                f'Platform {platform.getName()} does not support these '
                f'properties; ignoring them: {ignored}'
            )
        if not filtered_properties:
            filtered_properties = None
    return platform, filtered_properties


def is_state_xml_usable(p: Path) -> bool:
    """Whether this state.xml can be deserialized, so a gen can resume from it."""
    if not p.is_file() or p.stat().st_size == 0:
        return False
    try:
        mm.XmlSerializer.deserialize(p.read_text())
    except Exception as exc:
        print(f'is_state_xml_usable: deserialize failed for {p}: {exc}')
        return False
    return True


def state_xml_step_count(p: Path) -> int:
    """The step a state.xml stopped at. ValueError if it cannot be read."""
    import xml.etree.ElementTree as ET
    try:
        root = ET.parse(p).getroot()
    except ET.ParseError as exc:
        # ParseError is a SyntaxError, which the caller's ValueError guard
        # would miss, dropping the clone instead of cascading to an older gen.
        raise ValueError(f'state.xml at {p} is not parseable XML, so the step '
                         f'it stopped at cannot be read: {exc}') from exc
    sc = root.attrib.get('stepCount')
    if sc is None:
        raise ValueError(f'state.xml at {p} has no stepCount attribute '
                         '(pre-OpenMM-8 writeState?)')
    return int(sc)


def dcd_header_info(p: Path) -> dict:
    """Parse a DCD header. Returns nset, istart, nsavc, with_unitcell,
    n_atoms, and header_size (file offset where the first frame begins).
    """
    with open(p, 'rb') as f:
        bs1 = struct.unpack('<i', f.read(4))[0]
        if bs1 != 84:
            raise ValueError(f'DCD block-1 size {bs1} != 84 at {p}')
        magic = f.read(4)
        if magic != b'CORD':
            raise ValueError(f'DCD magic {magic!r} != b"CORD" at {p}')
        ints = struct.unpack('<20i', f.read(80))
        nset, istart, nsavc = ints[0], ints[1], ints[2]
        # ints[10], at byte 48, is 1 when frames carry the 6-double box record.
        with_unitcell = ints[10]
        be1 = struct.unpack('<i', f.read(4))[0]
        if be1 != 84:
            raise ValueError(f'DCD block-1 end marker {be1} != 84 at {p}')
        # title block
        bs2 = struct.unpack('<i', f.read(4))[0]
        f.read(bs2)
        be2 = struct.unpack('<i', f.read(4))[0]
        if be2 != bs2:
            raise ValueError(f'DCD title block markers disagree at {p}')
        # natoms block: always 4-byte payload
        bs3 = struct.unpack('<i', f.read(4))[0]
        if bs3 != 4:
            raise ValueError(f'DCD natoms block size {bs3} != 4 at {p}')
        n_atoms = struct.unpack('<i', f.read(4))[0]
        be3 = struct.unpack('<i', f.read(4))[0]
        if be3 != 4:
            raise ValueError(f'DCD natoms block end marker {be3} != 4 at {p}')
        header_size = f.tell()
    return {'nset': nset, 'istart': istart, 'nsavc': nsavc,
            'with_unitcell': bool(with_unitcell), 'n_atoms': n_atoms,
            'header_size': header_size}


def dcd_frame_size(with_unitcell: bool, n_atoms: int) -> int:
    """Bytes one DCD frame occupies: a 56-byte box record if the file has them,
    plus three Fortran coordinate records of 8 + 4*n_atoms each."""
    return (56 if with_unitcell else 0) + 24 + 12 * n_atoms


def truncate_dcd_to_nframes(p: Path, target_nframes: int) -> int:
    """Cut a DCD down to at most target_nframes, or however many whole
    frames its bytes actually hold if that is fewer. Rewrites nset and
    truncates trailing bytes so the header and the file length always
    agree afterward. Never grows the file. Idempotent.

    Returns the frame count actually achieved, which the caller must
    compare against target_nframes: they can disagree when the header's
    nset over- or under-states what is really on disk.
    """
    info = dcd_header_info(p)
    frame_size = dcd_frame_size(info['with_unitcell'], info['n_atoms'])
    header_size = info['header_size']
    # nset is not trusted: OpenMM bumps it before writing the frame it counts,
    # so a kill mid-write leaves it ahead of the data. Bytes are what happened.
    whole_frames = max(0, (p.stat().st_size - header_size) // frame_size)
    achievable = min(target_nframes, whole_frames)
    if achievable != target_nframes:
        print(f'{p}: asked to trim to {target_nframes} frames but only '
              f'{whole_frames} whole frames are actually on disk; '
              f'trimming to {achievable} instead.')
    with open(p, 'r+b') as f:
        f.seek(8)
        f.write(struct.pack('<i', achievable))
    os.truncate(str(p), header_size + achievable * frame_size)
    return achievable


def check_whole_frames(total_steps, write_interval, source='config_template'):
    """Refuse a step count that is not a whole number of write_intervals.

    The leftover steps write no frame and no checkpoint, so the run passes its
    last report and never registers as finished: calx_remaining_steps keeps
    asking for the remainder and every relaunch spends it again.

    source names where the numbers came from, since the Farmer checks a config
    template and a Clone checks the generation it is about to launch. Either
    value being absent or zero means there is nothing to check. Returns
    total_steps, so a caller can validate in place.
    """
    if not total_steps or not write_interval or total_steps % write_interval == 0:
        return total_steps
    raise ValueError(
        f'{source} steps={total_steps} is not a whole number of '
        f'write_interval={write_interval} steps. The remaining '
        f'{total_steps % write_interval} would write no frame and no '
        'checkpoint, so the generation would never finish.')


def calx_remaining_steps(traj_fn, top_fn, total_steps, write_interval):
    """Steps a generation still owes, from the frames already on disk.

    A negative answer means more frames than the generation should hold, usually
    duplicates from appending against a stale state.xml or a config whose
    write_interval or total_steps has changed since the run started. The caller
    treats that as complete, so it is warned about rather than passed silently.
    """
    traj_len = get_traj_len(traj_fn, top_fn)
    remaining = total_steps - traj_len * write_interval
    if remaining < 0:
        print(f'WARNING: calx_remaining_steps({traj_fn}) = {remaining}; '
              f'traj has {traj_len} frames at write_interval={write_interval} '
              f'but total_steps={total_steps}. Treating as complete; check '
              f'for over-appended trajectory.')
    return remaining


# Keys config.json carries for the run block rather than for a runner: the
# sim-block wrappers pop traj_list back out before calling one.
CONFIG_ONLY_KEYS = ('traj_list',)


def merge_args_defaults_dict(function, config_only_keys=CONFIG_ONLY_KEYS,
                             **kwargs):
    """A config dict recording the full call: every parameter and its value.

    Only as jsonizable as the values put in it.

    Two things stay out of the result, since both would carry the sentinel
    inspect._empty, a class, which json.dumps cannot write:

    * **kwargs-style catch-alls (gmx_generation has **_unused),
      which have no default because they collect leftovers;
    * parameters with no default that the caller did not supply. Those are
      required arguments, and are named in a TypeError instead.

    A keyword that is neither a parameter of the function nor one of
    config_only_keys is a typo: it is refused, not written into the config,
    where gmx_generation's **_unused would swallow it at run time and leave
    config.json describing a call that never happened.
    """
    sig = inspect.signature(function)
    variadic = (inspect.Parameter.VAR_KEYWORD, inspect.Parameter.VAR_POSITIONAL)
    config = {name: param.default
              for name, param in sig.parameters.items()
              if param.kind not in variadic}
    unknown = sorted(set(kwargs) - set(sig.parameters) - set(config_only_keys))
    if unknown:
        raise TypeError(
            f'{function.__name__} has no parameter {unknown}; check for a '
            f'typo. Only arguments {function.__name__} takes, plus '
            f'{list(config_only_keys)}, belong here, so config.json records '
            f'the call that really happens.')
    config.update(kwargs)
    missing = sorted(name for name, value in config.items()
                     if value is inspect.Parameter.empty)
    if missing:
        raise TypeError(
            f'{function.__name__} has no default for {missing}, and none was '
            f'supplied. Pass them here so config.json records the whole call.')
    return config


def fdir(dirpre, num, pad, sep='-', padchar='0'):
    return f'{dirpre}{sep}{num:{padchar}>{pad}}'


def dir_seeds_clones(top_lvl: Path, seed_index, clone_index, pad, sep='-',
                     padchar='0', mkdir=True):
    p = top_lvl / fdir('seed', seed_index, pad, sep=sep, padchar=padchar) / \
        fdir('clone', clone_index, pad, sep=sep, padchar=padchar)
    if mkdir:
        p.mkdir(exist_ok=True, parents=True)
    return p


# What every generation directory calls its run record. That file plus the job
# script is all a rerun of the generation needs.
CONFIG_NAME = 'config.json'


def earlier_gen_configs(top_lvl, seed_index, clone_index, gen_index, pad,
                        sep='-', config_name=CONFIG_NAME):
    """The config of every generation before this one, oldest first.

    Each generation records what it actually ran, so the ones before it are
    counted rather than assumed to match it. A seed may be given a different
    generation length between boots, and a wallclock-matched run ends its
    generations wherever the clock ran out.
    """
    for earlier in range(gen_index):
        gen_dir = dir_seeds_clones_gens(Path(top_lvl), seed_index, clone_index,
                                        earlier, pad, sep=sep, mkdir=False)
        config_p = gen_dir / config_name
        if not config_p.is_file():
            raise FileNotFoundError(
                f'{config_p} is missing, so there is no record of what '
                f'generation {earlier} ran. Every later generation is placed by '
                'counting the ones before it, and assuming they match this one '
                'would misplace every step and frame from here on. Restore that '
                "file, or write one recording that generation's steps_per_gen "
                'and write_interval.')
        yield json.loads(config_p.read_text())


def steps_before(top_lvl, seed_index, clone_index, gen_index, pad, sep='-',
                 config_name=CONFIG_NAME):
    """The absolute step this generation starts from."""
    return sum(c['steps_per_gen'] for c in earlier_gen_configs(
        top_lvl, seed_index, clone_index, gen_index, pad, sep=sep,
        config_name=config_name))


def dir_seeds_clones_gens(top_lvl: Path, seed_index, clone_index, gen_index, pad,
                          sep='-', padchar='0', mkdir=True):
    p = top_lvl / fdir('seed', seed_index, pad, sep=sep, padchar=padchar) / \
        fdir('clone', clone_index, pad, sep=sep, padchar=padchar) / \
        fdir('gen', gen_index, pad, sep=sep, padchar=padchar)
    if mkdir:
        p.mkdir(exist_ok=True, parents=True)
    return p


default_steps = int(2.5e7)  # 100 ns at a 0.004 ps timestep.
default_state_data_kwargs = dict(
    totalSteps=default_steps,
    step=True,
    speed=True,
    progress=True,
    potentialEnergy=True,
    temperature=True,
    separator=' '
)

# you should change these to match your own setup.
default_straight_sampling_config_template = dict(
    traj_dir_top_level='straight-sampling',
    append=True,
    integrator_xml='integrator.xml',
    dirname_pad=2,
    sep='-',
    traj_name='traj',
    traj_suffix='.xtc',
    restart_name='state.xml',
    # None takes the fastest platform that is not Reference; name one to force it.
    platform_name=None,
    # Filtered against what the platform supports, so CPU just drops this.
    platform_properties={'Precision': 'mixed'},
    steps=default_steps,
    state_data_kwargs=default_state_data_kwargs,
    eq_steps=None,
    write_interval=2500,
    minimize_first=False,
    temperature=300,
    new_velocities=False
)

# Replace all of these; they are here to make the slots obvious.
default_straight_sampling_init_config = dict(
    title='samplingX',  # This you should def overwrite for your own jobs!
    seeds=[
        # Expect that len(seeds) == len(top_fns) == len(system_fns)
        'state.xml'
    ],
    top_fns=[
        "my_system.pdb"
    ],
    system_fns=[
        'system.xml'
    ]
)


# Writes two trajectories, one dried and one downsampled but still solvated.
# harvest_generation picks its own backend, and a re-run of one does nothing.
default_harvest_shellscript = inspect.cleandoc("""#!/bin/bash
                #BSUB -J harvest
                #BSUB -o harvest.out
                #BSUB -q {queue_name}

                python -c 'from mdfarmer import harvest_generation; harvest_generation("config.json", "hconfig.json")'
                """)

default_harvest_shellscript_slurm = inspect.cleandoc("""#!/bin/bash
                #SBATCH -J harvest
                #SBATCH -o harvest.out
                #SBATCH -p {queue_name}
                #SBATCH --time={harvest_time}
                #SBATCH --cpus-per-task=1

                python -c 'from mdfarmer import harvest_generation; harvest_generation("config.json", "hconfig.json")'
                """)
