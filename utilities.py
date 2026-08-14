import inspect
import os
import struct
from pathlib import Path
import json
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
    try:
        top_p = Path(top_fn)
        top_ext = top_p.suffix
        topology = openmm_topology_readers[top_ext](top_fn).topology
    except KeyError:
        print('You seem to have used a topology format', top_ext,
              'for which we have not included a reader. Choices are:',
              *openmm_topology_readers.keys())
        raise
    return topology


# Frame counting. mdtraj.open(...) returns a format-specific file handle
# (DCD/XTC) whose __len__ reports the frame count without loading coordinates,
# and -- importantly -- without needing a topology at all. LOOS is used only as
# a fallback, because loos.createSystem() cannot read the one topology format a
# GROMACS run actually has:
#
#     >>> loos.createSystem('topol.top')
#     RuntimeError: Error- unknown system file type 'top'
#
# That RuntimeError is not a loos.LOOSError, so it escaped the except clause
# below, propagated out of calx_remaining_steps, and killed the Farmer process
# on its first tick of any GROMACS run in an environment where LOOS imports --
# which is every environment where the reimaging tools are installed.
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

if _mdtraj is None and loos is None:
    print('Neither mdtraj nor LOOS is importable; frame counting will fail.')


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

    `top_fn` is only consulted by the LOOS fallback; the mdtraj path does not
    need it, which is what makes this work for GROMACS runs whose `top_fn` is a
    `.top`.

    The fallback tries the harvester's solute-only topology as well, because
    after a harvest the trajectory name is a symlink to a *stripped* copy: LOOS
    built from the wet `top_fn` would hit an atom-count mismatch, get swallowed
    by the broad except below, and report 0 frames -- from which the tender
    concludes the generation never ran.
    """
    traj_p = Path(traj_fn)
    if not traj_p.is_file() or traj_p.stat().st_size == 0:
        return 0
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
                # Deliberately broad: createSystem raises RuntimeError for an
                # unsupported topology, LOOSError for an unreadable frame, and
                # the orchestrator must survive both.
                print(f'LOOS could not read {traj_fn} with topology '
                      f'{candidate}: {type(exc).__name__}: {exc}.')
        print(f'No usable topology for {traj_fn}; treating as empty.')
    return 0


"""
The harvest itself lives in `harvester.harvest_generation`; these two names are
kept because they are what existing submit scripts call. Both now pin a backend
and hand off, so an old script gets the frame-count guard, the sentinel and the
seam handling without being rewritten.

hconfig keys:
 - `'harvester_subset'`: selection string for the solute (LOOS syntax by
   default; set `'harvester_subset_syntax': 'mdtraj'` for the other dialect).
 - `'downsample_frq'`: keep every Nth frame in the solvated stream.
 - `'harvester_structure'`: structure file to build the model from. REQUIRED for
   GROMACS runs, whose `top_fn` is a force-field topology that neither LOOS nor
   mdtraj can build a model from. Defaults to `config['top_fn']`.
"""


def strip_and_downsample(config_fn, harvester_config_fn):
    """Backwards-compatible entry point pinning the LOOS backend."""
    from . import harvester
    return harvester.harvest_generation(
        config_fn, harvester_config_fn, backend=harvester.BACKEND_LOOS)


def strip_ds_mdtraj(config_fn, harvester_config_fn):
    """Backwards-compatible entry point pinning the mdtraj backend."""
    from . import harvester
    return harvester.harvest_generation(
        config_fn, harvester_config_fn, backend=harvester.BACKEND_MDTRAJ)



# These basic strings are useful in many cases on clusters using the scheduler named as the key.
# NOTE the format target '{job_name}' has to appear for the default queue parser to find the job.
# The 'NODE:' / 'GPU:' echoes are how BadNodeRegistry learns which host
# produced a failure and which device was on it -- keep them if you
# replace this fstring with your own and want bad-node blocking to work.
# {exclude_nodes} expands to a scheduler directive line excluding any
# nodes BadNodeRegistry has flagged (empty when none are blocked).
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
    # -J (not -j) is sbatch's job-name flag; -j is not an sbatch option at all,
    # so the old template was rejected outright. The shebang matters too: with
    # no interpreter line the job runs under the submitting user's login shell.
    # {exclude_nodes} expands to '' when nothing is blocked, and a bare blank
    # line is fine -- Slurm stops scanning #SBATCH directives at the first
    # non-comment, non-blank line, so keep all directives above the echoes.
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

# Preempt-aware variants: install a SIGTERM trap that touches the sentinel
# file SentinelReporter watches for, then background+wait the python
# invocation so the trap can fire (bash blocks signal delivery while a
# non-builtin foreground command runs). The trailing `sleep 70` keeps the
# script alive past the 60s preempt grace period so Slurm records the job
# as CANCELLED rather than FAILED. Pair with Farmer(handle_preempt=True).
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
    # --signal=B:TERM@120 is what makes this fire on a WALLTIME boundary as well
    # as on preemption: without it Slurm only signals at the very end of the
    # allocation, leaving no time to checkpoint. B: targets the batch shell, so
    # the trap below runs rather than the signal going straight to mdrun.
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

# MPS-packed variant: one job, one GPU, K replicas sharing the card through the
# CUDA Multi-Process Service. Pair with gmx_pack.gmx_pack_sim_block_json, which
# does the per-replica core pinning and the per-member outcome reporting.
#
# Two details that bite:
#   * The pipe and log directories are keyed on $SLURM_JOB_ID. Two packed jobs
#     landing on the same node otherwise share -- or clobber -- one daemon.
#   * If the daemon fails to start the mdruns still run, just time-sliced at
#     10-20% worse throughput. The runner checks and says so loudly; the script
#     echoes its own failure too, rather than tolerating it silently.
#
# --signal=B:TERM@120 is what makes the checkpoint handshake fire at a WALLTIME
# boundary as well as on preemption; B: sends it to the batch shell so the trap
# runs instead of the signal going straight to mdrun.
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

# Basic report to print _only_ a list of job ids associated to this runner.
# Should have 'title' fstring target somewhere to purify spurious jobids.
# update_jids calls:
#   self.scheduler_report_fstring.format(title=self.config_template['title'])
basic_scheduler_reports = {
    "lsf": "bjobs -o JOBID -noheader -J '{title}-*'",
    # -h -o '%i %j' prints JobID and untruncated JobName, two whitespace-separated columns.
    # The default -O Name truncates to 8 chars, which silently breaks title matching.
    # Curly braces must be escaped with curly braces when using awk via str.format.
    "slurm": "squeue --me -h -o '%i %j' | awk '/{title}/ {{print $1}}'"
}

# Basic report to print the name, and then the jobid, for each job with job title
# created by orchestrator. Allows scripts to associate currently running jobs to
# their seed, clone, and gen indexes. Output should be a string where each new line
# is a job, with the Job ID in the first field and the Job Name in the second.
# __init__ from Orchestrator calls:
#   self.scheduler_assoc_fstring.format(title=self.config_template['title'])
basic_scheduler_assoc_reports = {
    "lsf": "bjobs -o 'JOBID JOB_NAME' -noheader -J '{title}-*'",
    # awk (not grep) so a clean queue exits 0 instead of grep's exit-1-on-no-match,
    # which would crash the boot-time sp.check_output in Farmer.__init__.
    "slurm": "squeue --me -h -o '%i %j' | awk '/{title}/'"
}


# Per-scheduler post-mortem checks for whether a finished job was preempted.
# Used by Clone.check_start_gen to avoid charging preemptions against the
# per-gen restart_attempts budget. Returns False on any error (missing
# binary, timeout, accounting gap, unknown state) so a true failure still
# counts as a restart.
def slurm_was_preempted(jid):
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
    try:
        out = sp.check_output(
            ['bjobs', '-d', '-o', 'exit_reason', '-noheader', str(jid)],
            text=True, timeout=30
        ).strip()
    except (sp.CalledProcessError, sp.TimeoutExpired, FileNotFoundError):
        return False
    return 'TERM_PREEMPT' in out


preemption_checkers = {
    'sbatch': slurm_was_preempted,
    'bsub': lsf_was_preempted,
}


# Bad-node detection: when a clone's gen aborts at 0 steps because the
# *node* is broken (stale CUDA driver against a too-new PTX, libc/GLIBC
# ABI mismatch, missing GPU), retrying the same gen routes a new sbatch
# straight back to the same broken queue/partition and the scheduler
# tends to hand it to the same node. The result is a cascade: every
# clone the farmer tries to relaunch lands on the bad node and dies, and
# the farmer eventually exhausts each clone's restart budget.
#
# BadNodeRegistry breaks that cycle. When a 0-step abort is detected,
# the registry scans the gen's scheduler log for a known "this is the
# node, not the sim" pattern, harvests the NODE: line, and appends the
# node to a scheduler-directive line that's injected into the next
# submission via {exclude_nodes} in scheduler_fstring. It also writes a
# human-readable breadcrumb row to a persistence file (default
# bad_nodes.txt) so (a) a restarted farmer doesn't re-learn the same
# bad nodes and (b) the user has a paper trail showing which nodes
# failed how, when, and where to read the offending log.

default_bad_node_patterns = (
    'CUDA_ERROR_UNSUPPORTED_PTX_VERSION',
    'CUDA_ERROR_NO_DEVICE',
    'CUDA_ERROR_INVALID_DEVICE',
    'CUDA_ERROR_NOT_INITIALIZED',
    'CUDA driver version is insufficient',
    'No CUDA-capable device is detected',
    'Failed to initialize NVML',
    # GLIBC ABI mismatch ("version `GLIBC_2.34' not found"). Quoted
    # half is enough -- the rest of the line varies.
    "version `GLIBC_",
)

default_scheduler_log_names = {
    'sbatch': 'slurm.out',
    'slurm': 'slurm.out',
    'bsub': 'lsf.out',
    'lsf': 'lsf.out',
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
    scheduler log (slurm.out / lsf.out), looks for a `bad_node_patterns`
    hit, harvests `NODE:` (and `GPU:` if present), appends a breadcrumb
    row, and refreshes scheduler_kws['exclude_nodes'].

    Breadcrumb file is plain text, tab-separated. Comment-out (prefix
    `#`) or delete rows to clear entries; the farmer rereads the file
    on boot.
    """

    BREADCRUMB_HEADER = (
        '# mdfarmer bad-nodes blocklist\n'
        '# This file is appended to whenever a clone aborts at 0 steps\n'
        "# on a node whose log matches a known fatal-on-this-node pattern.\n"
        '# The farmer excludes these nodes on subsequent submissions via\n'
        "# the '{exclude_nodes}' placeholder in scheduler_fstring.\n"
        '#\n'
        '# To clear an entry: comment out (prefix `#`) or delete the row\n'
        '# and restart the farmer. Lines starting with `#` and blank lines\n'
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
            # New format: timestamp\tnode\tgpu\tpattern\tlog_path\tclone_tag
            # Tolerate older / hand-edited rows: pick the first field
            # that looks like a hostname (no spaces, not an ISO date).
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

    # Scan the gen's scheduler log for a known fatal-on-this-node
    # pattern. If matched, harvest the NODE: line (and GPU: if present),
    # write a breadcrumb row, add the node to the in-memory exclude set,
    # and refresh scheduler_kws['exclude_nodes']. Returns the node name
    # if recorded, else None.
    def scan_and_record(self, gen_dir, clone_tag):
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
                  'exclude. Add `echo "NODE: $SLURMD_NODENAME"` (or LSF '
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


# Resolve which OpenMM Platform to use and which platformProperties to apply.
# If platform_name is set, demand that exact platform (raise if it can't load).
# Otherwise pick the fastest available non-Reference platform; raise if only
# Reference is available, since AMOEBA / large-system MD on Reference is
# effectively a hang from the scheduler's perspective. Filter platform_properties
# to those the chosen platform actually exposes so e.g. {'Precision': 'mixed'}
# applies cleanly on CUDA/HIP/OpenCL but is silently dropped on CPU.
def select_platform(platform_name=None, platform_properties=None):
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


# Resume-correctness helpers: state.xml ↔ DCD alignment.
#
# On every resume we need the DCD's last frame's logical step to equal
# the state.xml's stepCount exactly, otherwise the next append lands at
# the wrong place in time and gen-to-gen concatenation drifts. The
# helpers below validate the state file, read the DCD header's frame
# accounting, and (if the kill happened between the DCD report and the
# checkpoint report) truncate the DCD to match.


def is_state_xml_usable(p: Path) -> bool:
    if not p.is_file() or p.stat().st_size == 0:
        return False
    try:
        mm.XmlSerializer.deserialize(p.read_text())
    except Exception as exc:
        print(f'is_state_xml_usable: deserialize failed for {p}: {exc}')
        return False
    return True


def state_xml_step_count(p: Path) -> int:
    import xml.etree.ElementTree as ET
    root = ET.parse(p).getroot()
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
        # ints[10] is at byte offset 48 from file start — the
        # with-unit-cell flag (1 if frames carry the 6-double box record).
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
    # PBC block: 4 + 6*8 + 4 = 56 bytes
    # Each coord record: 4 + 4*n_atoms + 4 = 8 + 4n
    # 3 coord records: 3*(8+4n) = 24 + 12n
    return (56 if with_unitcell else 0) + 24 + 12 * n_atoms


def truncate_dcd_to_nframes(p: Path, target_nframes: int) -> int:
    """Reduce a DCD's frame count to target_nframes by rewriting the
    header nset field and truncating trailing bytes. No-op if already
    at target. Refuses to grow. Idempotent on partial completion.
    """
    info = dcd_header_info(p)
    cur_nset = info['nset']
    if target_nframes > cur_nset:
        raise ValueError(f'truncate_dcd_to_nframes refuses to grow '
                         f'{p}: current nset={cur_nset}, target='
                         f'{target_nframes}')
    if target_nframes == cur_nset:
        return cur_nset
    frame_size = dcd_frame_size(info['with_unitcell'], info['n_atoms'])
    new_size = info['header_size'] + target_nframes * frame_size
    # Rewrite nset first, then truncate. If interrupted between, the
    # file's nset is below its byte length; the next call computes the
    # same target and is a no-op (trailing bytes stay as harmless
    # padding that mdtraj/LOOS ignore since they honor nset).
    with open(p, 'r+b') as f:
        f.seek(8)
        f.write(struct.pack('<i', target_nframes))
    os.truncate(str(p), new_size)
    return target_nframes


def calx_remaining_steps(traj_fn, top_fn, total_steps, write_interval):
    traj_len = get_traj_len(traj_fn, top_fn)
    remaining = total_steps - traj_len * write_interval
    # A negative result means the traj has more frames than the gen
    # should contain — usually duplicated frames from an append against
    # a stale state.xml, or a write_interval / total_steps mismatch
    # between the on-disk config and the current run. Surface it so it
    # doesn't masquerade as "gen complete."
    if remaining < 0:
        print(f'WARNING: calx_remaining_steps({traj_fn}) = {remaining}; '
              f'traj has {traj_len} frames at write_interval={write_interval} '
              f'but total_steps={total_steps}. Treating as complete; check '
              f'for over-appended trajectory.')
    return remaining


# This won't be nicely jsonizable unless all default and provided vals are.
def merge_args_defaults_dict(function, **kwargs):
    sig = inspect.signature(function)
    # create a dictionary of the parameters and their defaults.
    config = {p: sig.parameters[p].default for p in sig.parameters}
    # overwrite the defaults wherever an option was specified
    config.update(kwargs)
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


def dir_seeds_clones_gens(top_lvl: Path, seed_index, clone_index, gen_index, pad,
                          sep='-', padchar='0', mkdir=True):
    p = top_lvl / fdir('seed', seed_index, pad, sep=sep, padchar=padchar) / \
        fdir('clone', clone_index, pad, sep=sep, padchar=padchar) / \
        fdir('gen', gen_index, pad, sep=sep, padchar=padchar)
    if mkdir:
        p.mkdir(exist_ok=True, parents=True)
    return p


# All trajectory formats for which a reporter exists. Add grace in future.
traj_suffixes = ['.dcd',
                 '.xtc']


default_steps = int(2.5e7)  # Given 0.004 ps timestep,
# this is 100 ns of simulation.
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
    # None -> auto-select fastest available non-Reference platform.
    # Set explicitly (e.g. 'CUDA', 'HIP', 'OpenCL') if you want to force one.
    platform_name=None,
    # Precision is filtered against the chosen platform's supported properties,
    # so this works on CUDA/HIP/OpenCL and is silently dropped on CPU.
    platform_properties={'Precision': 'mixed'},
    steps=default_steps,
    state_data_kwargs=default_state_data_kwargs,
    eq_steps=None,
    write_interval=2500,
    minimize_first=False,
    temperature=300,
    new_velocities=False
)

# You'll need to replace all of these, but I wanted it to be more clear what the slots were.
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


#  make two trajs--one stripped of solvent, the _other_ downsampled by some integer factor but not dried.
# `harvest_generation` picks its backend from the box on the trajectory, so the
# same script is right for a rectangular cell (LOOS, streaming) and a triclinic
# one (mdtraj, chunked). It is idempotent: a requeued harvest job that already
# ran is a no-op, not a second pass over its own output.
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
