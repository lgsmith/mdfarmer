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


try:
    import loos
    from loos import pyloos

    def get_traj_len(traj_fn, top_fn):
        if Path(traj_fn).stat().st_size == 0:
            length = 0
        else:
            m = loos.createSystem(top_fn)
            try:
                t = pyloos.Trajectory(traj_fn, m)
                length = len(t)
            except loos.FileReadError:
                print('Assumming empty file; cannot read:', traj_fn)
                length = 0
        return length


except ImportError:
    print('LOOS not in import path; falling back to MDTraj.')
    import mdtraj

    # mdtraj.open(...) returns a format-specific file handle (DCD/XTC) whose
    # __len__ reports frame count without loading coordinates. Avoids the
    # mdtraj.load(...) round-trip that materialized the whole traj just to
    # measure it -- which was prohibitive for mature datasets.
    def get_traj_len(traj_fn, top_fn):
        if Path(traj_fn).stat().st_size == 0:
            return 0
        try:
            with mdtraj.open(traj_fn) as fh:
                return len(fh)
        except (OSError, ValueError) as exc:
            print(f'mdtraj.open could not read {traj_fn}: {exc}; '
                  'treating as empty.')
            return 0


"""
To run this one, hconfig needs to contain the following keys:
 - `'harvester_subset'`: a LOOS selection string that produces the desired system subsetting.
 - `'downsample_frq'`: An int---the number of frames to skip over before writing another solvated frame.
"""


def strip_and_downsample(config_fn, harvester_config_fn):
    import loos
    from loos import pyloos as pl
    config_fp = Path(config_fn)
    config = json.loads(config_fp.read_text())
    hconfig_fp = Path(harvester_config_fn)
    hconfig = json.loads(hconfig_fp.read_text())
    sep = config['sep']
    model = loos.createSystem(config['top_fn'])
    subset_selection = hconfig['harvester_subset']
    subset = loos.selectAtoms(model, subset_selection)

    traj_name = config['traj_name']
    traj_suffix = config['traj_suffix']
    traj_fn = f'{traj_name}{traj_suffix}'
    traj = pl.Trajectory(traj_fn, model)
    downsample_frq = hconfig['downsample_frq']
    dry_outfn = f'dry{sep}{traj_fn}'
    dry_outp = Path(dry_outfn)

    down_outfn = f'downsample{sep}{traj_fn}'
    down_outp = Path(down_outfn)
    if traj_suffix == ".xtc":
        dry_outtraj = loos.XTCWriter(dry_outfn)
        downsampe_outtraj = loos.XTCWriter(down_outfn)
    elif traj_suffix == '.dcd':
        dry_outtraj = loos.DCDWriter(dry_outfn)
        downsampe_outtraj = loos.DCDWriter(down_outfn)
    else:
        raise NotImplementedError(
            f'{traj_suffix}: not implemented for basic strip and downsample')
    print('Preparing to loop over trj in strip and downsample.')
    while next(traj, False):
        dry_outtraj.writeFrame(subset)
        if traj.index() % downsample_frq == 0:
            downsampe_outtraj.writeFrame(model)
    # dump to PDB for topology
    subset.pruneBonds()  # Need to do this to ensure connects are correct.
    pdb = loos.PDB.fromAtomicGroup(subset)
    Path('dry-top.pdb').write_text(str(pdb))

    # if we've subset and also dried the trajectories, remove the original.
    # Should raise a file not found error if the call to stat()
    # is applied to a file that was never created
    if dry_outp.stat().st_size > 0 and down_outp.stat().st_size > 0:
        traj_p = Path(traj_fn)
        traj_p.unlink()
        # leave a symlink to dry traj so that frame counting efforts don't go awry
        traj_p.symlink_to(dry_outp)
    else:
        print('either', dry_outp, 'or', down_outp,
              'are size zero, refusing to unlink')


def strip_ds_mdtraj(config_fn, harvester_config_fn, sep='-', image_molecules=True):
    import mdtraj as md
    config_fp = Path(config_fn)
    config = json.loads(config_fp.read_text())
    hconfig_fp = Path(harvester_config_fn)
    hconfig = json.loads(hconfig_fp.read_text())
    model_name = str(config['top_fn'])
    # subset = loos.selectAtoms(model, subset_selection)

    traj_name = config['traj_name']
    traj_suffix = config['traj_suffix']
    traj_fn = f'{traj_name}{traj_suffix}'
    traj = md.load(traj_fn, top=model_name)
    if image_molecules:
        traj.make_molecules_whole(inplace=True)
        traj.image_molecules(inplace=True)
    top = traj.top
    subset_selection = hconfig['harvester_subset']
    if subset_selection:
        subset_iis = top.select(subset_selection)
        dry_traj = traj.atom_slice(subset_iis)
    else:
        dry_traj = traj.remove_solvent()
    
    dry_outfn = f'dry{sep}{traj_fn}'
    dry_outp = Path(dry_outfn)
    dry_traj.save(dry_outfn)
    dry_topp = dry_outp.with_suffix('.pdb')
    # dump to PDB for topology
    dry_traj[-1].save(str(dry_topp))
    del dry_traj

    # make and save downsampled traj
    downsample_frq = hconfig['downsample_frq']
    downsample_traj = traj[::downsample_frq]

    down_outfn = f'downsample{sep}{traj_fn}'
    down_outp = Path(down_outfn)
    downsample_traj.save(down_outfn)

    # if we've subset and also dried the trajectories, remove the original.
    # Should raise a file not found error if the call to stat()
    # is applied to a file that was never created
    if dry_outp.stat().st_size > 0 and down_outp.stat().st_size > 0:
        traj_p = Path(traj_fn)
        traj_p.unlink()
        # leave a symlink to dry traj so that frame counting efforts don't go awry
        traj_p.symlink_to(dry_outp)
    else:
        print('either', dry_outp, 'or', down_outp,
              'are size zero, refusing to unlink')
        

# These basic strings are useful in many cases on clusters using the scheduler named as the key.
# NOTE the format target '{job_name}' has to appear for the default queue parser to find the job.
basic_scheduler_fstrings = {
    "lsf": inspect.cleandoc("""#!/bin/bash
                #BSUB -J {job_name}
                #BSUB -o lsf.out
                {gpu_line}
                #BSUB -q {queue_name}

                python {run_script_name}
                """),
    "slurm": inspect.cleandoc("""
                #SBATCH -j {job_name}
                #SBATCH -e slurm.out
                #SBATCH -o slurm.out
                {gpu_line}
                #SBATCH -p {queue_name}
                
                python {run_script_name}
                """)
}

basic_gpu_lines = {
    "lsf": "#BSUB -gpu 'num=1:j_exclusive=yes'",
    "slurm": "#SBATCH --gpus=1"
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
default_harvest_shellscript = inspect.cleandoc("""#!/bin/bash
                #BSUB -J harvest
                #BSUB -o harvest.out
                #BSUB -q {queue_name}

                python -c 'from mdfarmer.utilities import strip_ds_mdtraj; strip_ds_mdtraj("config.json", "hconfig.json")'
                """)
