import inspect
from pathlib import Path
import json
import openmm as mm
from openmm import app


openmm_topology_readers = {
    '.top': app.GromacsTopFile,
    '.prmtop': app.AmberPrmtopFile,
    '.psf': app.CharmmPsfFile,
    '.pdb': app.PDBFile
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
    print('Expect restarting a mature dataset to be slow.')
    import mdtraj

    def get_traj_len(traj_fn, top_fn):
        if Path(traj_fn).stat().st_size == 0:
            length = 0
        else:
            length = len(mdtraj.load(traj_fn, top=top_fn))
        return length


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
    # note when using pipes to awk, curly braces must be escaped with curly braces.
    "slurm": "squeue --me --noheader -O JobID,name | awk '/{title}/ {{print $1}}'"
}

# Basic report to print the name, and then the jobid, for each job with job title
# created by orchestrator. Allows scripts to associate currently running jobs to
# their seed, clone, and gen indexes. Output should be a string where each new line
# is a job, with the Job ID in the first field and the Job Name in the second.
# __init__ from Orchestrator calls:
#   self.scheduler_assoc_fstring.format(title=self.config_template['title'])
basic_scheduler_assoc_reports = {
    "lsf": "bjobs -o 'JOBID JOB_NAME' -noheader -J '{title}-*'",
    "slurm": "squeue --me --noheader -O JobID,name | grep '{title}'"
}


def calx_remaining_steps(traj_fn, top_fn, total_steps, write_interval):
    traj_len = get_traj_len(traj_fn, top_fn)
    return total_steps - traj_len * write_interval


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
    platform_name='CUDA',  # CHECK TO BE SURE!
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
