import inspect
from pathlib import Path

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



# These basic strings are useful in many cases on clusters using the scheduler named as the key.
# NOTE the format target '{job_name}' has to appear for the default queue parser to find the job.
basic_scheduler_fstrings = {
    "lsf": inspect.cleandoc("""#!/bin/bash
                #BSUB -J {job_name}
                #BSUB -e lsf.out
                #BSUB -o lsf.out
                {gpu_line}
                #BSUB -q bowman

                python -c 'from omm_rst_runner import omm_basic_sim_block_json as runner; runner("config.json")'
                """),
    "slurm": inspect.cleandoc("""
                #SBATCH -j {job_name}
                #SBATCH -e slurm.out
                #SBATCH -o slurm.out
                {gpu_line}
                #SBATCH -p all
                
                python -c 'from omm_rst_runner import omm_basic_sim_block_json as runner; runner("config.json")'
                """)
}

basic_scheduler_kws = {
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
# their run, clone, and gen indexes. Output should be a string where each new line
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


def dir_runs_clones(top_lvl: Path, run_index, clone_index, pad, sep='-',
                    padchar='0', mkdir=True):
    p = top_lvl / fdir('run', run_index, pad, sep=sep, padchar=padchar) / \
        fdir('clone', clone_index, pad, sep=sep, padchar=padchar)
    if mkdir:
        p.mkdir(exist_ok=True, parents=True)
    return p


def dir_runs_clones_gens(top_lvl: Path, run_index, clone_index, gen_index, pad,
                         sep='-', padchar='0', mkdir=True):
    p = top_lvl / fdir('run', run_index, pad, sep=sep, padchar=padchar) / \
        fdir('clone', clone_index, pad, sep=sep, padchar=padchar) / \
        fdir('gen', gen_index, pad, sep=sep, padchar=padchar)
    if mkdir:
        p.mkdir(exist_ok=True, parents=True)
    return p


# All trajectory formats for which a reporter exists. Add grace in future.
traj_suffixes = ['.dcd',
                 '.xtc']
