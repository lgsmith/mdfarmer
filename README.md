# mdfarmer

Farm out simulations, tend your servers, reap data. It's a combine seeder/harvester for Molecular Dynamics Simulation that is written in simple python and uses configurations written to disk to seed jobscripts that are then run by your cluster's scheduler (so long as it's LSF or SLURM).

## Install

Currently, you need the following libraries:

- `openmm` for running the simulations
- `mdtraj` or `loos` for reading trajectory files.
- This codebase needs to be somewhere on the `PYTHONPATH`.
  - If you're trying to use a stable version, once there are releases, make sure you check out the version you want.[^1]
  - Normally, I put software I install myself in `~/software` because most clusters tend to mount your home directory to the node that your jobs land on. This means that cloning this software into that directory, then adding that directory to the python path in `~/.bashrc`, pretty much works in a 'normally' configured academic cluster.

Here is a series of install commands that'll get you roughly this setup:

```bash
mkdir $HOME/software
echo 'export PYTHONPATH=$HOME/software:$PYTHONPATH' >> $HOME/.bashrc
cd $HOME/software
git clone git@github.com:lgsmith/mdfarmer.git
mamba create -n mdfarmer -c conda-forge loos mdtraj openmm
```

In keeping with our intention that this codebase stay as simple as possible, you might want it to do other things, or have custom scripts that do other stuff like cluster data iteratively at certain milestones. We're working on an Adaptive Farmer class, but it's not ready yet.

## Example

Following is an example farmer configuration using LSF and starting from one solvated configuration:

```python
import mdfarmer as mdf
from pathlib import Path

cfg_template = mdf.default_straight_sampling_config_template.copy()
cfg_template['title'] = 'mytitle'
cfg_template['integrator_xml'] = str(Path('lmi-4ps.xml').resolve())
steps = int(2.5e6)
# steps = 10000  # to get things started, try running really short gens 
                # as a smoke test for whether all the things your farmer
                # needs are available and working.
cfg_template['steps'] = steps
cfg_template['state_data_kwargs']['totalSteps'] = steps
cfg_template['write_interval'] = 2500 
cfg_template['traj_suffix'] = '.dcd'  # I have been having issues with openmm's xtc reporter for these scripts, but it should be an option too.
queue_name = 'tolbertgpu'
harvester = mdf.Harvester(  # The default harvester makes a dry trajectory at full sampling frequency, and a downsampled but still solvated trajectory
    harvester_template=mdf.default_harvest_shellscript.format(
        queue_name=queue_name),
    scheduler='bsub',
    run_config=dict(
        harvester_subset='resid < 199',  # assuming a 200 residue macromolecule
        downsample_frq=10  # how much to downsample by in frames
    )
)
# Create an instance of Farmer with the configs we need.
farmer = mdf.Farmer(
    n_seeds=1,
    n_clones=100,
    active_clone_threshold=50  # how many clones will we try to schedule simultaneously.
    n_gens=9,
    config_template=cfg_template,
    seed_structure_fns=['start-state.xml'],
    top_fns=['start-state.pdb'],
    system_fns=['start-sys.xml'],
    scheduler='bsub',
    scheduler_fstring=mdf.basic_scheduler_fstrings['lsf'],
    scheduler_report_cmd=mdf.basic_scheduler_reports['lsf'].format(title=cfg_template['title']),
    scheduler_assoc_rep_cmd=mdf.basic_scheduler_assoc_reports['lsf'].format(title=cfg_template['title']),
    scheduler_kws=dict(gpu_line='#BSUB -gpu num=1:mig=1', queue_name=queue_name),
    harvester=harvester,
    overwrite=True,
)
# Start launching jobs. Try to launch up to active_clones clones, 
# then wait for update_interval seconds before trying to (re) launch more.
farmer.start_tending_fields(update_interval=60)
```

The above should be saved in some python script--name it however you like. Here I've called mine `straight-sampling-farmer.py`.

Here we are asking for 50 clones to be run simultaneously, across a dataset of 100 clones. Each clone is going to do 2.5 million steps per generation, and each clone is going to do 9 sequential generations. As the comment says, shorten the generation time (and also the update interval) if you want to test whether things will be working properly with your system.

Right now, the actual tender process (i.e. the one that the instance of Farmer is being run by) should sit on your head-node and remain running even if you log off. If it checks up on its jobs every 60-120 seconds and sleeps the rest of the time, it is extremely unlikely to bog down your head node much. Choosing an extremely short `update_interval` (say `5`, or `1`) could be bad for a number of reasons, one being that schedulers take time to update and so the tender process could accidentally submit the same clone twice, because it believes the first time didn't work, which will sow chaos in your data fields. How short too short is could be different for different systems, but after you finish debugging I very much doubt you're gaining much by having that interval be on the long side. I recommend something like 20-30 seconds for debugging, and between 1 and 5 minutes for normal usage.

You launch the tender process by just calling the correct python on the script above. For debugging I recommend doing this in an interactive session, but normally these datasets take weeks or even months to collect so I often run them using the shell utility `nohup`, in a script such as the following. This allows me to call `tail` on `straight-sampling-farmer.out` to read what's going on.

```bash
#!/bin/bash

nohup python straight-sampling-farmer.py > straight-sampling-farmer.out 2>&1 &
```

This nohupped process can be annoying to stop. You can find it using `pgrep` and `pkill`, but if you want it to stop gracefully you can create a file in the directory you launched it from called `stop`. This file can be empty, it just needs to be present. At each update interval, if the script detects a file with that name, it'll exit and the process will return.

## Dataset structure

These tools generate datasets that look roughly like the following:

```text
my-sampling-project
└── seed-000
    ├── clone-000
    │   └── gen-000
    │   |    ├── bsub.sh
    │   |    ├── config.json
    │   |    ├── lsf.out
    │   |    ├── run.py
    │   |    ├── state.xml
    │   |    ├── traj.dcd
    │   |    └── traj.out
    │   └────gen-001
    │        ├── bsub.sh
    │        ├── config.json
    │        ├── lsf.out
    │        ├── run.py
    │        ├── traj.out
    │        └── traj.xtc
    │
    └── clone-001
        ├── gen-000
        |    ├── bsub.sh
        |    ├── config.json
        |    ├── lsf.out
        |    ├── run.py
        |    ├── traj.out
        |    └── traj.xtc
        └── gen-001
             ├── bsub.sh
             ├── config.json
             ├── lsf.out
             ├── run.py
             ├── traj.out
             └── traj.xtc

```

Here what you see is a folding-at-home like directory tree structure, where we've started two `Clone`s from one seed structure. Each clone has been run for two generations. Generations are always meant to be interpreted as junks of a contiguous trajectory. So the first frame of gen-001 would have been the next saved frame of `gen-000`'s simulation if it were extended for longer, given that most integrators are both chaotic and stochastic so it probably wouldn't be precisely that. Each seed ought to represent a starting configuration for the system---thus if I wanted to start from 10 different structures I'd have seeds 1-10 at the level right below the project name. Clones all start from the same configuration, but with different initial velocities (though you can change this to try to read velocities from starting configurations if you want).

The `bsub.sh` file is what got submitted to generate the `traj.out`, `traj.dcd`, and `state.xml` (which is the restart file made by OpenMM's state-based checkpoint reporter, not the binary `.chk` reporter). They're made from taking the configuration variables specified in `config.json` and reading them into a function that sets up and performs the simulation they specify. The output from the scheduler is in `lsf.out`. As such, if you're troubleshooting some issue where your simulations don't run, the `lsf.out` (or if you're using slurm the `slurm.out`) will be good places to start to figure out what's wrong, in addition to whatever python trace-back might be provided by the farmer process. Actually looking in the `config.json` with a text editor may also be helpful. A lot of the problems I've run into historically have been path errors, so making sure the file-paths in the configs are correct is a good first step.

### Analysis

There are two modes of analysis with this type of dataset. If you have fewer clones, but long length per clone, you could do a classical 'replicate' analysis of an observable across contiguous trajectories. If you have multiple seeds, or many clones with short generations, or some combination thereof, you're better off making some kind of transition-counting model from the data, such as a Markov State Model.

We're hoping to add some scripts for both modes of analysis--mostly these will be simple functions that just use the configurations you've given for the farmer and or the structure of the data-set tree to provide you with lists of trajectories that might be useful, such as a nested list of file-paths that follows the overall structure of the tree. If you're writing functions like this yourself, note that python's `glob` and `iterdir` functionalities provide sub-paths in no particular order. The reason the directory names are padded is so that the built-in `sorted` will 'just work' with a semantic sort on the file names, but you do have to bother to use sorted if you're writing your own iterator and you want the order to be 1. the same and 2. for the generations to be sequential each time you read the files. Note that the top level file titled `traj_list.txt` records the trajectory paths in the order they are produced, which could be good for some things like a function that surveys how much data has been collected thus far, but is probably not what you want for most analysis.

## Extensions

We'd love to get your PRs for any extensions you think could be helpful. Even if you have a hacky solution to some problem you face, it might be interesting if the problem seems general. Feel free to post these kinds of solutions, or requests for development, as issues. If we get to it we get to it! We'll at least try to talk with you about how to improve any hacky solutions you've made.

We are trying to use a [feature-branch workflow](https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow), so please branch off of main and make sure any feature you're making a pull-request for can be fastforward-merged with main at the time of merger.

### Adaptive sampling

As mentioned previously, the primary objective of structuring all of this stuff in this way is to set ourselves up to implement Adaptive Sampling (AS) schemes using this codebase. If you think you can help us add a cool AS method, or you think our framing of the `AdaptiveFarmer` class needs some change to have the flexibility that you need for some sampling scheme that excites you, raise it in discussion.

Most AS schemes have an iteration time or cycle that they go through, which looks roughly like this:

1. The walkers (clones) do some specified amount of unbiased sampling: A generation, in our framing.
2. Some analysis workflow gets run that selects N new configurations, where N is the number of walkers; this is called ranking.
   - It could also conceivably change other simulation control variables, though most methods probably don't do anything other than pick new configurations.
3. Instead of using the terminal configurations from the most recent generations of sampling, the selected N configurations are used to start new simulations of the walkers; this is called seeding.
4. This is repeated until some sort of exit condition is reached, which could either be based on some computed heuristic or a prescribed number of iterations.

The goal with our implementation will be to make the analysis and configuration selection processes occur within a scheduler submission that the farmer process waits on to start more generations. Our secondary goal is to keep things organized so that it would be easy to switch the project back to straight sampling, where generations are seeded off of one another sequentially, as described [above](#dataset-structure). This is worth maintaining because often adaptive sampling creates very thin sampling in parts of phase space that may actually be important, leading to low quality statistical estimates across those regions (methods favoring exploration and conformational diversity can do this very easily). Thus it is often common practice to turn off the adaptive portion of the sampling after some exploration has been done and 'burn in' the conformations that have been found.

Because it'll make the bookkeeping easier, and because switching back and forth between AS and straight sampling should be as frictionless as reasonably possible, the meaning of these directories is going to remain fixed. Seed directories will only contain trajectories started from the same configuration. Replicas starting from that configuration will always be called clones, and placed within those seed directories. Generations will always be sequential end-on-end extensions of the previous generation. What this means in practice is that if you run AS with 10 walkers, and you do 3 rounds or iterations of AS wherein new conformations are picked for each walker each time, you'll have 30 seeds each containing one clone and one generation. If your adaptive sampling scheme allows you to sample configurations with replacement, you might have 10 walkers that are restarted from some particularly high-ranking seed, in which case their data will be recorded as clones within that seed. If you stop doing rounds of adaptive sampling but chose to do more sampling run with these tools, any clones you have will produce subsequent generations within those clone directories.

### Rolling or non-seasonal Adaptive Sampling

One could imagine a different model where adaptive simulations are launched in a rolling fashion, with the walkers grabbing new configurations from a 'ranking' list whenever they finish, and with the analysis process being rerun to reorder and extend that ranking list on some other schedule--perhaps a fixed amount of wall-time. We hope to support this kind of strategy as well since it fits better with resources that are very asynchronus/heterogeneous, which tends to be true for larger clusters with heterogeneous hardware.

### Replica methods

The basic dataset structure discussed above may not make as much sense for replica based methods. We're interested to hear about how they might be considered, of course, but in our framing replicas--when they exchange--become discontinuous, which would require a new seed directory. Thus, implementing replica methods could be done with some of these tools but probably requires planning a different dataset structure alltogether.

### Lazy Farmer

I found it harder to implement, so I dropped it, but one could imagine a better structure where each time a clone finishes a generation, it starts up the farmer process and tells the farmer to take note of the state of the sampling more broadly and give it a new seed (whether that be a sequential generational seed, or an adaptively selected one). The farmer would do this by looking at a serialzed database (probably a json, but maybe something else?) that contains the relevant project state. One complexity here is that the farmer would need to work under some kind of lock condition for the database file. These are solvable problems, but they'd be a significant to-do at this stage. Having a lazy farmer may be a good goal for the first major version bump of this project, if it grows and develops a lot.
  
[^1]: As with most academic code distributed through GitHub, immediately post cloning you'll be on the most recent commit to `main`. We are trying to follow a development modality where [features](https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow) get built and trialed in other branches, then merged into main via a pull request once they seem to be working.
