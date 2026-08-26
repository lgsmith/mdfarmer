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

Launch the tender from the directory that holds your seed structures, system, integrator, and topology files. `Farmer.__init__` resolves any relative paths in the config to absolute paths once at boot, against the cwd it was started in. Absolute paths in your script work from anywhere; relative paths only work if you launch from the right directory.

```bash
#!/bin/bash

nohup python straight-sampling-farmer.py > straight-sampling-farmer.out 2>&1 &
```

This nohupped process can be annoying to stop. You can find it using `pgrep` and `pkill`, but if you want it to stop gracefully you can create a file in the directory you launched it from called `stop`. This file can be empty, it just needs to be present. At each update interval, if the script detects a file with that name, it'll exit and the process will return.

## Running GROMACS

Pass `runner=mdf.gmx_generation` and the Farmer picks the matching run script,
disk-recovery function and progress hook; supplying a mismatched set by hand is
refused rather than half-applied. Build the config template with
`gmx_config_template`, which records every `gmx_generation` parameter:

```python
cfg_template = mdf.gmx_config_template(
    traj_dir_top_level='trajectories',
    title='mysystem',
    structure_fn='start.gro',      # grompp -c for generation 0
    mdp_fn='prod.mdp',
    steps=12_500_000,
    steps_per_gen=12_500_000,
    write_interval=50_000,
    traj_suffix='.xtc',
    mdrun_args=['-nb', 'gpu', '-pme', 'gpu', '-update', 'gpu', '-pin', 'on'],
)
farmer = mdf.Farmer(..., runner=mdf.gmx_generation, config_template=cfg_template)
```

`gmx_config_template` fills placeholders for the five parameters a Clone
supplies per generation (`seed_index`, `clone_index`, `gen_index`, `seed_fn`,
`top_fn` — all overwritten before any job reads them) and omits the three
`gmx_pack` injects at runtime. Passing either group through to `gmx_generation`
raises, so building the template any other way means knowing both lists.

Generations chain with `gmx convert-tpr -nsteps` plus `mdrun -cpi`, not
`grompp -t`: the later generation inherits its predecessor's tpr, so the
integrator parameters and the Nose-Hoover / Parrinello-Rahman coupling state
carry across and the step and time counters stay globally continuous.

Only generation 0 reads the base `.mdp`, and the runner controls a handful of
its keys: `nsteps`, `init-step`, `nstxout-compressed`, `gen-vel`, `continuation`
and the seeds. `nsteps` is the campaign-absolute step a generation must reach,
so `init-step` is pinned to 0; a base `.mdp` that sets it gets a note on stdout
and is overridden. Everything else is inherited verbatim.

Every launch writes its own `prod.partNNNN.xtc`, merged with `gmx trjcat` when
the generation finishes. Where two parts cover the same time trjcat keeps the
**later** file's frames, so a part left behind by a relaunch that rewound to an
earlier checkpoint is renamed out of the way first — otherwise trjcat would
splice the abandoned branch into the one that actually continued.

Generation 0 is built with `grompp -c` from `seed_fn` when that is a structure
file, falling back to `structure_fn`. That is what lets seeds differ in
topology, and what lets an adaptive scheme reseed from a configuration it
picked rather than from the structure the campaign started with.

`n_gens` is a count, and generations are numbered from 0, so `n_gens=9` runs
generations 0 through 8. A clone retires once it has finished the last of them.

> Upgrading a campaign that ran before this was true: it may hold a generation
> numbered `n_gens` or higher, which the older code submitted by mistake. The
> tender now retires such a clone immediately and stops minding that job. It
> does not cancel it — mdfarmer never runs `scancel` — so check for those
> directories before you restart, and harvest or cancel them yourself.

### Packing replicas onto one GPU

Where Slurm exposes only a `gpu` gres — no `mps`, no `shard` — it cannot
co-schedule two jobs onto one card, so the packing has to happen inside one
job. `pack_size` members share one sbatch, one MPS daemon and one generation
step:

```python
farmer = mdf.Farmer(
    ...,
    pack_size=2,
    pack_cpus_per_task=16,
    # Optional: a policy other than consecutive runs of the priority order.
    pack_grouping=lambda clones: [...],
    # Optional: per-member core widths, for members with different core knees.
    pack_member_cores=[12, 4],
    # Optional: one dict per seed, applied over the shared template.
    seed_config_overrides=[dict(mdrun_args=[...]), dict(mdrun_args=[...])],
)
```

Each replica gets a private, contiguous block of cores (`-ntomp`, `-pinoffset`,
`-pinstride`), and `-ntmpi 1` on the thread-MPI builds that accept it — a
real-MPI build takes its rank count from `mpirun` and makes that flag fatal, so
it is probed for rather than assumed. Without the pinning, two replicas that
both say a bare `-pin on` start at core 0 and fight over the same cores, which
reads as node variance rather than a misconfiguration.

Failure is per-member: one replica raising does not abort the others, and the
tender fails exactly one clone. Preemption is the exception that must reach
everyone, so the sentinel is watched once and SIGTERM is fanned out to all K
mdruns — that handshake is what protects trajectory contiguity.

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

Another nice troubleshooting step can be simply trying to resubmit the job generated for a particular generation. The input files will expect the CWD to be the generation directory they reside in, but if you change to that directory and you have the correct env active you should simply be able to re-launch the scheduler script you can find there. Sometimes this can make issues with pathing more clear because you can edit the `config.json` and then relaunch just one gen until it works correctly, then go back and make those edits in the template to get the tender-process working with the corrected inputs.

### Analysis

There are two modes of analysis with this type of dataset. If you have fewer clones, but long length per clone, you could do a classical 'replicate' analysis of an observable across contiguous trajectories. If you have multiple seeds, or many clones with short generations, or some combination thereof, you're better off making some kind of transition-counting model from the data, such as a Markov State Model.

We're hoping to add some scripts for both modes of analysis--mostly these will be simple functions that just use the configurations you've given for the farmer and or the structure of the data-set tree to provide you with lists of trajectories that might be useful, such as a nested list of file-paths that follows the overall structure of the tree. If you're writing functions like this yourself, note that python's `glob` and `iterdir` functionalities provide sub-paths in no particular order. The reason the directory names are padded is so that the built-in `sorted` will 'just work' with a semantic sort on the file names, but you do have to bother to use sorted if you're writing your own iterator and you want the order to be 1. the same and 2. for the generations to be sequential each time you read the files. Note that the top level file titled `traj_list.txt` records the trajectory paths in the order they are produced, which could be good for some things like a function that surveys how much data has been collected thus far, but is probably not what you want for most analysis.

### Reimaging (making molecules whole again)

GROMACS writes whatever coordinates the integrator is holding, and it wraps
atoms into the box as it goes. Molecules that straddle a boundary come out
**split**, and there is no avoiding it — it is fundamental to how the engine
stores coordinates. On a real trp-cage system deliberately positioned across a
box face, the raw `.xtc` had 214 bonds longer than 2.5 Å, the worst of them
4.7 nm — a whole box length. Anything you compute per molecule on that
trajectory (R_g, RMSD, contacts, a picture) is wrong.

`mdfarmer.reimage` fixes this, with two backends chosen by the unit cell:

```python
from mdfarmer import reimage

# picks the backend from the box actually recorded in the trajectory
reimage.reimage_trajectory('prod.xtc', structure_fn='start.gro',
                           top_fn='topol.top', tpr_fn='prod.tpr')
# -> prod-whole.xtc   (the raw prod.xtc is never touched)
```

* **`'loos'` — orthorhombic ("box") cells only.** LOOS's periodic box is three
  numbers, and its readers keep only the diagonal of a triclinic box *without
  raising* — hand it a rhombic dodecahedron and it reports a rectangular cell
  and every minimum-image result downstream is quietly wrong. So this backend
  refuses a non-orthorhombic cell rather than producing plausible garbage.
* **`'trjconv'` — any cell, and the only option for triclinic.** Shells out to
  `gmx trjconv -pbc mol -ur compact`, which needs the run's `.tpr` because that
  is where molecule definitions live.

Two things worth knowing:

**Molecule membership never comes from the structure file.** A `.gro` carries no
bonds, and a bondless LOOS model makes `splitByMolecule()` return *one group
containing the whole system* — reimaging then degenerates into a single global
translation that looks like it worked. Bonds alone are not enough either: a
TIP4P-ice virtual site is bonded to nothing, so connected components over bonds
strand every `MW` in its own "molecule". Molecule blocks are therefore read from
the GROMACS `.top` through `openmm.app.GromacsTopFile`, whose chains reproduce
the `[ molecules ]` section exactly.

**Reimaging is checked, not trusted.** There are many ways imaging-by-atom goes
wrong quietly, so the LOOS backend verifies its own output against a physical
invariant — bond lengths, computed *without* the minimum-image convention,
against LOOS's `long-bond-finder` cutoff of 2.5 Å — and raises if any bond is
still overlong. You can run the same checks yourself:

```python
n_bad, violations = reimage.check_bond_lengths('prod-whole.xtc', top_fn='topol.top')
margin = reimage.check_anchor_distances('prod-whole.xtc', ranges)
```

`check_anchor_distances` reports the one thing that limits the LOOS backend:
`mergeImage()` minimum-images every atom against its molecule's *first* atom, so
it is only correct while no atom is more than half a box edge from that anchor.
Folded trp-cage in a 4.67 nm box already uses **89%** of that margin — an
extended conformation will exceed it, at which point the LOOS backend is outside
its safe regime and you want `backend='trjconv'`, which walks the bond graph
instead. The check reports the margin as a fraction so you can watch it.

### Harvesting (reducing a finished generation, safely)

When a generation completes, the `Harvester` submits a job that turns its raw
trajectory into the two streams you actually keep — a **dry** stream (every
frame, solute only) and a **downsampled** stream (every Nth frame, solvent
kept) — and then replaces the original with a symlink to the dry one. That last
step is irreversible, so the harvest is built to be interrupted:

```python
harvester = mdf.Harvester(
    harvester_template=mdf.default_harvest_shellscript_slurm.format(
        queue_name='ccb', harvest_time='02:00:00'),
    scheduler='sbatch',
    run_config=dict(
        harvester_subset='resid <= 20',      # LOOS syntax by default
        downsample_frq=10,
        # REQUIRED for GROMACS: `top_fn` is a force-field topology, and neither
        # LOOS nor mdtraj can build a model from a `.top`. Point this at the
        # structure the run actually started from.
        harvester_structure='reference/restarts/native_cluster000.gro',
        steps_per_gen=STEPS_PER_GEN,
    ),
)
```

Set `harvester_subset_syntax='mdtraj'` if you would rather write the selection
in mdtraj's language. Either way the selection is resolved to atom indices once
and both backends slice by index, so which backend runs cannot change which
atoms come out.

**The backend is chosen by looking at the box.** A rectangular cell goes to
LOOS, which streams frame by frame; anything triclinic goes to mdtraj's
`iterload`. This cannot be a `try`/`except` around the LOOS call, because LOOS
does not raise on a triclinic cell — it keeps the diagonal and carries on (see
the reimaging section) — so the handler would never fire on the case it exists
for. The off-diagonals are read off frame 0 and the choice is made before either
engine is touched. Neither path holds the trajectory in memory: a 1 µs
generation at 10 ps sampling is ~29 GB of coordinates.

**Frames are counted, not sized.** A harvest killed mid-write leaves both
outputs nonzero and truncated, so size cannot decide whether one finished. The
source is counted first, the output counts are *computed* from the frame plan,
and the written files are re-counted off disk. Any mismatch raises and leaves
the original alone.

**It is idempotent.** A `.harvested` sentinel carrying the counts is written
last; a re-run short-circuits on it. A generation whose original is already a
symlink but which has no sentinel is a harvest that died in that window — its
outputs are checked against the plan and the sentinel is written, rather than
the write loop running again with its own output as input. Under a requeueing
scheduler this is not optional.

**Two spacing rules, checked when a `Farmer` is built** — not when the harvest
runs, hours into a campaign:

```
steps_per_gen % write_interval == 0
(steps_per_gen / write_interval) % downsample_frq == 0
```

The first puts each generation's last frame exactly on the checkpoint the next
one restarts from; miss it and every seam silently drops the trajectory between
the two. The second keeps the downsampled stream evenly spaced across the
concatenation.

**The seam frame is dropped exactly once.** GROMACS writes a frame at the step
it restarts from, so generation N's last frame and generation N+1's first are
the same time point. Concatenating without handling that puts a duplicate at
every seam — which never fails loudly, it just biases lag times and kinetics.
Frame 0 of every generation after the first is dropped, and the downsample phase
is driven by a *global* frame index so it does not reset at each boundary. The
OpenMM reporters do not write that frame; the convention is detected from the
frame count rather than assumed.

**The time axis is carried through.** LOOS's `XTCWriter` numbers frames from
its own counters (`dt_ = 1.0`, `step_ = 0`, `steps_per_frame_ = 1`) and mdtraj
fills `step` with the frame index, so both backends read the source's step and
time and write them explicitly, verified against the source's last frame.

> Trajectories harvested before this was in place carry a fabricated time axis —
> 1 ps per frame, stamped with the frame index. Recompute from the generation's
> `write_interval` and `dt`; the coordinates are fine.

Two things to check on a campaign after the fact:

```python
# Generations that ran but have no sentinel. On a campaign that is still
# running, skip_newest=True leaves out the one each clone is mid-way through.
mdf.unharvested_gen_dirs('trajectories', skip_newest=True)
# Frame count, spacing and duplicates across a clone's harvested generations.
mdf.verify_dry_chain(sorted_gen_dirs)
```

The harvest does **not** reimage — it preserves the per-frame box (LOOS subsets
share the parent's `SharedPeriodicBox`, so the dry stream carries a live cell,
not a frozen one), which is what lets you run `mdfarmer.reimage` on the dry
stream afterwards.

## Tests

`tests/` holds plain scripts — no framework. `python tests/run_all.py` runs
them all and reports pass, fail or skip per suite; see `tests/README.md` for
what each one needs.

## AI assistance

Parts of this codebase have been developed with assistance from Anthropic's Claude (Opus 4.x family). Individual commits are not tagged with `Co-Authored-By` trailers; this section is the project-level attribution.

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
