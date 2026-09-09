# shakedown-gmx — 5 packed GROMACS pairs of trp-cage

The packed half of the shakedown, on the same molecular system as
`examples/openmm-trpcage`. Where Slurm exposes only a `gpu` gres — no `mps`, no
`shard` — it cannot co-schedule two jobs onto one card, so the packing happens
inside one job. This example exercises that path end to end on real MD that
takes seconds: one sbatch, one MPS daemon, one `pack.lock`, two pinned mdruns,
per-member failure and recovery, and a harvest per finished generation.

**10 clones in 5 packs of 2, 3 packs running at a time, 5 generations each.**
Three of five rather than all five, so the tender actually has to wait for a
slot and refill it, the way the OpenMM arm's 5-of-10 does.

## The system

Trp-cage (20 residues, 304 atoms) in 3287 four-site waters with 6 K⁺ and 7 Cl⁻,
13465 atoms, at 277 K under a C-rescale barostat, in **amber03 + TIP4P-ice**
(`sampling-trpcage/systems/fresh-omm/native-277`, whose `meta.json` marks the
amber03 as benchmark-only).

`inputs/` holds it gzipped — `gmx.gro.gz` and `gmx.top.gz`, 218 KB for 1.0 MB —
plus `prod-277.mdp` committed plain, because it is 2 KB and it is the file you
will actually want to read. `prepare_inputs()` inflates the two before the
Farmer is built. Unlike OpenMM, GROMACS has no choice about this: `grompp` takes
file names, not streams.

Inflation is a one-off, not a per-boot cost: `inflate()` returns immediately if
the inflated file is already in `prepared/`, and writes through a `.partial`
name so a boot killed halfway leaves nothing that a later boot could mistake for
a finished file. Delete `prepared/` to force it again. Nothing here compresses
trajectories or restarts *during* a run: `state.cpt` is GROMACS' own binary
checkpoint written plain, and the only compression in the stream is the one
built into the XTC format.

`prod-277.mdp` is `sampling-trpcage/mdp/prod-277-a19opc3.mdp` with only the
output strides and the step count shortened — `nstlog`, `nstenergy` and
`nstxout-compressed` from 50000 to 2000, `nsteps` from 5e8 to 10000. The physics
is untouched, so the shakedown runs the same integrator, coupling and
constraints a real campaign would.

### Which arm this is: `-update cpu`, not `-update gpu`

The water here is four-site (TIP4P-ice: one virtual site per molecule), and
GROMACS refuses GPU update outright with virtual sites:

```
Update task on the GPU was required, but the following condition(s) were not
satisfied: Virtual sites are not supported.
```

That is measured on this system, not inferred. So this is
`sampling-trpcage`'s **`vitrification`** arm rather than its `fulloffload` one,
and it takes that arm's core budget with it: the update runs on the CPU, so the
per-replica knee is 12 cores, not the 4 that suffices when the whole step is on
the card. Measured here, solo on an RTX A6000: **574 ns/day at 4 cores, 785 at
12.** Hence `CORES_PER_REPLICA = 12` and `PACK_CPUS = 24`.

Swap in a three-site system (amber19/OPC3) and `UPDATE_MODE = 'gpu'` with
`CORES_PER_REPLICA = 4` becomes the right setting; the mdp keeps v-rescale and
C-rescale precisely so that switch needs no mdp edit.

Those are the two arms, and they have to stay separate rather than be averaged
into one campaign: 12 cores per replica against 4 is a different allocation
shape and a slower replica, and because packing's speed loss is measured against
that baseline, it also changes how favourable packing looks. All of this was
settled in `sampling-trpcage`; here the system is only a framework for testing
the code, and no performance work belongs in this example.

## Setting it up

One conda environment for mdfarmer, and GROMACS from somewhere else entirely:

```bash
mamba create -n mdfarmer -c conda-forge python=3.12 loos mdtraj
mamba activate mdfarmer
pip install -e /path/to/mdfarmer      # editable: this checkout is what runs
```

No `gromacs` in that list: `gmx` reaches the compute node through `ENV_SETUP`,
which the job script sources there, and the tender never runs MD itself.

No `openmm` either. This arm never uses OpenMM to simulate, and it no longer
needs it to import: the handful of helpers that read a serialised state or
inspect a platform ask for OpenMM when they are called, so `import mdfarmer`
and a whole GROMACS campaign run without one installed. Call an OpenMM-only
entry point on this environment and it raises an ImportError naming the package
rather than failing somewhere further in.

`loos` and `mdtraj` are both needed, because the harvest picks between them per
trajectory: LOOS for a rectangular box, mdtraj for a triclinic one it cannot
represent.

Install **editable**. It puts a pointer to your checkout in the environment
rather than a copy, so the code you edit is the code the compute node runs. A
regular install or a `PYTHONPATH` entry work too, but then `import mdfarmer` can
quietly resolve to a different tree than the one you are reading — which is why
`--check` prints the file it landed on. Read that line.

Then name the environment when you launch:

```bash
./drive_gmx.sh --env mdfarmer --check
```

Nothing has to be installed on the compute node. `sbatch` exports the tender's
environment by default, and independently of that the job script pins the
absolute interpreter — `sys.executable` as the tender saw it — and refuses with
`MISSING INTERPRETER` rather than running the wrong python.

### What is site-specific

Set for Flatiron's `rusty`. These are the first things to change elsewhere, all
constants at the top of `farmer.py`:

| constant | here | what to check |
|---|---|---|
| `ENV_SETUP` | `module load modules/2.4-20250724 openmpi/cuda-4.1.8 gromacs/mpi-2024.4` | whatever puts a CUDA GROMACS on `PATH` for you; `--check` proves it in a subshell before submitting anything |
| `GMX_BIN` | `gmx_mpi` | `gmx` on a thread-MPI build. `gmx_pack` probes for `-ntmpi` rather than assuming, so either is fine |
| `PARTITION` | `gpu` | a partition your account can submit to |
| `GRES` | `gpu:rtx_pro_6000_blackwell:1` | `sinfo -o '%P %G'` for the names your cluster uses |
| `QOS` | unset | some sites require one |
| `EXTRA_SBATCH` | empty | account/reservation lines, if your site wants them |
| `PACK_CPUS` | 24 | cores per pack, split across its replicas |
| `HARVEST_PARTITION` | `ccb` | any CPU-only partition |
| `WALLTIME` | `00:20:00` | fine for seconds of MD; raise for a real campaign |

Packing two replicas onto one card only pays with an MPS daemon, which the job
script starts per job. Whether it pays *at all* is system-specific — see
`plans/NOTE-gpu-packing-scaling.md` and `tests/test_pack_scaling.py`.

## Run it

```bash
cd examples/gromacs-trpcage

./drive_gmx.sh --env mdfarmer --check      # run shape and input readiness; writes nothing
./drive_gmx.sh --env mdfarmer --dry-run    # every directory, config, pack.json and sbatch.sh; submits nothing
./drive_gmx.sh --env mdfarmer              # start the tender, detached, and return
./drive_gmx.sh --env mdfarmer --status     # up or down, its pid, the tail of its log
./drive_gmx.sh --env mdfarmer --stop       # brake it at its next tick
```

`drive_gmx.sh` activates the environment you name and runs the driver as a
direct child, appending to `shakedown-gmx.tend.out`, whose path it prints on the
way out. Direct rather than under `mamba run` on purpose: the loop reads the
driver's exit status to decide whether to re-enter, and you read its log live,
and a wrapper process sits in the way of both. `--env` names the environment;
`CONDA_ENV` in the environment does the same. `farmer.py` still runs perfectly
well by hand; the script is what makes it survivable.

**The tender is not a Slurm job.** It runs detached — `setsid nohup` — on the
login node or workstation you launch it from, and outlives the shell that
started it. It sleeps between ticks and needs nothing from the cluster but
`sbatch`, so an allocation of its own would idle for hours; worse, that
allocation's walltime or its preemption would end the campaign with it. Launch
it anywhere `sbatch` and your environment both work.

**One tender per campaign.** The loop holds `flock` on
`data/shakedown-gmx/tender.lock` for as long as it lives, and a second
`./drive_gmx.sh` refuses with the running one's pid rather than starting a rival
that would submit every pack a second time. (`pack.lock` inside the job is the
same idea one level down: it is what stops two mdruns entering one pack
directory. Neither substitutes for the other.) The kernel drops the lock when the
process dies, however it dies, so there is no stale pid file to reason about.

**It re-enters.** `Farmer.launch` drops a pack for good on a single transient
`sbatch` failure, and a fresh tender rebuilds every pack from disk and re-adopts
the job ids still running, so re-entering is the recovery. The loop does that
every 60 s (`GAP`) until the driver exits 0, which happens only when every clone
has finished. Three exits inside a minute in a row is a broken setup rather than
a scheduler hiccup, and the loop says so and gives up.

Before it detaches, the script checks that `sbatch` and an importable
`mdfarmer` are both there, prints which `mdfarmer` (in a git worktree that is not
necessarily the tree you are reading), and confirms in a **subshell** that
`module load modules/2.4-20250724 openmpi/cuda-4.1.8 gromacs/mpi-2024.4` really
yields a `gmx_mpi`: a campaign whose every generation would die at `mdrun` should
fail on the login node in one second, not across 25 pack jobs. The subshell is
deliberate — the tender never calls `gmx` itself, and leaving the module tree in
its environment would put module libraries ahead of the conda ones in every job
it submits. On the node it is `ENV_SETUP` at the top of `farmer.py` that loads
them, and the job script refuses to run if `gmx_mpi` or `mdfarmer` is missing.

`./drive_gmx.sh --stop` writes the `stop` brake file the driver watches for; the
tender exits at its next tick (20 s), the loop then exits too, and pack jobs
already submitted keep running. Everything a run writes lands in `data/` and
`prepared/`, both gitignored.

## What to expect

| | |
|---|---|
| steps per generation | 10 000 (20 ps at dt = 2 fs) |
| write interval | 2 000 steps (4 ps) → **5 new frames per generation** |
| frames on disk | 6 — GROMACS also writes the restart-step frame |
| downsample | every 5th frame → 1 wet frame per generation, 2 from generation 0 |
| generations | 5, so 100 ps and **26** dry frames per clone — see below |
| MD time per generation | **3.5–3.8 s** per replica, two sharing an RTX A6000 without MPS (454 and 498 ns/day, against 785 solo); 11.6 s for the whole pack job including grompp and start-up |
| whole campaign | 25 pack jobs of about a minute each; wall time is queue time |

The two step counts are not arbitrary. `utilities.check_whole_frames` requires
`10000 % 2000 == 0`, so the last frame of each generation lands exactly on the
checkpoint the next one restarts from; `harvester.check_commensurability`
additionally requires `(10000 / 2000) % 5 == 0`, so the downsampled stream keeps
its spacing across a generation boundary. Five is the *new*-frame count, which
is the one both rules are stated in.

**Twenty-six, not twenty-five.** Every generation holds 6 frames on disk, and
every generation after the first drops one as a duplicate of its predecessor's
last — but generation 0 has no predecessor, so its step-0 frame (the seed state)
is kept and belongs to no other generation. The total is `n_gens * 5 + 1`.
`first_global_index` counts frames *on disk*, duplicate included, so a global
index is always `step / write_interval` and the kept stream stays contiguous:
generation 0 supplies indices 0–5, generation 1 supplies 6–10, and its own
index 5 is the duplicate it drops. The OpenMM arm has no such frame — its
reporters write nothing at step 0 of a generation — so the two engines differ by
exactly one frame per clone, by design.

`maxh = 0.25` inside a 20-minute block: mdrun's own backstop, well clear of a
generation that takes seconds, so a generation ends at its step target and not
at the clock.

## How to tell it worked

**Frame counts.** Generation 0 writes 6 frames (0, 4, 8, 12, 16, 20 ps) and so
does every later generation, because GROMACS writes a frame at the step it
restarted from. The harvest drops that duplicate exactly once — the `.harvested`
sentinel says which way it went:

```
gen_00  n_orig 6  n_dry 6  n_down 2  first_global_index 0  skip_first false
gen_01  n_orig 6  n_dry 5  n_down 1  first_global_index 5  skip_first true
```

**A contiguous step axis, checked for you.** The dry stream is XTC and carries
per-frame times, so:

```python
gens = sorted((Path('data/shakedown-gmx/seed_00/clone_00')).glob('gen_*'))
mdf.verify_dry_chain(gens)
# {'n_frames': 26, 'expected': 26, 'contiguous': True,
#  'strictly_increasing': True, 'n_duplicate_times': 0,
#  'uniform_spacing': True, 'first_time': 0.0, 'last_time': 100.0}
```

Over two generations that reads `n_frames: 11, first_time 0.0, last_time 40.0` —
five new frames each plus generation 0's frame at step 0, no duplicate at the
seam. Generations chain through `gmx convert-tpr -nsteps` plus `mdrun -cpi`, so
the step and time counters are globally continuous by construction;
`verify_dry_chain` is what proves the harvest did not break that.

Per generation, `gen_status.json` records the same thing more cheaply:

```json
{"target_step": 10000, "reached_step": 10000, "complete": true}
```

with `target_step` climbing 10000, 20000, … across generations.

**The pack ran as a pack.** In `data/shakedown-gmx/packs/pack_s00c00_s00c01/`:
`pack.json` names both members' configs, `slurm-<jobid>.out` says
`MPS: daemon up`, and the log shows `Applying core pinning offset 12` for the
second replica — without that, both mdruns start at core 0 and fight, which
reads as node variance rather than a misconfiguration. `pack.lock` is what stops
a second job ever entering that directory.

**The harvest.** Each finished generation should end up with `dry_traj.xtc`
(304 atoms), `downsample_traj.xtc`, a `.harvested` sentinel, and `traj.xtc`
turned into a symlink to the dry stream. Anything the tender skipped:

```python
mdf.unharvested_gen_dirs('data/shakedown-gmx', skip_newest=True)
```

The harvest is triggered by the tender noticing a generation finish, so a
generation that completed while the tender was down is never harvested;
`unharvested_gen_dirs` is how you find it.

## When it does not work

A packed job writes one `slurm-<jobid>.out` in the pack directory, and the
cleanup trap copies it into each member's generation directory as `slurm.out`
so `BadNodeRegistry` can still scan it. Read that first, then `traj.part0001.log`
(mdrun's own log) and `config.json` in the generation directory. Failure is
per-member: one replica raising fails exactly one clone and leaves the other
running, so expect the tender to carry on with a short pack.

## A known blemish in the inputs

`gmx.gro` carries a 4.837 nm box with coordinates that were equilibrated in a
4.683 nm one, so this arm starts about 10% too low in density and the C-rescale
barostat spends the first few picoseconds compressing back. The OpenMM arm reads
its box from `state.xml` and does not. It is harmless for a pipeline shakedown —
nothing here is science — but it does mean the two arms are not bitwise the same
starting state, and it is worth fixing in the inputs before anyone reuses this
`.gro` for anything real.
