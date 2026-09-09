# shakedown-omm — 10 solo OpenMM clones of trp-cage

A pre-flight for the whole pipeline: seeding, generation chaining, the
preemption handshake, and a harvest per finished generation, on real MD that
takes seconds rather than weeks. Launch it before you commit a real campaign to
this cluster, and you find out about the missing interpreter, the wrong
`PYTHONPATH` or the partition that will not give you a card in about five
minutes instead of two days.

**10 clones from one seed, 5 running at a time, 5 generations each.** One card
per clone, no packing — the GROMACS example next door is the packed arm.

## The system

Trp-cage (20 residues, 304 atoms) in 3287 four-site waters with 6 K⁺ and 7 Cl⁻,
13465 atoms in a 4.68 nm cube, at 277 K under a Monte Carlo barostat, in
**amber03 + TIP4P-ice** (`sampling-trpcage/systems/fresh-omm/native-277`, whose
`meta.json` marks the amber03 as benchmark-only). It is the same system
`examples/gromacs-trpcage` runs, so the two datasets are comparable frame for
frame.

`inputs/` holds it gzipped — `system.xml.gz`, `state.xml.gz`,
`topology.pdb.gz`, 1.07 MB for 6.6 MB of XML. `prepare_inputs()` inflates them
into `prepared/` before the Farmer is built. It has to: `omm_generation` reads
the system with `Path(system_fn).read_text()` and hands the seed to
`Simulation.loadState`, and neither inflates. (`XmlSerializer.deserialize` will
happily take `gzip.open(p, 'rt').read()` — `tests/test_example_inputs.py`
checks that it round-trips — but nothing in the runner does that today.)

Inflation is a one-off, not a per-boot cost: `inflate()` returns immediately if
the inflated file is already in `prepared/`, and writes through a `.partial`
name so a boot killed halfway leaves nothing that a later boot could mistake for
a finished file. Delete `prepared/` to force it again. Nothing in either example
compresses trajectories or restarts *during* a run — the DCD and the `state.xml`
restarts are written plain, and the only lossy compression anywhere is the XTC
format the GROMACS arm writes natively.

`prepare_inputs()` also rewinds the seed state's `stepCount` from 150000 to 0.
That number is left over from the equilibration this system came from, and
`seeder._try_recover_gen` reads `state.xml`'s `stepCount` as a step counted from
the start of the campaign — so a seed carrying it would make every *resumed*
generation look finished before it started. A clean first run would not notice.

### Which arm this is

This is the **TIP4P-ice (4-site) arm**, the same one `examples/gromacs-trpcage`
runs, and the two physics arms have to stay separate because of what the water
does to the GROMACS side: one virtual site per water means GROMACS refuses
`mdrun -update gpu`, the update runs on the CPU, and the per-replica knee moves
from ~4 cores to 12. An **amber19 + OPC3** (3-site) build is the other arm:
`-update gpu`, ~4 cores per replica, faster per card — and because packing's
speed loss is measured against that faster baseline, it also changes how
favourable packing looks. Mixed into one campaign the two would need different
core budgets per replica and could not be compared.

OpenMM has no `-update gpu` switch to lose — the CUDA platform integrates and
places virtual sites on the card either way — so the four-site water costs this
arm nothing: `CPUS = 2`, one busy core and a spare, whichever water it is. That
asymmetry is why the GROMACS arm asks for 24 cores a pack and this one asks for
2 a clone. All of it was settled in `sampling-trpcage`; here the system is only
a framework for testing the code, and no performance work belongs in either
example.

## Setting it up

One conda environment, holding mdfarmer and everything it imports:

```bash
mamba create -n mdfarmer -c conda-forge python=3.12 openmm loos mdtraj
mamba activate mdfarmer
pip install -e /path/to/mdfarmer      # editable: this checkout is what runs
```

`openmm` runs the MD; `loos` and `mdtraj` are both needed because the harvest
picks between them per trajectory — LOOS for a rectangular box, mdtraj for a
triclinic one it cannot represent.

Install **editable**. It puts a pointer to your checkout in the environment
rather than a copy, so the code you edit is the code the compute node runs.
A regular install or a `PYTHONPATH` entry work too, but then `import mdfarmer`
can quietly resolve to a different tree than the one you are reading — which is
why `--check` prints the file it landed on. Read that line.

Then name the environment when you launch:

```bash
./drive_omm.sh --env mdfarmer --check
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
| `PARTITION` | `gpu` | a partition your account can submit to |
| `GRES` | `gpu:rtx_pro_6000_blackwell:1` | `sinfo -o '%P %G'` for the names your cluster uses |
| `QOS` | unset | some sites require one |
| `EXTRA_SBATCH` | empty | account/reservation lines, if your site wants them |
| `HARVEST_PARTITION` | `ccb` | any CPU-only partition |
| `WALLTIME` | `00:20:00` | fine for seconds of MD; raise for a real campaign |

`STEPS_PER_GEN` is calibrated so a generation is a few seconds on an RTX A6000.
On a slower card it is still short. The two spacing rules it satisfies are
arithmetic, not tuning — see the comment above them before changing it.

## Run it

```bash
cd examples/openmm-trpcage

./drive_omm.sh --env mdfarmer --check      # run shape and input readiness; writes nothing
./drive_omm.sh --env mdfarmer --dry-run    # every directory, config and sbatch.sh; submits nothing
./drive_omm.sh --env mdfarmer              # start the tender, detached, and return
./drive_omm.sh --env mdfarmer --status     # up or down, its pid, the tail of its log
./drive_omm.sh --env mdfarmer --stop       # brake it at its next tick
```

`drive_omm.sh` activates the environment you name and runs the driver as a
direct child, appending to `shakedown-omm.tend.out`, whose path it prints on the
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
`data/shakedown-omm/tender.lock` for as long as it lives, and a second
`./drive_omm.sh` refuses with the running one's pid rather than starting a rival
that would submit every clone a second time. The kernel drops the lock when the
process dies, however it dies, so there is no stale pid file to reason about.

**It re-enters.** `Farmer.launch` drops a clone for good on a single transient
`sbatch` failure, and a fresh tender rebuilds every clone from disk and re-adopts
the job ids still running, so re-entering is the recovery. The loop does that
every 60 s (`GAP`) until the driver exits 0, which happens only when every clone
has finished. It re-enters on exit 1 and not on exit 2: the driver reports
"asked to stop" and "died" as different statuses, so the loop reads the status
and forms no opinion about the brake file itself. Three exits inside a minute in a row is a broken setup rather than
a scheduler hiccup, and the loop says so and gives up.

Before it detaches, the script checks that `sbatch` and an importable
`mdfarmer` are both there, and prints which `mdfarmer` — `import mdfarmer`
resolves to whatever the environment installed, which in a git worktree is not
necessarily the tree you are reading. On the node, the submit script makes the
same check and refuses rather than failing mid-generation.

`./drive_omm.sh --stop` writes the `stop` brake file the driver watches for; the
tender exits at its next tick (20 s), the loop then exits too, and jobs already
submitted keep running. Everything a run writes lands in `data/` and `prepared/`,
both gitignored.

## What to expect

| | |
|---|---|
| steps per generation | 10 000 (20 ps at dt = 2 fs) |
| write interval | 2 000 steps (4 ps) → **5 frames per generation** |
| downsample | every 5th frame → 1 wet frame per generation |
| generations | 5, so 100 ps and 25 dry frames per clone |
| MD time per generation | **3.1 s** measured on an RTX A6000 (450 ns/day); 5.7 s including python and CUDA start-up |
| whole campaign | 50 jobs of about a minute each; wall time is queue time |

The two step counts are not arbitrary. `utilities.check_whole_frames` requires
`10000 % 2000 == 0`, so the last frame of each generation lands exactly on the
checkpoint the next one restarts from; `harvester.check_commensurability`
additionally requires `(10000 / 2000) % 5 == 0`, so the downsampled stream keeps
its spacing across a generation boundary. Five frames is the smallest count that
still shows a step axis marching across a seam, and 5 is the largest downsample
that divides it.

## How to tell it worked

**Frame counts.** Every generation's `traj.dcd` holds exactly 5 frames — OpenMM
does not write a frame at the step it restarts from, so there is no seam
duplicate to drop. After the harvest, the `.harvested` sentinel in each
generation directory records the plan and the counts it verified:

```json
{"n_orig": 5, "n_dry": 5, "n_down": 1, "first_global_index": 5,
 "frames_per_gen": 5, "skip_first": false, "n_subset_atoms": 304}
```

`first_global_index` climbing 0, 5, 10, 15, 20 across the five generations is
the thing to look at: it says the downsample is driven by a campaign-wide frame
index rather than resetting at each boundary.

**A contiguous step axis.** `traj.out` in each generation directory is the state
reporter's log. Read the `Step` column across a clone's generations and it must
march without a gap or a repeat:

```
gen_00   2000  4000  6000  8000 10000
gen_01  12000 14000 16000 18000 20000
gen_02  22000 ...
```

and `state.xml`'s `stepCount` at the end of generation *N* must be
`(N + 1) * 10000`.

`mdf.verify_dry_chain` does this check for you on the GROMACS arm, but **not
here**: it needs per-frame times, and DCD does not carry them. DCD is
deliberate — it is the only format `seeder._try_recover_gen` can trim when a
kill lands between the trajectory reporter and the checkpoint reporter, so an
interrupted generation resumes instead of being redone.

**The harvest.** Each finished generation should end up with `dry_traj.dcd`
(5 frames, 304 atoms), `downsample_traj.dcd` (1 frame, all 13465), a
`.harvested` sentinel, and `traj.dcd` turned into a symlink to the dry stream.
Anything the tender skipped:

```python
mdf.unharvested_gen_dirs('data/shakedown-omm', skip_newest=True)
```

Note that the harvest is triggered by the tender noticing a generation finish.
Restart the tender across a generation that completed while it was down and that
one is never harvested — `unharvested_gen_dirs` is how you find it.

## When it does not work

`slurm.out` in the generation directory first, then `config.json` beside it —
most failures here are path failures, and the config records every path the job
resolved. `node_history.tsv` accumulates one line per launch, so a clone that
keeps dying on one node is visible; `bad_nodes.txt` under `data/shakedown-omm`
is what the Farmer learned and will exclude.
