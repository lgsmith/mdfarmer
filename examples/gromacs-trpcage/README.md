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
13465 atoms, at 277 K under a C-rescale barostat.

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

### `-update cpu`, not `-update gpu`

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

Swap in a three-site system (a19/OPC3) and `UPDATE_MODE = 'gpu'` with
`CORES_PER_REPLICA = 4` becomes the right setting; the mdp keeps v-rescale and
C-rescale precisely so that switch needs no mdp edit.

## Run it

```bash
PY=/mnt/home/lsmith/miniforge3/envs/omm/bin/python
cd examples/gromacs-trpcage

$PY farmer.py --check      # run shape and input readiness; writes nothing
$PY farmer.py --dry-run    # every directory, config, pack.json and sbatch.sh; submits nothing

# for real:
nohup $PY -u farmer.py > shakedown-gmx.tend.out 2>&1 &
tail -f shakedown-gmx.tend.out
```

The job script loads GROMACS itself (`ENV_SETUP` at the top of `farmer.py`;
`modules/2.4-20250724 openmpi/cuda-4.1.8 gromacs/mpi-2024.4`, giving `gmx_mpi`)
and refuses to run if the binary or `mdfarmer` is missing, rather than failing
halfway through a generation. The tender itself needs `mdfarmer` importable by
`$PY`; Slurm carries that environment to the node.

`touch stop` in the launch directory to stop the tender at the next tick.
Everything it writes lands in `data/` and `prepared/`, both gitignored.

## What to expect

| | |
|---|---|
| steps per generation | 10 000 (20 ps at dt = 2 fs) |
| write interval | 2 000 steps (4 ps) → **5 new frames per generation** |
| frames on disk | 6 — GROMACS also writes the restart-step frame |
| downsample | every 5th frame → 1 wet frame per generation |
| generations | 5, so 100 ps and 25 dry frames per clone |
| MD time per generation | **3.5–3.8 s** per replica, two sharing an RTX A6000 without MPS (454 and 498 ns/day, against 785 solo); 11.6 s for the whole pack job including grompp and start-up |
| whole campaign | 25 pack jobs of about a minute each; wall time is queue time |

The two step counts are not arbitrary. `utilities.check_whole_frames` requires
`10000 % 2000 == 0`, so the last frame of each generation lands exactly on the
checkpoint the next one restarts from; `harvester.check_commensurability`
additionally requires `(10000 / 2000) % 5 == 0`, so the downsampled stream keeps
its spacing across a generation boundary. Five is the *new*-frame count, which
is the one both rules are stated in.

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
