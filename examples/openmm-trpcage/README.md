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
13465 atoms in a 4.68 nm cube, at 277 K under a Monte Carlo barostat. It is the
same system `examples/gromacs-trpcage` runs, so the two datasets are comparable
frame for frame.

`inputs/` holds it gzipped — `system.xml.gz`, `state.xml.gz`,
`topology.pdb.gz`, 1.07 MB for 6.6 MB of XML. `prepare_inputs()` inflates them
into `prepared/` before the Farmer is built. It has to: `omm_generation` reads
the system with `Path(system_fn).read_text()` and hands the seed to
`Simulation.loadState`, and neither inflates. (`XmlSerializer.deserialize` will
happily take `gzip.open(p, 'rt').read()` — `tests/test_example_inputs.py`
checks that it round-trips — but nothing in the runner does that today.)

`prepare_inputs()` also rewinds the seed state's `stepCount` from 150000 to 0.
That number is left over from the equilibration this system came from, and
`seeder._try_recover_gen` reads `state.xml`'s `stepCount` as a step counted from
the start of the campaign — so a seed carrying it would make every *resumed*
generation look finished before it started. A clean first run would not notice.

## Run it

```bash
PY=/mnt/home/lsmith/miniforge3/envs/omm/bin/python
cd examples/openmm-trpcage

$PY farmer.py --check      # run shape and input readiness; writes nothing
$PY farmer.py --dry-run    # every directory, config and sbatch.sh; submits nothing

# for real:
nohup $PY -u farmer.py > shakedown-omm.tend.out 2>&1 &
tail -f shakedown-omm.tend.out
```

`mdfarmer` must be importable by `$PY` — the submit script checks and refuses
rather than failing on the node. Slurm exports the tender's environment, so
launching from a shell where `python -c 'import mdfarmer'` works is enough.

To stop the tender gracefully, `touch stop` in the directory you launched it
from; it exits at the next tick. Everything it writes lands in `data/` and
`prepared/`, both gitignored.

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
