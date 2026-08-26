# tests

Plain scripts, no test framework. Run one:

```bash
python tests/test_pack_farmer.py
```

or all of them:

```bash
python tests/run_all.py
```

Each check prints a `[PASS]`/`[FAIL]` line naming what it checked, and the exit
status is the result. A suite exits `77` when its dependencies are missing,
which `run_all.py` reports as `skip` rather than a failure.

## What needs what

| suite | needs |
|---|---|
| `test_pack_farmer.py` | nothing beyond mdfarmer |
| `test_pack_wiring.py` | nothing beyond mdfarmer |
| `test_pack_findings.py` | a working `gmx` (runs two concurrent mdruns) |
| `test_harvest.py` | a working `gmx`, plus LOOS and mdtraj |

`harness.py` locates GROMACS through `gmx -version` and builds the test system
from the 216-molecule SPC water box and force field that ship with it, so
nothing is vendored here. Set `GMXBIN` if the binary is not called `gmx`:

```bash
GMXBIN=gmx_mpi python tests/run_all.py
```

Output goes to `$MDFARMER_TEST_DIR`, or a `mdfarmer-tests` directory under the
system temp directory. Set `MDFARMER_TEST_DIR` to keep it somewhere durable —
`/tmp` is periodically pruned, and the directories are worth reading after a
failure.

## What is not covered

Nothing here submits to a scheduler. The suites drive `check_start_gen` with
`dry_run=True`, so they cover everything up to and including the generated
submit script, and nothing after it. MPS itself is likewise unexercised: the
pack runs two real mdruns concurrently, but on whatever device the test machine
has, without an MPS daemon.
