# Conventions for agents working on mdfarmer

mdfarmer runs MD campaigns that take months and must not be restarted. A bug
here does not crash; it produces a trajectory that looks fine and is not. Write
accordingly.

## Comments

Docstrings say what a thing is for and what it promises. `#` comments elucidate
the line they sit on: one line, two at the outside.

A comment earns its place by recording **why**, where the why is not recoverable
from the code:

```python
# GROMACS writes a frame at the step it restarts from; OpenMM does not.
if n_orig == frames_per_gen + 1:
```

Do not narrate the change you are making. No "fixed the bug where...", no
"previously this did X", no dated notes, no ticket numbers. A reader a year from
now wants the constraint, not your afternoon. Bug rationale belongs in the
commit message and, if it is load-bearing, in a docstring stated as a property
of the code rather than as history.

Do not restate the code. `# increment i` is noise.

## Python

Plain, legible Python. Prefer the obvious construction over the clever one.
If numpy or scipy already does it, use that rather than hand-rolling.

Constants go at the top of the file **and** are passed in as named parameters,
with the default routed to the top-of-file value:

```python
WRITE_INTERVAL = 2000

def frames(steps, write_interval=WRITE_INTERVAL):
    return steps // write_interval
```

Never reach a module global implicitly from inside a function body. Renames then
break loudly instead of going silently stale.

Fail loudly and specifically. A function that cannot do the right thing raises
with a message naming the file, the number, and what it expected. Never guess
past missing data, never silently substitute a default for something the caller
needed to be true. Refusing is cheap; a wrong trajectory is not.

## Commits

One commit, one change, one line. The subject is a terse lowercase clause in the
repo's existing voice; read `git log` before writing one.

**If the subject needs a second line, the commit needs splitting.** A message
that wants to say "and also" is two commits. No body, no `Co-Authored-By`, no
AI-attribution trailer — the README's "AI assistance" section is the
project-level attribution and covers all of it.

Keep unrelated edits out. A drive-by fix noticed while doing something else is
its own commit.

## Branches

Feature work happens on a branch named for the feature, not for who or what
made it. Several small commits on it, then merge to `main` with `--no-ff` so the
branch reads as one unit of work in the log. Delete the branch after merging.

Group work by theme rather than by file: a branch is a swathe of commits that
belong to one idea, even when only one thing is being developed at a time.

Do not push, delete remote branches, or rewrite published history without being
asked.

## Tests

`tests/` holds plain scripts, no framework. `python tests/run_all.py` runs them;
each reports `pass`, `fail` or `skip`. Read two or three existing suites before
adding one — match the module docstring that explains *why the test exists*, the
`Suite`/`suite.check` shape, and the constants-at-top style.

A test must cover the **call site**, not only the helper. A helper tested in
isolation while its caller does something else is the standard way a bug ships
here. When a check matters, mutate the code to confirm the test actually fails.

Set `MDFARMER_TEST_DIR` when running suites concurrently; the scratch root is
shared. `GMXBIN` names the GROMACS binary if it is not `gmx`.

## Engines

Every feature works on both the OpenMM and GROMACS paths. Lift shared logic into
one function rather than writing it twice. Where the two genuinely differ,
say why in a comment at the point of difference.

Neither engine may be a hard requirement of the other: someone who runs only
GROMACS must be able to install and use mdfarmer without OpenMM.

## Trajectories

LOOS first, mdtraj second — mdtraj for what LOOS gets wrong, principally
triclinic boxes. A third backend, or shelling out to a CLI tool, needs asking
about first. `gmx trjcat` and `trjconv` are not used here and are not to be
reintroduced.

## Environments

`mamba run -n "$(cat .conda-env)" python -u ...`. `--no-capture-output` is
broken on this system; `-u` is what keeps output live.

Do not create conda environments — `/mnt/home` enforces a per-user file-count
quota and env trees dominate it. Never run `mamba clean`.
