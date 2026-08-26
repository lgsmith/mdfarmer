"""Shared scaffolding for the mdfarmer suites.

Each suite is a plain script: run it directly, or run them all with
`python tests/run_all.py`. Exit status is the result, and every check prints a
line, so a failure names itself without a framework.

The GROMACS-backed suites build their system from the water box and force field
that ship with GROMACS itself, located through `gmx -version`. Nothing is
vendored, and a machine without a `gmx` on PATH skips those suites rather than
failing them.
"""
import os
import shutil
import subprocess as sp
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
# Import mdfarmer as a package from its parent directory, so a suite runs
# against this checkout rather than whatever is installed.
sys.path.insert(0, str(REPO_ROOT.parent))

# Override to keep test output somewhere durable; /tmp is periodically pruned.
SCRATCH_ROOT = Path(os.environ.get(
    'MDFARMER_TEST_DIR', Path(tempfile.gettempdir()) / 'mdfarmer-tests'))

# Sites that build an MPI-only GROMACS have `gmx_mpi`.
GMX_BIN = os.environ.get('GMXBIN', 'gmx')

# GROMACS' own 216-molecule SPC water box, and the force field that types it.
WATER_BOX_NAME = 'spc216.gro'
WATER_FORCE_FIELD = 'oplsaa.ff'
WATER_MOLECULES = 216
WATER_ATOMS = WATER_MOLECULES * 3

# The box is 1.86 nm on a side, so the pair-list cutoff has to stay under half
# of that or grompp refuses to build a tpr.
WATER_CUTOFF_NM = 0.7
WATER_NSTLIST = 10

# Line `gmx -version` prints for its install prefix.
GMX_DATA_PREFIX_KEY = 'Data prefix:'


class Suite:
    """Counts checks and reports them, so a failure prints where it happened."""

    def __init__(self, name):
        self.name = name
        self.results = []

    def check(self, name, condition, detail=''):
        self.results.append(bool(condition))
        print(f'  [{"PASS" if condition else "FAIL"}] {name} {detail}',
              flush=True)
        return bool(condition)

    def section(self, title):
        print(f'\n=== {title} ===', flush=True)

    def report(self):
        passed = sum(self.results)
        total = len(self.results)
        print(f'\n{self.name}: {passed}/{total} checks passed', flush=True)
        return 0 if passed == total else 1


class Skip(Exception):
    """Raised when a suite's dependencies are absent, not when it fails."""


def gmx_data_prefix(gmx_bin=GMX_BIN, key=GMX_DATA_PREFIX_KEY):
    """The install prefix `gmx` reports, whose share/gromacs/top holds the
    standard structures and force fields."""
    try:
        result = sp.run([gmx_bin, '-version'], capture_output=True, text=True,
                        timeout=120)
    except (FileNotFoundError, sp.TimeoutExpired, OSError) as exc:
        raise Skip(f'{gmx_bin} is not runnable: {exc}')
    for line in (result.stdout + result.stderr).splitlines():
        if line.strip().startswith(key):
            return Path(line.split(key, 1)[1].strip())
    raise Skip(f'{gmx_bin} -version printed no {key!r} line')


def require_gmx(gmx_bin=GMX_BIN):
    """The GROMACS top directory, or Skip if this machine has no usable gmx."""
    top = gmx_data_prefix(gmx_bin=gmx_bin) / 'share' / 'gromacs' / 'top'
    if not (top / WATER_BOX_NAME).is_file():
        raise Skip(f'{top / WATER_BOX_NAME} is missing')
    return top


def workdir(name, scratch_root=SCRATCH_ROOT):
    """A fresh, empty directory for one suite."""
    path = Path(scratch_root) / name
    shutil.rmtree(path, ignore_errors=True)
    path.mkdir(parents=True)
    return path


def water_mdp(cutoff_nm=WATER_CUTOFF_NM, nstlist=WATER_NSTLIST):
    """An NPT mdp for the water box: Nose-Hoover and Parrinello-Rahman, which
    is the coupling whose state a generation chain has to carry across."""
    return '\n'.join((
        'integrator               = md',
        'dt                       = 0.002',
        'nsteps                   = 1000',
        'nstlog                   = 100',
        'nstenergy                = 100',
        'nstxout-compressed       = 100',
        'cutoff-scheme            = Verlet',
        f'nstlist                  = {nstlist}',
        'coulombtype              = PME',
        f'rcoulomb                 = {cutoff_nm}',
        'vdw-type                 = cut-off',
        f'rvdw                     = {cutoff_nm}',
        'tcoupl                   = nose-hoover',
        'tc-grps                  = System',
        'tau-t                    = 0.5',
        'ref-t                    = 300',
        'pcoupl                   = Parrinello-Rahman',
        'pcoupltype               = isotropic',
        'tau-p                    = 2.0',
        'ref-p                    = 1.0',
        'compressibility          = 4.5e-5',
        'constraints              = h-bonds',
        'constraint-algorithm     = lincs',
        'continuation             = no',
        'gen-vel                  = yes',
        'gen-temp                 = 300',
        'gen-seed                 = 42',
        '',
    ))


def water_top(force_field=WATER_FORCE_FIELD, molecules=WATER_MOLECULES):
    return '\n'.join((
        f'#include "{force_field}/forcefield.itp"',
        f'#include "{force_field}/spc.itp"',
        '',
        '[ system ]',
        'spc water box',
        '',
        '[ molecules ]',
        f'SOL   {molecules}',
        '',
    ))


def build_water_system(dest, gmx_bin=GMX_BIN, water_box_name=WATER_BOX_NAME):
    """Write conf.gro / topol.top / base.mdp into `dest`.

    Returns the three paths. Raises Skip when GROMACS is unavailable, so a
    suite can report "skipped" rather than "failed" on a machine without it.
    """
    top_dir = require_gmx(gmx_bin=gmx_bin)
    dest = Path(dest)
    structure = dest / 'conf.gro'
    topology = dest / 'topol.top'
    mdp = dest / 'base.mdp'
    shutil.copy(top_dir / water_box_name, structure)
    topology.write_text(water_top())
    mdp.write_text(water_mdp())
    return structure, topology, mdp


def run_suite(main):
    """Run a suite's `main(suite)`, turning Skip into exit status 77."""
    try:
        return main()
    except Skip as exc:
        print(f'SKIPPED: {exc}', flush=True)
        return 77
