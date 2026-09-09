"""A GROMACS-only site must be able to import the package with no OpenMM.

GROMACS is reached as a subprocess, so nothing on that path imports its engine,
and someone who runs only GROMACS campaigns should never have to install an
OpenMM to use this package. The asymmetry is real and one-directional: an
OpenMM-only site is already fine, because `gmx` is a binary, not an import.

Proving it needs an environment without OpenMM, and building one would cost a
whole conda environment. Instead this suite re-runs itself in a subprocess with
a meta_path finder that makes `import openmm` fail exactly the way an
uninstalled OpenMM fails, imports the checkout there, and drives a GROMACS
campaign up to its submit script. The OpenMM-only entry points are checked in
the same subprocess: each has to raise an ImportError naming what to install,
rather than an AttributeError about a name that quietly went missing.
"""
import importlib.util
import json
import subprocess as sp
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
# argv marker for the run of this file that has OpenMM blocked.
CHILD_FLAG = '--no-openmm-child'
# The child prints one JSON line behind this, so its own stdout stays readable.
RESULT_MARKER = 'RESULTS:'
BLOCKED_MODULE = 'openmm'
CHILD_TIMEOUT = 300

N_SEEDS = 1
N_CLONES = 1
N_GENS = 2
STEPS_PER_GEN = 1000
WRITE_INTERVAL = 100
CPUS = 8
# What plow_harrow_plant has to leave in a generation directory.
GEN_ARTEFACTS = ('config.json', 'run.py')
# An equilibrated OpenMM seed's step count, read without OpenMM.
SEED_STEP_COUNT = 150000


class OpenMMBlocker:
    """A meta_path finder that makes one package fail to import.

    Raising from find_spec rather than returning None, so the failure is a
    ModuleNotFoundError naming the module even if it is installed further down
    sys.path -- the point is to reproduce a site that never installed it.
    """

    def __init__(self, blocked=BLOCKED_MODULE):
        self.blocked = blocked

    def find_spec(self, fullname, path=None, target=None):
        if fullname == self.blocked or fullname.startswith(self.blocked + '.'):
            raise ModuleNotFoundError(f'No module named {fullname!r}',
                                      name=fullname)
        return None


def import_checkout(repo_root=REPO_ROOT):
    """Import the checkout this file lives in, under the name mdfarmer.

    By path rather than by sys.path, for the reason harness._import_checkout
    gives: in a git worktree the directory is not named mdfarmer and the import
    falls through to whatever is installed.
    """
    spec = importlib.util.spec_from_file_location(
        'mdfarmer', repo_root / '__init__.py',
        submodule_search_locations=[str(repo_root)])
    module = importlib.util.module_from_spec(spec)
    sys.modules['mdfarmer'] = module
    spec.loader.exec_module(module)
    return module


def raised(call):
    """The exception call() raised, or None."""
    try:
        call()
    except Exception as exc:
        return exc
    return None


def names_openmm(exc):
    """Whether an exception is an ImportError that says OpenMM is the problem."""
    return isinstance(exc, ImportError) and 'OpenMM' in str(exc)


def seed_inputs(work, n_seeds=N_SEEDS):
    """Placeholder structure, mdp and top per seed: a dry run reads none of them."""
    work.mkdir(parents=True, exist_ok=True)
    structures, systems, tops = [], [], []
    for seed_index in range(n_seeds):
        for name, fns in ((f'seed{seed_index}.gro', structures),
                          (f'seed{seed_index}.mdp', systems),
                          (f'seed{seed_index}.top', tops)):
            path = work / name
            path.write_text('placeholder\n')
            fns.append(str(path))
    return structures, systems, tops


def gromacs_farmer(mdfarmer, work, n_seeds=N_SEEDS, n_clones=N_CLONES,
                   n_gens=N_GENS, steps_per_gen=STEPS_PER_GEN,
                   write_interval=WRITE_INTERVAL, cpus=CPUS):
    """A dry-run GROMACS campaign over placeholder inputs: no gmx, no scheduler."""
    structures, systems, tops = seed_inputs(work / 'inputs', n_seeds=n_seeds)
    template = mdfarmer.gmx_config_template(
        traj_dir_top_level=str(work / 'farm'), title='nomm',
        structure_fn=structures[0], mdp_fn=systems[0],
        dirname_pad=2, sep='_', traj_name='prod', traj_suffix='.xtc',
        restart_name='state.cpt', steps=steps_per_gen,
        steps_per_gen=steps_per_gen, write_interval=write_interval,
        temperature=300, gen_seed_base=1, mdrun_args=['-nb', 'gpu'],
        traj_list=str(work / 'tl.txt'))
    return mdfarmer.Farmer(
        n_seeds=n_seeds, n_clones=n_clones, n_gens=n_gens,
        config_template=template,
        seed_structure_fns=structures, system_fns=systems, top_fns=tops,
        scheduler='sbatch',
        scheduler_fstring=mdfarmer.basic_scheduler_fstrings['slurm'],
        scheduler_kws=dict(gpu_line='', queue_name='gpu', exclude_nodes='',
                           cpus=cpus, run_script_name='run.py'),
        scheduler_report_cmd='true', scheduler_assoc_rep_cmd='true',
        sep='_', dirname_pad=2, runner=mdfarmer.gmx_generation, dry_run=True,
        overwrite=True, jids_file=work / 'jids.txt')


def state_xml(path, step_count=SEED_STEP_COUNT):
    """A state.xml with the one attribute the step reader looks at."""
    path.write_text(f'<?xml version="1.0" ?>\n<State openmmVersion="8.5" '
                    f'time="0.0" stepCount="{step_count}">\n</State>\n')
    return path


def probe_import(results):
    """The package imports at all, and exports what it promises to."""
    results.append(['import openmm really fails here',
                    isinstance(raised(lambda: __import__(BLOCKED_MODULE)),
                               ModuleNotFoundError), ''])
    mdfarmer = import_checkout()
    missing = [name for name in mdfarmer.__all__
               if not hasattr(mdfarmer, name)]
    results.append(['__all__ is honest without OpenMM', not missing,
                    f'-> {missing}'])
    results.append(['a star import works',
                    raised(lambda: exec('from mdfarmer import *', {})) is None,
                    ''])
    return mdfarmer


def probe_gromacs(results, mdfarmer, work, gen_artefacts=GEN_ARTEFACTS):
    """A GROMACS campaign boots and writes a generation's job."""
    exc = raised(lambda: gromacs_farmer(mdfarmer, work))
    results.append(['a GROMACS Farmer boots', exc is None, f'-> {exc!r}'[:90]])
    if exc is not None:
        return
    farmer = gromacs_farmer(mdfarmer, work)
    clone = list(farmer.priority_ordered_clones[0])[0]
    exc = raised(lambda: clone.check_start_gen(set()))
    results.append(['a generation is planted', exc is None,
                    f'-> {exc!r}'[:90]])
    if exc is not None:
        return
    gen_dir = clone.current_gen_dir
    for name in gen_artefacts:
        results.append([f'the job has its {name}', (gen_dir / name).is_file(),
                        f'-> {gen_dir / name}'])
    results.append(['the run script calls the GROMACS sim block',
                    'gmx_basic_sim_block_json'
                    in (gen_dir / 'run.py').read_text(), ''])


def probe_openmm_entry_points(results, mdfarmer):
    """Every OpenMM-only symbol is present and says what is missing."""
    util = mdfarmer.utilities
    cases = (('mdfarmer.omm_generation', lambda: mdfarmer.omm_generation()),
             ('mdfarmer.omm_basic_sim_block_json',
              lambda: mdfarmer.omm_basic_sim_block_json({})),
             ('mdfarmer.simulate', lambda: mdfarmer.simulate),
             ('select_platform', lambda: util.select_platform()),
             ('read_openmm_top', lambda: util.read_openmm_top('seed.pdb')),
             ('openmm_topology_readers',
              lambda: util.openmm_topology_readers))
    for name, call in cases:
        exc = raised(call)
        results.append([f'{name} raises ImportError naming OpenMM',
                        names_openmm(exc), f'-> {exc!r}'[:90]])
        results.append([f'{name} says how to install it',
                        exc is not None and 'conda-forge openmm' in str(exc),
                        f'-> {exc}'[:90]])


def probe_state_xml(results, mdfarmer, work, step_count=SEED_STEP_COUNT):
    """A state.xml is still counted, and honestly refused, without OpenMM."""
    util = mdfarmer.utilities
    path = state_xml(work / 'seed.xml', step_count=step_count)
    read = util.state_xml_step_count(path)
    results.append(["the seed's step count is read without OpenMM",
                    read == step_count, f'-> {read}'])
    results.append(['a seed carrying that count is the campaign origin',
                    util.state_xml_origin(path) == step_count, ''])
    results.append(['no OpenMM makes a state.xml unusable, not broken',
                    util.is_state_xml_usable(path) is False, ''])


def child_main(work, blocked=BLOCKED_MODULE, marker=RESULT_MARKER):
    """Run every probe with OpenMM blocked, and print the results as JSON."""
    sys.meta_path.insert(0, OpenMMBlocker(blocked=blocked))
    results = []
    mdfarmer = probe_import(results)
    probe_gromacs(results, mdfarmer, work)
    probe_openmm_entry_points(results, mdfarmer)
    probe_state_xml(results, mdfarmer, work)
    print(f'{marker}{json.dumps(results)}')
    return 0


# The child branch runs before harness is imported: importing harness imports
# the package, and the child has to install its blocker first.
if __name__ == '__main__' and CHILD_FLAG in sys.argv:
    sys.exit(child_main(Path(sys.argv[sys.argv.index(CHILD_FLAG) + 1])))

import harness  # noqa: E402
from harness import Suite  # noqa: E402


def run_child(work, marker=RESULT_MARKER, timeout=CHILD_TIMEOUT):
    """(results, stderr) from this file re-run with OpenMM blocked."""
    done = sp.run([sys.executable, '-u', str(Path(__file__).resolve()),
                   CHILD_FLAG, str(work)],
                  capture_output=True, text=True, timeout=timeout)
    for line in done.stdout.splitlines():
        if line.startswith(marker):
            return json.loads(line[len(marker):]), done.stderr
    return None, done.stdout + done.stderr


def check_without_openmm(suite, work):
    suite.section('a site that never installed OpenMM')
    results, stderr = run_child(work / 'child')
    if results is None:
        suite.check('the subprocess reported its results', False,
                    f'-> {stderr[-400:]}')
        return
    for name, ok, detail in results:
        suite.check(name, ok, detail)


def check_with_openmm(suite):
    """The surface OpenMM sites already had, which laziness must not shrink."""
    suite.section('the same package where OpenMM is installed')
    import mdfarmer
    from mdfarmer import utilities as util
    if not util.openmm_available():
        raise harness.Skip('OpenMM is absent, so there is nothing to compare')
    suite.check('the runner module is imported eagerly',
                mdfarmer.simulate.omm_generation is mdfarmer.omm_generation)
    suite.check('the reader table is the real openmm.app one',
                util.openmm_topology_readers['.pdb'].__name__ == 'PDBFile',
                f"-> {util.openmm_topology_readers['.pdb']}")
    suite.check('the table is built once and reused',
                util.openmm_topology_readers is util.openmm_topology_readers)
    suite.check('a Farmer still defaults to the OpenMM runner',
                mdfarmer.farmer.omm_generation_runner()
                is mdfarmer.omm_generation)
    suite.check('and still recognises it as the OpenMM one',
                mdfarmer.farmer.is_omm_generation(mdfarmer.omm_generation))
    suite.check('a GROMACS runner is not mistaken for it',
                not mdfarmer.farmer.is_omm_generation(mdfarmer.gmx_generation))


def main():
    suite = Suite('openmm_optional')
    work = harness.workdir('openmm_optional')
    check_without_openmm(suite, work)
    check_with_openmm(suite)
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
