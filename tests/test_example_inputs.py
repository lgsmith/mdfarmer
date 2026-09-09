"""The committed example inputs inflate, load, and add up to the run shape.

examples/ ships its trp-cage inputs gzipped so the two shakedown campaigns are
self-contained without carrying 7 MB of XML. Nothing else in the tree reads
them, so a corrupt blob or a driver whose step arithmetic drifted out of
commensurability would only show up when someone tried to launch a campaign.

The same goes for each example's farm-*.sh, which cannot be run here without
starting a real tender: what is checked is that it parses and that every path it
hardcodes -- campaign directory, brake file, log, lock, GROMACS modules -- still
agrees with the driver it launches.
"""
import gzip
import importlib.util
import os
import re
import subprocess as sp
import sys
from pathlib import Path

import openmm as mm
from openmm import unit

import harness
from harness import Suite

import mdfarmer

EXAMPLES = harness.REPO_ROOT / 'examples'
N_ATOMS = 13465
N_SOLUTE_ATOMS = 304
N_CONSTRAINTS = 10011
FRAMES_PER_GEN = 5


def load_driver(name, path):
    """Import an example's farmer.py under its own module name.

    By path, and with the example directory on sys.path, because both drivers
    are called farmer.py and neither is importable as a package.
    """
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def inflated_text(gz_path):
    with gzip.open(gz_path, 'rt') as fh:
        return fh.read()


def main():
    suite = Suite('example_inputs')
    work = harness.workdir('example_inputs')
    omm = load_driver('example_omm_farmer',
                      EXAMPLES / 'openmm-trpcage' / 'farmer.py')
    gmx = load_driver('example_gmx_farmer',
                      EXAMPLES / 'gromacs-trpcage' / 'farmer.py')

    suite.section('every committed input is there and is gzip')
    committed = [EXAMPLES / 'openmm-trpcage/inputs/system.xml.gz',
                 EXAMPLES / 'openmm-trpcage/inputs/state.xml.gz',
                 EXAMPLES / 'openmm-trpcage/inputs/topology.pdb.gz',
                 EXAMPLES / 'gromacs-trpcage/inputs/gmx.gro.gz',
                 EXAMPLES / 'gromacs-trpcage/inputs/gmx.top.gz']
    for path in committed:
        suite.check(f'{path.name} is present and gzip',
                    path.is_file()
                    and path.read_bytes()[:2] == b'\x1f\x8b',
                    f'-> {path.stat().st_size if path.is_file() else "missing"} bytes')
    suite.check('the .mdp is committed plain, being small and worth reading',
                (EXAMPLES / 'gromacs-trpcage/inputs/prod-277.mdp').is_file())

    suite.section('OpenMM inflates its inputs in process')
    system = mm.XmlSerializer.deserialize(
        inflated_text(EXAMPLES / 'openmm-trpcage/inputs/system.xml.gz'))
    suite.check('system.xml.gz deserializes to the whole system',
                system.getNumParticles() == N_ATOMS
                and system.getNumConstraints() == N_CONSTRAINTS,
                f'-> {system.getNumParticles()} particles, '
                f'{system.getNumConstraints()} constraints')
    suite.check('and it carries the barostat the campaign runs under',
                any(isinstance(f, mm.MonteCarloBarostat)
                    for f in system.getForces()))
    state = mm.XmlSerializer.deserialize(
        inflated_text(EXAMPLES / 'openmm-trpcage/inputs/state.xml.gz'))
    suite.check('state.xml.gz deserializes to positions for every particle',
                len(state.getPositions()) == N_ATOMS,
                f'-> {len(state.getPositions())}')

    suite.section('the driver prepares what the runner is given')
    paths = omm.prepare_inputs(prepared=work / 'omm-prepared')
    topology = mdfarmer.read_openmm_top(paths['top_fn'])
    suite.check('the inflated topology loads through the reader Farmer uses',
                topology.getNumAtoms() == N_ATOMS,
                f'-> {topology.getNumAtoms()}')
    suite.check('the inflated system loads the way omm_generation loads it',
                mm.XmlSerializer.deserialize(
                    Path(paths['system_fn']).read_text()).getNumParticles()
                == N_ATOMS)
    integrator = mm.XmlSerializer.deserialize(
        Path(paths['integrator_xml']).read_text())
    suite.check("the integrator it writes takes the mdp's 2 fs step",
                abs(integrator.getStepSize().value_in_unit(unit.picoseconds)
                    - omm.DT_PS) < 1e-12,
                f'-> {integrator.getStepSize()}')

    suite.section('the seed state starts the step axis at zero')
    # seeder._try_recover_gen reads state.xml's stepCount as a step counted
    # from the start of the campaign. The equilibration this system came from
    # left 150000 in there, which would make every resumed generation look
    # finished before it started.
    raw_step = mdfarmer.state_xml_step_count(
        work / 'omm-prepared' / 'state.xml')
    suite.check('prepare_inputs rewound stepCount', raw_step == 0,
                f'-> {raw_step}')
    suite.check('and the positions survived the rewind',
                len(mm.XmlSerializer.deserialize(
                    Path(paths['seed_fn']).read_text()).getPositions())
                == N_ATOMS)

    suite.section('GROMACS gets real files on disk, since grompp needs them')
    gmx_paths = gmx.prepare_inputs(prepared=work / 'gmx-prepared')
    gro_lines = Path(gmx_paths['structure_fn']).read_text().splitlines()
    suite.check('gmx.gro.gz inflates to a .gro of the same system',
                int(gro_lines[1]) == N_ATOMS and len(gro_lines) == N_ATOMS + 3,
                f'-> {gro_lines[1].strip()} atoms, {len(gro_lines)} lines')
    top_text = Path(gmx_paths['top_fn']).read_text()
    suite.check('gmx.top.gz inflates to a topology with a molecules block',
                '[ molecules ]' in top_text and '[ atomtypes ]' in top_text)
    suite.check('the .mdp the driver names is the committed one',
                Path(gmx_paths['mdp_fn']).is_file()
                and 'dt                       = 0.002' in
                Path(gmx_paths['mdp_fn']).read_text())

    suite.section('both engines are asked for the same run shape')
    for name, driver in (('openmm', omm), ('gromacs', gmx)):
        frames = mdfarmer.check_commensurability(
            driver.STEPS_PER_GEN, driver.WRITE_INTERVAL, driver.DOWNSAMPLE_FRQ)
        suite.check(f'{name}: {driver.STEPS_PER_GEN} steps at '
                    f'{driver.WRITE_INTERVAL} is {FRAMES_PER_GEN} frames, '
                    f'divisible by downsample {driver.DOWNSAMPLE_FRQ}',
                    frames == FRAMES_PER_GEN, f'-> {frames}')
        suite.check(f'{name}: check_whole_frames accepts the same pair',
                    mdfarmer.check_whole_frames(driver.STEPS_PER_GEN,
                                                driver.WRITE_INTERVAL)
                    == driver.STEPS_PER_GEN)
    suite.check('the two arms run the same number of clones and generations',
                (omm.N_CLONES, omm.N_GENS) == (gmx.N_CLONES, gmx.N_GENS)
                == (10, 5),
                f'-> omm {(omm.N_CLONES, omm.N_GENS)}, '
                f'gmx {(gmx.N_CLONES, gmx.N_GENS)}')
    suite.check('the GROMACS arm packs those clones two to a card',
                gmx.REPS_PER_PACK == 2
                and gmx.N_CLONES % gmx.REPS_PER_PACK == 0
                and gmx.PACK_CPUS == gmx.CORES_PER_REPLICA * gmx.REPS_PER_PACK,
                f'-> {gmx.N_CLONES // gmx.REPS_PER_PACK} packs of '
                f'{gmx.REPS_PER_PACK} at {gmx.PACK_CPUS} cores')
    # 4-site water: GROMACS refuses GPU update with virtual sites, measured as
    # "Update task on the GPU was required, but ... Virtual sites are not
    # supported." Flipping this on would make every generation fail at mdrun.
    suite.check('and leaves the update on the CPU, which the water forces',
                gmx.UPDATE_MODE == 'cpu'
                and list(gmx.BASE_MDRUN_ARGS)[
                    list(gmx.BASE_MDRUN_ARGS).index('-update') + 1] == 'cpu')

    suite.section('the harvest subset picks out the solute, not nothing')
    spec = mdfarmer.resolve_subset(gmx_paths['structure_fn'],
                                   gmx.HARVESTER_SUBSET)
    suite.check('the LOOS selection keeps only the peptide',
                len(spec['indices']) == N_SOLUTE_ATOMS,
                f'-> {len(spec["indices"])} atoms')
    suite.check('and both arms select with the same expression',
                omm.HARVESTER_SUBSET == gmx.HARVESTER_SUBSET)

    suite.section('inflating a second time does nothing, since a tender re-boots')
    # Deliberately last: the sentinel check overwrites an inflated file to prove
    # that a second prepare_inputs leaves what is already there alone.
    for name, driver, prepared, inflated in (
            ('openmm', omm, work / 'omm-prepared',
             ('system.xml', 'topology.pdb', 'state.xml')),
            ('gromacs', gmx, work / 'gmx-prepared', ('gmx.gro', 'gmx.top'))):
        before = {n: (prepared / n).stat().st_mtime_ns for n in inflated}
        driver.prepare_inputs(prepared=prepared)
        after = {n: (prepared / n).stat().st_mtime_ns for n in inflated}
        suite.check(f'{name}: a second prepare_inputs rewrites nothing',
                    before == after,
                    f'-> {sorted(n for n in inflated if before[n] != after[n])}')
        suite.check(f'{name}: and leaves no half-inflated file behind',
                    not list(prepared.glob('*.partial')))
        sentinel = prepared / inflated[0]
        sentinel.write_text('SENTINEL')
        driver.prepare_inputs(prepared=prepared)
        suite.check(f'{name}: what is on disk is what the runner gets',
                    sentinel.read_text() == 'SENTINEL')
    # state.xml is rewound while staged, so it is never on disk under its real
    # name carrying the equilibration's step count.
    suite.check('the seed is only ever written already rewound',
                mdfarmer.state_xml_step_count(
                    work / 'omm-prepared' / 'state.xml') == 0)

    suite.section('each example ships a drive script that agrees with its driver')
    gitignore = (EXAMPLES / '.gitignore').read_text()
    for driver, arm, script_name in ((omm, 'openmm-trpcage', 'farm-omm.sh'),
                                     (gmx, 'gromacs-trpcage', 'farm-gmx.sh')):
        script = EXAMPLES / arm / script_name
        suite.check(f'{script_name} is committed and executable',
                    script.is_file() and os.access(script, os.X_OK))
        parsed = sp.run(['bash', '-n', str(script)], capture_output=True,
                        text=True)
        suite.check(f'{script_name} parses', parsed.returncode == 0,
                    parsed.stderr.strip())
        text = script.read_text()
        suite.check(f'{script_name} stops at the first failed command',
                    'set -euo pipefail' in text)
        suite.check(f'{script_name} names the campaign the driver writes to',
                    f"PROJECT='{driver.PROJECT}'" in text
                    and 'CAMPAIGN="$HERE/data/$PROJECT"' in text
                    and Path(driver.TRAJ_TOP)
                    == Path(driver.HERE) / 'data' / driver.PROJECT)
        # Farmer.start_tending_fields brakes on Path('stop'), read from its cwd.
        suite.check(f'{script_name} writes the brake file where Farmer reads it',
                    'STOP="$HERE/stop"' in text and 'cd "$HERE"' in text)
        suite.check(f'{script_name} allows one tender per campaign, under flock',
                    'LOCK="$CAMPAIGN/tender.lock"' in text
                    and 'flock -n "$LOCK" true' in text
                    and 'flock -n 9' in text)
        suite.check(f'{script_name} logs where .gitignore already covers it',
                    'LOG="$HERE/$PROJECT.tend.out"' in text
                    and '*.tend.out' in gitignore)
        suite.check(f'{script_name} detaches, so the tender outlives the shell',
                    'setsid nohup bash "$SELF" --loop' in text)
        # Three outcomes, three statuses: 0 finished, BRAKED_EXIT asked to
        # stop, 1 died. The loop re-enters on 1 alone.
        farmer_text = (EXAMPLES / arm / 'farmer.py').read_text()
        suite.check(f'{arm}: the driver exits 0 only when every clone finished',
                    'raise SystemExit(0)' in farmer_text)
        suite.check(f'{arm}: and separates being braked from having died',
                    'SystemExit(BRAKED_EXIT if' in farmer_text)
        py_braked = re.search(r'^BRAKED_EXIT = (\d+)', farmer_text, re.M)
        sh_braked = re.search(r'^BRAKED_EXIT=(\d+)', text, re.M)
        # A disagreement here would re-enter a braked campaign forever.
        suite.check(f'{script_name} stops on the same status the driver exits',
                    py_braked and sh_braked
                    and py_braked.group(1) == sh_braked.group(1),
                    f'-> py={py_braked and py_braked.group(1)} '
                    f'sh={sh_braked and sh_braked.group(1)}')
        suite.check(f'{script_name} re-enters on a death and not on a brake',
                    '"$status" -eq "$BRAKED_EXIT"' in text
                    and 'not re-entering' in text)
    gmx_modules = re.search(
        r"GMX_MODULES='([^']*)'",
        (EXAMPLES / 'gromacs-trpcage' / 'farm-gmx.sh').read_text()).group(1)
    suite.check('farm-gmx.sh preflights the modules the job script loads',
                f'module load {gmx_modules}' in gmx.ENV_SETUP,
                f'-> {gmx_modules}')
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
