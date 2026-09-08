"""The committed example inputs inflate, load, and add up to the run shape.

examples/ ships its trp-cage inputs gzipped so the two shakedown campaigns are
self-contained without carrying 7 MB of XML. Nothing else in the tree reads
them, so a corrupt blob or a driver whose step arithmetic drifted out of
commensurability would only show up when someone tried to launch a campaign.
"""
import gzip
import importlib.util
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
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
