"""What a tandem velocity file really holds, per format and per reader.

The tandem reporters put velocities or forces in a trajectory's position slot,
so the only thing telling a scientist what scale the numbers are on is the
docstring. This suite pins that: it writes velocities it chose, reads them back
with mdtraj and with LOOS, and asserts the factor each one returns.
"""
import sys

import numpy as np

import harness
from harness import Skip, Suite

from mdfarmer import simulate as sim

N_ATOMS = 6
BOX_NM = 3.0

# Three decimals, so XTC's fixed 1000x precision stores them exactly.
VELOCITY_STEP = np.array([0.125, -0.250, 0.375])

# LOOS reports angstroms whatever the file's own length unit is.
LOOS_FACTOR = 10.0

# Float32 for DCD, 1e-3 nm quantization for XTC.
TOLERANCE = 2e-3


def _argon_box(n_atoms=N_ATOMS, box_nm=BOX_NM):
    """A periodic OpenMM (topology, system) of n_atoms uncharged argons."""
    import openmm as mm
    import openmm.app as app
    from openmm import unit

    vectors = (mm.Vec3(box_nm, 0, 0), mm.Vec3(0, box_nm, 0),
               mm.Vec3(0, 0, box_nm)) * unit.nanometer
    topology = app.Topology()
    residue = topology.addResidue('AR', topology.addChain())
    system = mm.System()
    force = mm.NonbondedForce()
    force.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
    force.setCutoffDistance(1.0 * unit.nanometer)
    for _ in range(n_atoms):
        topology.addAtom('AR', app.element.argon, residue)
        system.addParticle(40.0 * unit.amu)
        force.addParticle(0.0, 0.3, 0.5)
    topology.setPeriodicBoxVectors(vectors)
    system.setDefaultPeriodicBoxVectors(*vectors)
    system.addForce(force)
    return topology, system


def _seeded_simulation(velocities, work):
    """A Reference-platform Simulation holding exactly these velocities.

    Returns (simulation, pdb_path); the PDB is the model both readers need to
    open the tandem file.
    """
    import openmm as mm
    import openmm.app as app
    from openmm import unit

    topology, system = _argon_box(n_atoms=len(velocities))
    simulation = app.Simulation(
        topology, system, mm.VerletIntegrator(0.001 * unit.picosecond),
        platform=mm.Platform.getPlatformByName('Reference'))
    positions = np.array([[0.5 * (i + 1), 1.0, 1.5]
                          for i in range(len(velocities))])
    simulation.context.setPositions(positions * unit.nanometer)
    simulation.context.setVelocities(
        velocities * (unit.nanometer / unit.picosecond))
    pdb = work / 'model.pdb'
    with open(pdb, 'w') as handle:
        app.PDBFile.writeFile(topology, positions * unit.nanometer, handle)
    return simulation, pdb


def _write_tandem(simulation, path, n_frames=2):
    """Report the simulation's current velocities into path, n_frames times."""
    reporter = sim._TANDEM_REPORTER_CLS[path.suffix](str(path), 1, 'velocities')
    state = simulation.context.getState(getPositions=True, getVelocities=True)
    for _ in range(n_frames):
        reporter.report(simulation, state)
    del reporter
    return path


def _loos_first_frame(traj_path, pdb):
    """Atom coordinates of traj_path's first frame, as LOOS reports them."""
    import loos
    import loos.pyloos

    model = loos.createSystem(str(pdb))
    for frame in loos.pyloos.Trajectory(str(traj_path), model):
        return np.array([[a.coords().x(), a.coords().y(), a.coords().z()]
                         for a in frame])
    raise AssertionError(f'{traj_path} held no frames')


def _h5_node_paths(path):
    """Every node path in an mdtraj HDF5 trajectory file."""
    import tables

    with tables.open_file(str(path), 'r') as handle:
        return [node._v_pathname for node in handle.walk_nodes('/')]


def main(tolerance=TOLERANCE, loos_factor=LOOS_FACTOR):
    try:
        import mdtraj as md
        import tables  # noqa: F401  -- mdtraj's HDF5 backend
    except ImportError as exc:
        raise Skip(f'mdtraj with HDF5 support is required: {exc}')

    suite = Suite('tandem units')
    work = harness.workdir('tandem-units')
    written = np.array([VELOCITY_STEP * (i + 1) for i in range(N_ATOMS)])
    simulation, pdb = _seeded_simulation(written, work)

    suite.section('mdtraj hands back native nm/ps for every tandem format')
    for suffix in ('.dcd', '.xtc', '.h5'):
        path = _write_tandem(simulation, work / f'velocities{suffix}')
        traj = (md.load(str(path)) if suffix == '.h5'
                else md.load(str(path), top=str(pdb)))
        error = np.abs(traj.xyz[0] - written).max()
        suite.check(f'mdtraj reads {suffix} back as written',
                    error < tolerance, f'-> max error {error:.2e} nm/ps')

    suite.section('the DCD bytes are ten times native, the XTC bytes are not')
    with md.formats.DCDTrajectoryFile(str(work / 'velocities.dcd')) as handle:
        dcd_raw = handle.read()[0][0]
    suite.check('DCD stores 10*x, as DCDFile.writeModel does for positions',
                np.abs(dcd_raw - loos_factor * written).max() < tolerance,
                f'-> atom 0 {dcd_raw[0]} for {written[0]}')
    with md.formats.XTCTrajectoryFile(str(work / 'velocities.xtc')) as handle:
        xtc_raw = handle.read()[0][0]
    suite.check('XTC stores the raw number, so the two formats differ on disk',
                np.abs(xtc_raw - written).max() < tolerance,
                f'-> atom 0 {xtc_raw[0]} for {written[0]}')

    suite.section('LOOS hands back ten times native for both formats')
    for suffix in ('.dcd', '.xtc'):
        coords = _loos_first_frame(work / f'velocities{suffix}', pdb)
        error = np.abs(coords - loos_factor * written).max()
        suite.check(f'LOOS reads {suffix} as angstroms, {loos_factor:g}x native',
                    error < tolerance * loos_factor,
                    f'-> atom 0 {coords[0]} for {written[0]}')

    suite.section('the tandem HDF5 file has no velocities node')
    nodes = _h5_node_paths(work / 'velocities.h5')
    suite.check('velocities land in /coordinates', '/coordinates' in nodes)
    suite.check('there is no /velocities node to read them from',
                '/velocities' not in nodes, f'-> {nodes}')

    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
