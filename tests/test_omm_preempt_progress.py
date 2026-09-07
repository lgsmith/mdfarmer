"""A preempted OpenMM generation still has to say how far it got.

The GROMACS runner needs an explicit write_gen_status on its preempt path,
because there the orchestrator reads gen_status.json rather than mdrun's own
checkpoint, and a preempt that left only a checkpoint read as a fresh start.
The OpenMM path keeps no second record: the trajectory and state.xml that
calx_remaining_steps and _try_recover_gen read ARE what the reporters wrote,
and SentinelReporter is appended last, so every writer for that cycle has
already fired before Preempted is raised. This suite pins that down.
"""
import json
import os
import sys
from pathlib import Path

import harness
from harness import Suite

import openmm as mm
import openmm.app as app
from openmm import unit

import mdfarmer
from mdfarmer import simulate as sim
from mdfarmer import seeder
from mdfarmer import utilities as util

STEPS_PER_GEN = 400
WRITE_INTERVAL = 100
PREEMPT_AFTER_CYCLES = 2        # cycles the launch gets before the sentinel
SEED_INDEX = 0
CLONE_INDEX = 0
DIRNAME_PAD = 2
SEP = '-'

# A lattice of argon, big enough to integrate honestly and small enough that
# 400 CPU steps cost nothing. Spacing sits past the LJ minimum, so it settles.
LATTICE_SIDE = 4
LATTICE_SPACING_NM = 0.5
BOX_NM = 2.4
CUTOFF_NM = 1.0
ARGON_SIGMA_NM = 0.34
ARGON_EPSILON_KJ = 0.996
ARGON_MASS_AMU = 39.948
TEMPERATURE_K = 300
FRICTION_PER_PS = 1.0
TIMESTEP_PS = 0.002
PLATFORM_NAME = 'CPU'           # never take a GPU from a running campaign


def build_argon_box(dest, lattice_side=LATTICE_SIDE,
                    spacing=LATTICE_SPACING_NM, box_nm=BOX_NM):
    """(system_fn, integrator_fn, pdb_fn, seed_fn) for a small periodic argon box."""
    dest = Path(dest)
    box = [mm.Vec3(box_nm, 0, 0), mm.Vec3(0, box_nm, 0), mm.Vec3(0, 0, box_nm)]
    argon = app.Element.getBySymbol('Ar')

    topology = app.Topology()
    chain = topology.addChain()
    system = mm.System()
    system.setDefaultPeriodicBoxVectors(*box)
    nonbonded = mm.NonbondedForce()
    nonbonded.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
    nonbonded.setCutoffDistance(CUTOFF_NM * unit.nanometer)
    positions = []
    for i in range(lattice_side):
        for j in range(lattice_side):
            for k in range(lattice_side):
                residue = topology.addResidue('AR', chain)
                topology.addAtom('AR', argon, residue)
                system.addParticle(ARGON_MASS_AMU * unit.amu)
                nonbonded.addParticle(
                    0.0, ARGON_SIGMA_NM * unit.nanometer,
                    ARGON_EPSILON_KJ * unit.kilojoule_per_mole)
                positions.append(mm.Vec3(i, j, k) * spacing)
    system.addForce(nonbonded)
    topology.setPeriodicBoxVectors(box)
    positions = positions * unit.nanometer

    pdb_fn = dest / 'argon.pdb'
    with pdb_fn.open('w') as fh:
        app.PDBFile.writeFile(topology, positions, fh)

    integrator = mm.LangevinMiddleIntegrator(
        TEMPERATURE_K * unit.kelvin, FRICTION_PER_PS / unit.picosecond,
        TIMESTEP_PS * unit.picosecond)
    system_fn = dest / 'system.xml'
    integrator_fn = dest / 'integrator.xml'
    system_fn.write_text(mm.XmlSerializer.serialize(system))
    integrator_fn.write_text(mm.XmlSerializer.serialize(integrator))

    seed_fn = dest / 'seed.xml'
    simulation = app.Simulation(
        topology, system, integrator,
        platform=mm.Platform.getPlatformByName(PLATFORM_NAME))
    simulation.context.setPositions(positions)
    simulation.context.setVelocitiesToTemperature(TEMPERATURE_K * unit.kelvin)
    simulation.saveState(str(seed_fn))
    return str(system_fn), str(integrator_fn), str(pdb_fn), str(seed_fn)


class SentinelAfter(sim.SentinelReporter):
    """Drop the preempt sentinel on the Nth cycle, then act as the real one.

    Stands in for the batch script's SIGTERM trap without a race: the touch and
    the check happen in one report, so the cycle a preempt lands on is fixed.
    """

    def __init__(self, reportInterval, sentinel_path=None,
                 after=PREEMPT_AFTER_CYCLES):
        super().__init__(reportInterval, sentinel_path=sentinel_path)
        self._after = after
        self._cycles = 0

    def report(self, simulation, state):
        self._cycles += 1
        if self._cycles >= self._after:
            self._sentinel.touch()
        super().report(simulation, state)


def write_gen_config(gen_dir, config):
    """The config.json _try_recover_gen reads to place a generation."""
    (Path(gen_dir) / util.CONFIG_NAME).write_text(json.dumps(config, indent=2))


def main():
    suite = Suite('omm_preempt_progress')
    work = harness.workdir('omm_preempt_progress')
    system_fn, integrator_fn, pdb_fn, seed_fn = build_argon_box(work)
    farm = work / 'farm'
    gen_dir = mdfarmer.dir_seeds_clones_gens(
        farm, SEED_INDEX, CLONE_INDEX, 0, DIRNAME_PAD, sep=SEP)

    config = dict(
        traj_dir_top_level=str(farm), system_fn=system_fn, top_fn=pdb_fn,
        seed_index=SEED_INDEX, clone_index=CLONE_INDEX, gen_index=0,
        title='omm-preempt', integrator_xml=integrator_fn, seed_fn=seed_fn,
        append=False, dirname_pad=DIRNAME_PAD, sep=SEP, traj_name='positions',
        traj_suffix='.dcd', restart_name='state.xml',
        platform_name=PLATFORM_NAME, steps=STEPS_PER_GEN,
        write_interval=WRITE_INTERVAL, handle_preempt=True)
    write_gen_config(gen_dir, dict(config, steps_per_gen=STEPS_PER_GEN))

    suite.section('a launch is preempted partway through')
    real_sentinel_cls = sim.SentinelReporter
    sim.SentinelReporter = SentinelAfter
    cwd = Path.cwd()
    try:
        os.chdir(gen_dir)       # state.xml and the sentinel are cwd-relative
        sim.omm_generation(**config)
        preempted = False
    except sim.Preempted:
        preempted = True
    finally:
        sim.SentinelReporter = real_sentinel_cls
        os.chdir(cwd)
    suite.check('the preempt reaches the orchestrator', preempted)

    reached = PREEMPT_AFTER_CYCLES * WRITE_INTERVAL
    traj_p = gen_dir / 'positions.dcd'
    frames = util.get_traj_len(str(traj_p), pdb_fn)
    suite.check('the trajectory holds every frame written before it raised',
                frames == PREEMPT_AFTER_CYCLES,
                f'-> {frames} vs {PREEMPT_AFTER_CYCLES}')

    restart_p = gen_dir / 'state.xml'
    suite.check('and a loadable state.xml sits beside it',
                util.is_state_xml_usable(restart_p))
    step_count = util.state_xml_step_count(restart_p)
    suite.check('recording the same step the last frame was written at',
                step_count == reached, f'-> {step_count} vs {reached}')

    suite.section('so the orchestrator reads progress, not a fresh start')
    remaining = util.calx_remaining_steps(
        str(traj_p), pdb_fn, STEPS_PER_GEN, WRITE_INTERVAL)
    suite.check('the steps still owed are what the launch did not reach',
                remaining == STEPS_PER_GEN - reached,
                f'-> {remaining} vs {STEPS_PER_GEN - reached}')
    suite.check('which is under total_steps, so no restart is charged',
                remaining < STEPS_PER_GEN, f'-> {remaining}')

    recovered = seeder._try_recover_gen(
        gen_dir, append_mode=True, restart_name='state.xml',
        traj_name='positions', traj_suffix='.dcd',
        write_interval=WRITE_INTERVAL, total_steps=STEPS_PER_GEN,
        top_fn=pdb_fn)
    suite.check('and recovery resumes this generation rather than redoing it',
                recovered is not None and recovered[0] == 0
                and recovered[2] == STEPS_PER_GEN - reached
                and recovered[3] is True,
                f'-> {recovered}')
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
