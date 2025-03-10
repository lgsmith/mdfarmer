import json
import openmm.app as app
import openmm as mm
from openmm import unit
from . import utilities as util
from pathlib import Path

from os import environ

# This function is written so that you could use jug's 'Task' class to uplift
# instances of calls. It returns the path to the trajectory written.


def omm_generation(traj_dir_top_level: str,
                   system_fn: str,
                   top_fn: str,
                   seed_index: int,
                   clone_index: int,
                   gen_index: int,
                   title: str,
                   integrator_xml: str,
                   # Begin from these coordinates/velocities. Expects a path to an OpenMM State.
                   seed_fn: str,
                   # if restarting, Do we start fresh or do we append?
                   append=False,
                   # dir-name indexes zero padded by this value. If 'None', then no padding.
                   dirname_pad=2,
                   # job and dir name separator
                   sep='-',
                   # mimics fah directory structure for processed trjs.
                   traj_name='traj',
                   traj_suffix='.xtc',
                   restart_name='state.xml',
                   # If None, mm picks whichever platform it thinks is fastest at runtime.
                   platform_name=None,
                   # passed to PlatformProperties
                   platform_properties=None,
                   # State data reporter's kwargs, true or false options, see omm docs.
                   state_data_kwargs=None,
                   # 2e7 is 100 ns, given 0.004 ps dt. Needs to be an int.
                   steps=25 * 10 ** 6,
                   # If provided, take this many steps without writing any output before starting to report.
                   eq_steps=None,
                   # Given 0.004 ps dt, 10 ps write freq.
                   write_interval=2500,
                   # If true, run minimizeEnergy on simulation before taking steps.
                   # simulation parameters below here; standard values for normal solvated protein inserted.
                   # Note units in comments
                   minimize_first=False,
                   # Integrator parameters
                   temperature=None,  # kelvin
                   new_velocities=False
                   ):

    # make reporter by extension
    reporters = {
        '.dcd': app.DCDReporter,
        '.xtc': app.XTCReporter,
    }

    if not state_data_kwargs:
        state_data_kwargs = dict(
            totalSteps=steps,
            step=True,
            speed=True,
            progress=True,
            potentialEnergy=True,
            temperature=True,
            separator='\t'
        )

    print('starting', title, seed_index, clone_index, gen_index)
    traj_dir = util.dir_seeds_clones_gens(Path(traj_dir_top_level), seed_index,
                                          clone_index,
                                          gen_index, dirname_pad, sep=sep)
    traj_path = (traj_dir / traj_name).with_suffix(traj_suffix)

    try:
        if traj_path.is_file():
            traj_reporter = reporters[traj_suffix](str(traj_path),
                                                   write_interval,
                                                   append=append)
        else:
            traj_reporter = reporters[traj_suffix](str(traj_path),
                                                   write_interval)
    except KeyError:
        print('You seem to have used a trajectory extension,',
              traj_suffix, 'for which no reporter is implemented yet.\n',
              'Your current choices are:', *reporters.keys())
        raise

    data_reporter_p = traj_path.with_suffix('.out')
    data_reporter = app.StateDataReporter(
        str(data_reporter_p),
        write_interval,
        append=append,
        **state_data_kwargs)

    # This will write xmls with system velocities in them
    restart_reporter = app.CheckpointReporter(
        restart_name,
        write_interval,
        writeState=True)

    system = mm.XmlSerializer.deserialize(Path(system_fn).read_text())
    topology = app.PDBFile(top_fn).topology

    integrator = mm.XmlSerializer.deserialize(Path(integrator_xml).read_text())
    if platform_name:
        platform = mm.Platform.getPlatformByName(platform_name)
    else:
        platform = None
    if platform_properties:
        simulation = app.Simulation(topology, system, integrator,
                                    platform=platform,
                                    platformProperties=platform_properties)
    else:
        simulation = app.Simulation(topology, system, integrator,
                                    platform=platform)
    simulation.loadState(seed_fn)

    if new_velocities:
        simulation.context.setVelocitiesToTemperature(
            temperature * unit.kelvin)
    # Because of openmm issue with xtcwriter, check if seed is current gen's restart
    if traj_path.suffix == '.xtc':
        seed_p = Path(seed_fn)
        if traj_path.parent != seed_p.parent or not append:
            simulation.currentStep = 0
            print('Not appending or restarting, so setting steps to 0.',
                  flush=True)
    if minimize_first:
        print('Performing energy minimization...')
        simulation.minimizeEnergy()
    if eq_steps:
        print('Equilibrating...')
        simulation.step(eq_steps)
        simulation.currentStep = simulation.currentStep - eq_steps
        print(f'Ran {eq_steps} of equilibration.')
    # run simulation here
    print('Simulating...')
    simulation.reporters.append(traj_reporter)
    simulation.reporters.append(data_reporter)
    simulation.reporters.append(restart_reporter)
    simulation.step(steps)

    # To cover the case where the restart reporter's interval isn't an
    # integer multiple of the number of steps taken, make sure the last frame
    # gets reported.
    state = simulation.context.getState(getPositions=True, getVelocities=True,
                                        enforcePeriodicBox=True)
    restart_reporter.report(simulation, state)

    return traj_path.resolve()


def omm_basic_sim_block_json(config):
    with open(config, 'r') as f:
        conf_dict = json.load(f)

    try:
        device_idx = environ['OMM_DeviceIndex']
    except KeyError:
        device_idx = None
        print('OMM_DeviceIndex not set; it defines what device(s) to use.',
              'Continuing with default behavior.')
    if device_idx:
        if conf_dict['platform_properties']:
            conf_dict['platform_properties']['DeviceIndex'] = device_idx
        else:
            conf_dict['platform_properties'] = {'DeviceIndex': device_idx}

    traj_list_path = Path(conf_dict['traj_list'])
    del conf_dict['traj_list']
    new_traj_path = omm_generation(**conf_dict)
    with traj_list_path.open('a') as tl:
        tl.write(str(new_traj_path) + '\n')
