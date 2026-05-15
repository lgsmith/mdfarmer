import json
import openmm.app as app
import openmm as mm
from openmm import unit
from . import utilities as util
from pathlib import Path

from os import environ


# Subclass that pushes Python's user-space buffer to the kernel after each
# frame write. A SIGKILL from Slurm preempt otherwise loses bytes still
# sitting in the BufferedWriter; flushing promotes them to the kernel page
# cache, which survives the process death (node stays up). flush() is
# microseconds and never blocks on disk — the kernel writes to disk
# asynchronously, on its own schedule, independent of the simulation loop.
#
# Don't be tempted to substitute buffering=0 on the underlying open():
# DCDFile.writeModel emits many small struct.pack writes per frame, each
# of which would become its own syscall under buffering=0, slowing the
# write loop dramatically. Buffered writes + explicit flush is the
# correct combination.
class FlushingDCDReporter(app.DCDReporter):
    def report(self, simulation, state):
        super().report(simulation, state)
        try:
            self._dcd._file.flush()
        except AttributeError:
            raise RuntimeError(
                'FlushingDCDReporter could not reach self._dcd._file; '
                'OpenMM internals likely changed. Update the wrapper.'
            )


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
        '.dcd': FlushingDCDReporter,
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
    # set up trajectory reporter
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
    
    topology = util.read_openmm_top(top_fn)


    # Set up reporters
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

    integrator = mm.XmlSerializer.deserialize(Path(integrator_xml).read_text())
    platform, platform_properties = util.select_platform(
        platform_name=platform_name, platform_properties=platform_properties)
    print(f'Using OpenMM platform: {platform.getName()}')
    if platform_properties is not None:
        simulation = app.Simulation(topology, system, integrator,
                                    platform=platform,
                                    platformProperties=platform_properties)
    else:
        simulation = app.Simulation(topology, system, integrator,
                                    platform=platform)
    simulation.loadState(seed_fn)

    # CUDA context warmup. When several GPU jobs initialize concurrently
    # on a shared node, the first getState(getPositions=True) readback
    # can race with kernel completion and return uninitialized device
    # memory — observed as a single garbage frame at frame 0 in
    # clone-{028,030,031,032}/gen-000 and clone-{046..049}/gen-001 of
    # the antifreeze dataset (four H200 jobs per node, same SLURM job
    # ID block). Pulling positions back here forces the device→host
    # sync before any reporter fires, so DCDReporter's first frame is
    # from a settled context.
    _ = simulation.context.getState(getPositions=True)

    if new_velocities:
        simulation.context.setVelocitiesToTemperature(temperature * unit.kelvin)

    if minimize_first:
        print('Performing energy minimization...')
        simulation.minimizeEnergy()
    # Equilibration runs only when starting a generation from scratch. On a
    # resume (append=True) the eq has already been done; re-running it would
    # also wind simulation.currentStep backwards, desyncing reporter triggers
    # from the loaded state.
    if eq_steps and not append:
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
    # The CheckpointReporter writes at every write_interval, so the most
    # recent state.xml on disk already aligns with the trajectory's last
    # frame. Writing one final state at a non-write_interval boundary
    # desyncs state.xml from the traj frame count, breaking the next
    # resume's calx_remaining_steps math.
    print('Done!')
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
