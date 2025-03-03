import json
import openmm.app as app
import openmm as mm

# This function is written so that you could use jug's 'Task' class to uplift
# instances of calls. It returns the path to the trajectory written.
def omm_generation(traj_dir_top_level: str,
                   prmtop_fn: str,
                   run_index: int,
                   clone_index: int,
                   gen_index: int,
                   title: str,
                   # Normally the project name. Note this is obligatory.
                   # options to write output, simulation options.
                   # If provided, begin from these coordinates/velocities. Checks for pdb or rst7.
                   # If None, use naming conventions below to look for restart from previous gen.
                   initial=None,
                   # if restarting, Do we start fresh or do we append?
                   append=False,
                   # dir-name indexes zero padded by this value. If 'None', then no padding.
                   dirname_pad=2,
                   # job and dir name separator
                   sep='-',
                   # mimics fah directory structure for processed trjs.
                   traj_name='traj',
                   # options are only .dcd at present
                   traj_suffix='.dcd',
                   restart_name='restart.rst7',
                   # If None, mm picks whichever platform it thinks is fastest at runtime.
                   platform_name=None,
                   # passed to PlatformProperties
                   platform_properties=None,
                   # State data reporter's kwargs, true or false options, see omm docs.
                   state_data_kwargs=None,
                   # 2e7 is 100 ns, given 0.005 ps dt. Needs to be an int.
                   steps=2 * 10 ** 7,
                   # If provided, take this many steps without writing any output before starting to report.
                   eq_steps=None,
                   # Given 0.005 ps dt, 10 ps write freq.
                   write_interval=2000,
                   # If true, run minimizeEnergy on simulation before taking steps.
                   # simulation parameters below here; standard values for normal solvated protein inserted.
                   # Note units in comments
                   minimize_first=False,
                   # Integrator parameters
                   temperature=300,  # kelvin
                   new_velocities=False,
                   # If true, get new velocities at temperature when initializing simulation.
                   pressure=1.0,  # atmospheres
                   # for LangevinMiddleIntegrator; 1/picosecond
                   friction=1.0,
                   barostat_interval=100,  # for MC barostat
                   # Default for Langevin; for mixed precision should be ok
                   constraint_tolerance=1e-8,
                   integrator_name='LangevinMiddleIntegrator',
                   dt=0.004,
                   # Picoseconds; chosen to work with selected hydrogen mass.
                   # System configuration
                   hydrogen_mass=2,
                   # AMU; chosen to work with selected dt. Passed to createSystem.
                   constraints='HBonds',
                   # One of constraints.keys() in this file, or None for none.
                   nonbond_method='PME',
                   # One of nonbonded_methods.keys() in this file.
                   nonbond_cutoff=1.0,
                   # Nanometers, 1 is a decent number for a normal water + protein system.
                   rigid_water=True,  # Passed to createSystem
                   ):
    if not platform_properties:
        platform_properties = {'Precision': 'mixed'}

    nonbonded_methods = {
        'NoCutoff': app.NoCutoff,
        'CutoffNonPeriodic': app.CutoffNonPeriodic,
        'CutoffPeriodic': app.CutoffPeriodic,
        'Ewald': app.Ewald,
        # note, you probably want PME unless you're doing something special
        'PME': app.PME,
        # Include attractive LJ forces in PME calculation. Useful for membranes.
        'LJPME': app.LJPME
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

    constraint_types = {
        'HBonds': app.HBonds,
        'AllBonds': app.AllBonds,
        'HAngles': app.HAngles
    }
    print('starting', title, run_index, clone_index, gen_index)
    traj_dir = dir_runs_clones_gens(Path(traj_dir_top_level), run_index,
                                    clone_index,
                                    gen_index, dirname_pad, sep=sep)
    traj_path = (traj_dir / traj_name).with_suffix(traj_suffix)

    # make reporter by extension
    # TODO: add more reporters and integrators
    reporters = {
        '.dcd': app.DCDReporter,
        '.rst7': parmed.openmm.reporters.RestartReporter,
        '.out': app.StateDataReporter
    }

    integrators = {
        "LangevinMiddleIntegrator": mm.LangevinMiddleIntegrator
    }
    try:
        if traj_path.is_file():
            traj_reporter = reporters[traj_suffix](str(traj_path),
                                                   write_interval,
                                                   append=append)
        else:
            traj_reporter = reporters[traj_suffix](str(traj_path),
                                                   write_interval)
    except KeyError:
        print(
            f'You seem to have used a trajectory extension, {traj_suffix}, '
            f'for which no reporter is implemented yet.')
        print('Your current choices are:', *reporters.keys())
        raise

    data_reporter_p = traj_path.with_suffix('.out')
    if data_reporter_p.is_file():
        data_reporter = reporters['.out'](str(data_reporter_p),
                                          write_interval,
                                          append=append,
                                          **state_data_kwargs)
    else:
        data_reporter = reporters['.out'](str(data_reporter_p),
                                          write_interval,
                                          **state_data_kwargs)

    restart_reporter = reporters['.rst7'](restart_name,
                                          write_interval,
                                          write_velocities=True)

    if initial:
        prmtop_crds = parmed.load_file(prmtop_fn, initial)
    elif append and Path(restart_name).is_file():
        prmtop_crds = parmed.load_file(prmtop_fn, restart_name)
    else:
        prev_gen_dir = dir_runs_clones_gens(Path(traj_dir_top_level), run_index,
                                            clone_index, gen_index - 1,
                                            dirname_pad, sep=sep)
        prev_rst_p = prev_gen_dir / restart_name
        if not prev_rst_p.is_file():
            print(prev_rst_p,
                  'Does not exist, and no restart was provided,',
                  'so no crds to run with. Exiting...')
            raise
        prmtop_crds = parmed.load_file(prmtop_fn, str(prev_rst_p))
    system = prmtop_crds.createSystem(
        nonbondedMethod=nonbonded_methods[nonbond_method],
        nonbondedCutoff=nonbond_cutoff,
        constraints=constraint_types[constraints],
        rigidWater=rigid_water,
        hydrogenMass=hydrogen_mass * u.amu
    )
    topology = prmtop_crds.topology
    positions = prmtop_crds.positions
    if barostat_interval > 0:
        system.addForce(
            mm.MonteCarloBarostat(pressure * u.atmospheres,
                                  temperature * u.kelvin, barostat_interval))

    integrator = integrators[integrator_name](temperature * u.kelvin,
                                              friction / u.picoseconds,
                                              dt * u.picoseconds)
    integrator.setConstraintTolerance(constraint_tolerance)
    if platform_name:
        platform = mm.Platform_getPlatformByName(platform_name)
    else:
        platform = None
    simulation = app.Simulation(topology, system, integrator,
                                platform=platform,
                                platformProperties=platform_properties)
    simulation.context.setPositions(positions)
    if prmtop_crds.box_vectors:
        simulation.context.setPeriodicBoxVectors(*prmtop_crds.box_vectors)
    if new_velocities or prmtop_crds.velocities is None:
        simulation.context.setVelocitiesToTemperature(temperature * u.kelvin)
    else:
        simulation.context.setVelocities(
            prmtop_crds.velocities * u.angstrom / u.picosecond)

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


