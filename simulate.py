import json
import openmm.app as app
import openmm as mm
from openmm import unit
from openmm.app.dcdfile import DCDFile
from openmm.app.xtcfile import XTCFile
from . import utilities as util
from pathlib import Path

from os import environ


class FlushingDCDReporter(app.DCDReporter):
    """DCDReporter that pushes Python's buffer to the kernel after each frame.

    A SIGKILL from Slurm preempt loses bytes still sitting in the BufferedWriter;
    the kernel page cache they are flushed into survives the process death.
    flush() costs microseconds and never blocks on disk. Do not reach for
    buffering=0 instead: DCDFile.writeModel emits many small struct.pack writes
    per frame, each of which would then become its own syscall.
    """

    def report(self, simulation, state):
        super().report(simulation, state)
        try:
            self._dcd._file.flush()
        except AttributeError:
            raise RuntimeError(
                'FlushingDCDReporter could not reach self._dcd._file; '
                'OpenMM internals likely changed. Update the wrapper.'
            )


def _flush_dcd_file(reporter):
    """Flush the underlying file handle of a reporter that owns a DCDFile at self._dcd."""
    try:
        reporter._dcd._file.flush()
    except AttributeError:
        raise RuntimeError(
            f'{type(reporter).__name__} could not reach self._dcd._file; '
            'OpenMM internals likely changed. Update the wrapper.'
        )


class _TandemDCDReporter:
    """Write velocities or forces into the position slot of a DCD file.

    DCD carries no velocity or force field, so the position slot is repurposed
    and the file rides alongside the position trajectory frame-for-frame. The
    numbers on disk are in the quantity's native units — nm/ps for velocities,
    kJ/(mol·nm) for forces — but labelled nanometers, so a reader has to know
    which quantity a given file holds.
    """

    def __init__(self, file, reportInterval, quantity, append=False,
                 enforcePeriodicBox=None):
        """
        Parameters
        ----------
        quantity : 'velocities' or 'forces'
        """
        if quantity not in ('velocities', 'forces'):
            raise ValueError(f"quantity must be 'velocities' or 'forces', got {quantity!r}")
        self._quantity = quantity
        self._reportInterval = reportInterval
        self._enforcePeriodicBox = enforcePeriodicBox
        self._append = append
        self._dcd = None
        self._out = open(file, 'r+b' if append else 'wb')

    def describeNextReport(self, simulation):
        # OpenMM 8.x dict format; 'include' limits what getState() populates.
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return {'steps': steps, 'periodic': self._enforcePeriodicBox,
                'include': [self._quantity]}

    def report(self, simulation, state):
        if self._dcd is None:
            self._dcd = DCDFile(
                self._out, simulation.topology,
                simulation.integrator.getStepSize(),
                self._reportInterval, self._reportInterval, self._append
            )
        # writeModel dimension-checks via .value_in_unit(nanometers), which dies
        # on velocity/force Quantities, so re-tag the raw numbers as nanometers.
        if self._quantity == 'velocities':
            raw = state.getVelocities(asNumpy=True).value_in_unit(
                unit.nanometer / unit.picosecond)
        else:
            raw = state.getForces(asNumpy=True).value_in_unit(
                unit.kilojoule_per_mole / unit.nanometer)
        vectors = raw * unit.nanometer
        self._dcd.writeModel(vectors, periodicBoxVectors=state.getPeriodicBoxVectors())
        _flush_dcd_file(self)

    def __del__(self):
        self._out.close()


class _TandemXTCReporter:
    """Write velocities or forces into the position slot of an XTC file.

    Same repurposing and same unit labelling as _TandemDCDReporter, but XTC's
    writer additionally requires abs(value) * 1000 to fit in int32 (~2.1e6 nm
    after scaling). Velocities (a few nm/ps) and bonded forces (up to ~1e5
    kJ/mol/nm) both clear that bound, but XTC compression is lossy — use .dcd
    if you need full precision on saved velocities or forces.
    """

    def __init__(self, file, reportInterval, quantity, append=False,
                 enforcePeriodicBox=None):
        if quantity not in ('velocities', 'forces'):
            raise ValueError(f"quantity must be 'velocities' or 'forces', got {quantity!r}")
        self._quantity = quantity
        self._reportInterval = reportInterval
        self._enforcePeriodicBox = enforcePeriodicBox
        self._append = append
        self._xtc = None
        self._fileName = file
        if not append:
            open(file, 'wb').close()

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return {'steps': steps, 'periodic': self._enforcePeriodicBox,
                'include': [self._quantity]}

    def report(self, simulation, state):
        if self._xtc is None:
            self._xtc = XTCFile(
                self._fileName, simulation.topology,
                simulation.integrator.getStepSize(),
                self._reportInterval, self._reportInterval, self._append
            )
        # writeModel calls .value_in_unit(nanometers) here too; re-tag as above.
        if self._quantity == 'velocities':
            raw = state.getVelocities(asNumpy=True).value_in_unit(
                unit.nanometer / unit.picosecond)
        else:
            raw = state.getForces(asNumpy=True).value_in_unit(
                unit.kilojoule_per_mole / unit.nanometer)
        vectors = raw * unit.nanometer
        self._xtc.writeModel(vectors, periodicBoxVectors=state.getPeriodicBoxVectors())


class _TandemHDF5Reporter:
    """Write velocities or forces to an HDF5 trajectory file.

    Velocities go into the native velocities field; forces are shoehorned into
    the coordinates field, since HDF5TrajectoryFile.write has no forces kwarg.
    """

    def __init__(self, file, reportInterval, quantity, append=False,
                 enforcePeriodicBox=None):
        if quantity not in ('velocities', 'forces'):
            raise ValueError(f"quantity must be 'velocities' or 'forces', got {quantity!r}")
        try:
            from mdtraj.formats import HDF5TrajectoryFile
        except ImportError as e:
            raise ImportError('mdtraj is required for HDF5 reporter support') from e
        self._quantity = quantity
        self._reportInterval = reportInterval
        self._enforcePeriodicBox = enforcePeriodicBox
        self._is_initialized = False
        mode = 'a' if append else 'w'
        self._traj_file = HDF5TrajectoryFile(file, mode)

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return {'steps': steps, 'periodic': self._enforcePeriodicBox,
                'include': [self._quantity]}

    def _initialize(self, simulation):
        from mdtraj import Topology as MDTopology
        self._traj_file.topology = MDTopology.from_openmm(simulation.topology)
        self._is_initialized = True

    def report(self, simulation, state):
        if not self._is_initialized:
            self._initialize(simulation)
        if self._quantity == 'velocities':
            vels = state.getVelocities(asNumpy=True)
            # nm/ps — matches HDF5TrajectoryFile's expected velocity units
            vels_nm_ps = vels.value_in_unit(unit.nanometer / unit.picosecond)
            self._traj_file.write(coordinates=vels_nm_ps, velocities=None)
        else:
            forces = state.getForces(asNumpy=True)
            # kJ/(mol·nm) numbers, written into the coordinate slot.
            forces_kj = forces.value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
            self._traj_file.write(coordinates=forces_kj)
        if hasattr(self._traj_file, 'flush'):
            self._traj_file.flush()

    def __del__(self):
        self._traj_file.close()


_TANDEM_REPORTER_CLS = {
    '.dcd': _TandemDCDReporter,
    '.xtc': _TandemXTCReporter,
    '.h5':  _TandemHDF5Reporter,
    # No .trr: mdtraj's TRRTrajectoryFile.write() accepts only xyz (positions).
}

_SUPPORTED_TANDEM_SUFFIXES = set(_TANDEM_REPORTER_CLS)


# Touched by the batch script's SIGTERM trap when Slurm preempts the job.
PREEMPT_SENTINEL_NAME = 'PREEMPT_SIGTERM'


class Preempted(Exception):
    pass


class SentinelReporter:
    """Raise Preempted once the batch script's trap handler drops the sentinel.

    Checked on each write_interval cycle, and appended LAST in the reporter list
    so the position / state / data writers for that cycle have already fired and
    produced aligned on-disk output before the simulation loop unwinds.
    """

    def __init__(self, reportInterval, sentinel_path=None):
        self._reportInterval = reportInterval
        self._sentinel = Path(sentinel_path or PREEMPT_SENTINEL_NAME)

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return {'steps': steps, 'periodic': None, 'include': []}

    def report(self, simulation, state):
        if self._sentinel.is_file():
            raise Preempted(f'preempt sentinel detected at {self._sentinel.resolve()}')


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
                   traj_name='positions',
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
                   new_velocities=False,
                   # Basename for the parallel velocity trajectory (no suffix).
                   velocity_name='velocities',
                   # Basename for the parallel force trajectory (no suffix).
                   force_name='forces',
                   # Write velocities into the main trajectory file.
                   # Only supported for traj_suffix='.h5' (mdtraj HDF5Reporter).
                   embed_velocities=False,
                   # Writing forces into the main trajectory is not supported by any
                   # mdtraj reporter; always raises ValueError if True.
                   embed_forces=False,
                   # Extension for a parallel velocity file. One of '.dcd', '.xtc', '.h5'.
                   # None means no parallel velocity file is written.
                   velocity_traj_suffix=None,
                   # Extension for a parallel force file. One of '.dcd', '.xtc', '.h5'.
                   # None means no parallel force file is written.
                   force_traj_suffix=None,
                   # If True, install a SentinelReporter that watches for a
                   # PREEMPT_SIGTERM file in cwd (touched by the batch
                   # script's SIGTERM trap on Slurm preempt) and raises
                   # Preempted at the next write_interval cycle. Requires the
                   # batch script to install the trap and background+wait
                   # the python invocation; see basic_scheduler_fstrings_preempt.
                   handle_preempt=False,
                   ):

    # Validate embedded-output requests up front.
    if embed_forces:
        raise ValueError(
            'embed_forces=True is not supported: no mdtraj reporter exposes a forces '
            'field in its write() API. Use force_traj_suffix for a parallel force file.'
        )
    if embed_velocities and traj_suffix != '.h5':
        raise ValueError(
            f'embed_velocities=True requires traj_suffix=".h5" (got {traj_suffix!r}). '
            'mdtraj HDF5Reporter is the only reporter that can embed velocities in the '
            'main trajectory file.'
        )
    if velocity_traj_suffix is not None and velocity_traj_suffix not in _SUPPORTED_TANDEM_SUFFIXES:
        raise ValueError(
            f'velocity_traj_suffix={velocity_traj_suffix!r} is not supported. '
            f'Choose from: {sorted(_SUPPORTED_TANDEM_SUFFIXES)}'
        )
    if force_traj_suffix is not None and force_traj_suffix not in _SUPPORTED_TANDEM_SUFFIXES:
        raise ValueError(
            f'force_traj_suffix={force_traj_suffix!r} is not supported. '
            f'Choose from: {sorted(_SUPPORTED_TANDEM_SUFFIXES)}'
        )

    # make reporter by extension
    reporters = {
        '.dcd': FlushingDCDReporter,
        '.xtc': app.XTCReporter,
    }
    try:
        import mdtraj.reporters as _mdt_reporters
        reporters['.h5'] = _mdt_reporters.HDF5Reporter
    except ImportError:
        pass

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
        if traj_suffix == '.h5':
            h5_cls = reporters['.h5']
            # HDF5Reporter does not accept append via init kwarg — open in append
            # mode by passing an already-open HDF5TrajectoryFile if appending.
            if traj_path.is_file() and append:
                from mdtraj.formats import HDF5TrajectoryFile
                h5_file = HDF5TrajectoryFile(str(traj_path), 'a')
            else:
                h5_file = str(traj_path)
            traj_reporter = h5_cls(
                h5_file,
                write_interval,
                coordinates=True,
                time=True,
                cell=True,
                potentialEnergy=False,
                kineticEnergy=False,
                temperature=False,
                velocities=embed_velocities,
            )
        elif traj_path.is_file():
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

    # Build tandem velocity/force reporters if requested. On a resume
    # (append=True), only activate the tandem reporter if its file is
    # frame-aligned with the position trajectory — same frame count.
    # Skip cases:
    #   - file doesn't exist (pre-velocity gen).
    #   - file exists but has fewer frames than the position traj
    #     (crashed mid-gen — e.g. the 0-frame `velocities.dcd` left by
    #     the nm/ps-vs-nm units bug). Appending would write frame N of
    #     velocity while position writes frame N+k, permanently offset
    #     for the rest of the gen.
    # Skipping preserves the user's invariant that velocity frame N
    # corresponds to position frame N. Next fresh gen creates a clean
    # from-frame-0 tandem file.
    extra_reporters = []

    def _tandem_aligned(tandem_path):
        if not tandem_path.is_file():
            return False
        if not traj_path.is_file():
            return False
        try:
            return util.get_traj_len(str(tandem_path), top_fn) == \
                util.get_traj_len(str(traj_path), top_fn)
        except Exception as exc:
            print(f'Could not measure frame count of {tandem_path}: '
                  f'{exc}; treating as not aligned and skipping.')
            return False

    if velocity_traj_suffix is not None:
        vel_path = (traj_dir / velocity_name).with_suffix(velocity_traj_suffix)
        if append and not _tandem_aligned(vel_path):
            print(f'Skipping velocity reporter at {vel_path}: file '
                  f'absent or not frame-aligned with {traj_path}. The '
                  'next fresh gen will start a clean velocity trajectory.')
        else:
            cls = _TANDEM_REPORTER_CLS[velocity_traj_suffix]
            extra_reporters.append(cls(str(vel_path), write_interval, 'velocities', append=append))
    if force_traj_suffix is not None:
        force_path = (traj_dir / force_name).with_suffix(force_traj_suffix)
        if append and not _tandem_aligned(force_path):
            print(f'Skipping force reporter at {force_path}: file '
                  f'absent or not frame-aligned with {traj_path}. The '
                  'next fresh gen will start a clean force trajectory.')
        else:
            cls = _TANDEM_REPORTER_CLS[force_traj_suffix]
            extra_reporters.append(cls(str(force_path), write_interval, 'forces', append=append))

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
    for r in extra_reporters:
        simulation.reporters.append(r)
    # SentinelReporter must be appended LAST so DCD / state.xml / .out
    # writes for the current cycle have already landed on disk before it
    # raises. Stale sentinel from a previous preempted run in this gen
    # dir would fire immediately, so clear it first.
    if handle_preempt:
        sentinel_p = Path(PREEMPT_SENTINEL_NAME)
        if sentinel_p.exists():
            sentinel_p.unlink()
        simulation.reporters.append(SentinelReporter(write_interval))
    try:
        simulation.step(steps)
    except Preempted as exc:
        print(f'Preempt received: {exc}; exiting cleanly at last reporter '
              f'cycle. The partial gen will resume via append on next launch.')
        raise
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
    try:
        new_traj_path = omm_generation(**conf_dict)
    except Preempted:
        # Skip the traj_list append: the gen is incomplete. The orchestrator
        # will detect the partial DCD on the next boot and resume in append
        # mode. Exit 0 so Slurm records the job as cancelled, not failed.
        return
    with traj_list_path.open('a') as tl:
        tl.write(str(new_traj_path) + '\n')
