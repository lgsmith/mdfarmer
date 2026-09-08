#!/usr/bin/env python
"""shakedown-omm: 10 solo trp-cage clones, 5 generations of a few seconds each.

A pre-flight for the whole mdfarmer pipeline on Slurm: seeding, generation
chaining, preemption handling, harvesting. Nothing here is science -- the
generations are deliberately far too short for that. What it proves is that
the tender, the scheduler scripts, the interpreter on the compute node and the
harvest all work before a real campaign is committed to them.

drive_omm.sh is how it is launched: it holds the campaign's tender lock, logs
somewhere findable, and detaches. The driver runs on its own just as well.

    ./drive_omm.sh --check      # readiness table, no writes
    ./drive_omm.sh --dry-run    # dirs, configs, scripts; submit nothing
    ./drive_omm.sh              # start the tender, detached
    ./drive_omm.sh --stop       # graceful stop, at the next tick
"""
import argparse
import gzip
import re
import shutil
import sys
from pathlib import Path

import openmm as mm
from openmm import unit

import mdfarmer as mdf

PROJECT = 'shakedown-omm'
HERE = Path(__file__).resolve().parent
INPUTS = HERE / 'inputs'
PREPARED = HERE / 'prepared'
TRAJ_TOP = HERE / 'data' / PROJECT
# The interpreter a compute node runs a generation with. The tender is already
# in the right environment, so its own python is the honest default; --python
# overrides it when the node needs a different one.
PYTHON_CMD = sys.executable

# ---- run shape --------------------------------------------------------------
# The whole point of the shakedown is that a generation is seconds long but
# still holds a handful of frames, so the two spacing rules mdfarmer enforces
# are exercised for real rather than trivially:
#
#   utilities.check_whole_frames:    STEPS_PER_GEN % WRITE_INTERVAL == 0
#                                    10000 % 2000 == 0
#   harvester.check_commensurability: (STEPS_PER_GEN // WRITE_INTERVAL)
#                                         % DOWNSAMPLE_FRQ == 0
#                                    5 % 5 == 0
#
# 5 frames per generation is the smallest count that still shows a step axis
# marching across a seam, and 5 is the largest downsample that divides it, so
# the wet stream keeps one frame per generation -- enough to prove the global
# frame index does not reset at a boundary. At dt=0.002 ps a generation is
# 20 ps of trajectory at 4 ps/frame; five of them is 100 ps per clone.
DT_PS = 0.002
STEPS_PER_GEN = 10_000           # measured: 4.4 s of MD on an RTX A6000
WRITE_INTERVAL = 2_000           # 4 ps/frame -> 5 frames per generation
DOWNSAMPLE_FRQ = 5               # wet stream keeps 1 frame per generation
N_SEEDS = 1
N_CLONES = 10
ACTIVE_CLONES = 5                # half the clones run at once, so the tender
                                 # has to wait for a slot and fill it again
N_GENS = 5
TEMPERATURE = 277                # K, the temperature the system was built at
FRICTION_PER_PS = 1.0
GEN_SEED_BASE = 42
UPDATE_INTERVAL = 20             # tender tick (s); short, because gens are too
RESTARTS_PER_GEN = 3
SEP = '_'
TRAJ_NAME = 'traj'
TRAJ_SUFFIX = '.dcd'             # the only format seeder can trim on recovery
RESTART_NAME = 'state.xml'

# ---- slurm ------------------------------------------------------------------
# Same shape as sampling-trpcage/farmer.py, with the clock wound down: seconds
# of MD do not need 48 h blocks, and a shakedown that sits in the queue behind
# a 48 h request is not a fast pre-flight.
MEM = '32G'
CPUS = 2                         # OpenMM wants one busy core, plus one spare
WALLTIME = '00:20:00'
PARTITION = 'gpu'
QOS = ''
GRES = 'gpu:rtx_pro_6000_blackwell:1'
EXTRA_SBATCH = ''
HARVEST_PARTITION = 'ccb'        # harvesting is CPU-only
HARVEST_TIME = '00:30:00'

# LOOS syntax, against topology.pdb. resolve_subset refuses an empty selection,
# so a wrong residue name fails loudly rather than harvesting nothing.
HARVESTER_SUBSET = '!(resname == "HOH" || resname == "K" || resname == "CL")'

SCHEDULER_FSTRING = """#!/bin/bash
#SBATCH -J {job_name}
#SBATCH -o slurm.out
#SBATCH -e slurm.out
#SBATCH -p {partition}
#SBATCH --gres={gres}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --time={walltime}
#SBATCH --signal=B:TERM@120
{qos_line}{extra_sbatch}{exclude_nodes}

echo "JOB_NAME: {job_name}"
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo "NODE: $SLURMD_NODENAME"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | paste -sd, -)"
echo "DATE: $(date -Is)"
printf '%s\\t%s\\t%s\\n' "$(date -Is)" "$SLURM_JOB_ID" "$SLURMD_NODENAME" >> node_history.tsv

# The trap is what turns a preemption into a clean stop on a whole frame:
# SentinelReporter watches for PREEMPT_SIGTERM, and the sleep outlives the
# grace period so Slurm records CANCELLED rather than killing python outright.
preempt_handler() {{ touch PREEMPT_SIGTERM; sleep 70; }}
trap preempt_handler SIGTERM

[ -x {python} ] || {{ echo "MISSING INTERPRETER {python}"; exit 1; }}
{python} -c 'import mdfarmer' || {{ echo "mdfarmer is not importable by {python}; check PYTHONPATH"; exit 1; }}
{python} -u {run_script_name} &
wait
"""

HARVEST_FSTRING = """#!/bin/bash
#SBATCH -J harvest
#SBATCH -o harvest-%j.out
#SBATCH -p {harvest_partition}
#SBATCH --time={harvest_time}
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

{python} -u -c 'from mdfarmer import harvest_generation; harvest_generation("config.json", "hconfig.json")'
"""

# stepCount/time on the seed State. OpenMM writes the campaign's absolute step
# into every restart, and seeder._try_recover_gen reads it back as a step
# counted from zero, so a seed carrying its equilibration's 150000 would make
# every resumed generation look finished before it started.
STATE_COUNTER_RE = re.compile(r'stepCount="\d+" time="[^"]*"')
STATE_COUNTER_ZERO = 'stepCount="0" time="0.0"'


def inflate(src_gz, dest, finish=None):
    """Decompress one committed .gz input to `dest` once, and return `dest`.

    The inputs are committed gzipped so the example is self-contained without
    carrying 7 MB of XML, but omm_generation reads system_fn with
    `Path(system_fn).read_text()` and hands seed_fn to `Simulation.loadState`,
    neither of which inflates. So they are inflated once, here, rather than in
    the runner.

    A `dest` already on disk is left alone: every tender boot calls this, and
    re-inflating 7 MB each time would be waste, not safety. The inflated file
    only appears under its real name once it is whole, so a boot killed
    mid-inflate leaves nothing a later boot can mistake for done. `finish` is
    called on the staged file first, for edits that must happen exactly once.
    """
    dest = Path(dest)
    if dest.is_file():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    staged = dest.with_name(dest.name + '.partial')
    with gzip.open(src_gz, 'rb') as fin, staged.open('wb') as fout:
        shutil.copyfileobj(fin, fout)
    if finish is not None:
        finish(staged)
    return staged.replace(dest)


def zero_state_counters(state_p, pattern=STATE_COUNTER_RE,
                        replacement=STATE_COUNTER_ZERO):
    """Rewind the seed State's step and time counters to the origin.

    Edited as text rather than round-tripped through XmlSerializer: a State's
    stepCount is not settable from python, and the positions and velocities
    must come through bit-identical.
    """
    text = state_p.read_text()
    rewound, n = pattern.subn(replacement, text, count=1)
    if n != 1:
        raise SystemExit(f'{state_p} has no stepCount/time attributes to zero')
    state_p.write_text(rewound)
    return state_p


def write_integrator(dest, temperature=TEMPERATURE, dt_ps=DT_PS,
                     friction_per_ps=FRICTION_PER_PS):
    """The integrator the campaign runs, serialised where the config can name it.

    Written here rather than committed so the three numbers that set the step
    arithmetic are visible in the driver next to STEPS_PER_GEN. dt matches the
    GROMACS arm's mdp, and h-bonds are constrained in system.xml, so 2 fs is
    the same step both engines take.
    """
    integrator = mm.LangevinMiddleIntegrator(temperature * unit.kelvin,
                                             friction_per_ps / unit.picosecond,
                                             dt_ps * unit.picoseconds)
    Path(dest).write_text(mm.XmlSerializer.serialize(integrator))
    return Path(dest)


def prepare_inputs(inputs=INPUTS, prepared=PREPARED):
    """Inflate the committed inputs and write the integrator beside them.

    Cheap enough to call on every tender boot: the inflation is skipped once
    the inflated file is there, and the integrator is three numbers.
    """
    prepared = Path(prepared)
    paths = dict(
        system_fn=inflate(inputs / 'system.xml.gz', prepared / 'system.xml'),
        top_fn=inflate(inputs / 'topology.pdb.gz', prepared / 'topology.pdb'),
        # Rewound while staged, so the seed is never on disk under its real
        # name still carrying its equilibration's step count.
        seed_fn=inflate(inputs / 'state.xml.gz', prepared / 'state.xml',
                        finish=zero_state_counters),
    )
    paths['integrator_xml'] = write_integrator(prepared / 'integrator.xml')
    return {key: str(p.resolve()) for key, p in paths.items()}


def config_template(paths, traj_top=TRAJ_TOP, project=PROJECT,
                    steps_per_gen=STEPS_PER_GEN,
                    write_interval=WRITE_INTERVAL, temperature=TEMPERATURE,
                    traj_name=TRAJ_NAME, traj_suffix=TRAJ_SUFFIX,
                    restart_name=RESTART_NAME):
    """The shared per-generation config; Clone fills the per-clone keys in."""
    template = mdf.default_straight_sampling_config_template.copy()
    template.update(
        title=project,
        traj_dir_top_level=str(Path(traj_top).resolve()),
        # Kept beside the data rather than in whatever cwd the tender booted in.
        traj_list=str((Path(traj_top) / 'traj_list.txt').resolve()),
        integrator_xml=paths['integrator_xml'],
        traj_name=traj_name,
        traj_suffix=traj_suffix,
        restart_name=restart_name,
        steps=steps_per_gen,
        # Recorded so the harvest knows the full length even when a resumed
        # generation's own `steps` is only the remainder it still owed.
        steps_per_gen=steps_per_gen,
        write_interval=write_interval,
        temperature=temperature,
        state_data_kwargs=dict(step=True, time=True, potentialEnergy=True,
                               temperature=True, volume=True, speed=True,
                               elapsedTime=True, totalSteps=steps_per_gen),
    )
    return template


def build_harvester(paths, steps_per_gen=STEPS_PER_GEN,
                    write_interval=WRITE_INTERVAL,
                    downsample_frq=DOWNSAMPLE_FRQ, subset=HARVESTER_SUBSET,
                    harvest_partition=HARVEST_PARTITION,
                    harvest_time=HARVEST_TIME, python=PYTHON_CMD):
    """Reduce each finished generation to a dry stream and a downsampled one."""
    mdf.check_commensurability(steps_per_gen, write_interval, downsample_frq)
    return mdf.Harvester(
        harvester_template=HARVEST_FSTRING,
        scheduler='sbatch',
        run_config=dict(harvester_structure=paths['top_fn'],
                        harvester_subset=subset,
                        downsample_frq=downsample_frq,
                        steps_per_gen=steps_per_gen,
                        harvest_partition=harvest_partition,
                        harvest_time=harvest_time,
                        python=python))


def scheduler_kws(partition=PARTITION, qos=QOS, gres=GRES, mem=MEM, cpus=CPUS,
                  walltime=WALLTIME, extra_sbatch=EXTRA_SBATCH,
                  python=PYTHON_CMD):
    return dict(
        partition=partition, gres=gres, cpus=cpus, mem=mem, walltime=walltime,
        qos_line=(f'#SBATCH -q {qos}\n' if qos else ''),
        extra_sbatch=(extra_sbatch + '\n' if extra_sbatch else ''),
        exclude_nodes='', python=python, run_script_name='run.py')


def build_farmer(paths, n_clones=N_CLONES, n_gens=N_GENS,
                 active_clones=ACTIVE_CLONES, traj_top=TRAJ_TOP,
                 project=PROJECT, steps_per_gen=STEPS_PER_GEN,
                 write_interval=WRITE_INTERVAL,
                 downsample_frq=DOWNSAMPLE_FRQ, harvest=True,
                 report_cmd=None, assoc_cmd=None, dry_run=False,
                 restarts_per_gen=RESTARTS_PER_GEN, python=PYTHON_CMD):
    harvester = build_harvester(
        paths, steps_per_gen=steps_per_gen, write_interval=write_interval,
        downsample_frq=downsample_frq, python=python) if harvest else None
    return mdf.Farmer(
        n_seeds=N_SEEDS, n_clones=n_clones, n_gens=n_gens,
        config_template=config_template(
            paths, traj_top=traj_top, project=project,
            steps_per_gen=steps_per_gen, write_interval=write_interval),
        seed_structure_fns=[paths['seed_fn']],
        system_fns=[paths['system_fn']],
        top_fns=[paths['top_fn']],
        scheduler='sbatch',
        scheduler_fstring=SCHEDULER_FSTRING,
        scheduler_kws=scheduler_kws(python=python),
        scheduler_report_cmd=(
            report_cmd or mdf.basic_scheduler_reports['slurm']),
        scheduler_assoc_rep_cmd=(
            assoc_cmd or mdf.basic_scheduler_assoc_reports['slurm']),
        runner=mdf.omm_generation,
        harvester=harvester,
        handle_preempt=True,
        active_clone_threshold=active_clones,
        dirname_pad=2,
        sep=SEP,
        bad_node_persist=str(Path(traj_top) / 'bad_nodes.txt'),
        seed_labels=['native-277'],
        restarts_per_gen=restarts_per_gen,
        jids_file=Path(traj_top) / f'{project}-jids.txt',
        dry_run=dry_run,
    )


def report_readiness(steps_per_gen=STEPS_PER_GEN,
                     write_interval=WRITE_INTERVAL,
                     downsample_frq=DOWNSAMPLE_FRQ, dt_ps=DT_PS,
                     n_clones=N_CLONES, n_gens=N_GENS,
                     active_clones=ACTIVE_CLONES, inputs=INPUTS,
                     project=PROJECT):
    frames = mdf.check_commensurability(steps_per_gen, write_interval,
                                        downsample_frq)
    print(f'{project}: {n_clones} clones x {n_gens} gens, '
          f'{active_clones} running at once')
    print(f'  {steps_per_gen:,} steps/gen ({steps_per_gen * dt_ps:g} ps), '
          f'{frames} frames at {write_interval * dt_ps:g} ps, '
          f'wet every {downsample_frq} -> {write_interval * dt_ps * downsample_frq:g} ps')
    print(f'  {n_gens * frames} dry frames per clone over '
          f'{n_gens * steps_per_gen * dt_ps:g} ps')
    missing = []
    for name in ('system.xml.gz', 'state.xml.gz', 'topology.pdb.gz'):
        path = Path(inputs) / name
        print(f'  {name:18s} ' + ('READY' if path.is_file() else 'MISSING'))
        if not path.is_file():
            missing.append(str(path))
    return missing


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--check', action='store_true',
                    help='print the run shape and input readiness, write nothing')
    ap.add_argument('--dry-run', action='store_true',
                    help='write every directory, config and job script; submit nothing')
    ap.add_argument('--no-harvest', action='store_true')
    ap.add_argument('--n-clones', type=int, default=N_CLONES)
    ap.add_argument('--n-gens', type=int, default=N_GENS)
    ap.add_argument('--steps', type=int, default=STEPS_PER_GEN)
    ap.add_argument('--write-interval', type=int, default=WRITE_INTERVAL)
    ap.add_argument('--traj-top', default=str(TRAJ_TOP))
    ap.add_argument('--update-interval', type=int, default=UPDATE_INTERVAL)
    ap.add_argument('--python', default=PYTHON_CMD,
                    help='interpreter a compute node runs a generation with '
                         f'(default: this tender\'s own, {PYTHON_CMD})')
    args = ap.parse_args()

    missing = report_readiness(steps_per_gen=args.steps,
                               write_interval=args.write_interval,
                               n_clones=args.n_clones, n_gens=args.n_gens)
    if missing:
        raise SystemExit(f'missing committed inputs: {missing}')
    if args.check:
        return

    paths = prepare_inputs()
    print('prepared inputs:')
    for key, value in sorted(paths.items()):
        print(f'  {key:15s} {value}')

    farmer = build_farmer(paths, n_clones=args.n_clones, n_gens=args.n_gens,
                          traj_top=Path(args.traj_top),
                          steps_per_gen=args.steps,
                          write_interval=args.write_interval,
                          harvest=not args.no_harvest, dry_run=args.dry_run,
                          python=args.python)
    finished = farmer.start_tending_fields(update_interval=args.update_interval)
    # The exit status a re-entering tender loop reads: 0 only when every clone
    # finished, so a braked or failed run is re-entered rather than called done.
    raise SystemExit(0 if finished else 1)


if __name__ == '__main__':
    main()
