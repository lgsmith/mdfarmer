#!/usr/bin/env python
"""shakedown-gmx: 5 packs of 2 trp-cage clones, 5 generations of a few seconds.

The GROMACS half of the shakedown, on the same molecular system as
examples/openmm-trpcage but with two replicas sharing one card under CUDA MPS.
It proves the packed path end to end -- one sbatch, one MPS daemon, one pack
lock, two pinned mdruns, per-member failure and recovery, and a harvest per
generation -- before a real campaign is committed to it.

drive_gmx.sh is how it is launched: it holds the campaign's tender lock, logs
somewhere findable, and detaches. The driver runs on its own just as well.

    ./drive_gmx.sh --check      # readiness table, no writes
    ./drive_gmx.sh --dry-run    # dirs, configs, scripts; submit nothing
    ./drive_gmx.sh              # start the tender, detached
    ./drive_gmx.sh --stop       # graceful stop, at the next tick
"""
import argparse
import gzip
import shutil
import sys
from pathlib import Path

import mdfarmer as mdf

PROJECT = 'shakedown-gmx'
HERE = Path(__file__).resolve().parent
INPUTS = HERE / 'inputs'
PREPARED = HERE / 'prepared'
TRAJ_TOP = HERE / 'data' / PROJECT
# The interpreter a compute node runs a generation with. The tender is already
# in the right environment, so its own python is the honest default; --python
# overrides it when the node needs a different one.
PYTHON_CMD = sys.executable

# Whose presence stops the tender, and the status it exits with when it does.
# Distinct from 1 so the launching shell can tell "asked to stop" from "died".
BRAKE_FILE = mdf.farmer.BRAKE_FILE
BRAKED_EXIT = 2
# Every CUDA GROMACS in this module tree is an MPI build, and such a binary
# calls MPI_Init even for grompp -- unconditionally, from main(), before any
# subcommand is dispatched. Inside a Slurm step that init finds SLURM_STEP_ID,
# looks for a PMIx server a plain batch step does not run, and aborts before
# any MD. Measured here: SLURM_STEP_ID alone is the trigger, and removing it
# alone is the cure.
#
# So the binary runs with those variables taken out of its environment. `env`
# execs rather than forks, so the process mdfarmer waits on IS gmx_mpi and the
# preemption SIGTERM reaches GROMACS itself -- which mpirun -n 1, the other way
# to satisfy this MPI, would put a launcher in front of. The PMIX_/PMI_ entries
# are insurance for a step launched by another PMI plugin.
#
# A site with a thread-MPI gmx wants the bare string instead:
#
#     GMX_BIN = 'gmx'
#
# which is the default, and is why nothing in mdfarmer knows what any of this
# is: gmx_bin is a command vector, and what goes in it is the site's business.
# The README's "Naming the gmx binary" covers both forms.
GMX_SCRUB = ['env',
             '-u', 'SLURM_STEP_ID', '-u', 'SLURM_STEPID',
             '-u', 'PMIX_NAMESPACE', '-u', 'PMIX_RANK',
             '-u', 'PMIX_SERVER_URI41', '-u', 'PMIX_SERVER_URI3',
             '-u', 'PMI_FD', '-u', 'PMI_RANK', '-u', 'PMI_SIZE']
GMX_BINARY = 'gmx_mpi'           # the binary itself, for the job's guard
GMX_BIN = [*GMX_SCRUB, GMX_BINARY]

# Loaded inside the job, since a compute node inherits no module environment
# worth relying on. gmx_pack probes for -ntmpi rather than assuming it, so an
# MPI build is fine here.
ENV_SETUP = """source /etc/profile.d/modules.sh
module load modules/2.4-20250724 openmpi/cuda-4.1.8 gromacs/mpi-2024.4"""

# ---- run shape --------------------------------------------------------------
# Identical arithmetic to the OpenMM arm, so the two datasets are comparable
# frame for frame. Both spacing rules mdfarmer enforces are satisfied by
# construction rather than by luck:
#
#   utilities.check_whole_frames:     10000 % 2000 == 0
#   harvester.check_commensurability: (10000 // 2000) % 5 == 0
#
# 5 frames per generation is the smallest count that still shows a step axis
# marching across a seam, and 5 is the largest downsample that divides it. At
# dt=0.002 ps a generation is 20 ps at 4 ps/frame; five of them is 100 ps.
#
# GROMACS writes a frame at the step it restarts from, so generation 0 holds 5
# frames and every later generation holds 6 on disk; the harvest drops that
# duplicate seam frame exactly once. Both counts are what
# harvester.resolve_seam expects, which is why 5 has to be the NEW-frame count.
DT_PS = 0.002
STEPS_PER_GEN = 10_000           # measured: 6.6 s wall on an RTX A6000
WRITE_INTERVAL = 2_000           # 4 ps/frame -> 5 new frames per generation
DOWNSAMPLE_FRQ = 5               # wet stream keeps 1 frame per generation
N_SEEDS = 1
N_CLONES = 10
N_GENS = 5
TEMPERATURE = 277                # K, the temperature the system was built at
GEN_SEED_BASE = 42
UPDATE_INTERVAL = 20             # tender tick (s); short, because gens are too
RESTARTS_PER_GEN = 3
SEP = '_'
TRAJ_NAME = 'traj'
TRAJ_SUFFIX = '.xtc'
RESTART_NAME = 'state.cpt'

# ---- slurm ------------------------------------------------------------------
# Same shape as sampling-trpcage/farmer.py, with the clock wound down: seconds
# of MD do not need 48 h blocks, and maxh is the backstop that stops mdrun
# inside the allocation rather than the thing that ends a generation.
REPS_PER_PACK = 2
ACTIVE_PACKS = 3                 # 3 of 5 packs at once, so the tender has to
                                 # wait for a slot and refill it
MEM = '32G'
WALLTIME = '00:20:00'
MAXH = 0.25                      # 15 min, comfortably inside the 20 min block
PARTITION = 'gpu'
QOS = ''
# NOT the Blackwell card the OpenMM arm uses. This GROMACS module is built
# --generate-code code=sm_70;sm_80;sm_90 with no code=compute_XX among them,
# so it embeds no PTX and cannot JIT for an architecture it was not built for.
# A cubin runs only within its own major arch, and RTX PRO 6000 Blackwell is
# sm_120, so mdrun would fail there with "no kernel image is available" -- and
# only once it reached the GPU, after grompp had already succeeded. The A100 is
# sm_80 and covered; h100_pcie (sm_90) is the other option on this partition.
# OpenMM is unaffected and keeps the Blackwell card: it compiles its kernels at
# runtime rather than shipping cubins.
GRES = 'gpu:a100-sxm4-80gb:1'
EXTRA_SBATCH = ''
HARVEST_PARTITION = 'ccb'        # harvesting is CPU-only
HARVEST_TIME = '00:30:00'

# The water here is 4-site (TIP4P-ice: one virtual site per molecule), and
# GROMACS refuses `-update gpu` outright with virtual sites -- measured, not
# assumed: "Update task on the GPU was required, but ... Virtual sites are not
# supported." So this is sampling-trpcage's `vitrification` arm, not its
# `fulloffload` one, and it takes that arm's core budget with it: the update
# runs on the CPU, so the per-replica knee is 12 cores, not 4. Measured on this
# system: 574 ns/day at 4 cores, 785 ns/day at 12.
UPDATE_MODE = 'cpu'
CORES_PER_REPLICA = 12
PACK_CPUS = CORES_PER_REPLICA * REPS_PER_PACK

# -notunepme keeps the real/reciprocal split where grompp put it. Thread and
# pin flags are added per replica by gmx_pack, which is what stops two members
# both starting at core 0.
BASE_MDRUN_ARGS = ('-nb', 'gpu', '-bonded', 'gpu', '-pme', 'gpu',
                   '-update', UPDATE_MODE, '-nstlist', '200', '-notunepme')

# LOOS syntax, against the .gro. resolve_subset refuses an empty selection, so
# a wrong residue name fails loudly rather than harvesting nothing.
HARVESTER_SUBSET = '!(resname == "HOH" || resname == "K" || resname == "CL")'

PACK_FSTRING = """#!/bin/bash
#SBATCH -J {job_name}
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.out
#SBATCH -p {partition}
#SBATCH --gres={gres}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --time={walltime}
#SBATCH --signal=B:TERM@120
{qos_line}{extra_sbatch}{exclude_nodes}

# One mdrun set per pack dir: two mdruns sharing a checkpoint destroy the
# trajectory, and no checkpoint can undo it. Held on fd 9 until the job ends.
exec 9>pack.lock
flock -n 9 || {{ echo "REFUSING TO RUN: pack.lock is held by another job"; exit 0; }}

echo "JOB_NAME: {job_name}"
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo "NODE: $SLURMD_NODENAME"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | paste -sd, -)"
echo "DATE: $(date -Is)"
printf '%s\\t%s\\t%s\\n' "$(date -Is)" "$SLURM_JOB_ID" "$SLURMD_NODENAME" >> node_history.tsv

{env_setup}
for cmd in {gmx_check}; do
  command -v "$cmd" >/dev/null || {{ echo "MISSING $cmd after env setup"; exit 1; }}
done

# Per-job MPS daemon, keyed on the job id so two packed jobs on one node never
# share or clobber each other's. Without it the replicas time-slice the card.
export CUDA_MPS_PIPE_DIRECTORY="/tmp/mps-$USER-$SLURM_JOB_ID/pipe"
export CUDA_MPS_LOG_DIRECTORY="/tmp/mps-$USER-$SLURM_JOB_ID/log"
mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
nvidia-cuda-mps-control -d && echo "MPS: daemon up" \\
  || echo "MPS: WARNING daemon FAILED to start -- replicas will time-slice the GPU"

cleanup() {{
  echo quit | nvidia-cuda-mps-control 2>/dev/null || true
  rm -rf "/tmp/mps-$USER-$SLURM_JOB_ID"
  # BadNodeRegistry looks for slurm.out in each gen dir, which a packed job
  # never writes there.
  for cfg in $(grep -o '"/[^"]*config\\.json"' pack.json | tr -d '"'); do
    cp -f "slurm-$SLURM_JOB_ID.out" "$(dirname "$cfg")/slurm.out" 2>/dev/null || true
  done
}}
preempt_handler() {{ touch PREEMPT_SIGTERM; sleep 70; }}
trap preempt_handler SIGTERM
trap cleanup EXIT

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


def inflate(src_gz, dest):
    """Decompress one committed .gz input to `dest` once, and return `dest`.

    The inputs are committed gzipped so the example is self-contained without
    carrying a megabyte of .gro and .top. GROMACS cannot read them that way --
    grompp takes file names, not streams -- so they are inflated once, here,
    before the Farmer is built.

    A `dest` already on disk is left alone: every tender boot calls this, and
    re-inflating each time would be waste, not safety. The inflated file only
    appears under its real name once it is whole, so a boot killed mid-inflate
    leaves nothing a later boot can mistake for done.
    """
    dest = Path(dest)
    if dest.is_file():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    staged = dest.with_name(dest.name + '.partial')
    with gzip.open(src_gz, 'rb') as fin, staged.open('wb') as fout:
        shutil.copyfileobj(fin, fout)
    return staged.replace(dest)


def prepare_inputs(inputs=INPUTS, prepared=PREPARED):
    """Inflate the committed inputs; the .mdp is small enough to commit plain.

    Cheap enough to call on every tender boot: the inflation is skipped once
    the inflated file is there.
    """
    paths = dict(
        structure_fn=inflate(inputs / 'gmx.gro.gz', Path(prepared) / 'gmx.gro'),
        top_fn=inflate(inputs / 'gmx.top.gz', Path(prepared) / 'gmx.top'),
        mdp_fn=Path(inputs) / 'prod-277.mdp',
    )
    return {key: str(p.resolve()) for key, p in paths.items()}


def config_template(paths, traj_top=TRAJ_TOP, project=PROJECT,
                    steps_per_gen=STEPS_PER_GEN,
                    write_interval=WRITE_INTERVAL, temperature=TEMPERATURE,
                    maxh=MAXH, gmx_bin=GMX_BIN, traj_name=TRAJ_NAME,
                    traj_suffix=TRAJ_SUFFIX, restart_name=RESTART_NAME,
                    base_mdrun_args=BASE_MDRUN_ARGS,
                    gen_seed_base=GEN_SEED_BASE):
    """The shared per-generation config; Clone fills the per-clone keys in."""
    return mdf.gmx_config_template(
        title=project,
        traj_dir_top_level=str(Path(traj_top).resolve()),
        # Kept beside the data rather than in whatever cwd the tender booted in.
        traj_list=str((Path(traj_top) / 'traj_list.txt').resolve()),
        structure_fn=paths['structure_fn'],
        mdp_fn=paths['mdp_fn'],
        restart_name=restart_name,
        traj_name=traj_name,
        traj_suffix=traj_suffix,
        deffnm=traj_name,        # so the parts are traj.partNNNN.* too
        steps=steps_per_gen,
        steps_per_gen=steps_per_gen,
        write_interval=write_interval,
        temperature=temperature,
        gen_seed_base=gen_seed_base,
        maxh=maxh,
        gmx_bin=gmx_bin,
        grompp_maxwarn=3,
        mdrun_args=list(base_mdrun_args),
        append=False,
    )


def build_harvester(paths, steps_per_gen=STEPS_PER_GEN,
                    write_interval=WRITE_INTERVAL,
                    downsample_frq=DOWNSAMPLE_FRQ, subset=HARVESTER_SUBSET,
                    harvest_partition=HARVEST_PARTITION,
                    harvest_time=HARVEST_TIME, python=PYTHON_CMD):
    """Reduce each finished generation to a dry stream and a downsampled one.

    harvester_structure is REQUIRED here: top_fn is a force field topology and
    neither LOOS nor mdtraj builds a model from one.
    """
    mdf.check_commensurability(steps_per_gen, write_interval, downsample_frq)
    return mdf.Harvester(
        harvester_template=HARVEST_FSTRING,
        scheduler='sbatch',
        run_config=dict(harvester_structure=paths['structure_fn'],
                        harvester_subset=subset,
                        downsample_frq=downsample_frq,
                        steps_per_gen=steps_per_gen,
                        harvest_partition=harvest_partition,
                        harvest_time=harvest_time,
                        python=python))


def scheduler_kws(partition=PARTITION, qos=QOS, gres=GRES, mem=MEM,
                  walltime=WALLTIME, cpus=PACK_CPUS, extra_sbatch=EXTRA_SBATCH,
                  python=PYTHON_CMD, gmx_bin=GMX_BIN, env_setup=ENV_SETUP,
                  gmx_binary=GMX_BINARY):
    # The guard checks the binary by name, so a command vector cannot reach the
    # shell as a python list -- which would fail every job at the guard.
    checks = [gmx_binary]
    return dict(
        partition=partition, gres=gres, cpus=cpus, mem=mem, walltime=walltime,
        qos_line=(f'#SBATCH -q {qos}\n' if qos else ''),
        extra_sbatch=(extra_sbatch + '\n' if extra_sbatch else ''),
        exclude_nodes='', python=python, gmx_bin=gmx_bin,
        gmx_check=' '.join(checks),
        env_setup=env_setup, run_script_name='run.py')


def build_farmer(paths, n_clones=N_CLONES, n_gens=N_GENS,
                 active_packs=ACTIVE_PACKS, reps_per_pack=REPS_PER_PACK,
                 pack_cpus=PACK_CPUS, traj_top=TRAJ_TOP, project=PROJECT,
                 steps_per_gen=STEPS_PER_GEN, write_interval=WRITE_INTERVAL,
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
        seed_structure_fns=[paths['structure_fn']],
        system_fns=[paths['mdp_fn']],      # repurposed as the .mdp
        top_fns=[paths['top_fn']],
        scheduler='sbatch',
        # Packing submits the pack template and never the solo one; this is
        # what a member's own gen dir gets written with.
        scheduler_fstring=PACK_FSTRING,
        scheduler_kws=scheduler_kws(cpus=pack_cpus, python=python),
        scheduler_report_cmd=(
            report_cmd or mdf.basic_scheduler_reports['slurm']),
        scheduler_assoc_rep_cmd=(
            assoc_cmd or mdf.basic_scheduler_assoc_reports['slurm']),
        runner=mdf.gmx_generation,
        harvester=harvester,
        handle_preempt=True,
        active_clone_threshold=active_packs,   # counts packs, not clones
        dirname_pad=2,
        sep=SEP,
        bad_node_persist=str(Path(traj_top) / 'bad_nodes.txt'),
        seed_labels=['native-277'],
        restarts_per_gen=restarts_per_gen,
        jids_file=Path(traj_top) / f'{project}-jids.txt',
        pack_size=reps_per_pack,
        pack_cpus_per_task=pack_cpus,
        pack_scheduler_fstring=PACK_FSTRING,
        dry_run=dry_run,
    )


def report_readiness(steps_per_gen=STEPS_PER_GEN,
                     write_interval=WRITE_INTERVAL,
                     downsample_frq=DOWNSAMPLE_FRQ, dt_ps=DT_PS,
                     n_clones=N_CLONES, n_gens=N_GENS,
                     reps_per_pack=REPS_PER_PACK, active_packs=ACTIVE_PACKS,
                     pack_cpus=PACK_CPUS, inputs=INPUTS, project=PROJECT):
    frames = mdf.check_commensurability(steps_per_gen, write_interval,
                                        downsample_frq)
    n_packs = -(-n_clones // reps_per_pack)
    print(f'{project}: {n_clones} clones in {n_packs} packs of '
          f'{reps_per_pack} x {n_gens} gens, {active_packs} packs at once')
    print(f'  {steps_per_gen:,} steps/gen ({steps_per_gen * dt_ps:g} ps), '
          f'{frames} new frames at {write_interval * dt_ps:g} ps '
          f'({frames + 1} on disk every gen, the seam included), '
          f'wet every {downsample_frq} -> {write_interval * dt_ps * downsample_frq:g} ps')
    # One more than the generations contribute between them: GROMACS writes a
    # frame at the step it restarts from, so generation 0's own first frame --
    # the seed state, at step 0 -- is a frame no later generation repeats.
    print(f'  {n_gens * frames + 1} dry frames per clone over '
          f'{n_gens * steps_per_gen * dt_ps:g} ps')
    print(f'  {pack_cpus} cores/pack, {pack_cpus // reps_per_pack} per '
          f'replica, -update {UPDATE_MODE}')
    missing = []
    for name in ('gmx.gro.gz', 'gmx.top.gz', 'prod-277.mdp'):
        path = Path(inputs) / name
        print(f'  {name:16s} ' + ('READY' if path.is_file() else 'MISSING'))
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
    ap.add_argument('--brake-file', default=BRAKE_FILE,
                    help='stop at the next tick once this file exists '
                         f'(default: {BRAKE_FILE}, relative to the working '
                         'directory the tender was started in)')
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
    finished = farmer.start_tending_fields(update_interval=args.update_interval,
                                           brake_file=args.brake_file)
    # Three outcomes, three statuses, so the shell that launched this needs no
    # opinion about brake files: 0 every clone finished, BRAKED_EXIT someone
    # asked it to stop, 1 something went wrong and re-entering is the recovery.
    if finished:
        raise SystemExit(0)
    raise SystemExit(BRAKED_EXIT if Path(args.brake_file).is_file() else 1)


if __name__ == '__main__':
    main()
