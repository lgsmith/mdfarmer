"""A tender restart must not copy the campaign's opening structure over a
generation's own restart file.

After an orchestrator restart a clone that recovers nothing falls back to
initial_seed_fn, which lives outside the generation directory. The launch
preparation then copies that seed in under restart_name -- the same name the
generation's own checkpoint has. Every step already run is thrown away, the
generation starts again from zero, and nothing says so.
"""
import shutil
import sys
from pathlib import Path

import harness
from harness import Suite

import mdfarmer
from mdfarmer import gmx_simulate as gs, seeder

STEPS_PER_GEN = 400
WRITE_INTERVAL = 100
PARTIAL_STEPS = 200             # steps the launch gets through before it stops
SEED_INDEX = 0
CLONE_INDEX = 0
DIRNAME_PAD = 2
SEP = '-'
CPU_ARGS = ('-nb', 'cpu', '-pme', 'cpu', '-ntomp', '2')
JOB_NAME = '{title}-{seed_index}-{clone_index}-{gen_index}'


def rebuilt_clone(farm, structure, topology, mdp, template):
    """The Clone a tender rebuilds from disk, as Farmer builds it."""
    return seeder.Clone.from_disk(
        farm, SEED_INDEX, CLONE_INDEX, initial_seed_fn=str(structure),
        top_fn=str(topology), system_fn=str(mdp), structure_fn=str(structure),
        config_overrides={}, config_template=template,
        scheduler='slurm', scheduler_fstring='echo {job_name}',
        scheduler_kws={}, dirname_pad=DIRNAME_PAD, sep=SEP,
        job_number_re='[1-9][0-9]*', job_name_fstring=JOB_NAME,
        recover_fn=gs.gmx_try_recover_gen, progress_fn=gs.gmx_gen_progress,
        run_script=gs.default_gmx_run_script)


def main(gmx_bin=harness.GMX_BIN):
    suite = Suite('restart_seed')
    work = harness.workdir('restart_seed')
    structure, topology, mdp = harness.build_water_system(work, gmx_bin=gmx_bin)
    farm = work / 'farm'

    template = gs.gmx_config_template(
        traj_dir_top_level=str(farm), title='restart-seed',
        structure_fn=str(structure), mdp_fn=str(mdp), dirname_pad=DIRNAME_PAD,
        sep=SEP, traj_name='prod', traj_suffix='.xtc', restart_name='state.cpt',
        steps=STEPS_PER_GEN, steps_per_gen=STEPS_PER_GEN,
        write_interval=WRITE_INTERVAL, temperature=300, gen_seed_base=1,
        gmx_bin=gmx_bin, grompp_maxwarn=3, mdrun_args=CPU_ARGS, maxh=23.5)

    suite.section('a generation is left partway through')
    try:
        gs.gmx_generation(**dict(
            template, seed_index=SEED_INDEX, clone_index=CLONE_INDEX,
            gen_index=0, seed_fn=str(structure), top_fn=str(topology),
            new_velocities=True,
            mdrun_args=CPU_ARGS + ('-nsteps', str(PARTIAL_STEPS))))
    except gs.GenIncomplete:
        pass
    gen_dir = mdfarmer.dir_seeds_clones_gens(
        farm, SEED_INDEX, CLONE_INDEX, 0, DIRNAME_PAD, sep=SEP)
    own_cpt = gen_dir / template['restart_name']
    reached = gs.checkpoint_step(own_cpt, gmx_bin=gmx_bin)
    suite.check('its checkpoint holds the steps it ran',
                gs.is_checkpoint(own_cpt) and reached == PARTIAL_STEPS,
                f'-> step {reached}')
    before = own_cpt.read_bytes()

    suite.section('the tender restarts and prepares the next launch')
    clone = rebuilt_clone(farm, structure, topology, mdp, template)
    outside = clone.current_seed.parent != gen_dir
    suite.check('the rebuilt clone seeds from outside the generation directory',
                outside, f'-> {clone.current_seed}')
    clone.plow_harrow_plant()

    suite.check('the generation keeps its own checkpoint',
                gs.is_checkpoint(own_cpt))
    suite.check('byte for byte', own_cpt.read_bytes() == before)
    suite.check('so the steps it ran are still there',
                gs.is_checkpoint(own_cpt)
                and gs.checkpoint_step(own_cpt, gmx_bin=gmx_bin) == reached,
                f'-> step {reached}')
    suite.check('and the config names the file the job will resume from',
                Path(clone.config['seed_fn']) == own_cpt,
                f"-> {clone.config['seed_fn']}")

    suite.section('a fresh generation directory still gets its seed copied')
    fresh = mdfarmer.dir_seeds_clones_gens(
        farm, SEED_INDEX, CLONE_INDEX, 1, DIRNAME_PAD, sep=SEP, mkdir=True)
    clone.config['gen_index'] = 1
    clone.current_gen_dir = fresh
    clone.set_seed(own_cpt)
    clone.check_copy_set_restart_seed()
    copied = fresh / template['restart_name']
    suite.check('the seed lands in the empty directory',
                copied.is_file() and copied.read_bytes() == before)

    suite.section('an older restart file is replaced, not kept')
    stale = mdfarmer.dir_seeds_clones_gens(
        farm, SEED_INDEX, CLONE_INDEX, 2, DIRNAME_PAD, sep=SEP, mkdir=True)
    old = stale / template['restart_name']
    old.write_bytes(b'stale')
    import os
    os.utime(old, (0, 0))                # older than any real seed
    clone.config['gen_index'] = 2
    clone.current_gen_dir = stale
    clone.set_seed(own_cpt)
    clone.check_copy_set_restart_seed()
    suite.check('a restart file older than the seed is overwritten',
                old.read_bytes() == before)
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
