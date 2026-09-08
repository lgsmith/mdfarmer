"""An OpenMM generation must write its artefacts under its own gen directory.

A solo clone's job runs with cwd set to its generation directory, so a bare
'state.xml' landed in the right place by accident. A packed job gives every
member the one pack directory as cwd, so K replicas would have written one
shared state.xml and clobbered each other's restart state on every write
interval. This suite runs two generations from a cwd that is neither of their
generation directories and pins down that each keeps its own restart file,
its own .out and its own preempt sentinel, and that a caller may still point
several members at one shared sentinel.
"""
import os
import sys
from pathlib import Path

import harness
from harness import Suite

import mdfarmer
from mdfarmer import simulate as sim
from mdfarmer import utilities as util

from test_omm_preempt_progress import (
    build_argon_box, SentinelAfter, PLATFORM_NAME, WRITE_INTERVAL,
)

SEED_INDEX = 0
CLONE_INDEX = 0
DIRNAME_PAD = 2
SEP = '-'
# Different lengths per gen, so a clobbered state.xml cannot pass by accident.
STEPS_BY_GEN = (200, 400)
PREEMPT_STEPS = 400


def gen_config(farm, system_fn, pdb_fn, integrator_fn, seed_fn, gen_index,
               steps, **overrides):
    """The omm_generation kwargs for one generation of the argon box."""
    config = dict(
        traj_dir_top_level=str(farm), system_fn=system_fn, top_fn=pdb_fn,
        seed_index=SEED_INDEX, clone_index=CLONE_INDEX, gen_index=gen_index,
        title='omm-gen-dir', integrator_xml=integrator_fn, seed_fn=seed_fn,
        append=False, dirname_pad=DIRNAME_PAD, sep=SEP, traj_name='positions',
        traj_suffix='.dcd', restart_name='state.xml',
        platform_name=PLATFORM_NAME, steps=steps,
        write_interval=WRITE_INTERVAL)
    config.update(overrides)
    return config


def run_from(cwd, config):
    """Run one generation with the process parked in cwd, returning Preempted."""
    here = Path.cwd()
    try:
        os.chdir(cwd)
        try:
            sim.omm_generation(**config)
        except sim.Preempted:
            return True
        return False
    finally:
        os.chdir(here)


def main():
    suite = Suite('omm_gen_dir_artefacts')
    work = harness.workdir('omm_gen_dir_artefacts')
    system_fn, integrator_fn, pdb_fn, seed_fn = build_argon_box(work)
    farm = work / 'farm'
    # Stands in for a pack directory: one cwd shared by every member.
    pack_dir = work / 'pack'
    pack_dir.mkdir()

    gen_dirs = [mdfarmer.dir_seeds_clones_gens(
        farm, SEED_INDEX, CLONE_INDEX, g, DIRNAME_PAD, sep=SEP)
        for g in range(len(STEPS_BY_GEN))]

    suite.section('two generations run from a cwd that is neither gen dir')
    for gen_index, steps in enumerate(STEPS_BY_GEN):
        run_from(pack_dir, gen_config(farm, system_fn, pdb_fn, integrator_fn,
                                      seed_fn, gen_index, steps))

    suite.check('nothing was written into the shared cwd',
                not list(pack_dir.iterdir()),
                f'-> {sorted(p.name for p in pack_dir.iterdir())}')

    for gen_index, steps in enumerate(STEPS_BY_GEN):
        gen_dir = gen_dirs[gen_index]
        restart_p = gen_dir / 'state.xml'
        suite.check(f'gen {gen_index} kept its own loadable state.xml',
                    util.is_state_xml_usable(restart_p), f'-> {restart_p}')
        step_count = (util.state_xml_step_count(restart_p)
                      if restart_p.is_file() else None)
        suite.check(f'gen {gen_index} state.xml records its own step count',
                    step_count == steps, f'-> {step_count} vs {steps}')
        suite.check(f'gen {gen_index} .out landed beside it',
                    (gen_dir / 'positions.out').is_file())

    suite.section('the default preempt sentinel is the gen dir, not cwd')
    real_sentinel_cls = sim.SentinelReporter
    sim.SentinelReporter = SentinelAfter
    own_dir = mdfarmer.dir_seeds_clones_gens(
        farm, SEED_INDEX, CLONE_INDEX, 2, DIRNAME_PAD, sep=SEP)
    shared_dir = mdfarmer.dir_seeds_clones_gens(
        farm, SEED_INDEX, CLONE_INDEX, 3, DIRNAME_PAD, sep=SEP)
    shared_sentinel = pack_dir / sim.PREEMPT_SENTINEL_NAME
    try:
        own = run_from(pack_dir, gen_config(
            farm, system_fn, pdb_fn, integrator_fn, seed_fn, 2, PREEMPT_STEPS,
            handle_preempt=True))
        suite.check('the preempt reaches the caller', own)
        suite.check('the sentinel was dropped in the gen dir',
                    (own_dir / sim.PREEMPT_SENTINEL_NAME).is_file())
        suite.check('and not in the shared cwd', not shared_sentinel.is_file())

        suite.section('a caller may point a member at one shared sentinel')
        shared = run_from(pack_dir, gen_config(
            farm, system_fn, pdb_fn, integrator_fn, seed_fn, 3, PREEMPT_STEPS,
            handle_preempt=True, sentinel_path=str(shared_sentinel)))
        suite.check('the preempt reaches the caller', shared)
        suite.check('the shared sentinel is the one that was watched',
                    shared_sentinel.is_file())
        suite.check('and the gen dir has none of its own',
                    not (shared_dir / sim.PREEMPT_SENTINEL_NAME).is_file())
        suite.check('the shared gen still kept its own state.xml',
                    util.is_state_xml_usable(shared_dir / 'state.xml'))
    finally:
        sim.SentinelReporter = real_sentinel_cls
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
