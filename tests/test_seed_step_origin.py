"""A seed that arrives already carrying a step count must not skew recovery.

An equilibrated OpenMM seed's state.xml records the step its equilibration
reached -- 150000 is typical -- and every later stepCount counts on from there.
Recovery compares that number against the step this generation started at,
which it counts from the campaign's own zero. Without subtracting what the seed
brought, the difference is too large by the whole equilibration, so a barely
started generation looks finished and recovery advances past unrun work.
"""
import json
import sys
from pathlib import Path

import harness
from harness import Suite

import mdfarmer
from mdfarmer import seeder, utilities as util

STEPS_PER_GEN = 1000
WRITE_INTERVAL = 100
# What an equilibration left on the seed before this campaign began.
SEED_ORIGIN = 150000
SEED_INDEX = 0
CLONE_INDEX = 0
DIRNAME_PAD = 2
SEP = '-'


def write_state_xml(path, step_count, time_ps=0.0):
    """A state.xml with only the attributes the step reader looks at."""
    path.write_text(
        f'<?xml version="1.0" ?>\n'
        f'<State openmmVersion="8.5" time="{time_ps}" '
        f'stepCount="{step_count}">\n</State>\n')
    return path


def main(steps_per_gen=STEPS_PER_GEN, write_interval=WRITE_INTERVAL,
         seed_origin=SEED_ORIGIN):
    suite = Suite('seed_step_origin')
    work = harness.workdir('seed_step_origin')

    suite.section('reading what a seed already carries')
    equilibrated = write_state_xml(work / 'seed.xml', seed_origin)
    suite.check('an equilibrated seed reports its own step count',
                util.state_xml_origin(equilibrated) == seed_origin,
                f'-> {util.state_xml_origin(equilibrated)}')
    fresh = write_state_xml(work / 'fresh.xml', 0)
    suite.check('a seed that starts at zero reports zero',
                util.state_xml_origin(fresh) == 0)
    gro = work / 'conf.gro'
    gro.write_text('not xml\n')
    suite.check('a structure file that is not a state.xml reports zero',
                util.state_xml_origin(gro) == 0)
    suite.check('a seed that is not there reports zero',
                util.state_xml_origin(work / 'absent.xml') == 0)

    suite.section('a clone records the origin once, at construction')
    farm = work / 'farm'
    config = dict(
        traj_dir_top_level=str(farm), seed_index=SEED_INDEX,
        clone_index=CLONE_INDEX, gen_index=0, dirname_pad=DIRNAME_PAD,
        sep=SEP, steps=steps_per_gen, steps_per_gen=steps_per_gen,
        write_interval=write_interval, restart_name='state.xml',
        traj_name='traj', traj_suffix='.dcd', top_fn=str(work / 'top.pdb'),
        seed_fn=str(equilibrated), title='origin')
    clone = seeder.Clone(
        dict(config), scheduler='slurm', scheduler_fstring='echo {job_name}',
        scheduler_kws={}, seed_fn=str(equilibrated), dry_run=True,
        job_name_fstring='{title}-{seed_index}-{clone_index}-{gen_index}')
    suite.check('the clone records the seed origin in its config',
                clone.config.get('step_origin') == seed_origin,
                f"-> {clone.config.get('step_origin')}")

    suite.section('recovery measures this campaign, not the equilibration')
    # Generation 1, with generation 0 recorded as a full-length generation.
    for gen in (0, 1):
        gen_dir = mdfarmer.dir_seeds_clones_gens(
            farm, SEED_INDEX, CLONE_INDEX, gen, DIRNAME_PAD, sep=SEP)
        util.write_json_atomic(gen_dir / util.CONFIG_NAME,
                               dict(config, gen_index=gen,
                                    step_origin=seed_origin), indent=4)
    gen_one = mdfarmer.dir_seeds_clones_gens(
        farm, SEED_INDEX, CLONE_INDEX, 1, DIRNAME_PAD, sep=SEP)
    # Generation 1 has run half its budget, on top of the seed's own history.
    half = steps_per_gen // 2
    write_state_xml(gen_one / 'state.xml', seed_origin + steps_per_gen + half)

    start = (json.loads((gen_one / util.CONFIG_NAME).read_text())
             .get('step_origin', 0)
             + util.steps_before(str(farm), SEED_INDEX, CLONE_INDEX, 1,
                                 DIRNAME_PAD, sep=SEP))
    suite.check('the generation is taken to start where the campaign put it',
                start == seed_origin + steps_per_gen,
                f'-> {start} vs {seed_origin + steps_per_gen}')
    offset = (seed_origin + steps_per_gen + half) - start
    suite.check('so the steps it has run are its own, not the seed history',
                offset == half, f'-> {offset} vs {half}')
    suite.check('and that is a sane number of frames, not hundreds',
                offset // write_interval == half // write_interval,
                f'-> {offset // write_interval} frames')

    ignoring_origin = (seed_origin + steps_per_gen + half) - steps_per_gen
    suite.check('ignoring the origin would have claimed far more',
                ignoring_origin // write_interval > 100,
                f'-> {ignoring_origin // write_interval} frames')
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
