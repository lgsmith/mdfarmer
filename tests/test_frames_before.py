"""Where a generation sits in the whole trajectory is counted, not multiplied.

Each generation's config.json records what that generation actually ran, and
together with its job script is all a rerun of it needs. The harvester places a
generation's frames by adding up what the earlier ones recorded, so a seed whose
generation length was changed between boots -- or a wallclock-matched scheme
whose generations end wherever the clock ran out -- still lands every frame in
the right place. Multiplying this generation's length by its index would assume
every earlier generation matched it, and silently misplace them all.
"""
import json
import sys

import harness
from harness import Suite

import mdfarmer
from mdfarmer import gmx_simulate as gs, harvester as hv

WRITE_INTERVAL = 100
DOWNSAMPLE_FRQ = 2
SEED_INDEX = 0
CLONE_INDEX = 0
DIRNAME_PAD = 2
SEP = '-'
# Generation lengths, in order: the third is twice the others.
GEN_STEPS = (1000, 1000, 2000, 1000)


def write_gen_config(farm, gen_index, steps_per_gen, write_interval, sep=SEP,
                     dirname_pad=DIRNAME_PAD):
    """A generation directory holding just the run record the harvest reads."""
    gen_dir = mdfarmer.dir_seeds_clones_gens(
        farm, SEED_INDEX, CLONE_INDEX, gen_index, dirname_pad, sep=sep)
    config = dict(
        traj_dir_top_level=str(farm), seed_index=SEED_INDEX,
        clone_index=CLONE_INDEX, gen_index=gen_index, dirname_pad=dirname_pad,
        sep=sep, steps_per_gen=steps_per_gen, write_interval=write_interval)
    (gen_dir / hv.CONFIG_NAME).write_text(json.dumps(config, indent=4))
    return gen_dir, config


def main(gen_steps=GEN_STEPS, write_interval=WRITE_INTERVAL,
         downsample_frq=DOWNSAMPLE_FRQ):
    suite = Suite('frames_before')
    work = harness.workdir('frames_before')
    farm = work / 'farm'
    configs = [write_gen_config(farm, i, steps, write_interval)[1]
               for i, steps in enumerate(gen_steps)]
    per_gen = [s // write_interval for s in gen_steps]

    suite.section('a seed whose generation length changed part way through')
    running = 0
    for gen_index, config in enumerate(configs):
        got = hv.frames_before(config, downsample_frq)
        suite.check(f'generation {gen_index} starts at frame {running}',
                    got == running, f'-> {got}')
        running += per_gen[gen_index]

    suite.section('and the multiplying shortcut would have been wrong')
    longest = configs[-1]
    naive = longest['gen_index'] * (longest['steps_per_gen'] // write_interval)
    counted = hv.frames_before(longest, downsample_frq)
    suite.check('the two disagree once a generation differs',
                naive != counted, f'-> multiplied {naive}, counted {counted}')

    suite.section('and the harvest itself uses the counted offset')
    real = harness.workdir('frames_before_harvest')
    structure, topology, mdp = harness.build_water_system(
        real, gmx_bin=harness.GMX_BIN)
    real_farm = real / 'farm'
    gen_zero, _ = write_gen_config(real_farm, 0, gen_steps[2], write_interval)
    run = dict(
        traj_dir_top_level=str(real_farm), top_fn=str(topology),
        seed_index=SEED_INDEX, clone_index=CLONE_INDEX, title='fb',
        structure_fn=str(structure), mdp_fn=str(mdp), dirname_pad=DIRNAME_PAD,
        sep=SEP, traj_name='prod', traj_suffix='.xtc',
        restart_name='state.cpt', write_interval=write_interval,
        temperature=300, gen_seed_base=1, gmx_bin=harness.GMX_BIN,
        grompp_maxwarn=5, append=False,
        mdrun_args=('-nb', 'cpu', '-pme', 'cpu', '-ntomp', '2'))

    # Generation 0 runs the longer budget its own config records.
    gs.gmx_generation(gen_index=0, seed_fn=str(structure), steps=gen_steps[2],
                      steps_per_gen=gen_steps[2], target_step=gen_steps[2],
                      new_velocities=True, **run)
    gen_one = mdfarmer.dir_seeds_clones_gens(
        real_farm, SEED_INDEX, CLONE_INDEX, 1, DIRNAME_PAD, sep=SEP)
    seed = gen_one / run['restart_name']
    seed.write_bytes((gen_zero / run['restart_name']).read_bytes())
    # This generation adds its own length to where generation 0 ended, which is
    # further along than twice its own would put it.
    config = dict(run, gen_index=1, seed_fn=str(seed), steps=gen_steps[0],
                  steps_per_gen=gen_steps[0], new_velocities=False,
                  target_step=gen_steps[2] + gen_steps[0])
    (gen_one / hv.CONFIG_NAME).write_text(json.dumps(config, indent=4))
    gs.gmx_generation(**config)
    (gen_one / 'hconfig.json').write_text(json.dumps(dict(
        downsample_frq=downsample_frq, harvester_structure=str(structure),
        harvester_subset='all')))
    record = hv.harvest_generation(gen_one / hv.CONFIG_NAME,
                                   gen_one / 'hconfig.json')
    counted = gen_steps[2] // write_interval
    multiplied = gen_steps[0] // write_interval
    suite.check('the sentinel records where generation 0 really ended',
                record['first_global_index'] == counted,
                f"-> {record['first_global_index']}, counted {counted}, "
                f'multiplied would be {multiplied}')

    suite.section('a missing record is refused, not guessed past')
    gen_dir = mdfarmer.dir_seeds_clones_gens(
        farm, SEED_INDEX, CLONE_INDEX, 1, DIRNAME_PAD, sep=SEP)
    (gen_dir / hv.CONFIG_NAME).unlink()
    try:
        hv.frames_before(configs[-1], downsample_frq)
        suite.check('a generation with no config.json raises', False,
                    '-> no exception')
    except FileNotFoundError as exc:
        suite.check('a generation with no config.json raises', True,
                    f'-> {str(exc)[:60]}')
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
