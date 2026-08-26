"""Which generation a clone is on has one home: config['gen_index'].

Adaptive-sampling drivers redirect a clone by moving that field, which the
comment on plow_harrow_plant tells them to do. Anything that tracked the
generation separately would go stale the moment they did, and a stale value is
what decides whether a clone retires.
"""
import sys

import harness
from harness import Suite

from mdfarmer import utilities as util
from mdfarmer.seeder import Clone

STEPS_PER_GEN = 1000
WRITE_INTERVAL = 100
LAST_GEN_INDEX = 2
MOVED_GEN = 5


def base_config(work, steps_per_gen=STEPS_PER_GEN,
                write_interval=WRITE_INTERVAL):
    return dict(
        traj_dir_top_level=str(work / 'farm'), seed_index=0, clone_index=0,
        gen_index=0, title='track', structure_fn=str(work / 'seed.gro'),
        dirname_pad=2, sep='_', traj_name='prod', traj_suffix='.xtc',
        restart_name='state.cpt', steps=steps_per_gen,
        steps_per_gen=steps_per_gen, write_interval=write_interval,
        new_velocities=True, append=False, mdrun_args=[])


def main(steps_per_gen=STEPS_PER_GEN, last_gen_index=LAST_GEN_INDEX,
         moved_gen=MOVED_GEN):
    suite = Suite('gen_tracking')
    work = harness.workdir('gen_tracking')
    (work / 'seed.gro').write_text('a seed\n')

    clone = Clone(
        base_config(work), 'sbatch', util.basic_scheduler_fstrings['slurm'],
        dict(gpu_line='', queue_name='gpu', exclude_nodes='',
             run_script_name='run.py'),
        seed_fn=str(work / 'seed.gro'), sep='_', dirname_pad=2,
        steps_per_gen=steps_per_gen, dry_run=True,
        last_gen_index=last_gen_index)

    suite.section('a fresh clone')
    suite.check('current_gen is the generation in the config',
                clone.current_gen == clone.config['gen_index'],
                f'-> {clone.current_gen}')
    suite.check('it is not done yet', clone.is_done is False)

    suite.section('start_next moves both together')
    (clone.current_gen_dir / clone.config['restart_name']).write_text('ckpt\n')
    clone.start_next()
    suite.check('current_gen still matches the config',
                clone.current_gen == clone.config['gen_index'] == 1,
                f'-> {clone.current_gen} vs {clone.config["gen_index"]}')

    suite.section('a caller that moves the config, as the comment invites')
    clone.config['gen_index'] = moved_gen
    suite.check('current_gen follows it', clone.current_gen == moved_gen,
                f'-> {clone.current_gen}')
    suite.check('and so does is_done, past the last generation asked for',
                clone.is_done is True)

    suite.section('there is nothing else to assign')
    try:
        clone.current_gen = 0
        suite.check('assigning current_gen is refused', False,
                    '-> no exception')
    except AttributeError:
        suite.check('assigning current_gen is refused', True)
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
