"""A generation's config.json, as plow_harrow_plant leaves it on disk.

It is the one file the orchestrator and a running job both read, and the step
arithmetic reads earlier generations' copies, so it goes through the shared
atomic writer and keeps the indent it has always had. Reshaping a file already
on disk would make a stored config differ from itself on the next boot.
"""
import json
import sys

import harness
from harness import Suite

from mdfarmer import utilities as util
from mdfarmer.seeder import Clone

STEPS_PER_GEN = 1000
WRITE_INTERVAL = 100
# What config.json has always been written at, and what the harvester expects.
CONFIG_INDENT = 4


def base_config(work, steps_per_gen=STEPS_PER_GEN,
                write_interval=WRITE_INTERVAL):
    return dict(
        traj_dir_top_level=str(work / 'farm'), seed_index=0, clone_index=0,
        gen_index=0, title='cfg', structure_fn=str(work / 'seed.gro'),
        dirname_pad=2, sep='_', traj_name='prod', traj_suffix='.xtc',
        restart_name='state.cpt', steps=steps_per_gen,
        steps_per_gen=steps_per_gen, write_interval=write_interval,
        new_velocities=True, append=False, mdrun_args=[])


def main(config_indent=CONFIG_INDENT, steps_per_gen=STEPS_PER_GEN):
    suite = Suite('gen_config_write')
    work = harness.workdir('gen_config_write')
    (work / 'seed.gro').write_text('a seed\n')

    clone = Clone(
        base_config(work), 'sbatch', util.basic_scheduler_fstrings['slurm'],
        dict(gpu_line='', queue_name='gpu', exclude_nodes='',
             run_script_name='run.py'),
        seed_fn=str(work / 'seed.gro'), sep='_', dirname_pad=2,
        steps_per_gen=steps_per_gen, dry_run=True)
    clone.plow_harrow_plant()
    config_p = clone.current_gen_dir / 'config.json'

    suite.section('what plow_harrow_plant leaves beside the job')
    suite.check('config.json is written', config_p.is_file(), f'-> {config_p}')
    suite.check('it round-trips as the clone config',
                json.loads(config_p.read_text()) == clone.config)
    suite.check(f'it keeps its indent of {config_indent}',
                config_p.read_text() == json.dumps(clone.config,
                                                   indent=config_indent),
                f'-> {config_p.read_text()[:40]!r}')
    suite.check('the atomic write leaves no temp file behind',
                sorted(p.name for p in clone.current_gen_dir.iterdir()
                       if p.name.endswith('.tmp')) == [],
                f'-> {sorted(p.name for p in clone.current_gen_dir.iterdir())}')

    suite.section('a stale config is replaced even without overwrite')
    config_p.write_text('{"steps": 1}')
    clone.plow_harrow_plant(overwrite=False)
    suite.check('the rewrite happens anyway, so no half-finished gen restarts '
                'from scratch',
                json.loads(config_p.read_text()) == clone.config)
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
