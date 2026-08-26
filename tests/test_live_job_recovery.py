"""Disk recovery must not touch a generation whose job is still running.

Recovery trims and renames files, and at boot the scheduler query has already
told us which generations are alive. Running recovery over one of those rewrites
a trajectory the live job still holds open.
"""
import json
import sys

import harness
from harness import Suite

from mdfarmer import utilities as util
from mdfarmer import gmx_simulate as gs
from mdfarmer.seeder import Clone

SEED_INDEX = 0
CLONE_INDEX = 0
GEN_INDEX = 0
LIVE_JOB = 987654
STEPS_PER_GEN = 1000
WRITE_INTERVAL = 100


def main(seed_index=SEED_INDEX, clone_index=CLONE_INDEX, gen_index=GEN_INDEX,
         live_job=LIVE_JOB, steps_per_gen=STEPS_PER_GEN,
         write_interval=WRITE_INTERVAL):
    suite = Suite('live_job_recovery')
    work = harness.workdir('live_job_recovery')
    for name in ('a.gro', 'topol.top', 'base.mdp'):
        (work / name).write_text('placeholder\n')

    template = gs.gmx_config_template(
        traj_dir_top_level=str(work / 'farm'), title='live',
        structure_fn=str(work / 'a.gro'), mdp_fn=str(work / 'base.mdp'),
        dirname_pad=2, sep='_', traj_name='prod', traj_suffix='.xtc',
        restart_name='state.cpt', steps=steps_per_gen,
        steps_per_gen=steps_per_gen, write_interval=write_interval,
        traj_list=str(work / 'tl.txt'))

    tdir = work / 'farm'
    gen_dir = util.dir_seeds_clones_gens(tdir, seed_index, clone_index,
                                         gen_index, 2, sep='_', mkdir=True)
    (gen_dir / 'config.json').write_text(json.dumps(dict(template)))
    (gen_dir / 'prod.xtc').write_text('a trajectory the live job is writing\n')
    (gen_dir / 'state.cpt').write_text('a checkpoint the live job is writing\n')

    touched = []

    def spy(gen_path, **context):
        touched.append(gen_path)
        return None

    def build(rep_dict):
        del touched[:]
        return Clone.from_disk(
            tdir, seed_index, clone_index,
            initial_seed_fn=str(work / 'a.gro'), top_fn=str(work / 'topol.top'),
            system_fn=str(work / 'base.mdp'), config_template=dict(template),
            scheduler='sbatch',
            scheduler_fstring=util.basic_scheduler_fstrings['slurm'],
            scheduler_kws=dict(gpu_line='', queue_name='gpu', exclude_nodes='',
                               run_script_name='run.py'),
            dirname_pad=2, sep='_', job_number_re='[1-9][0-9]*',
            job_name_fstring='{title}_{seed_index}_{clone_index}_{gen_index}',
            recover_fn=spy, rep_dict=rep_dict, dry_run=True)

    suite.section('no live job: recovery runs as usual')
    clone = build({})
    suite.check('the generation on disk is examined', touched == [gen_dir],
                f'-> {[p.name for p in touched]}')
    suite.check('no job is bound', clone.job_number is None)

    suite.section('a live job: its files are left alone')
    clone = build({(seed_index, clone_index, gen_index): live_job})
    suite.check('recovery never looks at the live generation', touched == [],
                f'-> {[p.name for p in touched]}')
    suite.check('the clone binds to the running job',
                clone.job_number == live_job, f'-> {clone.job_number}')
    suite.check('it stays on that generation',
                clone.config['gen_index'] == gen_index)
    suite.check('it is treated as a resume, not a fresh start',
                clone.config['new_velocities'] is False
                and clone.config['append'] is True)
    suite.check('it seeds from the checkpoint in that directory',
                clone.config['seed_fn'] == str((gen_dir / 'state.cpt').resolve()),
                f"-> {clone.config['seed_fn']}")

    suite.section('a live job whose generation has no checkpoint yet')
    (gen_dir / 'state.cpt').unlink()
    clone = build({(seed_index, clone_index, gen_index): live_job})
    suite.check('recovery still leaves it alone', touched == [])
    suite.check('it falls back to the seed it was given',
                clone.config['seed_fn'] == str(work / 'a.gro'))
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
