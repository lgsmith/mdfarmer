"""A failed scheduler query must not read as an empty queue.

Every launch decision rests on which jobs the scheduler says are alive. A
squeue that dies and prints nothing looks exactly like a quiet queue, and read
that way it relaunches every live clone on top of itself. The shipped reports
are single commands whose own exit status settles that; pipefail is still set
for the pipelines a site may substitute, where the last stage would otherwise
exit 0 over a dead first stage.
"""
import sys

import harness
from harness import Suite

from mdfarmer import utilities as util
from mdfarmer import farmer as fm
from mdfarmer import gmx_simulate as gs

# What LSF prints when the queue holds none of our jobs. It exits non-zero.
LSF_EMPTY = "echo 'No unfinished job found' >&2; exit 255"
# A pipeline whose first stage fails, which is the shape of the real commands.
FAILED_PIPELINE = "sh -c 'exit 1' | awk '/anything/'"


def main(lsf_empty=LSF_EMPTY, failed_pipeline=FAILED_PIPELINE):
    suite = Suite('scheduler_query')
    work = harness.workdir('scheduler_query')

    suite.section('a query that worked')
    suite.check('an empty queue is trusted and empty',
                util.scheduler_query('true') == (True, ''))
    suite.check('output comes back stripped',
                util.scheduler_query('echo " 123 job "') == (True, '123 job'))

    suite.section('a query that failed')
    trusted, _ = util.scheduler_query(failed_pipeline)
    suite.check('a broken pipeline is not trusted', trusted is False)
    suite.check('a plain non-zero exit is not trusted',
                util.scheduler_query('false')[0] is False)
    suite.check('a query that cannot run at all is not trusted',
                util.scheduler_query('exit 127')[0] is False)
    suite.check('a query that never comes back is not trusted',
                util.scheduler_query('sleep 30', timeout=1)[0] is False)

    suite.section('a scheduler that reports an empty queue by exiting non-zero')
    suite.check('LSF\'s empty queue is trusted, not treated as a failure',
                util.scheduler_query(lsf_empty) == (True, ''))

    suite.section('the Farmer refuses to boot on an untrusted query')
    for name in ('a.gro', 'topol.top', 'base.mdp'):
        (work / name).write_text('placeholder\n')
    template = gs.gmx_config_template(
        traj_dir_top_level=str(work / 'farm'), title='sq',
        structure_fn=str(work / 'a.gro'), mdp_fn=str(work / 'base.mdp'),
        dirname_pad=2, sep='_', steps=1000, steps_per_gen=1000,
        write_interval=100, traj_list=str(work / 'tl.txt'))

    def boot(assoc_cmd):
        return fm.Farmer(
            n_seeds=1, n_clones=1, n_gens=1, config_template=dict(template),
            seed_structure_fns=[str(work / 'a.gro')],
            system_fns=[str(work / 'base.mdp')],
            top_fns=[str(work / 'topol.top')],
            scheduler='sbatch',
            scheduler_fstring=util.basic_scheduler_fstrings['slurm'],
            scheduler_kws=dict(gpu_line='', queue_name='gpu',
                               exclude_nodes='', run_script_name='run.py'),
            scheduler_report_cmd='true', scheduler_assoc_rep_cmd=assoc_cmd,
            sep='_', dirname_pad=2, runner=gs.gmx_generation, dry_run=True,
            overwrite=True, jids_file=work / 'jids.txt')

    try:
        boot(failed_pipeline)
        suite.check('a failed boot query stops the Farmer', False,
                    '-> booted anyway')
    except RuntimeError as exc:
        suite.check('a failed boot query stops the Farmer', True,
                    f'-> {str(exc)[:55]}')
    try:
        boot('true')
        suite.check('a genuinely empty queue still boots', True)
    except RuntimeError as exc:
        suite.check('a genuinely empty queue still boots', False, f'-> {exc}')
    try:
        boot(lsf_empty)
        suite.check('an empty LSF queue still boots', True)
    except RuntimeError as exc:
        suite.check('an empty LSF queue still boots', False, f'-> {exc}')
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
