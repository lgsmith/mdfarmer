"""A submission whose job number cannot be read must not read as a failure.

The job is already running when the id is parsed, so mistaking an unparseable
id for a failed submission leaves that job writing into a generation directory
the tender has stopped minding. What the operator needs to be told is that the
submission SUCCEEDED, so the orphan can be found and cancelled.
"""
import sys

import harness
from harness import Suite

from mdfarmer import utilities as util
from mdfarmer import seeder as sd
from mdfarmer.seeder import Clone, ClonePack

STEPS_PER_GEN = 1000
WRITE_INTERVAL = 100
# What a wrapper that swallows sbatch's own output looks like.
NO_ID_STDOUT = 'submitted, have a nice day\n'
NO_ID_STDERR = ''


class FakeResult:
    def __init__(self, stdout, stderr='', returncode=0):
        self.stdout, self.stderr, self.returncode = stdout, stderr, returncode


class FakeSp:
    """Stands in for the subprocess module, so nothing is ever submitted."""

    def __init__(self, result):
        self.result = result
        self.calls = 0

    def run(self, *args, **kwargs):
        self.calls += 1
        return self.result


def base_config(work, clone_index=0, steps_per_gen=STEPS_PER_GEN,
                write_interval=WRITE_INTERVAL):
    return dict(
        traj_dir_top_level=str(work / 'farm'), seed_index=0,
        clone_index=clone_index, gen_index=0, title='submit',
        structure_fn=str(work / 'seed.gro'), dirname_pad=2, sep='_',
        traj_name='prod', traj_suffix='.xtc', restart_name='state.cpt',
        steps=steps_per_gen, steps_per_gen=steps_per_gen,
        write_interval=write_interval, new_velocities=True, append=False,
        mdrun_args=[])


def make_clone(work, clone_index=0, steps_per_gen=STEPS_PER_GEN):
    return Clone(
        base_config(work, clone_index=clone_index), 'sbatch',
        util.basic_scheduler_fstrings['slurm'],
        dict(gpu_line='', queue_name='gpu', exclude_nodes='',
             run_script_name='run.py'),
        seed_fn=str(work / 'seed.gro'), sep='_', dirname_pad=2,
        steps_per_gen=steps_per_gen, dry_run=False)


def submit_with(sp_stub, call):
    """Run call() with seeder's subprocess replaced, capturing what it printed."""
    import io
    import contextlib
    saved = sd.sp
    sd.sp = sp_stub
    printed = io.StringIO()
    try:
        with contextlib.redirect_stdout(printed):
            result = call()
    finally:
        sd.sp = saved
    print(printed.getvalue(), end='', flush=True)
    return result, printed.getvalue()


def main(no_id_stdout=NO_ID_STDOUT, no_id_stderr=NO_ID_STDERR):
    suite = Suite('submission_guards')
    work = harness.workdir('submission_guards')
    (work / 'seed.gro').write_text('a seed\n')

    suite.section('a solo clone whose sbatch prints no job number')
    clone = make_clone(work)
    stub = FakeSp(FakeResult(no_id_stdout, no_id_stderr))
    started, output = submit_with(stub, clone.start_current)
    suite.check('the submission is not retried', stub.calls == 1,
                f'-> {stub.calls} calls')
    suite.check('it does not raise, and reports no launch', started is False,
                f'-> {started}')
    suite.check('the message says the job was submitted',
                'SUBMITTED' in output)
    suite.check('the message says the job is untracked',
                'UNTRACKED' in output)
    suite.check('the message names the directory to look in',
                str(clone.current_gen_dir) in output)
    suite.check('both streams are printed',
                'stdout:' in output and 'stderr:' in output)
    suite.check('no job number is invented', clone.job_number is None,
                f'-> {clone.job_number}')

    suite.section('a pack whose sbatch prints no job number')
    pack = ClonePack([make_clone(work, clone_index=1)], work / 'pack',
                     'sbatch', util.basic_scheduler_fstrings['slurm'],
                     dict(gpu_line='', queue_name='gpu', exclude_nodes=''),
                     run_script='print("hi")', cpus_per_task=4, sep='_')
    stub = FakeSp(FakeResult(no_id_stdout, no_id_stderr))
    started, output = submit_with(stub, lambda: pack.check_start_gen(set()))
    suite.check('the pack reports no launch', started is False, f'-> {started}')
    suite.check('the pack message says the job was submitted',
                'SUBMITTED' in output and 'UNTRACKED' in output)
    suite.check('the pack names its own directory',
                str(pack.pack_dir) in output)
    suite.check('no job number is invented for the pack',
                pack.job_number is None, f'-> {pack.job_number}')
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
