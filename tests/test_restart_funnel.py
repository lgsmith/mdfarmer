"""Every failure to advance is counted against one budget: restart_attempts.

Clone.spend_restart is the only place a clone is written off, so a refused
sbatch, an unreadable job id and a raise inside the tender all cost a restart
rather than the whole campaign. A clone is given up on only when that one
counter runs past restarts_per_gen.

Hand-built Clones with subprocess replaced, so nothing is ever submitted and
neither GROMACS nor a scheduler is needed.
"""
import contextlib
import io
import sys

import harness
from harness import Suite

from mdfarmer import utilities as util
from mdfarmer import seeder as sd
from mdfarmer.seeder import Clone, ClonePack

STEPS_PER_GEN = 1000
WRITE_INTERVAL = 100
RESTARTS = 2
# An sbatch that refuses: bad QOS, bad account, a scheduler hiccup.
REFUSED_STDERR = 'sbatch: error: Batch job submission failed\n'


class FakeResult:
    def __init__(self, stdout='', stderr='', returncode=0):
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
        clone_index=clone_index, gen_index=0, title='funnel',
        structure_fn=str(work / 'seed.gro'), dirname_pad=2, sep='_',
        traj_name='prod', traj_suffix='.xtc', restart_name='state.cpt',
        steps=steps_per_gen, steps_per_gen=steps_per_gen,
        write_interval=write_interval, new_velocities=True, append=False,
        mdrun_args=[])


def make_clone(work, clone_index=0, restarts_per_gen=RESTARTS, cls=Clone,
               steps_per_gen=STEPS_PER_GEN):
    return cls(
        base_config(work, clone_index=clone_index), 'sbatch',
        util.basic_scheduler_fstrings['slurm'],
        dict(gpu_line='', queue_name='gpu', exclude_nodes='',
             run_script_name='run.py'),
        seed_fn=str(work / 'seed.gro'), sep='_', dirname_pad=2,
        steps_per_gen=steps_per_gen, restarts_per_gen=restarts_per_gen,
        dry_run=False)


class HalfDoneClone(Clone):
    """Reports this generation as half run, so a check-in sees progress."""

    def gen_remaining_steps(self):
        return self.total_steps // 2

    def was_preempted(self):
        return False


def submit_with(sp_stub, call):
    """Run call() with seeder's subprocess replaced, discarding what it prints."""
    saved = sd.sp
    sd.sp = sp_stub
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            return call()
    finally:
        sd.sp = saved


def main(restarts=RESTARTS, steps_per_gen=STEPS_PER_GEN):
    suite = Suite('restart_funnel')
    work = harness.workdir('restart_funnel')
    (work / 'seed.gro').write_text('a seed\n')

    suite.section('spend_restart is what decides a clone is finished')
    clone = make_clone(work)
    spent = [clone.spend_restart() for _ in range(restarts)]
    suite.check(f'the first {restarts} calls are granted', all(spent),
                f'-> {spent}')
    suite.check('each granted call cost exactly one restart',
                clone.restart_attempts == restarts,
                f'-> {clone.restart_attempts}')
    with contextlib.redirect_stdout(io.StringIO()):
        refused = [clone.spend_restart() for _ in range(2)]
    suite.check('the call past the budget is refused, and stays refused',
                refused == [False, False], f'-> {refused}')
    suite.check('a refused call does not run the counter up further',
                clone.restart_attempts == restarts,
                f'-> {clone.restart_attempts}')

    suite.section('a refused sbatch costs a restart, not the clone')
    clone = make_clone(work, clone_index=1)
    stub = FakeSp(FakeResult(stderr=REFUSED_STDERR, returncode=1))
    for attempt in range(1, restarts + 1):
        started = submit_with(stub, clone.start_current)
        suite.check(f'refusal {attempt} of {restarts} keeps the clone',
                    started is True and clone.restart_attempts == attempt,
                    f'-> {started}, {clone.restart_attempts}')
    started = submit_with(stub, clone.start_current)
    suite.check('the attempt past the budget gives the clone up',
                started is False, f'-> {started}')
    suite.check('and it was never submitted again', stub.calls == restarts,
                f'-> {stub.calls} calls')

    suite.section('a preempted restart that cannot be submitted still counts')
    clone = make_clone(work, clone_index=2)
    stub = FakeSp(FakeResult(stderr=REFUSED_STDERR, returncode=1))
    started = submit_with(
        stub, lambda: clone.start_current(count_as_restart=False))
    suite.check('a waived restart is charged once the submission fails',
                started is True and clone.restart_attempts == 1,
                f'-> {started}, {clone.restart_attempts}')

    suite.section('a launch that gets somewhere clears the budget')
    clone = make_clone(work, clone_index=3, cls=HalfDoneClone)
    clone.dry_run = True
    clone.restart_attempts = restarts
    with contextlib.redirect_stdout(io.StringIO()):
        clone.check_start_gen(set(), overwrite=True)
    suite.check('progress resets the counter the failures were building',
                clone.restart_attempts == 0, f'-> {clone.restart_attempts}')

    suite.section('a pack spends its members budgets, then retires them')
    members = [make_clone(work, clone_index=4 + i) for i in range(2)]
    pack = ClonePack(members, work / 'pack', 'sbatch',
                     util.basic_scheduler_fstrings_mps['slurm'],
                     dict(gpu_line='', queue_name='gpu', exclude_nodes=''),
                     run_script='', cpus_per_task=4, sep='_', dry_run=True)
    with contextlib.redirect_stdout(io.StringIO()):
        outcomes = [pack.spend_restart() for _ in range(restarts + 1)]
    suite.check('the pack keeps going while any member has budget left',
                outcomes[:restarts] == [True] * restarts, f'-> {outcomes}')
    suite.check('every member was charged the same number of restarts',
                all(c.restart_attempts == restarts for c in members),
                f'-> {[c.restart_attempts for c in members]}')
    suite.check('the pack reports failure once the last member is retired',
                outcomes[-1] is False, f'-> {outcomes}')
    suite.check('both members are recorded as retired',
                pack.retired == {0, 1}, f'-> {pack.retired}')
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
