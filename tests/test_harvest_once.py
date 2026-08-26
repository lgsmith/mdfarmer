"""One harvest job per generation, however many times it is checked in on.

check_start_gen reaps a finished generation and then calls start_next, which
raises if the checkpoint it needs is missing. A ClonePack catches that per
member and keeps the member, so the same finished generation comes back round
next tick. Two harvest jobs in one directory read and write the same files at
once, which is the thing the sentinel exists to prevent.
"""
import sys

import harness
from harness import Suite

from mdfarmer import utilities as util
from mdfarmer import harvester as hv
from mdfarmer.seeder import Clone

STEPS_PER_GEN = 1000
WRITE_INTERVAL = 100


class FakeSp:
    """Stands in for the subprocess module, so nothing is ever submitted."""

    def __init__(self):
        self.calls = 0

    def run(self, *args, **kwargs):
        self.calls += 1
        raise AssertionError('a harvest was submitted into a harvested '
                             'directory')


class ReapCounter:
    """Stands in for a Harvester, so no scheduler job is ever submitted."""

    def __init__(self):
        self.reaped = []

    def reap(self, gen_dir, dry_run=False):
        self.reaped.append(gen_dir)


class FinishedClone(Clone):
    """Every check-in reports this generation as already finished."""

    def gen_remaining_steps(self):
        return 0

    def was_preempted(self):
        return False


def base_config(work, steps_per_gen=STEPS_PER_GEN,
                write_interval=WRITE_INTERVAL):
    return dict(
        traj_dir_top_level=str(work / 'farm'), seed_index=0, clone_index=0,
        gen_index=0, title='once', structure_fn=str(work / 'seed.gro'),
        dirname_pad=2, sep='_', traj_name='prod', traj_suffix='.xtc',
        restart_name='state.cpt', steps=steps_per_gen,
        steps_per_gen=steps_per_gen, write_interval=write_interval,
        new_velocities=True, append=False, mdrun_args=[])


def main(steps_per_gen=STEPS_PER_GEN):
    suite = Suite('harvest_once')
    work = harness.workdir('harvest_once')
    (work / 'seed.gro').write_text('a seed\n')

    suite.section('a generation whose start_next cannot proceed')
    counter = ReapCounter()
    clone = FinishedClone(
        base_config(work), 'sbatch', util.basic_scheduler_fstrings['slurm'],
        dict(gpu_line='', queue_name='gpu', exclude_nodes='',
             run_script_name='run.py'),
        seed_fn=str(work / 'seed.gro'), sep='_', dirname_pad=2,
        steps_per_gen=steps_per_gen, dry_run=True, harvester=counter)
    # No checkpoint in the gen dir, so start_next's set_seed raises, exactly as
    # it does for the pack member ClonePack catches and keeps.
    for _ in range(3):
        try:
            clone.check_start_gen(set())
        except (IOError, FileNotFoundError):
            pass
    suite.check('three check-ins reap the generation once',
                len(counter.reaped) == 1, f'-> {len(counter.reaped)} reaps')

    suite.section('an already harvested directory')
    gen_dir = work / 'harvested'
    gen_dir.mkdir()
    (gen_dir / hv.SENTINEL_NAME).write_text('{}')
    harvester = hv.Harvester('#!/bin/bash\necho hi\n', 'sbatch')
    # The scheduler is replaced rather than trusted: a regression here would
    # otherwise put a real job on a real queue every time the suite is run.
    stub, saved = FakeSp(), hv.sp
    hv.sp = stub
    try:
        submitted = harvester.reap(gen_dir)
    finally:
        hv.sp = saved
    suite.check('reap submits nothing and writes no script',
                submitted is None and stub.calls == 0
                and not (gen_dir / harvester.scriptname).is_file())

    suite.section('a directory with no sentinel still gets its script')
    fresh = work / 'fresh'
    fresh.mkdir()
    suite.check('a dry-run reap writes the script it would submit',
                harvester.reap(fresh, dry_run=True) is None
                and (fresh / harvester.scriptname).is_file())
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
