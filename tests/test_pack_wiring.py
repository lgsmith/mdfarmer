"""Node blocking, the restart budget, and velocity seeds, for packed jobs.

Needs no GROMACS: every check here is about orchestration bookkeeping, and the
files only have to exist.
"""
import sys

import harness
from harness import Suite

from mdfarmer import utilities as util
from mdfarmer import gmx_pack as gp
from mdfarmer import gmx_simulate as gs
from mdfarmer.seeder import Clone, ClonePack

CPUS = 8
STEPS_PER_GEN = 1000
WRITE_INTERVAL = 100
N_MEMBERS = 2
FAILING_NODE = 'workstation099'
# One of default_bad_node_patterns, so the registry recognises it.
NODE_FAILURE = 'CUDA_ERROR_NO_DEVICE'


def clone_config(work, clone_index, steps_per_gen=STEPS_PER_GEN,
                 write_interval=WRITE_INTERVAL):
    return dict(
        traj_dir_top_level=str(work / 'farm'), top_fn=str(work / 't.top'),
        seed_index=0, clone_index=clone_index, gen_index=0, title='w',
        structure_fn=str(work / 'c.gro'), mdp_fn=str(work / 'b.mdp'),
        dirname_pad=2, sep='_', traj_name='prod', traj_suffix='.xtc',
        restart_name='state.cpt', steps=steps_per_gen,
        steps_per_gen=steps_per_gen, write_interval=write_interval,
        temperature=300, gen_seed_base=1, new_velocities=True, append=False,
        mdrun_args=[], traj_list=str(work / 'tl.txt'))


def make_clone(work, clone_index, scheduler_kws, restarts_per_gen=3,
               cls=Clone):
    return cls(clone_config(work, clone_index), 'sbatch',
               util.basic_scheduler_fstrings_mps['slurm'], scheduler_kws,
               seed_fn=str(work / 'c.gro'), sep='_', dirname_pad=2,
               steps_per_gen=STEPS_PER_GEN, dry_run=True,
               restarts_per_gen=restarts_per_gen)


class ScriptedClone(Clone):
    """Clone whose reported progress is set by the test. Clone uses __slots__,
    so a method cannot be patched onto an instance."""

    owed = 0

    def gen_remaining_steps(self):
        return type(self).owed

    def was_preempted(self):
        return False


def main(cpus=CPUS, n_members=N_MEMBERS, failing_node=FAILING_NODE,
         node_failure=NODE_FAILURE, steps_per_gen=STEPS_PER_GEN):
    suite = Suite('pack_wiring')
    work = harness.workdir('pack_wiring')
    for name in ('c.gro', 't.top', 'b.mdp'):
        (work / name).write_text('placeholder\n')
    scheduler_kws = dict(gpu_line='', queue_name='gpu', exclude_nodes='',
                         cpus=cpus, run_script_name='run.py')

    suite.section('a pack member\'s node scan reads the pack log')
    clones = [make_clone(work, i, scheduler_kws) for i in range(n_members)]
    pack = ClonePack(clones, work / 'pack-00', 'sbatch',
                     util.basic_scheduler_fstrings_mps['slurm'], scheduler_kws,
                     run_script=gp.default_gmx_pack_run_script,
                     cpus_per_task=cpus, dry_run=True)
    suite.check('every member points at the pack directory',
                all(c.scheduler_log_dir == pack.pack_dir for c in clones))

    registry = util.BadNodeRegistry(work / 'bad_nodes.csv', 'slurm',
                                    scheduler_kws)
    (pack.pack_dir / 'slurm.out').write_text(
        f'JOB_NAME: w_0_0_0\nNODE: {failing_node}\nGPU: RTX\n'
        f'Fatal error: {node_failure}\n')
    hit = registry.scan_and_record(
        clones[0].scheduler_log_dir or clones[0].current_gen_dir,
        clones[0].get_tag())
    suite.check('a node-local failure in the pack log is detected', bool(hit),
                f'-> {hit}')
    suite.check('the failing node reaches exclude_nodes',
                failing_node in scheduler_kws.get('exclude_nodes', ''),
                f'-> {scheduler_kws.get("exclude_nodes", "")!r}')

    suite.section('exclusions learned after construction reach the pack')
    pack.check_start_gen(set(), overwrite=True)
    script = (pack.pack_dir / 'sbatch.sh').read_text()
    suite.check('the submit script carries the live exclusion',
                failing_node in script)
    suite.check('the pack still overlays its own cpus-per-task',
                f'--cpus-per-task={cpus}' in script)
    suite.check('the pack does not mutate the shared scheduler_kws',
                scheduler_kws['cpus'] == cpus)

    suite.section('the restart budget counts consecutive dead launches')
    clone = make_clone(work, 0, dict(scheduler_kws), cls=ScriptedClone)
    clone.remaining_steps = steps_per_gen
    clone.restart_attempts = 2
    ScriptedClone.owed = steps_per_gen // 2          # a launch that progressed
    clone.check_start_gen(set(), overwrite=True)
    suite.check('progress clears the accumulated budget',
                clone.restart_attempts == 0, f'-> {clone.restart_attempts}')
    clone.restart_attempts = 2
    clone.remaining_steps = steps_per_gen // 2
    ScriptedClone.owed = steps_per_gen // 2          # a launch that did nothing
    clone.check_start_gen(set(), overwrite=True)
    suite.check('a dead launch still charges the budget',
                clone.restart_attempts == 3, f'-> {clone.restart_attempts}')

    suite.section('gen-seed differs across seeds as well as clones')
    seeds = {}
    for seed_index in (0, 1):
        for clone_index in (0, 1):
            gen_dir = work / f's{seed_index}c{clone_index}'
            gen_dir.mkdir(parents=True, exist_ok=True)
            base = gen_dir / 'base.mdp'
            base.write_text('integrator = md\ngen-seed = 0\nld-seed = -1\n')
            gs.write_gen_mdp(
                str(base), str(gen_dir / 'gen.mdp'), nsteps=10,
                nstxout_compressed=5, gen_vel=True, continuation=False,
                gen_seed=1 + gs.GEN_SEED_STRIDE * seed_index + clone_index,
                gen_temp=300)
            for line in (gen_dir / 'gen.mdp').read_text().splitlines():
                if line.startswith('gen-seed'):
                    seeds[(seed_index, clone_index)] = int(line.split('=')[1])
    suite.check('all four replicas draw distinct velocity seeds',
                len(set(seeds.values())) == 4, f'-> {sorted(seeds.values())}')
    suite.check('the same clone index under two seeds differs',
                seeds[(0, 0)] != seeds[(1, 0)])
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
