"""Per-replica mdrun arguments, config templates, and a real packed job.

Runs two concurrent mdruns, so it needs a working `gmx` (set GMXBIN if the
binary is not called `gmx`). Everything else about packing is covered by
test_pack_farmer, which needs no GROMACS.
"""
import json
import sys

import harness
from harness import Suite

from mdfarmer import utilities as util
from mdfarmer import gmx_pack as gp
from mdfarmer import gmx_simulate as gs

CPUS = 8
PACK_CPUS = 16
N_REPLICAS = 2
STEPS_PER_GEN = 500
WRITE_INTERVAL = 100
BASE_MDRUN_ARGS = ['-nb', 'gpu', '-pme', 'gpu', '-ntmpi', '4', '-nt', '32',
                   '-nstlist', '200']
CPU_MDRUN_ARGS = ['-nb', 'cpu', '-pme', 'cpu']


def main(cpus=CPUS, pack_cpus=PACK_CPUS, n_replicas=N_REPLICAS,
         steps_per_gen=STEPS_PER_GEN, write_interval=WRITE_INTERVAL,
         base_mdrun_args=BASE_MDRUN_ARGS, cpu_mdrun_args=CPU_MDRUN_ARGS,
         gmx_bin=harness.GMX_BIN):
    suite = Suite('pack_findings')
    work = harness.workdir('pack_findings')
    structure, topology, mdp = harness.build_water_system(work,
                                                          gmx_bin=gmx_bin)

    suite.section('per-replica thread and pinning flags')
    first = gp.replica_mdrun_args(base_mdrun_args, 0, n_replicas, pack_cpus)
    second = gp.replica_mdrun_args(base_mdrun_args, 1, n_replicas, pack_cpus)
    print('   replica 0:', ' '.join(first), flush=True)
    suite.check('an inherited -ntmpi is replaced, not merely stripped',
                first.count('-ntmpi') == 1
                and first[first.index('-ntmpi') + 1] == '1')
    without = gp.replica_mdrun_args(base_mdrun_args, 0, n_replicas, pack_cpus,
                                    ntmpi=None)
    suite.check('ntmpi=None emits none, for a build that rejects the flag',
                '-ntmpi' not in without)
    suite.check('-nt is stripped, since it fixes total threads',
                '-nt' not in first)
    suite.check('replicas get distinct, contiguous core blocks',
                first[first.index('-pinoffset') + 1] == '0'
                and second[second.index('-pinoffset') + 1] == '8')
    suite.check('physics flags are preserved',
                '-nb' in first and '-nstlist' in first)
    print(f'   {gmx_bin} is a thread-MPI build:',
          gp.gmx_supports_ntmpi(gmx_bin), flush=True)

    suite.section('config templates')
    raw = util.merge_args_defaults_dict(
        gs.gmx_generation, traj_dir_top_level=str(work), title='t',
        seed_index=0, clone_index=0, gen_index=0, seed_fn=str(structure),
        top_fn=str(topology))
    suite.check('a raw template does record the runtime-only keys',
                sorted(k for k in gp.RUNTIME_ONLY_KEYS if k in raw)
                == sorted(gp.RUNTIME_ONLY_KEYS))
    suite.check('a raw template omits the **kwargs catch-all',
                '_unused' not in raw)
    template = gs.gmx_config_template(
        traj_dir_top_level=str(work / 'farm'), top_fn=str(topology),
        title='pack', structure_fn=str(structure), mdp_fn=str(mdp),
        dirname_pad=2, sep='_', traj_name='prod', traj_suffix='.xtc',
        restart_name='state.cpt', steps=steps_per_gen,
        steps_per_gen=steps_per_gen, write_interval=write_interval,
        temperature=300, gen_seed_base=42, gmx_bin=gmx_bin, grompp_maxwarn=3,
        mdrun_args=cpu_mdrun_args, traj_list=str(work / 'traj_list.txt'))
    suite.check('gmx_config_template drops the runtime-only keys',
                not any(k in template for k in gp.RUNTIME_ONLY_KEYS))
    suite.check('gmx_config_template is JSON-serializable',
                bool(json.dumps(template)))

    suite.section(f'{n_replicas} real concurrent mdruns through the pack')
    configs = []
    for clone_index in range(n_replicas):
        gen_dir = (work / 'farm' / 'seed_00' / f'clone_{clone_index:02d}'
                   / 'gen_00')
        gen_dir.mkdir(parents=True)
        config = dict(template, seed_index=0, clone_index=clone_index,
                      gen_index=0, seed_fn=str(structure),
                      new_velocities=True, append=False)
        config_p = gen_dir / 'config.json'
        config_p.write_text(json.dumps(config, indent=2))
        configs.append(config_p)
    gp.write_pack_manifest(work, configs, cpus_per_task=cpus,
                           reps_per_card=n_replicas)
    status = gp.gmx_pack_sim_block_json(work / 'pack.json')
    states = [member['status'] for member in status['members']]
    suite.check('every replica of a template-built pack completes',
                states == ['complete'] * n_replicas, f'-> {states}')
    for member in status['members']:
        if member['status'] != 'complete':
            print('      ', member.get('detail', '')[:200], flush=True)

    import mdtraj as md
    for member in status['members']:
        if member['status'] == 'complete':
            frames = md.load(member['traj'], top=str(structure)).n_frames
            suite.check(f'replica {member["replica"]} wrote a full generation',
                        frames == steps_per_gen // write_interval + 1,
                        f'-> {frames} frames')
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
