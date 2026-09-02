"""The pin mode and stride are configurable, and default to today's values.

'-pin on -pinstride 1' is right only where a logical core is a physical core
and the offsets are cpuset-relative. Nothing in GROMACS says that by
specification, so it has to be a parameter; these checks pin both the default
and the override reaching the command line.
"""
import json
import sys

import harness
from harness import Suite

from mdfarmer import gmx_pack as gp

CPUS = 16
N_REPLICAS = 2
INHERITED_ARGS = ['-nb', 'gpu', '-pin', 'off', '-pinstride', '4',
                  '-pinoffset', '9', '-ntomp', '32']

# Flags replica_mdrun_args derives itself, so an inherited copy must be dropped.
EMITTED_FLAGS = ('-ntmpi', '-ntomp', '-pin', '-pinoffset', '-pinstride')


def _flag_value(args, flag):
    """The single value following flag, or None if it is absent."""
    return args[args.index(flag) + 1] if flag in args else None


def _pack_member_args(work, configs, cpus=CPUS, **pin_kwargs):
    """mdrun_args a packed job hands each member, with gmx_generation stubbed.

    The stub raises GenIncomplete, which the pack treats as an ordinary
    unfinished generation, so nothing is launched and nothing is logged as a
    crash; the flags are the only thing under test.
    """
    seen = []

    def fake_generation(**conf):
        seen.append(conf['mdrun_args'])
        raise gp.gmx.GenIncomplete('stubbed out')

    gp.write_pack_manifest(work, configs, cpus_per_task=cpus)
    real = gp.gmx.gmx_generation
    gp.gmx.gmx_generation = fake_generation
    try:
        gp.gmx_pack_sim_block_json(work / 'pack.json', **pin_kwargs)
    finally:
        gp.gmx.gmx_generation = real
    return seen


def _member_configs(work, n_replicas=N_REPLICAS):
    """One config.json per replica, complete enough for the pack to dispatch."""
    paths = []
    for index in range(n_replicas):
        config = dict(
            traj_dir_top_level=str(work / 'farm'), top_fn=str(work / 't.top'),
            seed_index=0, clone_index=index, gen_index=0, title='p',
            structure_fn=str(work / 'c.gro'), mdp_fn=str(work / 'b.mdp'),
            steps=100, steps_per_gen=100, write_interval=10, temperature=300,
            gen_seed_base=1, new_velocities=True, append=False,
            mdrun_args=[], traj_list=str(work / f'tl{index}.txt'))
        path = work / f'config_{index}.json'
        path.write_text(json.dumps(config))
        paths.append(path)
    return paths


def main(cpus=CPUS, n_replicas=N_REPLICAS, inherited=INHERITED_ARGS,
         emitted_flags=EMITTED_FLAGS):
    suite = Suite('pack_pinning')
    work = harness.workdir('pack_pinning')

    suite.section("the defaults are today's behaviour")
    suite.check('PIN_MODE is on', gp.PIN_MODE == 'on', f'-> {gp.PIN_MODE!r}')
    suite.check('PIN_STRIDE is 1', gp.PIN_STRIDE == 1, f'-> {gp.PIN_STRIDE!r}')
    default = gp.replica_mdrun_args([], 0, n_replicas, cpus)
    print('   default:', ' '.join(default), flush=True)
    suite.check('an unconfigured replica still asks for -pin on',
                _flag_value(default, '-pin') == 'on')
    suite.check('an unconfigured replica still asks for -pinstride 1',
                _flag_value(default, '-pinstride') == '1')

    suite.section('a non-default mode and stride reach the command line')
    tuned = gp.replica_mdrun_args([], 1, n_replicas, cpus,
                                  pin_mode='off', pin_stride=2)
    print('   tuned:  ', ' '.join(tuned), flush=True)
    suite.check('pin_mode is what -pin becomes',
                _flag_value(tuned, '-pin') == 'off')
    suite.check('pin_stride is what -pinstride becomes',
                _flag_value(tuned, '-pinstride') == '2')
    suite.check('the core block is still derived per replica',
                _flag_value(tuned, '-pinoffset') == '8'
                and _flag_value(tuned, '-ntomp') == '8')

    suite.section('PER_REPLICA_MDRUN_FLAGS still covers every derived flag')
    for flag in emitted_flags:
        suite.check(f'{flag} is stripped from an inherited template',
                    flag in gp.PER_REPLICA_MDRUN_FLAGS)
    from_template = gp.replica_mdrun_args(inherited, 0, n_replicas, cpus)
    print('   inherited:', ' '.join(from_template), flush=True)
    suite.check('an inherited -pin off does not survive alongside the new one',
                from_template.count('-pin') == 1
                and _flag_value(from_template, '-pin') == 'on')
    suite.check('an inherited -pinstride 4 does not survive either',
                from_template.count('-pinstride') == 1
                and _flag_value(from_template, '-pinstride') == '1')
    suite.check('feeding the result back in is a fixed point',
                gp.replica_mdrun_args(from_template, 0, n_replicas, cpus)
                == from_template)
    suite.check('physics flags are preserved', '-nb' in from_template)

    suite.section('a packed job threads both settings through to its members')
    configs = _member_configs(work, n_replicas=n_replicas)
    tuned_args = _pack_member_args(work, configs, cpus=cpus,
                                   pin_mode='off', pin_stride=2)
    print('   member 0:', ' '.join(tuned_args[0]), flush=True)
    suite.check('every member was dispatched', len(tuned_args) == n_replicas,
                f'-> {len(tuned_args)}')
    suite.check('every member got the overridden mode and stride',
                all(_flag_value(a, '-pin') == 'off'
                    and _flag_value(a, '-pinstride') == '2'
                    for a in tuned_args))
    default_args = _pack_member_args(work, configs, cpus=cpus)
    suite.check('an unconfigured pack still pins on with stride 1',
                all(_flag_value(a, '-pin') == 'on'
                    and _flag_value(a, '-pinstride') == '1'
                    for a in default_args))

    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
