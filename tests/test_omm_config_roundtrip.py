"""The config a Clone writes must be runnable by the engine it names.

Every suite that drives a runner hands it a hand-built argument list, which is
not how a job starts. A real job reads the config.json Clone wrote and splats
it, so a key Clone records but the runner does not accept is a TypeError before
any MD happens -- invisible to a test that never goes through the file.

Clone records the whole config, including entries that address the chain rather
than one run, so each runner has to tolerate the other's.
"""
import inspect
import json
import sys

import harness
from harness import Suite

from mdfarmer import gmx_simulate as gs, simulate, utilities as util

# Keys Clone writes into every config.json regardless of engine.
CLONE_WRITES = ('steps_per_gen', 'target_step', 'structure_fn')
RUNNERS = (('omm_generation', simulate.omm_generation),
           ('gmx_generation', gs.gmx_generation))


def takes_everything(runner, keys):
    """(accepted, refused) of these keys, counting **kwargs as accepting all."""
    params = inspect.signature(runner).parameters
    if any(p.kind is p.VAR_KEYWORD for p in params.values()):
        return list(keys), []
    return ([k for k in keys if k in params],
            [k for k in keys if k not in params])


def main():
    suite = Suite('omm_config_roundtrip')
    work = harness.workdir('omm_config_roundtrip')

    suite.section('every runner tolerates what a Clone records')
    for name, runner in RUNNERS:
        _, refused = takes_everything(runner, CLONE_WRITES)
        suite.check(f'{name} accepts the chain-level keys',
                    not refused, f'-> would refuse {refused}')

    suite.section('and tolerates the other engine\'s config keys')
    for name, runner in RUNNERS:
        other = next(r for n, r in RUNNERS if n != name)
        other_keys = [k for k, p in inspect.signature(other).parameters.items()
                      if p.kind is not p.VAR_KEYWORD]
        _, refused = takes_everything(runner, other_keys)
        suite.check(f'{name} survives a config written for the other engine',
                    not refused, f'-> would refuse {sorted(refused)[:4]}')

    suite.section('a config.json on disk splats without a TypeError')
    # Every argument either runner requires, plus the chain-level keys Clone
    # adds on top -- which is what a real config.json holds.
    config = dict(
        traj_dir_top_level=str(work / 'farm'), system_fn=str(work / 'sys.xml'),
        top_fn=str(work / 'top.pdb'), seed_index=0, clone_index=0, gen_index=0,
        title='rt', seed_fn=str(work / 'seed.xml'), steps=10,
        integrator_xml=str(work / 'integrator.xml'),
        steps_per_gen=10, target_step=10, structure_fn=str(work / 'a.gro'),
        write_interval=2, traj_list=str(work / 'tl.txt'))
    config_p = work / util.CONFIG_NAME
    util.write_json_atomic(config_p, config, indent=4)
    loaded = json.loads(config_p.read_text())
    suite.check('the config round-trips through the file',
                loaded == config)

    # Bind the arguments without running: a missing or refused key raises here,
    # which is exactly where a real job would die.
    for name, runner in RUNNERS:
        try:
            inspect.signature(runner).bind(**loaded)
            bound, detail = True, ''
        except TypeError as exc:
            bound, detail = False, f'-> {str(exc)[:60]}'
        suite.check(f'{name} binds every key in that config', bound, detail)
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
