"""How the harvest job's own config file reaches the generation directory.

prep_and_write_inputs writes it just before reap submits the job that reads it,
so the two can race: it goes through the atomic helper like every other file a
runner and the orchestrator share. It is also indented now, which is safe
because harvest_generation parses it with json.loads and nothing reads it as
text -- but that is worth a check, since the only reader is a job on a queue
whose failure surfaces hours later.
"""
import json
import pathlib
import sys

import harness
from harness import Suite

from mdfarmer import harvester as hv

TEMPLATE = '#!/bin/bash\necho {queue_name}\n'
RUN_CONFIG = dict(queue_name='ccb', downsample_frq=5, harvester_subset='all',
                  harvester_structure='system.pdb')


def prep(gen_dir, run_config=RUN_CONFIG, template=TEMPLATE):
    """Lay down one generation's harvest script and config, submitting nothing.

    Returns (config path, renames onto it, direct writes to its final name).
    """
    harvester = hv.Harvester(template, 'sbatch', run_config=run_config)
    config_p = gen_dir / harvester.run_config_name
    renames, direct = [], []
    real_replace, real_write_text = pathlib.Path.replace, pathlib.Path.write_text

    def watching_replace(self, target):
        if pathlib.Path(target) == config_p:
            renames.append(self.read_text() if self.is_file() else '')
        return real_replace(self, target)

    def watching_write_text(self, data, *args, **kwargs):
        if self == config_p:
            direct.append(data)
        return real_write_text(self, data, *args, **kwargs)

    pathlib.Path.replace = watching_replace
    pathlib.Path.write_text = watching_write_text
    try:
        harvester.prep_and_write_inputs(gen_dir)
    finally:
        pathlib.Path.replace = real_replace
        pathlib.Path.write_text = real_write_text
    return config_p, renames, direct


def main(run_config=RUN_CONFIG):
    suite = Suite('hconfig_write')
    work = harness.workdir('hconfig_write')
    gen_dir = work / 'gen-00'
    gen_dir.mkdir(parents=True)
    config_p, renames, direct = prep(gen_dir)

    suite.section('the config lands whole or not at all')
    suite.check('nothing writes the config name directly',
                not direct, f'-> {direct}')
    suite.check('it arrives by exactly one rename',
                len(renames) == 1, f'-> {len(renames)} renames')
    suite.check('the whole config is on disk before the rename',
                json.loads(renames[0] if renames else 'null') == run_config,
                f'-> {(renames[0] if renames else None)!r}')

    suite.section('the shape the harvest job parses')
    suite.check('it round-trips through json.loads, which is its only reader',
                json.loads(config_p.read_text()) == run_config)
    suite.check('and it is indented, so a human can read it too',
                config_p.read_text() == json.dumps(run_config,
                                                   indent=hv.util.JSON_INDENT),
                f'-> {config_p.read_text()!r}')

    suite.section('a harvester with no config writes only its script')
    bare = work / 'bare'
    bare.mkdir()
    script_p = hv.Harvester('#!/bin/bash\necho hi\n',
                            'sbatch').prep_and_write_inputs(bare)
    suite.check('no config file appears beside it',
                [p.name for p in bare.iterdir()] == [script_p.name],
                f'-> {sorted(p.name for p in bare.iterdir())}')
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
