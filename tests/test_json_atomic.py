"""write_json_atomic: what lands on disk, and what a racing reader can see.

Every config, manifest and status file the orchestrator and a running job share
goes through this one helper, so the properties that make it safe -- a temp name
in the target's own directory, an atomic rename, nothing left behind -- are
checked here rather than once per call site.
"""
import json
import pathlib
import sys

import harness
from harness import Suite

from mdfarmer import utilities as util


def check_written_content(suite, work):
    suite.section('what the file holds afterwards')
    path = work / 'status.json'
    obj = {'target_step': 5000, 'complete': False, 'members': [1, 2]}
    returned = util.write_json_atomic(path, obj)
    suite.check('the path is returned, so a caller can chain on it',
                returned == path, f'-> {returned}')
    suite.check('the object round-trips', json.loads(path.read_text()) == obj)

    path.write_text('{"stale": true}')
    util.write_json_atomic(path, obj)
    suite.check('an existing file is replaced, not appended to',
                json.loads(path.read_text()) == obj)

    util.write_json_atomic(str(path), obj)
    suite.check('a str path works as well as a Path',
                json.loads(path.read_text()) == obj)


def check_indent(suite, work):
    """config.json is written at 4 and everything else at 2, so the indent has
    to travel with the call rather than being fixed in the helper."""
    suite.section('the indent a caller asks for')
    two = work / 'two.json'
    four = work / 'four.json'
    obj = {'a': 1}
    util.write_json_atomic(two, obj)
    util.write_json_atomic(four, obj, indent=4)
    suite.check('the default is 2, matching manifests and statuses',
                two.read_text() == json.dumps(obj, indent=2),
                f'-> {two.read_text()!r}')
    suite.check('indent=4 is honoured, so config.json keeps its shape',
                four.read_text() == json.dumps(obj, indent=4),
                f'-> {four.read_text()!r}')


def check_atomicity(suite, work):
    """A rename is only atomic within one filesystem, so the temp file has to be
    a sibling of the target rather than in a temp directory."""
    suite.section('what makes the replace atomic')
    sub = work / 'gen-000'
    sub.mkdir()
    util.write_json_atomic(sub / 'config.json', {'steps': 10})
    leftovers = sorted(p.name for p in sub.iterdir())
    suite.check('no temp file survives the write',
                leftovers == ['config.json'], f'-> {leftovers}')

    seen = []
    real_replace = pathlib.Path.replace

    def watching_replace(self, target):
        seen.append({'tmp_dir': self.parent,
                     'tmp_written': self.is_file(),
                     'tmp_text': self.read_text() if self.is_file() else '',
                     'target_exists': pathlib.Path(target).exists()})
        return real_replace(self, target)

    pathlib.Path.replace = watching_replace
    try:
        util.write_json_atomic(sub / 'manifest.json', {'members': [3]})
    finally:
        pathlib.Path.replace = real_replace

    state = seen[0] if seen else {}
    suite.check('the write goes through a rename at all', len(seen) == 1,
                f'-> {len(seen)} renames')
    suite.check('the temp file is a sibling of the target, so one filesystem',
                state.get('tmp_dir') == sub, f'-> {state.get("tmp_dir")}')
    suite.check('the whole object is on disk before the rename',
                json.loads(state.get('tmp_text') or 'null') == {'members': [3]},
                f'-> {state.get("tmp_text")!r}')
    suite.check('the target does not exist until the rename',
                state.get('target_exists') is False,
                f'-> {state.get("target_exists")}')


def check_seed_map_uses_it(suite, work):
    suite.section('the seed map goes through the same helper')
    seed_map_p = work / 'nested' / util.SEED_MAP_NAME
    util.write_seed_map(seed_map_p, {0: 'apo', 1: 'holo'})
    suite.check('a missing parent directory is still created',
                seed_map_p.is_file(), f'-> {seed_map_p}')
    suite.check('it reads back as it was written',
                util.read_seed_map(seed_map_p) == {0: 'apo', 1: 'holo'})
    suite.check('and is indented at 2 like the other status files',
                seed_map_p.read_text().splitlines()[1].startswith('  "0"'),
                f'-> {seed_map_p.read_text()!r}')


def main():
    suite = Suite('json_atomic')
    work = harness.workdir('json_atomic')
    check_written_content(suite, work)
    check_indent(suite, work)
    check_atomicity(suite, work)
    check_seed_map_uses_it(suite, work)
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
