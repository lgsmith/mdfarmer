"""Where a topology's #include files live has to be sayable.

A .top that names a force field resolves its #include lines against the process
cwd or, failing that, against whatever openmm guesses from the environment. A
module-installed GROMACS sets GMXBIN to a binary name rather than a directory,
which makes that guess a path that does not exist, so the campaign has to be
able to say the directory itself.
"""
import json
import sys

import harness
from harness import Suite

from mdfarmer import reimage

INCLUDE_DIR = '/some/force/fields/share/gromacs/top'


def main(include_dir=INCLUDE_DIR):
    suite = Suite('reimage_includes')
    work = harness.workdir('reimage_includes')

    suite.section('a .top whose force field is nowhere to be found')
    top_p = work / 'topol.top'
    top_p.write_text('#include "not-a-real.ff/forcefield.itp"\n')
    try:
        reimage.gromacs_topology(top_p)
        suite.check('the failure is explained, not just raised', False,
                    '-> no exception')
    except ValueError as exc:
        message = str(exc)
        print(f'   {message}', flush=True)
        suite.check('the message names both environment variables',
                    'GMXDATA' in message and 'GMXBIN' in message)
        suite.check('the message says to pass include_dir',
                    'include_dir' in message)

    suite.section('a generation directory that records its include_dir')
    gen_dir = work / 'gen-00'
    gen_dir.mkdir()
    (gen_dir / 'prod.xtc').write_text('a trajectory\n')
    (gen_dir / 'config.json').write_text(json.dumps(dict(
        traj_name='prod', traj_suffix='.xtc', gen_index=0,
        top_fn=str(top_p), structure_fn=str(work / 'conf.gro'),
        include_dir=include_dir)))

    seen = {}

    def spy(traj_fn, **kwargs):
        seen.update(kwargs)
        return traj_fn

    saved = reimage.reimage_trajectory
    reimage.reimage_trajectory = spy
    try:
        reimage.reimage_gen_dir(gen_dir)
    finally:
        reimage.reimage_trajectory = saved
    suite.check('the config supplies include_dir when the caller does not',
                seen.get('include_dir') == include_dir,
                f'-> {seen.get("include_dir")}')

    reimage.reimage_trajectory = spy
    try:
        reimage.reimage_gen_dir(gen_dir, include_dir='/asked/for/this/one')
    finally:
        reimage.reimage_trajectory = saved
    suite.check('an explicit include_dir still wins',
                seen.get('include_dir') == '/asked/for/this/one',
                f'-> {seen.get("include_dir")}')
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
