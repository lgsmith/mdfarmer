"""Everything the package promises to export actually exists.

Deleting a dead function is easy; forgetting it was named in __all__ is easier,
and `from mdfarmer import *` then fails with a message about the wrong thing.
"""
import sys

import harness
from harness import Suite

import mdfarmer


def main():
    suite = Suite('package_surface')

    suite.section('__all__ is honest')
    missing = [name for name in mdfarmer.__all__ if not hasattr(mdfarmer, name)]
    suite.check('every exported name exists', not missing, f'-> {missing}')
    duplicated = sorted({n for n in mdfarmer.__all__
                         if mdfarmer.__all__.count(n) > 1})
    suite.check('no name is listed twice', not duplicated, f'-> {duplicated}')
    suite.check('a star import works',
                exec('from mdfarmer import *', {}) is None)

    suite.section('the pieces a driver is told to use')
    for name in ('Farmer', 'Clone', 'ClonePack', 'Harvester',
                 'gmx_generation', 'gmx_config_template', 'harvest_generation',
                 'reimage', 'gmx_pack', 'frame_timing', 'get_traj_len'):
        suite.check(f'mdfarmer.{name}', hasattr(mdfarmer, name))
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
