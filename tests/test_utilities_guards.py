"""The guards in utilities: what each refuses, and how it says so.

Every check here is about an error a scientist reads at 3am, so what matters is
that the type matches what the caller catches and that the message names the
thing that is actually wrong.
"""
import contextlib
import io
import sys

import harness
from harness import Suite

from mdfarmer import harvester
from mdfarmer import utilities as util


def raises(fn, *args, **kwargs):
    """The exception fn raised, or None."""
    try:
        fn(*args, **kwargs)
    except Exception as exc:
        return exc
    return None


def check_topology_reader(suite, readers=None):
    readers = util.openmm_topology_readers if readers is None else readers

    suite.section('a topology format with no reader')
    exc = raises(util.read_openmm_top, 'system.xyz')
    suite.check('refusing to read it is a ValueError',
                isinstance(exc, ValueError), f'-> {type(exc).__name__}')
    suite.check('the message names the format and the alternatives',
                exc is not None and '.xyz' in str(exc) and '.pdb' in str(exc),
                f'-> {exc}')

    suite.section("a broken topology in a format we do read")

    def reader_that_blames_the_file(fn):
        raise KeyError('HWX')

    readers['.fake'] = reader_that_blames_the_file
    said = io.StringIO()
    try:
        with contextlib.redirect_stdout(said):
            exc = raises(util.read_openmm_top, 'system.fake')
    finally:
        del readers['.fake']
    suite.check("the reader's own KeyError reaches the caller",
                isinstance(exc, KeyError) and 'HWX' in str(exc),
                f'-> {type(exc).__name__}: {exc}')
    suite.check('a broken topology is not blamed on a missing reader',
                'reader' not in said.getvalue(),
                f'-> {said.getvalue().strip()!r}')


def check_state_xml_step_count(suite, work):
    suite.section('the step a state.xml stopped at')
    good = work / 'good-state.xml'
    good.write_text('<State stepCount="1200" time="4.8"></State>')
    suite.check('a well-formed state gives its step count',
                util.state_xml_step_count(good) == 1200)

    # A state.xml written by a job the scheduler killed mid-write.
    truncated = work / 'truncated-state.xml'
    truncated.write_text('<State stepCount="1200" tim')
    exc = raises(util.state_xml_step_count, truncated)
    suite.check('an unparseable state.xml raises the ValueError the caller '
                'cascades on', isinstance(exc, ValueError),
                f'-> {type(exc).__name__}')

    old = work / 'old-state.xml'
    old.write_text('<State time="4.8"></State>')
    suite.check('a state with no stepCount also raises ValueError',
                isinstance(raises(util.state_xml_step_count, old), ValueError))


def check_harvest_entry_points(suite):
    """The old shims must not pin LOOS: it keeps only the diagonal of a
    triclinic cell, and only the auto choice looks at the box at all."""
    suite.section('the harvest entry points scripts on disk still call')
    calls = []
    real = harvester.harvest_generation
    harvester.harvest_generation = lambda *args, **kws: calls.append(kws)
    try:
        util.strip_and_downsample('config.json', 'hconfig.json')
        util.strip_ds_mdtraj('config.json', 'hconfig.json')
    finally:
        harvester.harvest_generation = real
    suite.check('strip_and_downsample lets the box pick the backend',
                calls[0].get('backend', harvester.BACKEND_AUTO)
                == harvester.BACKEND_AUTO, f'-> {calls[0]}')
    suite.check('strip_ds_mdtraj still pins mdtraj, right for either box',
                calls[1].get('backend') == harvester.BACKEND_MDTRAJ,
                f'-> {calls[1]}')
    suite.check('only the auto choice consults the box at all',
                harvester.select_backend('traj.dcd',
                                         backend=harvester.BACKEND_LOOS)
                == harvester.BACKEND_LOOS)


def main():
    suite = Suite('utilities_guards')
    work = harness.workdir('utilities_guards')
    check_topology_reader(suite)
    check_state_xml_step_count(suite, work)
    check_harvest_entry_points(suite)
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
