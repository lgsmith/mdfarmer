"""The guards in utilities: what each refuses, and how it says so.

Every check here is about an error a scientist reads at 3am, so what matters is
that the type matches what the caller catches and that the message names the
thing that is actually wrong.
"""
import contextlib
import io
import subprocess as sp
import sys

import harness
from harness import Suite

from mdfarmer import harvester
from mdfarmer import utilities as util

# The module file itself, loaded standalone by the missing-backend check.
UTILITIES = str(harness.REPO_ROOT / 'utilities.py')


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
    suite.section('the harvest entry point scripts on disk still call')
    calls = []
    real = harvester.harvest_generation
    harvester.harvest_generation = lambda *args, **kws: calls.append(kws)
    try:
        util.strip_and_downsample('config.json', 'hconfig.json')
    finally:
        harvester.harvest_generation = real
    suite.check('strip_and_downsample lets the box pick the backend',
                calls[0].get('backend', harvester.BACKEND_AUTO)
                == harvester.BACKEND_AUTO, f'-> {calls[0]}')
    suite.check('it is the only old name left, so no pinned twin to pick wrong',
                not hasattr(util, 'strip_ds_mdtraj'),
                f'-> {[n for n in dir(util) if n.startswith("strip")]}')
    suite.check('only the auto choice consults the box at all',
                harvester.select_backend('traj.dcd',
                                         backend=harvester.BACKEND_LOOS)
                == harvester.BACKEND_LOOS)


def check_config_from_signature(suite):
    def a_runner(traj_name='traj', write_interval=100, steps=None, **_unused):
        return traj_name, write_interval, steps

    suite.section('the config dict recording a call')
    config = util.merge_args_defaults_dict(a_runner, write_interval=250)
    suite.check('defaults and overrides are both recorded',
                config == dict(traj_name='traj', write_interval=250,
                               steps=None), f'-> {config}')

    # traj_list is the run block's, not the runner's, and rides in the config.
    with_list = util.merge_args_defaults_dict(a_runner, traj_list='tl.txt')
    suite.check('a key the run block reads is carried, not refused',
                with_list.get('traj_list') == 'tl.txt', f'-> {with_list}')

    exc = raises(util.merge_args_defaults_dict, a_runner, write_intervall=250)
    suite.check('a misspelled keyword is refused, not written to the config',
                isinstance(exc, TypeError) and 'write_intervall' in str(exc),
                f'-> {type(exc).__name__}: {exc}')

    def needs_a_seed(seed_fn, steps=10):
        return seed_fn, steps

    exc = raises(util.merge_args_defaults_dict, needs_a_seed)
    suite.check('a required argument nobody supplied is still refused',
                isinstance(exc, TypeError) and 'seed_fn' in str(exc),
                f'-> {type(exc).__name__}: {exc}')


def check_whole_frames(suite):
    """One implementation of "steps has to land on a write_interval", so the
    Farmer's template check and an engine's per-generation check cannot drift.
    """
    suite.section('a step count off the write_interval grid')
    suite.check('a whole number of intervals passes through unchanged',
                util.check_whole_frames(10000, 500) == 10000)
    suite.check('one interval exactly is still whole',
                util.check_whole_frames(500, 500) == 500)

    exc = raises(util.check_whole_frames, 10001, 500)
    suite.check('a remainder is a ValueError, which the Farmer cascades on',
                isinstance(exc, ValueError), f'-> {type(exc).__name__}')
    suite.check('the message names both numbers and the leftover',
                exc is not None and all(s in str(exc)
                                        for s in ('10001', '500', '1')),
                f'-> {exc}')
    suite.check('fewer steps than one interval is a remainder too',
                isinstance(raises(util.check_whole_frames, 250, 500),
                           ValueError))

    exc = raises(util.check_whole_frames, 10001, 500, source='gen-03 config')
    suite.check('the caller names where the numbers came from',
                exc is not None and 'gen-03 config' in str(exc), f'-> {exc}')

    # Neither engine's config template is required to carry both keys.
    suite.check('a missing steps is nothing to check',
                util.check_whole_frames(None, 500) is None)
    suite.check('a missing write_interval is nothing to check',
                util.check_whole_frames(10001, None) == 10001)
    suite.check('a zero write_interval does not raise ZeroDivisionError',
                util.check_whole_frames(10001, 0) == 10001)


def check_frame_counting_backend(suite, work):
    """With neither mdtraj nor LOOS every trajectory would measure as empty and
    the orchestrator would delete it, so get_traj_len refuses to answer.

    The refusal is at count time rather than import time: a scheduler-only
    install never counts a frame, and failing its `import mdfarmer` would cost
    it the whole package for a backend it does not use.
    """
    suite.section('an install with no way to count frames')
    traj_p = work / 'has-frames.dcd'
    traj_p.write_bytes(b'not a real dcd, but it is a file with bytes in it')
    script = '\n'.join((
        'import importlib.util, sys',
        # None in sys.modules is what makes an import raise ImportError.
        "sys.modules['mdtraj'] = None",
        "sys.modules['loos'] = None",
        f"spec = importlib.util.spec_from_file_location('u', {UTILITIES!r})",
        'module = importlib.util.module_from_spec(spec)',
        'spec.loader.exec_module(module)',
        "print('IMPORTED')",
        f"print('MISSING:', module.get_traj_len({str(work / 'gone.dcd')!r}, None))",
        'try:',
        f'    module.get_traj_len({str(traj_p)!r}, None)',
        'except ImportError as exc:',
        "    print('REFUSED:', exc)",
        '    sys.exit(0)',
        'sys.exit(1)',
    ))
    result = sp.run([sys.executable, '-c', script], capture_output=True,
                    text=True)
    output = (result.stdout + result.stderr).strip()
    suite.check('importing utilities still works, so the scheduler side runs',
                'IMPORTED' in result.stdout, f'-> {output[-90:]}')
    suite.check('a trajectory that is not there is still 0, not an error',
                'MISSING: 0' in result.stdout, f'-> {output[-90:]}')
    suite.check('counting a trajectory that exists raises instead of saying 0',
                result.returncode == 0 and 'REFUSED:' in result.stdout,
                f'-> {output[-90:]}')
    suite.check('the message names the trajectory it was asked about',
                'has-frames.dcd' in result.stdout, f'-> {output[-120:]}')


def main():
    suite = Suite('utilities_guards')
    work = harness.workdir('utilities_guards')
    check_topology_reader(suite)
    check_state_xml_step_count(suite, work)
    check_harvest_entry_points(suite)
    check_config_from_signature(suite)
    check_whole_frames(suite)
    check_frame_counting_backend(suite, work)
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
