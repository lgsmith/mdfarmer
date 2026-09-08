"""A DCD rewritten by the harvest keeps the step axis the source carried.

Neither writer that rewrites a trajectory here fills a DCD's header in. LOOS
hardcodes istart and nsavc to 1 in DCDWriter::writeHeader and takes no per-frame
timing at all; mdtraj's DCD writer takes none either. So a harvested DCD used to
claim frame k sat at step k+1 -- a clock no reporter ever wrote, and one that
disagreed with the .xtc arm of the same campaign.

The axis is put back after the writer closes, which is where a DCD keeps it. The
two output streams keep different frames, so each has to state its own spacing
rather than the source's: the dry stream every frame, the downsampled one every
Nth. Both are checked against the source's own header, and against the frames
that actually landed in each file.
"""
import sys
from pathlib import Path

import numpy as np

import harness
from harness import Suite

from mdfarmer import harvester as hv, reimage, simulate as sim, utilities as util

from test_omm_preempt_progress import (
    build_argon_box, PLATFORM_NAME, WRITE_INTERVAL, TIMESTEP_PS,
)
from test_dcd_frame_timing import run_gen

SEED_INDEX = 0
CLONE_INDEX = 0
DIRNAME_PAD = 2
SEP = '-'
DOWNSAMPLE_FRQ = 2


def frame_steps(timing, n_frames):
    """The step each frame sits at, under one linear rule."""
    step0, steps_per_frame, _, _ = timing
    return [step0 + k * steps_per_frame for k in range(n_frames)]


def main(downsample_frq=DOWNSAMPLE_FRQ, write_interval=WRITE_INTERVAL):
    suite = Suite('dcd_axis_carry')
    work = harness.workdir('dcd_axis_carry')
    inputs = build_argon_box(work)
    farm = work / 'farm'
    source = run_gen(farm, inputs, 0, '.dcd')
    gen_dir = source.parent

    suite.section('what the reporter wrote')
    timing = util.dcd_frame_timing(source)
    suite.check('the source DCD carries an axis', timing is not None,
                f'-> {timing}')
    n_source = util.dcd_header_info(source)['nset']
    suite.check('frame_timing answers for a .dcd now',
                util.frame_timing(source) == timing,
                f'-> {util.frame_timing(source)}')

    suite.section('a stream that keeps every frame keeps the source spacing')
    dry_out = gen_dir / 'dry.dcd'
    axis = reimage.WrittenAxis()
    for index in range(n_source):
        axis.took(index)
    suite.check('an unbroken run of indices is uniform', axis.uniform())
    # Write a real DCD through the LOOS path, then carry the axis onto it.
    import loos
    model = loos.createSystem(str(inputs[2]))
    writer = reimage._loos_writer(dry_out)
    for _ in range(n_source):
        writer.writeFrame(model)
    del writer
    # Not None: LOOS writes istart 1, nsavc 1 and a 0.001 AKMA delta, all
    # positive, so the reconstruction answers with a fabricated axis rather
    # than refusing. That is what makes stamping necessary and not merely tidy.
    left_behind = util.dcd_frame_timing(dry_out)
    suite.check('LOOS leaves an axis that is fiction, not an absence',
                left_behind is not None and left_behind[:2] == (1, 1),
                f'-> {left_behind}')
    suite.check('and it disagrees with what the reporter wrote',
                left_behind[:2] != timing[:2], f'-> {left_behind} vs {timing}')
    stamped = reimage.stamp_written_axis(dry_out, timing, axis)
    suite.check('stamping reports that it wrote a header', stamped)
    carried = util.dcd_frame_timing(dry_out)
    suite.check('the rewritten DCD reports the source axis',
                carried is not None and carried[:2] == timing[:2],
                f'-> {carried} vs {timing}')
    suite.check('so its frames sit at the steps the reporter used',
                frame_steps(carried, n_source) == frame_steps(timing, n_source),
                f'-> {frame_steps(carried, n_source)}')

    suite.section('a downsampled stream states its own spacing, not the source\'s')
    down_out = gen_dir / 'down.dcd'
    kept = list(range(0, n_source, downsample_frq))
    down_axis = reimage.WrittenAxis()
    writer = reimage._loos_writer(down_out)
    for index in kept:
        writer.writeFrame(model)
        down_axis.took(index)
    del writer
    reimage.stamp_written_axis(down_out, timing, down_axis)
    down_timing = util.dcd_frame_timing(down_out)
    suite.check(f'its spacing is {downsample_frq}x the source\'s',
                down_timing[1] == timing[1] * downsample_frq,
                f'-> {down_timing[1]} vs {timing[1] * downsample_frq}')
    source_steps = frame_steps(timing, n_source)
    suite.check('and its frames land on the source steps it actually kept',
                frame_steps(down_timing, len(kept))
                == [source_steps[k] for k in kept],
                f'-> {frame_steps(down_timing, len(kept))} vs '
                f'{[source_steps[k] for k in kept]}')

    suite.section('an axis it cannot state is refused, not rounded')
    gappy = reimage.WrittenAxis()
    for index in (0, 1, 5):
        gappy.took(index)
    suite.check('an uneven set of kept frames is not uniform',
                not gappy.uniform())
    try:
        reimage.stamp_written_axis(down_out, timing, gappy)
        suite.check('stamping an uneven stream raises', False, '-> no exception')
    except ValueError as exc:
        suite.check('stamping an uneven stream raises', True,
                    f'-> {str(exc)[:60]}')

    suite.section('nothing is stamped when there is nothing to carry')
    suite.check('no timing means no header written',
                reimage.stamp_written_axis(down_out, None, down_axis) is False)
    suite.check('an .xtc output is left to its own per-frame stamps',
                reimage.stamp_written_axis(
                    gen_dir / 'x.xtc', timing, down_axis) is False)
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
