"""Every format a generation can be written in has to be trimmable.

A run killed between the trajectory report and the checkpoint report leaves
more frames than state.xml accounts for, and the campaign wants exactly
contiguous trajectories, so the extra frames have to come off before the
generation resumes. Only a .dcd could be trimmed; an .xtc or an .h5 sent the
whole generation back to be redone.

Each format gets the same three launches of a real argon box: run to a short
generation, keep a copy of the checkpoint that matches it, run further so the
trajectory overshoots that checkpoint, then trim back and resume from it. The
frames that survive the trim have to be the ones that were there, and the
finished trajectory has to have one unbroken, evenly spaced clock across the
seam -- which is the whole point of trimming rather than accepting the overshoot.
"""
import shutil
import sys
from pathlib import Path

import numpy as np

import harness
from harness import Suite

import mdfarmer
from mdfarmer import simulate as sim
from mdfarmer import utilities as util

from test_omm_preempt_progress import (
    build_argon_box, PLATFORM_NAME, WRITE_INTERVAL, TIMESTEP_PS,
)

SEED_INDEX = 0
CLONE_INDEX = 0
DIRNAME_PAD = 2
SEP = '-'
SUFFIXES = ('.dcd', '.xtc', '.h5')
# The checkpoint sits at KEEP frames while the trajectory reaches KEEP + EXTRA.
KEEP_FRAMES = 3
EXTRA_FRAMES = 2
TOTAL_FRAMES = KEEP_FRAMES + EXTRA_FRAMES
OVER_ASK = 9
PS_PER_FRAME = WRITE_INTERVAL * TIMESTEP_PS


def run_gen(farm, inputs, gen_index, suffix, steps, seed_fn, append):
    """One launch of the argon box into gen `gen_index`. Returns its path."""
    system_fn, integrator_fn, pdb_fn, _ = inputs
    return sim.omm_generation(
        traj_dir_top_level=str(farm), system_fn=system_fn, top_fn=pdb_fn,
        seed_index=SEED_INDEX, clone_index=CLONE_INDEX, gen_index=gen_index,
        title='traj-trunc', integrator_xml=integrator_fn, seed_fn=str(seed_fn),
        append=append, dirname_pad=DIRNAME_PAD, sep=SEP, traj_name='positions',
        traj_suffix=suffix, restart_name='state.xml',
        platform_name=PLATFORM_NAME, steps=steps,
        write_interval=WRITE_INTERVAL)


def coordinates(traj_p, pdb_fn):
    """Every frame's coordinates in nm, whichever format holds them."""
    import mdtraj as md
    if traj_p.suffix.lower() == '.h5':
        return md.load(str(traj_p)).xyz
    return md.load(str(traj_p), top=str(pdb_fn)).xyz


def frame_times(traj_p, pdb_fn):
    """Every frame's time in ps, from wherever the format keeps it.

    A DCD keeps one linear rule in its header instead of a stamp per frame, so
    its axis comes from dcd_frame_timing and the other two from the frames.
    """
    import mdtraj as md
    if traj_p.suffix.lower() == '.dcd':
        step0, per_step, time0, per_time = util.dcd_frame_timing(traj_p)
        n_frames = util.get_traj_len(str(traj_p), str(pdb_fn))
        return time0 + np.arange(n_frames) * per_time
    if traj_p.suffix.lower() == '.h5':
        return md.load(str(traj_p)).time
    return md.load(str(traj_p), top=str(pdb_fn)).time


def tear_last_frame(traj_p):
    """Leave the file ending partway through its final frame.

    An .h5 is torn a node at a time rather than a byte at a time: mdtraj writes
    coordinates, then time, then the cell, so a kill between two of those
    appends leaves one array a row ahead of another.
    """
    if traj_p.suffix.lower() == '.h5':
        import tables
        with tables.open_file(str(traj_p), 'a') as handle:
            handle.root.time.truncate(handle.root.time.shape[0] - 1)
        return
    data = traj_p.read_bytes()
    if traj_p.suffix.lower() == '.xtc':
        last_start = util.xtc_frame_offsets(traj_p)[-2]
    else:
        info = util.dcd_header_info(traj_p)
        frame = util.dcd_frame_size(info['with_unitcell'], info['n_atoms'])
        last_start = len(data) - frame
    traj_p.write_bytes(data[:last_start + (len(data) - last_start) // 2])


def check_format(suite, work, inputs, suffix):
    """Trim one format's generation back to its checkpoint and resume it."""
    system_fn, integrator_fn, pdb_fn, seed_fn = inputs
    farm = work / f'farm{suffix.replace(".", "_")}'
    traj_p = Path(run_gen(farm, inputs, 0, suffix,
                          KEEP_FRAMES * WRITE_INTERVAL, seed_fn, False))
    gen_dir = traj_p.parent
    kept_state = work / f'kept_state{suffix.replace(".", "_")}.xml'
    shutil.copy(gen_dir / 'state.xml', kept_state)
    kept_coords = coordinates(traj_p, pdb_fn)
    kept_times = frame_times(traj_p, pdb_fn)

    run_gen(farm, inputs, 0, suffix, EXTRA_FRAMES * WRITE_INTERVAL,
            gen_dir / 'state.xml', True)
    suite.check(f'{suffix}: the launch overshoots the kept checkpoint',
                util.get_traj_len(str(traj_p), pdb_fn) == TOTAL_FRAMES,
                f'-> {util.get_traj_len(str(traj_p), pdb_fn)}')

    reached = util.truncate_traj_to_nframes(traj_p, KEEP_FRAMES)
    suite.check(f'{suffix}: trimming reports the frame count asked for',
                reached == KEEP_FRAMES, f'-> {reached}')
    suite.check(f'{suffix}: and the file holds exactly that many',
                util.get_traj_len(str(traj_p), pdb_fn) == KEEP_FRAMES,
                f'-> {util.get_traj_len(str(traj_p), pdb_fn)}')
    suite.check(f'{suffix}: the frames that stayed are the ones that were there',
                np.array_equal(coordinates(traj_p, pdb_fn), kept_coords))
    suite.check(f'{suffix}: with their own clock untouched',
                np.allclose(frame_times(traj_p, pdb_fn), kept_times,
                            rtol=1e-6))

    size_before = traj_p.stat().st_size
    suite.check(f'{suffix}: trimming again changes nothing',
                util.truncate_traj_to_nframes(traj_p, KEEP_FRAMES) == KEEP_FRAMES
                and traj_p.stat().st_size == size_before)
    suite.check(f'{suffix}: asking for more frames than exist never grows it',
                util.truncate_traj_to_nframes(traj_p, KEEP_FRAMES + OVER_ASK)
                == KEEP_FRAMES and traj_p.stat().st_size == size_before)

    shutil.copy(kept_state, gen_dir / 'state.xml')
    run_gen(farm, inputs, 0, suffix, EXTRA_FRAMES * WRITE_INTERVAL,
            gen_dir / 'state.xml', True)
    suite.check(f'{suffix}: resuming on the trim fills the generation back up',
                util.get_traj_len(str(traj_p), pdb_fn) == TOTAL_FRAMES,
                f'-> {util.get_traj_len(str(traj_p), pdb_fn)}')
    resumed = coordinates(traj_p, pdb_fn)
    suite.check(f'{suffix}: leaving the trimmed frames still in place',
                np.array_equal(resumed[:KEEP_FRAMES], kept_coords))
    spacing = np.diff(frame_times(traj_p, pdb_fn))
    suite.check(f'{suffix}: on one unbroken clock across the seam',
                len(spacing) == TOTAL_FRAMES - 1
                and np.allclose(spacing, PS_PER_FRAME, rtol=1e-5),
                f'-> {spacing}')

    torn = gen_dir / f'torn{suffix}'
    shutil.copy(traj_p, torn)
    tear_last_frame(torn)
    whole = util.truncate_traj_to_nframes(torn, TOTAL_FRAMES)
    suite.check(f'{suffix}: a torn final frame is not counted as whole',
                whole == TOTAL_FRAMES - 1, f'-> {whole}')
    suite.check(f'{suffix}: and is gone, so what is left reads',
                len(coordinates(torn, pdb_fn)) == TOTAL_FRAMES - 1)
    return traj_p


def main():
    suite = Suite('traj_truncation')
    work = harness.workdir('traj_truncation')
    inputs = build_argon_box(work)

    suite.section('the dispatch covers what a generation can be written in')
    suite.check('every writable format can be trimmed',
                util.TRUNCATABLE_SUFFIXES == set(SUFFIXES),
                f'-> {sorted(util.TRUNCATABLE_SUFFIXES)}')
    unsupported = None
    try:
        util.truncate_traj_to_nframes(work / 'nothing.trr', 1)
    except ValueError as exc:
        unsupported = str(exc)
    suite.check('a format with no truncator says so rather than half-trimming',
                unsupported is not None and '.trr' in unsupported,
                f'-> {unsupported}')

    trimmed = {}
    for suffix in SUFFIXES:
        suite.section(f'a {suffix} generation is trimmed back and resumed')
        trimmed[suffix] = check_format(suite, work, inputs, suffix)

    suite.section('an xtc is moved, never re-encoded')
    precise = work / 'precision.xtc'
    shutil.copy(trimmed['.xtc'], precise)
    original = precise.read_bytes()
    first_frame_end = util.xtc_frame_offsets(precise)[1]
    before = mdfarmer.gmx_simulate.xtc_precision(precise)
    util.truncate_traj_to_nframes(precise, 1)
    suite.check('the frame left behind is byte for byte the one that was there',
                precise.read_bytes() == original[:first_frame_end])
    suite.check('so no re-encode can have coarsened its precision',
                before is not None
                and mdfarmer.gmx_simulate.xtc_precision(precise) == before,
                f'-> {before}')

    suite.section('the suffix is matched however it is spelled')
    shouty = work / 'SHOUTY.XTC'
    shutil.copy(trimmed['.xtc'], shouty)
    reached = util.truncate_traj_to_nframes(shouty, 1)
    frames = len(util.xtc_frame_offsets(shouty)) - 1
    suite.check('an upper-case suffix routes to the same truncator',
                reached == 1 and frames == 1, f'-> {reached}, {frames} frames')
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
