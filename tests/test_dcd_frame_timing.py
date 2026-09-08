"""A DCD's time axis lives in its header, and reading it must not invent one.

An .xtc stamps a step and a time on every frame; a DCD stamps neither, and its
header's istart, nsavc and delta describe the whole file at once. For anything
OpenMM's DCDReporter wrote that map is exact, so this suite runs one generation
of argon twice from the same seed, once to .dcd and once to .xtc, and checks
the reconstruction agrees frame for frame with the clock the XTC carries.

The other half is refusal. A DCD whose header nothing filled in -- mdtraj's --
would reconstruct to a frame 0 at step 0, a frame no reporter writes. That has
to come back as None rather than as a plausible-looking axis, so each field the
reconstruction rests on is knocked out in turn and the answer checked.

A DCD written by LOOS is not caught that way and cannot be: istart 1, nsavc 1
and a 0.001 delta are all positive, so the reconstruction answers with fiction.
That is why a rewritten DCD has its real axis stamped back into the header
instead -- see tests/test_dcd_axis_carry.py.
"""
import struct
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
N_FRAMES = 4
STEPS = N_FRAMES * WRITE_INTERVAL

# Header fields the reconstruction rests on, and the byte each one starts at.
HEADER_FIELD_OFFSETS = {'istart': 12, 'nsavc': 16, 'delta': 44}


def run_gen(farm, inputs, gen_index, suffix):
    """One generation of the argon box, written in `suffix`. Returns its path."""
    system_fn, integrator_fn, pdb_fn, seed_fn = inputs
    return sim.omm_generation(
        traj_dir_top_level=str(farm), system_fn=system_fn, top_fn=pdb_fn,
        seed_index=SEED_INDEX, clone_index=CLONE_INDEX, gen_index=gen_index,
        title='dcd-timing', integrator_xml=integrator_fn, seed_fn=seed_fn,
        append=False, dirname_pad=DIRNAME_PAD, sep=SEP, traj_name='positions',
        traj_suffix=suffix, restart_name='state.xml',
        platform_name=PLATFORM_NAME, steps=STEPS,
        write_interval=WRITE_INTERVAL)


def xtc_frame_clock(traj_p):
    """(steps, times) as the .xtc itself stamped them."""
    import mdtraj as md
    with md.open(str(traj_p)) as fh:
        _, time, step, _ = fh.read()
    return np.asarray(step), np.asarray(time)


def zeroed_copy(src_p, dest_p, field, offsets=HEADER_FIELD_OFFSETS):
    """A copy of a DCD with one header field set to zero."""
    dest_p.write_bytes(src_p.read_bytes())
    packer = '<f' if field == 'delta' else '<i'
    with open(dest_p, 'r+b') as fh:
        fh.seek(offsets[field])
        fh.write(struct.pack(packer, 0))
    return dest_p


def mdtraj_dcd(dest_p, n_frames=N_FRAMES):
    """A DCD written by mdtraj, which stamps no timing at all."""
    import mdtraj as md
    top = md.Topology()
    chain = top.add_chain()
    residue = top.add_residue('AR', chain)
    for i in range(2):
        top.add_atom('AR', md.element.argon, residue)
    xyz = np.zeros((n_frames, 2, 3), dtype=np.float32)
    lengths = np.tile([2.0] * 3, (n_frames, 1)).astype(np.float32)
    angles = np.tile([90.0] * 3, (n_frames, 1)).astype(np.float32)
    md.Trajectory(xyz, top, unitcell_lengths=lengths,
                  unitcell_angles=angles).save_dcd(str(dest_p))
    return dest_p


def main():
    suite = Suite('dcd_frame_timing')
    work = harness.workdir('dcd_frame_timing')
    inputs = build_argon_box(work)
    farm = work / 'farm'
    dcd_p = run_gen(farm, inputs, 0, '.dcd')
    xtc_p = run_gen(farm, inputs, 1, '.xtc')

    suite.section('the header carries the integrator that wrote it')
    info = util.dcd_header_info(dcd_p)
    dt_ps = info['delta'] * util.DCD_AKMA_PICOSECONDS
    suite.check('delta reads back as the timestep in picoseconds',
                np.isclose(dt_ps, TIMESTEP_PS, rtol=1e-6),
                f'-> {dt_ps} vs {TIMESTEP_PS}')
    suite.check('istart and nsavc are both the write interval',
                info['istart'] == WRITE_INTERVAL == info['nsavc'],
                f'-> istart={info["istart"]} nsavc={info["nsavc"]}')

    suite.section('the reconstruction agrees with the clock an xtc carries')
    timing = util.dcd_frame_timing(dcd_p)
    suite.check('a DCD written by a reporter has a timing', timing is not None)
    steps, times = xtc_frame_clock(xtc_p)
    suite.check('the two runs wrote the same number of frames',
                len(steps) == util.get_traj_len(str(dcd_p), inputs[2]) == N_FRAMES,
                f'-> {len(steps)} xtc, '
                f'{util.get_traj_len(str(dcd_p), inputs[2])} dcd')
    step0, steps_per_frame, time0, time_per_frame = timing
    frames = np.arange(len(steps))
    suite.check('every frame lands on the step the xtc stamped',
                np.array_equal(step0 + frames * steps_per_frame, steps),
                f'-> {step0 + frames * steps_per_frame} vs {steps}')
    suite.check('and on the time the xtc stamped',
                np.allclose(time0 + frames * time_per_frame, times, rtol=1e-5),
                f'-> {time0 + frames * time_per_frame} vs {times}')
    suite.check('which is the write interval of real timesteps apart',
                np.isclose(time_per_frame, WRITE_INTERVAL * TIMESTEP_PS,
                           rtol=1e-6),
                f'-> {time_per_frame}')

    suite.section('a header nothing filled in yields no axis at all')
    placeholder = mdtraj_dcd(work / 'mdtraj.dcd')
    placeholder_info = util.dcd_header_info(placeholder)
    suite.check("mdtraj's writer really does leave the placeholder header",
                placeholder_info['istart'] == 0 and placeholder_info['nsavc'] == 1,
                f'-> istart={placeholder_info["istart"]} '
                f'nsavc={placeholder_info["nsavc"]}')
    suite.check('so a DCD it wrote reports no timing',
                util.dcd_frame_timing(placeholder) is None)
    for field in HEADER_FIELD_OFFSETS:
        knocked = zeroed_copy(dcd_p, work / f'no_{field}.dcd', field)
        suite.check(f'a DCD with {field} zeroed reports no timing',
                    util.dcd_frame_timing(knocked) is None,
                    f'-> {util.dcd_frame_timing(knocked)}')

    suite.section('frame_timing answers for a DCD, from the same header')
    suite.check('the general entry point agrees with the DCD-specific one',
                util.frame_timing(dcd_p) == util.dcd_frame_timing(dcd_p),
                f'-> {util.frame_timing(dcd_p)}')
    suite.check('while an xtc answers as it always did',
                util.frame_timing(xtc_p) is not None)
    suite.check('and an .h5 answers for neither, having no step field',
                util.frame_timing(work / 'nothing.h5') is None)
    # A DCD with nothing in its header still declines, which is the case that
    # kept frame_timing out of this format: an answer there would be invented.
    suite.check('a DCD written by mdtraj still declines',
                util.frame_timing(mdtraj_dcd(work / 'unstamped.dcd')) is None)
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
