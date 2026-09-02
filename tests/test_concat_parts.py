"""concat_parts merges a generation's parts into one contiguous trajectory.

The merge is what the whole campaign's contiguity rests on, so the checks are
against first principles rather than against another tool: real mdrun parts are
merged and the result must carry exactly the step and time axis the run's own
budget and write interval predict, with every step present once.

Where two parts cover the same steps the later part's frames win. That is
checked directly, by giving the two parts recognisably different coordinates.

Which writer does the merging follows from the parts' own compressed-x-precision:
mdtraj at the default 1000, LOOS above it. Both are checked by breaking one
writer and merging, and the grid each merge lands on is measured rather than
assumed.
"""
import contextlib
import io
import shutil
import struct
import subprocess as sp
import sys
from pathlib import Path

import numpy as np
import mdtraj as md
from mdtraj.formats import XTCTrajectoryFile

import harness
from harness import Suite

from mdfarmer import gmx_simulate as gs

# Long enough that a deliberate overlap covers several frames.
STEPS = 2000
WRITE_INTERVAL = 100

# Where the checkpoint that splits the run into two real mdrun parts falls.
SPLIT_STEP = 800

# Frames of the deliberate overlap: part0002 restamped over part0001's tail.
OVERLAP_FRAMES = 6

# compressed-x-precision a campaign might raise to, and which mdtraj's xtc
# writer -- fixed at 1000 -- cannot carry.
FINE_PRECISION = 10000

# How far a measured coordinate grid may sit from 1 / precision and still count.
QUANTUM_RTOL = 0.01

# Byte offset and length of the 3x3 box in an .xtc frame: it follows the magic
# number, the atom count, the step and the time.
XTC_BOX_OFFSET = 16
XTC_BOX_BYTES = 36

# Off-diagonal box element, nm, that makes a copied cell unmistakably triclinic.
BOX_SHEAR_NM = 0.5

# Atoms dropped to build a structure that is not the system the parts hold.
SHORT_STRUCTURE_DROP = 3

# Timestep in harness.water_mdp, ps; frame time is step * dt.
DT_PS = 0.002


def gmx(cmd, cwd, gmx_bin=harness.GMX_BIN, stdin=None):
    """Run a gmx subcommand that builds the fixture, or Skip if it fails."""
    result = sp.run([str(gmx_bin), *[str(c) for c in cmd]], cwd=str(cwd),
                    input=stdin, text=True, capture_output=True)
    if result.returncode != 0:
        raise harness.Skip(
            f'gmx {cmd[0]} failed, so there are no parts to merge:'
            f'\n{result.stderr[-2000:]}')
    return result


def read_xtc(path):
    """(xyz, time, step, box) for a whole .xtc, as arrays."""
    with md.open(str(path)) as handle:
        return tuple(np.asarray(a) for a in handle.read())


def write_xtc(path, xyz, time, step, box):
    with XTCTrajectoryFile(str(path), 'w') as handle:
        handle.write(xyz, time=time, step=step, box=box)


def build_real_parts(work, precision=None, gmx_bin=harness.GMX_BIN,
                     steps=STEPS, split_step=SPLIT_STEP,
                     write_interval=WRITE_INTERVAL):
    """A gen dir holding two parts mdrun really wrote, split at split_step."""
    work.mkdir(parents=True, exist_ok=True)
    structure, topology, mdp = harness.build_water_system(work, gmx_bin=gmx_bin)
    gs.write_gen_mdp(mdp, mdp, nsteps=steps, nstxout_compressed=write_interval,
                     gen_vel=True, continuation=False)
    if precision is not None:
        mdp.write_text(f'{mdp.read_text()}'
                       f'compressed-x-precision = {precision}\n')
    gmx(['grompp', '-f', mdp.name, '-c', structure.name, '-p', topology.name,
         '-o', 'prod.tpr', '-maxwarn', '3'], work, gmx_bin=gmx_bin)
    common = ['mdrun', '-s', 'prod.tpr', '-deffnm', 'prod', '-noappend',
              '-cpo', 'prod.cpt', '-nb', 'cpu', '-pme', 'cpu', '-ntomp', '2']
    gmx(common + ['-nsteps', str(split_step)], work, gmx_bin=gmx_bin)
    gmx(common + ['-cpi', 'prod.cpt'], work, gmx_bin=gmx_bin)
    return structure, topology


def coordinate_quantum(xyz):
    """The smallest gap between distinct coordinate values, nm.

    An .xtc stores coordinates as integers divided by its precision, so with
    this many atoms the occupied values land densely enough on that grid for
    the closest pair to be one quantum apart.
    """
    values = np.unique(np.asarray(xyz, dtype=np.float64).ravel())
    gaps = np.diff(values)
    return float(gaps[gaps > 0].min())


def expected_frames(early, later):
    """The frames a merge of two overlapping parts should carry, and no others."""
    keep = early[2] < later[2][0]
    return tuple(np.concatenate([e[keep], l]) for e, l in zip(early, later))


def shear_boxes(src, dst, shear_nm=BOX_SHEAR_NM, box_offset=XTC_BOX_OFFSET,
                box_bytes=XTC_BOX_BYTES):
    """Copy an .xtc frame by frame, tilting every cell and touching nothing else.

    The compressed coordinate block of each frame is carried over byte for byte,
    so the copy is triclinic at exactly the precision the source was written at.
    """
    with XTCTrajectoryFile(str(src)) as handle:
        offsets = [int(o) for o in handle.offsets]
    raw = src.read_bytes()
    out = bytearray()
    for start, stop in zip(offsets, offsets[1:] + [len(raw)]):
        frame = bytearray(raw[start:stop])
        box = list(struct.unpack('>9f', frame[box_offset:box_offset + box_bytes]))
        box[3] = shear_nm            # yx, which no rectangular cell carries
        frame[box_offset:box_offset + box_bytes] = struct.pack('>9f', *box)
        out += frame
    dst.write_bytes(bytes(out))


def short_structure(src, dst, drop=SHORT_STRUCTURE_DROP):
    """Copy a .gro without its last few atoms, so it is a different system."""
    lines = src.read_text().splitlines()
    kept = lines[2:-1][:-drop]
    dst.write_text('\n'.join([lines[0], str(len(kept)), *kept, lines[-1]]) + '\n')


def refuses(suite, label, call, wanted, exc_type=Exception):
    """Check a call raises, and that the message says which constraint bit."""
    try:
        call()
    except exc_type as exc:
        message = str(exc)
        print(f'   {message}', flush=True)
        suite.check(label, all(w in message for w in wanted),
                    f'-> {message[:80]}')
        return
    suite.check(label, False, '-> no exception')


@contextlib.contextmanager
def broken(module, name):
    """Swap a module attribute for a raiser, to prove it is never reached."""
    def refuse(*args, **kwargs):
        raise AssertionError(f'{name} was called')
    original = getattr(module, name)
    setattr(module, name, refuse)
    try:
        yield
    finally:
        setattr(module, name, original)


def merged_without(suite, label, name, call, module=gs):
    """Merge with one writer broken, checking it was not the one that ran."""
    try:
        with broken(module, name):
            result = call()
    except AssertionError:
        suite.check(label, False, f'-> {name} did the merge')
        return call()
    suite.check(label, True)
    return result


def contiguous_frames(suite, label, merged, last_step,
                      write_interval=WRITE_INTERVAL, dt_ps=DT_PS):
    """Check a merge carries exactly the clock its run's budget predicts."""
    xyz, time, step, box = merged
    wanted = np.arange(0, last_step + 1, write_interval)
    suite.check(f'{label}: one frame per write interval, none missing',
                len(step) == len(wanted), f'-> {len(step)} vs {len(wanted)}')
    suite.check(f'{label}: the MD steps are exactly the expected sequence',
                np.array_equal(step, wanted), f'-> {step}')
    suite.check(f'{label}: no step is written twice',
                len(set(step.tolist())) == len(step))
    suite.check(f'{label}: frame times follow the steps at dt',
                np.allclose(time, wanted * dt_ps, atol=1e-6),
                f'-> {time[:3]} ...')
    suite.check(f'{label}: every frame carries a box',
                bool(np.all(np.abs(box).sum(axis=(1, 2)) > 0)))


def main(gmx_bin=harness.GMX_BIN, overlap_frames=OVERLAP_FRAMES,
         fine_precision=FINE_PRECISION):
    suite = Suite('concat_parts')
    work = harness.workdir('concat_parts')
    harness.require_gmx(gmx_bin=gmx_bin)

    suite.section('two parts mdrun really wrote, merged both ways')
    gen = work / 'real'
    build_real_parts(gen, gmx_bin=gmx_bin)
    parts = gs.part_files(gen)
    starts = [gs._first_frame_step(p) for p in parts]
    suite.check('mdrun left two parts to merge', len(parts) == 2,
                f'-> {[p.name for p in parts]} starting at {starts}')
    merged = gs.concat_parts(gen, gen / 'prod.xtc')
    contiguous_frames(suite, 'real parts', read_xtc(merged), STEPS)
    suite.check('the parts are left on disk', all(p.is_file() for p in parts))

    suite.section('a multi-frame overlap: the later part wins')
    # part0002 is restamped onto part0001's last frames while carrying the
    # coordinates of the run's opening frames, so a merge that kept the earlier
    # part's frames is unmistakable rather than a rounding argument.
    over = work / 'overlap'
    over.mkdir()
    xyz, time, step, box = read_xtc(gen / 'prod.part0001.xtc')
    n = len(step)
    keep = slice(0, n)
    write_xtc(over / 'prod.part0001.xtc', xyz[keep], time[keep], step[keep],
              box[keep])
    tail = slice(n - overlap_frames, n)
    write_xtc(over / 'prod.part0002.xtc', xyz[0:overlap_frames], time[tail],
              step[tail], box[tail])
    over_merged = gs.concat_parts(over, over / 'prod.xtc')
    contiguous_frames(suite, 'overlap', read_xtc(over_merged), SPLIT_STEP)

    mine = read_xtc(over_merged)
    early = read_xtc(over / 'prod.part0001.xtc')
    later = read_xtc(over / 'prod.part0002.xtc')
    overlap_steps = np.intersect1d(early[2], later[2])
    suite.check('the two parts really do overlap',
                len(overlap_steps) == overlap_frames, f'-> {overlap_steps}')
    in_merge = np.isin(mine[2], overlap_steps)
    from_later = np.isin(later[2], overlap_steps)
    from_early = np.isin(early[2], overlap_steps)
    to_later = float(np.abs(mine[0][in_merge] - later[0][from_later]).max())
    to_early = float(np.abs(mine[0][in_merge] - early[0][from_early]).max())
    suite.check('the overlapping frames are the later part\'s, exactly',
                to_later == 0.0, f'-> {to_later:g} nm from part0002')
    suite.check('and are not the earlier part\'s', to_early > 0.1,
                f'-> {to_early:.3f} nm from part0001')
    suite.section('the verification pass catches a merge that kept the wrong frames')
    # Built by hand out of the earlier part's overlap frames, which is exactly
    # what an inverted overlap rule would produce.
    wrong = over / 'wrong.xtc'
    after = later[2] > early[2][-1]
    write_xtc(wrong, *[np.concatenate([e, l[after]])
                       for e, l in zip(early, later)])
    try:
        gs.verify_merge([over / 'prod.part0001.xtc',
                         over / 'prod.part0002.xtc'], wrong)
        suite.check('a merge built the wrong way round is refused', False,
                    '-> no exception')
    except RuntimeError as exc:
        suite.check('a merge built the wrong way round is refused', True,
                    f'-> {str(exc)[:60]}')

    suite.section('a single part is copied, not rewritten')
    solo = work / 'solo'
    solo.mkdir()
    shutil.copy(gen / 'prod.part0001.xtc', solo / 'prod.part0001.xtc')
    solo_out = gs.concat_parts(solo, solo / 'prod.xtc')
    suite.check('the one part comes through byte for byte',
                solo_out.read_bytes()
                == (gen / 'prod.part0001.xtc').read_bytes())
    suite.check('and is still on disk to be merged again',
                (solo / 'prod.part0001.xtc').is_file())
    suite.check('no temporary file is left behind',
                not list(solo.glob('*concat-tmp*')))

    suite.section('parts out of order are refused, not spliced')
    backwards = work / 'backwards'
    backwards.mkdir()
    write_xtc(backwards / 'prod.part0001.xtc', later[0], later[1], later[2],
              later[3])
    write_xtc(backwards / 'prod.part0002.xtc', early[0], early[1], early[2],
              early[3])
    try:
        gs.concat_parts(backwards, backwards / 'prod.xtc')
        suite.check('a part starting before its predecessor is refused', False,
                    '-> no exception')
    except ValueError as exc:
        suite.check('a part starting before its predecessor is refused', True,
                    f'-> {str(exc)[:60]}')

    suite.section('a format that cannot survive the round trip is refused')
    try:
        gs.concat_parts(gen, gen / 'prod.trr', traj_suffix='.trr')
        suite.check('a .trr is refused rather than silently stripped', False,
                    '-> no exception')
    except ValueError as exc:
        suite.check('a .trr is refused rather than silently stripped',
                    'velocities' in str(exc), f'-> {str(exc)[:60]}')

    suite.section('the default merge is mdtraj\'s, and bit exact')
    default_want = expected_frames(read_xtc(parts[0]), read_xtc(parts[1]))
    suite.check('every axis of the merge is the parts\' own bytes',
                all(np.array_equal(g, w)
                    for g, w in zip(read_xtc(merged), default_want)))
    suite.check(f'the parts read back at precision {gs.MDTRAJ_XTC_PRECISION:g}',
                gs.xtc_precision(parts[0]) == gs.MDTRAJ_XTC_PRECISION)
    merged_without(suite, 'and it never reaches the LOOS writer',
                   '_write_merged_loos',
                   lambda: gs.concat_parts(gen, gen / 'again.xtc'))

    suite.section(f'a run at compressed-x-precision {fine_precision} keeps it')
    fine = work / 'fine'
    structure, _ = build_real_parts(fine, precision=fine_precision,
                                    gmx_bin=gmx_bin)
    fine_parts = gs.part_files(fine)
    suite.check('the precision is read out of the part header',
                gs.xtc_precision(fine_parts[0]) == fine_precision,
                f'-> {gs.xtc_precision(fine_parts[0])}')
    source_quantum = coordinate_quantum(read_xtc(fine_parts[0])[0])
    suite.check('and the parts really do sit on that finer grid',
                np.isclose(source_quantum, 1 / fine_precision,
                           rtol=QUANTUM_RTOL), f'-> {source_quantum:.3e} nm')

    gs._warned_precisions.clear()
    said = io.StringIO()
    with contextlib.redirect_stdout(said):
        fine_merged = gs.concat_parts(fine, fine / 'prod.xtc',
                                      structure_fn=structure)
    warning = said.getvalue()
    print(f'   {warning.strip()}', flush=True)
    suite.check('the merge asks whether the finer grid is really wanted',
                'Are you sure' in warning)
    suite.check('and says what the accuracy is being spent on',
                'force field' in warning and 'noise' in warning)
    said = io.StringIO()
    with contextlib.redirect_stdout(said):
        gs.concat_parts(fine, fine / 'twice.xtc', structure_fn=structure)
    suite.check('the warning is said once, not once per generation',
                'Are you sure' not in said.getvalue())
    merged_without(suite, 'and the merge never reaches the mdtraj writer',
                   '_write_merged',
                   lambda: gs.concat_parts(fine, fine / 'loos-only.xtc',
                                           structure_fn=structure))

    merged_quantum = coordinate_quantum(read_xtc(fine_merged)[0])
    suite.check('the merge lands on the same grid the parts did',
                np.isclose(merged_quantum, 1 / fine_precision,
                           rtol=QUANTUM_RTOL), f'-> {merged_quantum:.3e} nm')
    contiguous_frames(suite, 'fine parts', read_xtc(fine_merged), STEPS)
    fine_want = expected_frames(read_xtc(fine_parts[0]), read_xtc(fine_parts[1]))
    fine_got = read_xtc(fine_merged)
    for axis, name in enumerate(('coordinates', 'frame times', 'MD steps',
                                 'boxes')):
        suite.check(f'the merged {name} are the parts\' own, exactly',
                    np.array_equal(fine_got[axis], fine_want[axis]))

    suite.section('what the fine-precision writer cannot do, it refuses')
    refuses(suite, 'a fine merge with no structure to write frames out of',
            lambda: gs.concat_parts(fine, fine / 'nostructure.xtc'),
            ['structure_fn', 'compressed-x-precision'], ValueError)
    suite.check('and published nothing',
                not (fine / 'nostructure.xtc').is_file()
                and not list(fine.glob('nostructure*concat-tmp*')))

    short = fine / 'short.gro'
    short_structure(structure, short)
    refuses(suite, 'a structure that is not the system the parts hold',
            lambda: gs.concat_parts(fine, fine / 'mismatch.xtc',
                                    structure_fn=short),
            ['atoms', 'compressed-x-grps'], ValueError)

    tric = work / 'triclinic'
    tric.mkdir()
    for part in fine_parts:
        shear_boxes(part, tric / part.name)
    suite.check('the sheared copies kept the precision they were written at',
                gs.xtc_precision(tric / fine_parts[0].name) == fine_precision)
    suite.check('and are triclinic',
                not gs.is_orthorhombic(read_xtc(tric / fine_parts[0].name)[3][0]))
    refuses(suite, 'a triclinic cell at a precision no writer holds',
            lambda: gs.concat_parts(tric, tric / 'prod.xtc',
                                    structure_fn=structure),
            ['non-orthorhombic', 'compressed-x-precision', 'three numbers'],
            gs.BoxTypeError)
    suite.check('and published nothing',
                not (tric / 'prod.xtc').is_file()
                and not list(tric.glob('*concat-tmp*')))

    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
