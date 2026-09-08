"""truncate_dcd_to_nframes must trust the file's bytes, not the header.

A DCD's nset field can disagree with what is actually on disk: a node crash
can let the header page reach disk while the tail data pages do not, or a
kill can land between the nset bump and the frame write it announces. Every
case here builds an honest DCD with mdtraj, corrupts it by hand the same way,
and checks that truncation ends with the header and the byte length agreeing
on a number of frames the file can actually back up.
"""
import os
import struct
import sys

import numpy as np

import harness
from harness import Suite

from mdfarmer import utilities as util

N_ATOMS = 4
N_FRAMES = 10
BOX_NM = 2.0
# md.open().read() is the raw DCD reader and hands back its native unit,
# Angstroms, unlike Trajectory.xyz which mdtraj always keeps in nanometres.
ANGSTROM_PER_NM = 10.0


def _make_topology(n_atoms):
    import mdtraj as md
    top = md.Topology()
    chain = top.add_chain()
    residue = top.add_residue('RES', chain)
    for i in range(n_atoms):
        top.add_atom(f'A{i}', md.element.carbon, residue)
    return top


def _write_dcd(path, n_frames=N_FRAMES, n_atoms=N_ATOMS, box_nm=BOX_NM):
    """An honest DCD: n_frames of distinct, reproducible coordinates."""
    top = _make_topology(n_atoms)
    xyz = np.random.RandomState(0).rand(n_frames, n_atoms, 3).astype(np.float32)
    lengths = np.tile([box_nm] * 3, (n_frames, 1)).astype(np.float32)
    angles = np.tile([90.0] * 3, (n_frames, 1)).astype(np.float32)
    import mdtraj as md
    traj = md.Trajectory(xyz, top, unitcell_lengths=lengths, unitcell_angles=angles)
    traj.save_dcd(str(path))
    return xyz


def _set_nset(path, nset):
    with open(path, 'r+b') as f:
        f.seek(8)
        f.write(struct.pack('<i', nset))


def _load_xyz(path, angstrom_per_nm=ANGSTROM_PER_NM):
    import mdtraj as md
    return md.open(str(path)).read()[0] / angstrom_per_nm


def main(n_frames=N_FRAMES, n_atoms=N_ATOMS):
    suite = Suite('dcd_truncation')
    work = harness.workdir('dcd_truncation')

    xyz = _write_dcd(work / 'src.dcd', n_frames=n_frames, n_atoms=n_atoms)
    info = util.dcd_header_info(work / 'src.dcd')
    header_size = info['header_size']
    frame_size = util.dcd_frame_size(info['with_unitcell'], info['n_atoms'])

    suite.section('sanity: our frame-size math matches what mdtraj wrote')
    suite.check('nset matches the frames we asked for', info['nset'] == n_frames)
    suite.check('file length matches header_size + nset*frame_size',
                (work / 'src.dcd').stat().st_size
                == header_size + n_frames * frame_size)

    def fresh(name):
        """A private copy of the honest source DCD to corrupt."""
        import shutil
        dest = work / name
        shutil.copy(work / 'src.dcd', dest)
        return dest

    suite.section('a healthy DCD truncates exactly as before')
    p = fresh('healthy.dcd')
    before_size = p.stat().st_size
    result = util.truncate_dcd_to_nframes(p, 6)
    after_info = util.dcd_header_info(p)
    suite.check('returns the requested count', result == 6, f'-> {result}')
    suite.check('nset matches the request', after_info['nset'] == 6)
    suite.check('file length matches header_size + 6*frame_size',
                p.stat().st_size == header_size + 6 * frame_size)
    suite.check('never grows the file', p.stat().st_size <= before_size)
    kept = _load_xyz(p)
    suite.check('retained frames are byte-identical to the source',
                np.allclose(kept, xyz[:6], atol=1e-6))

    suite.section('header claims more frames than the bytes hold')
    # A node crash whose header page reached disk but whose tail data pages
    # did not: nset still says n_frames, but the file was cut short.
    p = fresh('overreport.dcd')
    os.truncate(str(p), header_size + 6 * frame_size)
    before_size = p.stat().st_size
    result = util.truncate_dcd_to_nframes(p, 8)  # between 6 real and 10 claimed
    after_info = util.dcd_header_info(p)
    suite.check('returns the honest count, not the request',
                result == 6, f'-> {result}')
    suite.check('does not grow the file to manufacture frames',
                p.stat().st_size <= before_size, f'-> {p.stat().st_size}')
    suite.check('header nset now matches the bytes', after_info['nset'] == 6)
    suite.check('file length matches header_size + 6*frame_size',
                p.stat().st_size == header_size + 6 * frame_size)
    kept = _load_xyz(p)
    suite.check('the 6 real frames are untouched', np.allclose(kept, xyz[:6], atol=1e-6))

    p2 = fresh('overreport_at_nset.dcd')
    os.truncate(str(p2), header_size + 6 * frame_size)
    result2 = util.truncate_dcd_to_nframes(p2, n_frames)  # request == old nset
    suite.check('asking for the stale nset still gets the honest count',
                result2 == 6, f'-> {result2}')
    suite.check('bytes and header agree when the target equals nset',
                p2.stat().st_size == header_size + 6 * frame_size
                and util.dcd_header_info(p2)['nset'] == 6)

    suite.section('header claims fewer frames than the bytes hold')
    # A stale nset: the real data is all there, but the field lags behind it.
    p = fresh('underreport.dcd')
    _set_nset(p, 6)
    before_size = p.stat().st_size
    result = util.truncate_dcd_to_nframes(p, 8)  # between stale 6 and real 10
    after_info = util.dcd_header_info(p)
    suite.check('does not refuse or raise on a request above the stale nset',
                result == 8, f'-> {result}')
    suite.check('never grows the file', p.stat().st_size <= before_size)
    suite.check('header nset matches the achieved count', after_info['nset'] == 8)
    suite.check('file length matches header_size + 8*frame_size',
                p.stat().st_size == header_size + 8 * frame_size)
    kept = _load_xyz(p)
    suite.check('the kept frames are the real first 8', np.allclose(kept, xyz[:8], atol=1e-6))

    p2 = fresh('underreport_at_nset.dcd')
    _set_nset(p2, 6)
    result2 = util.truncate_dcd_to_nframes(p2, 6)  # request == the stale nset
    suite.check('a request matching the stale nset still trims real bytes',
                result2 == 6, f'-> {result2}')
    suite.check('the file actually shrinks to 6 real frames, not left at 10',
                p2.stat().st_size == header_size + 6 * frame_size)

    suite.section('file length is not a whole number of frames')
    # nset bumped for a frame whose bytes never fully landed.
    p = fresh('torn.dcd')
    with open(p, 'ab') as f:
        f.write(b'\x00' * (frame_size // 2))
    _set_nset(p, n_frames + 1)
    before_size = p.stat().st_size
    result = util.truncate_dcd_to_nframes(p, n_frames + 1)  # request == nset
    after_info = util.dcd_header_info(p)
    suite.check('the partial trailing frame is not counted',
                result == n_frames, f'-> {result}')
    suite.check('never grows the file', p.stat().st_size <= before_size)
    suite.check('the torn tail bytes are discarded',
                p.stat().st_size == header_size + n_frames * frame_size)
    suite.check('header nset matches the whole frames left',
                after_info['nset'] == n_frames)
    kept = _load_xyz(p)
    suite.check('the whole frames are untouched', np.allclose(kept, xyz, atol=1e-6))

    suite.section('the two frame-counting paths agree once truncated')
    for name in ('overreport.dcd', 'underreport.dcd', 'torn.dcd'):
        p = work / name
        via_header = util.dcd_header_info(p)['nset']
        via_bytes = util.get_traj_len(str(p), None)
        suite.check(f'{name}: header nset matches mdtraj byte-derived length',
                    via_header == via_bytes, f'-> {via_header} vs {via_bytes}')

    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
