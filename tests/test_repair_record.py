"""What a repaired harvest is allowed to claim about the original it deleted.

_repair_from_symlink reconstructs the sentinel from the two outputs. After
generation 0, a generation of frames_per_gen frames and one of frames_per_gen+1
produce exactly the same two counts, so the original's frame count cannot be
recovered from them. The record has to say so rather than assert a guess.
"""
import json
import sys

import numpy as np

import harness
from harness import Suite

from mdfarmer import harvester as hv

N_ATOMS = 4
FRAMES_PER_GEN = 10
DOWNSAMPLE_FRQ = 5
FRAME_DT = 0.2


def write_traj(path, n_frames, n_atoms=N_ATOMS, frame_dt=FRAME_DT):
    """An xtc holding n_frames evenly spaced frames."""
    import mdtraj as md
    xyz = np.zeros((n_frames, n_atoms, 3), dtype=np.float32)
    box = np.tile(np.eye(3, dtype=np.float32) * 2.0, (n_frames, 1, 1))
    with md.formats.XTCTrajectoryFile(str(path), 'w') as fh:
        fh.write(xyz, time=np.arange(n_frames, dtype=np.float32) * frame_dt,
                 step=np.arange(n_frames), box=box)
    return path


def repair(gen_dir, gen_index, n_orig, frames_per_gen=FRAMES_PER_GEN,
           downsample_frq=DOWNSAMPLE_FRQ):
    """Lay out a generation whose harvest deleted the original, then repair it."""
    gen_dir.mkdir(parents=True, exist_ok=True)
    skip_first = hv.resolve_seam(n_orig, frames_per_gen, gen_index)
    first_global_index = gen_index * frames_per_gen
    n_dry, n_down = hv.expected_counts(n_orig, first_global_index,
                                       downsample_frq, skip_first)
    dry_p = write_traj(gen_dir / 'dry-prod.xtc', n_dry)
    down_p = write_traj(gen_dir / 'downsample-prod.xtc', n_down)
    traj_p = gen_dir / 'prod.xtc'
    traj_p.symlink_to(dry_p.name)
    return hv._repair_from_symlink(
        traj_p, dry_p, down_p, gen_dir / hv.SENTINEL_NAME,
        gen_dir / hv.DRY_TOPOLOGY_NAME, str(dry_p),
        frames_per_gen=frames_per_gen, gen_index=gen_index,
        first_global_index=first_global_index,
        downsample_frq=downsample_frq)


def main(frames_per_gen=FRAMES_PER_GEN, downsample_frq=DOWNSAMPLE_FRQ):
    suite = Suite('repair_record')
    work = harness.workdir('repair_record')

    suite.section('generation 0, where the two candidates differ')
    record = repair(work / 'g0', 0, frames_per_gen + 1)
    suite.check('the repair is not marked ambiguous',
                record['ambiguous'] is False, f'-> {record["ambiguous"]}')
    suite.check('it recovers the original frame count',
                record['n_orig'] == frames_per_gen + 1,
                f'-> {record["n_orig"]}')
    suite.check('the sentinel on disk carries the same record',
                json.loads((work / 'g0' / hv.SENTINEL_NAME).read_text())
                == record)

    suite.section('a later generation, where they do not')
    record = repair(work / 'g1', 1, frames_per_gen + 1)
    suite.check('the repair is marked ambiguous',
                record['ambiguous'] is True, f'-> {record["ambiguous"]}')
    suite.check('the counts it can see are still right',
                (record['n_dry'], record['n_down'])
                == (frames_per_gen, frames_per_gen // downsample_frq),
                f'-> {record["n_dry"]}, {record["n_down"]}')

    suite.section('verify_dry_chain reads n_orig only where it is trustworthy')
    # A chain is measured against generation 0's extra step-0 frame, and that
    # is the one generation whose repair is unambiguous.
    chain = hv.verify_dry_chain([work / 'g0', work / 'g1'])
    suite.check('a repaired chain is contiguous', chain['contiguous'],
                f"-> {chain['n_frames']} vs {chain['expected']}")
    ambiguous_p = work / 'g1' / hv.SENTINEL_NAME
    ambiguous_p.write_text(json.dumps(
        dict(json.loads(ambiguous_p.read_text()), n_orig=frames_per_gen)))
    flipped = hv.verify_dry_chain([work / 'g0', work / 'g1'])
    suite.check('flipping the guessed n_orig changes no verdict',
                flipped == chain, f'-> {flipped}')
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
