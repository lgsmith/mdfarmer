"""What the harvest refuses, and what it counts.

Everything here is bookkeeping, so none of it needs GROMACS.
"""
import json
import sys

import numpy as np

import harness
from harness import Suite

from mdfarmer import harvester as hv

N_ATOMS = 4
FRAMES_PER_GEN = 10
WRITE_INTERVAL = 100
DOWNSAMPLE_FRQ = 5
FRAME_DT = 0.2


def write_traj(path, times, n_atoms=N_ATOMS):
    """A trajectory whose frames carry the given times."""
    import mdtraj as md
    xyz = np.zeros((len(times), n_atoms, 3), dtype=np.float32)
    box = np.tile(np.eye(3, dtype=np.float32) * 2.0, (len(times), 1, 1))
    if str(path).endswith('.xtc'):
        with md.formats.XTCTrajectoryFile(str(path), 'w') as fh:
            fh.write(xyz, time=np.asarray(times, dtype=np.float32),
                     step=np.arange(len(times)) * WRITE_INTERVAL, box=box)
    else:
        with md.formats.DCDTrajectoryFile(str(path), 'w') as fh:
            fh.write(xyz * 10.0)
    return path


def harvested_gen(gen_dir, gen_index, times, frames_per_gen=FRAMES_PER_GEN,
                  downsample_frq=DOWNSAMPLE_FRQ):
    """A generation directory that looks already harvested."""
    gen_dir.mkdir(parents=True, exist_ok=True)
    write_traj(gen_dir / 'dry-prod.xtc', times)
    n_orig = len(times) + (0 if gen_index == 0 else 1)
    (gen_dir / hv.SENTINEL_NAME).write_text(json.dumps(dict(
        status='harvested', dry='dry-prod.xtc', downsample='downsample-prod.xtc',
        gen_index=gen_index, frames_per_gen=frames_per_gen, n_orig=n_orig,
        n_dry=len(times), n_down=frames_per_gen // downsample_frq,
        downsample_frq=downsample_frq, skip_first=gen_index > 0)))
    return gen_dir


def main(frames_per_gen=FRAMES_PER_GEN, downsample_frq=DOWNSAMPLE_FRQ,
         frame_dt=FRAME_DT):
    suite = Suite('harvest_guards')
    work = harness.workdir('harvest_guards')

    suite.section('per-frame times come only from a format that has them')
    xtc = write_traj(work / 'a.xtc', [0.0, 0.2, 0.4])
    dcd = write_traj(work / 'a.dcd', [0.0, 0.2, 0.4])
    suite.check('an xtc gives its real times',
                np.allclose(hv._frame_times(xtc), [0.0, 0.2, 0.4]))
    suite.check('a dcd gives None, not its cell lengths',
                hv._frame_times(dcd) is None)

    suite.section('a chain that starts partway through a campaign')
    # Generation 0 keeps its step-0 frame; later ones drop the seam.
    first = [i * frame_dt for i in range(frames_per_gen + 1)]
    later = [(frames_per_gen + i) * frame_dt for i in range(1, frames_per_gen + 1)]
    from_zero = [harvested_gen(work / 'z0', 0, first),
                 harvested_gen(work / 'z1', 1, later)]
    report = hv.verify_dry_chain(from_zero)
    suite.check('a chain from generation 0 expects the extra frame',
                report['contiguous'] and report['expected'] == 2 * frames_per_gen + 1,
                f"-> {report['n_frames']} vs {report['expected']}")

    mid = [harvested_gen(work / 'm1', 1, later),
           harvested_gen(work / 'm2', 2,
                         [(2 * frames_per_gen + i) * frame_dt
                          for i in range(1, frames_per_gen + 1)])]
    report = hv.verify_dry_chain(mid)
    suite.check('a chain from generation 1 expects no extra frame',
                report['contiguous'] and report['expected'] == 2 * frames_per_gen,
                f"-> {report['n_frames']} vs {report['expected']}")
    try:
        hv.verify_dry_chain([])
        suite.check('no generations at all is refused', False, '-> no exception')
    except hv.HarvestError:
        suite.check('no generations at all is refused', True)

    suite.section('the two configs must agree on the generation length')
    suite.check('agreeing configs are fine',
                hv._steps_per_gen({'steps_per_gen': 1000, 'steps': 1000},
                                  {'steps_per_gen': 1000}) == 1000)
    suite.check('the run config wins when only it says',
                hv._steps_per_gen({'steps_per_gen': 1000, 'steps': 40},
                                  {}) == 1000)
    try:
        hv._steps_per_gen({'steps_per_gen': 1000, 'steps': 1000},
                          {'steps_per_gen': 500})
        suite.check('disagreeing configs are refused', False, '-> no exception')
    except hv.HarvestError as exc:
        suite.check('disagreeing configs are refused', True, f'-> {str(exc)[:50]}')

    suite.section('a downsample frequency below one')
    gen_dir = work / 'bad'
    gen_dir.mkdir()
    (gen_dir / 'config.json').write_text(json.dumps(dict(
        sep='-', traj_name='prod', traj_suffix='.xtc', gen_index=0,
        write_interval=100, steps=1000, steps_per_gen=1000,
        top_fn=str(work / 'a.xtc'))))
    for bad in (0, -1):
        (gen_dir / 'hconfig.json').write_text(json.dumps(
            dict(downsample_frq=bad, harvester_structure=str(work / 'a.xtc'))))
        try:
            hv.harvest_generation(gen_dir / 'config.json',
                                  gen_dir / 'hconfig.json')
            suite.check(f'downsample_frq={bad} is refused', False,
                        '-> no exception')
        except hv.HarvestError as exc:
            suite.check(f'downsample_frq={bad} is refused',
                        'downsample_frq' in str(exc), f'-> {str(exc)[:45]}')

    suite.section('a dry-run harvest submits nothing')
    harvester = hv.Harvester('#!/bin/bash\necho hi\n', 'sbatch')
    suite.check('reap returns nothing and writes only its script',
                harvester.reap(work, dry_run=True) is None
                and (work / 'harvest.sh').is_file())
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
