"""The sentinel is the commit point of a harvest, so it lands whole or not at all.

Everything a harvest does before the sentinel is redoable; the sentinel is what
says not to redo it. A torn sentinel fails json.loads, the generation reads as
unharvested, and the harvest runs again over a generation whose original has
already been unlinked and replaced by a symlink. Both writers -- the normal
harvest and the repair that finishes an interrupted one -- are checked here.
"""
import json
import pathlib
import sys

import numpy as np

import harness
from harness import Suite

import mdfarmer
from mdfarmer import harvester as hv

N_WATERS = 3
N_FRAMES = 11
WRITE_INTERVAL = 100
STEPS_PER_GEN = 1000
DOWNSAMPLE_FRQ = 5
SUBSET = 'name == "O"'
FRAME_DT = 0.2
SEED_INDEX = 0
CLONE_INDEX = 0
DIRNAME_PAD = 2
SEP = '-'


def write_structure(path, n_waters=N_WATERS):
    """A PDB of n_waters rigid waters, enough for a subset selection to bite."""
    import mdtraj as md
    topology = md.Topology()
    chain = topology.add_chain()
    for _ in range(n_waters):
        residue = topology.add_residue('HOH', chain)
        topology.add_atom('O', md.element.oxygen, residue)
        topology.add_atom('H1', md.element.hydrogen, residue)
        topology.add_atom('H2', md.element.hydrogen, residue)
    xyz = np.zeros((1, topology.n_atoms, 3), dtype=np.float32)
    frame = md.Trajectory(xyz, topology)
    frame.unitcell_vectors = np.eye(3, dtype=np.float32).reshape(1, 3, 3) * 2.0
    frame.save_pdb(str(path))
    return path


def write_traj(path, n_atoms, n_frames=N_FRAMES, frame_dt=FRAME_DT,
               write_interval=WRITE_INTERVAL):
    """An xtc of n_frames evenly spaced frames of n_atoms."""
    import mdtraj as md
    xyz = np.zeros((n_frames, n_atoms, 3), dtype=np.float32)
    box = np.tile(np.eye(3, dtype=np.float32) * 2.0, (n_frames, 1, 1))
    with md.formats.XTCTrajectoryFile(str(path), 'w') as fh:
        fh.write(xyz, time=np.arange(n_frames, dtype=np.float32) * frame_dt,
                 step=np.arange(n_frames) * write_interval, box=box)
    return path


def build_generation(farm, n_waters=N_WATERS, n_frames=N_FRAMES,
                     write_interval=WRITE_INTERVAL,
                     steps_per_gen=STEPS_PER_GEN,
                     downsample_frq=DOWNSAMPLE_FRQ, subset=SUBSET):
    """Generation 0 of a campaign, ready to harvest, without running an engine."""
    gen_dir = mdfarmer.dir_seeds_clones_gens(
        farm, SEED_INDEX, CLONE_INDEX, 0, DIRNAME_PAD, sep=SEP)
    structure = write_structure(gen_dir / 'system.pdb', n_waters=n_waters)
    write_traj(gen_dir / 'prod.xtc', 3 * n_waters, n_frames=n_frames,
               write_interval=write_interval)
    (gen_dir / hv.CONFIG_NAME).write_text(json.dumps(dict(
        traj_dir_top_level=str(farm), seed_index=SEED_INDEX,
        clone_index=CLONE_INDEX, gen_index=0, dirname_pad=DIRNAME_PAD, sep=SEP,
        traj_name='prod', traj_suffix='.xtc', write_interval=write_interval,
        steps=steps_per_gen, steps_per_gen=steps_per_gen,
        top_fn=str(structure)), indent=4))
    (gen_dir / 'hconfig.json').write_text(json.dumps(dict(
        downsample_frq=downsample_frq, harvester_structure=str(structure),
        harvester_subset=subset)))
    return gen_dir


class ReplaceProbe:
    """Records what each os-level rename saw, and whether anything bypassed one.

    Path.write_text is watched alongside Path.replace, because a direct write to
    the final name is exactly the regression this suite exists to catch: it
    leaves a window in which the sentinel exists but is not yet valid JSON.
    """

    def __init__(self, watched):
        self.watched = pathlib.Path(watched)
        self.renames = []
        self.direct_writes = []

    def __enter__(self):
        self._real_replace = pathlib.Path.replace
        self._real_write_text = pathlib.Path.write_text
        probe = self

        def watching_replace(self, target):
            if pathlib.Path(target) == probe.watched:
                probe.renames.append(dict(
                    tmp_dir=self.parent,
                    tmp_text=self.read_text() if self.is_file() else '',
                    target_exists=pathlib.Path(target).exists()))
            return probe._real_replace(self, target)

        def watching_write_text(self, data, *args, **kwargs):
            if self == probe.watched:
                probe.direct_writes.append(data)
            return probe._real_write_text(self, data, *args, **kwargs)

        pathlib.Path.replace = watching_replace
        pathlib.Path.write_text = watching_write_text
        return self

    def __exit__(self, *exc):
        pathlib.Path.replace = self._real_replace
        pathlib.Path.write_text = self._real_write_text
        return False


def check_one_write(suite, probe, sentinel_p, expected):
    """The checks both sentinel writers have to pass."""
    suite.check('nothing writes the sentinel name directly',
                not probe.direct_writes, f'-> {probe.direct_writes}')
    suite.check('the sentinel arrives by exactly one rename',
                len(probe.renames) == 1, f'-> {len(probe.renames)} renames')
    seen = probe.renames[0] if probe.renames else {}
    suite.check('the temp file is a sibling, so the rename is within one '
                'filesystem',
                seen.get('tmp_dir') == sentinel_p.parent,
                f'-> {seen.get("tmp_dir")}')
    suite.check('the whole record is on disk before the rename',
                json.loads(seen.get('tmp_text') or 'null') == expected,
                f'-> {seen.get("tmp_text")!r}')
    suite.check('the sentinel does not exist until the rename',
                seen.get('target_exists') is False,
                f'-> {seen.get("target_exists")}')
    suite.check('and it reads back as the record that was returned',
                json.loads(sentinel_p.read_text()) == expected)


def main():
    suite = Suite('sentinel_atomic')
    work = harness.workdir('sentinel_atomic')
    gen_dir = build_generation(work / 'farm')
    sentinel_p = gen_dir / hv.SENTINEL_NAME
    config_args = (gen_dir / hv.CONFIG_NAME, gen_dir / 'hconfig.json')

    suite.section('the sentinel a completed harvest writes')
    with ReplaceProbe(sentinel_p) as probe:
        record = hv.harvest_generation(*config_args)
    suite.check('the harvest reports it harvested',
                record['status'] == 'harvested', f'-> {record["status"]}')
    check_one_write(suite, probe, sentinel_p, record)

    suite.section('the sentinel the repair path writes')
    # The harvest unlinked the original, so dropping the sentinel leaves the
    # symlink-with-no-sentinel state that _repair_from_symlink finishes.
    sentinel_p.unlink()
    suite.check('the original really is a symlink now',
                (gen_dir / 'prod.xtc').is_symlink())
    with ReplaceProbe(sentinel_p) as probe:
        repaired = hv.harvest_generation(*config_args)
    suite.check('the repair reports it repaired',
                repaired['status'] == 'repaired', f'-> {repaired["status"]}')
    check_one_write(suite, probe, sentinel_p, repaired)

    suite.section('and a sentinel on disk still stops a second harvest')
    again = hv.harvest_generation(*config_args)
    suite.check('a harvested generation is left alone',
                again['status'] == 'already-harvested', f'-> {again["status"]}')
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
