"""Sorting the generations a campaign never harvested, and acting on the few
that are provably safe.

A harvest deletes the only copy of the original, so a generation the tender
missed is judged on evidence before anything is touched: the engine's own
record must say it reached the step the chain implies, and the trajectory must
hold exactly the frames the configs predict. This builds one generation of
every category the classifier knows, checks each verdict, and then checks that
acting on the report touches the safe ones and nothing else.

The last section is the one that matters. A generation whose trajectory is full
but whose checkpoint is behind harvests perfectly happily -- the frames are all
there -- even though they are a branch the next generation rewound past. The
classifier is the only thing that stops it, so that is asserted directly: with
the guard in place the harvest never runs, and with the guard widened it does.
"""
import json
import shutil
import sys
from pathlib import Path

import numpy as np

import harness
from harness import Suite

import mdfarmer
from mdfarmer import harvester as hv
from mdfarmer import gmx_simulate

N_ATOMS = 6
BOX_NM = 2.0
TIMESTEP_PS = 0.002
WRITE_INTERVAL = 100
STEPS_PER_GEN = 1000
FRAMES_PER_GEN = STEPS_PER_GEN // WRITE_INTERVAL
DOWNSAMPLE_FRQ = 5
SEED_INDEX = 0
DIRNAME_PAD = 2
SEP = '-'
TRAJ_NAME = 'prod'
TRAJ_SUFFIX = '.xtc'
GMX_RESTART_NAME = 'state.cpt'
OMM_RESTART_NAME = 'state.xml'
STRUCTURE_NAME = 'conf.pdb'

# Frames a harvest of generation 1 writes, worked out by hand so the suite does
# not check the plan against itself.
GEN1_DRY, GEN1_DOWN = 10, 2


# ---------------------------------------------------------------------------
# Building generation directories
# ---------------------------------------------------------------------------

def gen_dir(top_level, clone_index, gen_index):
    """The directory a Clone would give this seed/clone/generation."""
    return mdfarmer.dir_seeds_clones_gens(
        Path(top_level), SEED_INDEX, clone_index, gen_index, DIRNAME_PAD,
        sep=SEP)


def write_traj(path, n_frames, first_step=0, n_atoms=N_ATOMS, box_nm=BOX_NM,
               write_interval=WRITE_INTERVAL, timestep_ps=TIMESTEP_PS):
    """An evenly spaced XTC carrying real steps and times."""
    import mdtraj as md
    step = first_step + np.arange(n_frames) * write_interval
    with md.formats.XTCTrajectoryFile(str(path), 'w') as fh:
        fh.write(np.zeros((n_frames, n_atoms, 3), dtype=np.float32),
                 time=(step * timestep_ps).astype(np.float32), step=step,
                 box=np.tile(np.eye(3, dtype=np.float32) * box_nm,
                             (n_frames, 1, 1)))
    return Path(path)


def write_structure(path, n_atoms=N_ATOMS, box_nm=BOX_NM):
    """A PDB both backends can build a model from, matching write_traj."""
    import mdtraj as md
    top = md.Topology()
    chain = top.add_chain()
    for _ in range(n_atoms):
        top.add_atom('O', md.element.oxygen, top.add_residue('HOH', chain))
    traj = md.Trajectory(np.zeros((1, n_atoms, 3), dtype=np.float32), top)
    traj.unitcell_vectors = np.eye(3, dtype=np.float32).reshape(1, 3, 3) * box_nm
    traj.save_pdb(str(path))
    return Path(path)


def write_config(gen_p, top_level, clone_index, gen_index, structure,
                 restart_name=GMX_RESTART_NAME, **overrides):
    """The run record a generation directory carries."""
    config = dict(traj_dir_top_level=str(top_level), seed_index=SEED_INDEX,
                  clone_index=clone_index, gen_index=gen_index,
                  dirname_pad=DIRNAME_PAD, sep=SEP, traj_name=TRAJ_NAME,
                  traj_suffix=TRAJ_SUFFIX, write_interval=WRITE_INTERVAL,
                  steps=STEPS_PER_GEN, steps_per_gen=STEPS_PER_GEN,
                  top_fn=str(structure), restart_name=restart_name)
    config.update(overrides)
    (gen_p / hv.CONFIG_NAME).write_text(json.dumps(config))
    return config


def write_hconfig(gen_p, structure, **overrides):
    """The harvester config Harvester.reap would have written."""
    hconfig = dict(downsample_frq=DOWNSAMPLE_FRQ,
                   harvester_structure=str(structure),
                   steps_per_gen=STEPS_PER_GEN)
    hconfig.update(overrides)
    (gen_p / hv.HARVESTER_CONFIG_NAME).write_text(json.dumps(hconfig))
    return hconfig


def write_status(gen_p, target_step, reached_step, complete):
    """What the GROMACS runner leaves behind as its progress record."""
    (gen_p / hv.GEN_STATUS_NAME).write_text(json.dumps(dict(
        target_step=target_step, reached_step=reached_step,
        complete=complete)))


def write_checkpoint(gen_p, restart_name=GMX_RESTART_NAME):
    """A file carrying the GROMACS checkpoint magic, which is all is_checkpoint
    looks at."""
    (gen_p / restart_name).write_bytes(
        gmx_simulate.CHECKPOINT_MAGIC + b'\x00' * 16)


def write_state_xml(gen_p, step_count, restart_name=OMM_RESTART_NAME):
    """A real serialized OpenMM State stopped at step_count."""
    import openmm as mm
    from openmm import unit
    system = mm.System()
    system.addParticle(1.0)
    context = mm.Context(system, mm.VerletIntegrator(0.002 * unit.picoseconds),
                         mm.Platform.getPlatformByName('Reference'))
    context.setPositions([mm.Vec3(0, 0, 0) * unit.nanometer])
    context.setStepCount(step_count)
    (gen_p / restart_name).write_text(mm.XmlSerializer.serialize(
        context.getState(getPositions=True, getVelocities=True)))


def finished_gen(top_level, clone_index, gen_index, structure):
    """An earlier generation that was harvested, so a chain has somewhere to
    start and a sibling harvester config exists to borrow."""
    gen_p = gen_dir(top_level, clone_index, gen_index)
    write_config(gen_p, top_level, clone_index, gen_index, structure)
    write_hconfig(gen_p, structure)
    (gen_p / hv.SENTINEL_NAME).write_text(json.dumps(dict(
        status='harvested', gen_index=gen_index,
        frames_per_gen=FRAMES_PER_GEN)))
    return gen_p


# ---------------------------------------------------------------------------
# One generation per category
# ---------------------------------------------------------------------------

def case_harvested(top_level, clone_index, structure):
    """A generation with its sentinel; nothing to recover."""
    gen_p = finished_gen(top_level, clone_index, 0, structure)
    (gen_p / hv.SENTINEL_NAME).write_text(json.dumps(dict(
        status='harvested', gen_index=0, frames_per_gen=FRAMES_PER_GEN)))
    return gen_p


def case_gromacs_complete(top_level, clone_index, structure):
    """Generation 1 of a chain: complete, checkpointed, seam frame written."""
    finished_gen(top_level, clone_index, 0, structure)
    gen_p = gen_dir(top_level, clone_index, 1)
    write_config(gen_p, top_level, clone_index, 1, structure)
    write_hconfig(gen_p, structure)
    write_traj(gen_p / f'{TRAJ_NAME}{TRAJ_SUFFIX}', FRAMES_PER_GEN + 1,
               first_step=STEPS_PER_GEN)
    write_status(gen_p, target_step=2 * STEPS_PER_GEN,
                 reached_step=2 * STEPS_PER_GEN, complete=True)
    write_checkpoint(gen_p)
    return gen_p


def case_openmm_complete(top_level, clone_index, structure):
    """The same, on the other engine, and with no harvester config of its own."""
    finished_gen(top_level, clone_index, 0, structure)
    gen_p = gen_dir(top_level, clone_index, 1)
    write_config(gen_p, top_level, clone_index, 1, structure,
                 restart_name=OMM_RESTART_NAME)
    write_traj(gen_p / f'{TRAJ_NAME}{TRAJ_SUFFIX}', FRAMES_PER_GEN,
               first_step=STEPS_PER_GEN + WRITE_INTERVAL)
    write_state_xml(gen_p, 2 * STEPS_PER_GEN)
    return gen_p


def case_repairable(top_level, clone_index, structure):
    """A harvest that swapped the original for a symlink and then died."""
    gen_p = gen_dir(top_level, clone_index, 0)
    write_config(gen_p, top_level, clone_index, 0, structure)
    write_hconfig(gen_p, structure)
    dry_name = f'{hv.DRY_PREFIX}{SEP}{TRAJ_NAME}{TRAJ_SUFFIX}'
    write_traj(gen_p / dry_name, FRAMES_PER_GEN + 1)
    write_traj(gen_p / f'{hv.DOWNSAMPLE_PREFIX}{SEP}{TRAJ_NAME}{TRAJ_SUFFIX}', 3)
    write_structure(gen_p / hv.DRY_TOPOLOGY_NAME)
    (gen_p / f'{TRAJ_NAME}{TRAJ_SUFFIX}').symlink_to(dry_name)
    return gen_p


def case_partly_harvested(top_level, clone_index, structure):
    """A harvest that died while writing, before it swapped anything."""
    gen_p = gen_dir(top_level, clone_index, 0)
    write_config(gen_p, top_level, clone_index, 0, structure)
    write_hconfig(gen_p, structure)
    write_traj(gen_p / f'{TRAJ_NAME}{TRAJ_SUFFIX}', FRAMES_PER_GEN)
    write_traj(gen_p / f'{hv.DRY_PREFIX}{SEP}{TRAJ_NAME}{TRAJ_SUFFIX}', 4)
    write_checkpoint(gen_p)
    write_status(gen_p, target_step=STEPS_PER_GEN,
                 reached_step=STEPS_PER_GEN, complete=True)
    return gen_p


def case_gromacs_unfinished(top_level, clone_index, structure):
    """A full trajectory whose status record still owes steps."""
    gen_p = gen_dir(top_level, clone_index, 0)
    write_config(gen_p, top_level, clone_index, 0, structure)
    write_hconfig(gen_p, structure)
    write_traj(gen_p / f'{TRAJ_NAME}{TRAJ_SUFFIX}', FRAMES_PER_GEN)
    write_status(gen_p, target_step=STEPS_PER_GEN, reached_step=500,
                 complete=False)
    write_checkpoint(gen_p)
    return gen_p


def case_openmm_unfinished(top_level, clone_index, structure):
    """A full trajectory whose checkpoint is behind it: the tail is a branch
    the next generation rewound past."""
    gen_p = gen_dir(top_level, clone_index, 0)
    write_config(gen_p, top_level, clone_index, 0, structure,
                 restart_name=OMM_RESTART_NAME)
    write_hconfig(gen_p, structure)
    write_traj(gen_p / f'{TRAJ_NAME}{TRAJ_SUFFIX}', FRAMES_PER_GEN)
    write_state_xml(gen_p, 500)
    return gen_p


def case_frame_count(top_level, clone_index, structure):
    """A trajectory holding a number of frames no plan for this length allows."""
    gen_p = gen_dir(top_level, clone_index, 0)
    write_config(gen_p, top_level, clone_index, 0, structure)
    write_hconfig(gen_p, structure)
    write_traj(gen_p / f'{TRAJ_NAME}{TRAJ_SUFFIX}', FRAMES_PER_GEN - 3)
    write_status(gen_p, target_step=STEPS_PER_GEN,
                 reached_step=STEPS_PER_GEN, complete=True)
    write_checkpoint(gen_p)
    return gen_p


def case_off_chain_target(top_level, clone_index, structure):
    """Complete, but against a step target the chain does not agree with."""
    gen_p = gen_dir(top_level, clone_index, 0)
    write_config(gen_p, top_level, clone_index, 0, structure)
    write_hconfig(gen_p, structure)
    write_traj(gen_p / f'{TRAJ_NAME}{TRAJ_SUFFIX}', FRAMES_PER_GEN)
    write_status(gen_p, target_step=999999, reached_step=999999, complete=True)
    write_checkpoint(gen_p)
    return gen_p


def case_no_witness(top_level, clone_index, structure):
    """A trajectory of the right length and no record of how it got there."""
    gen_p = gen_dir(top_level, clone_index, 0)
    write_config(gen_p, top_level, clone_index, 0, structure)
    write_hconfig(gen_p, structure)
    write_traj(gen_p / f'{TRAJ_NAME}{TRAJ_SUFFIX}', FRAMES_PER_GEN)
    return gen_p


def case_no_steps_per_gen(top_level, clone_index, structure):
    """Neither config records the generation length, so it would be a guess."""
    gen_p = gen_dir(top_level, clone_index, 0)
    write_config(gen_p, top_level, clone_index, 0, structure,
                 steps_per_gen=None)
    write_hconfig(gen_p, structure, steps_per_gen=None)
    write_traj(gen_p / f'{TRAJ_NAME}{TRAJ_SUFFIX}', FRAMES_PER_GEN)
    write_status(gen_p, target_step=STEPS_PER_GEN,
                 reached_step=STEPS_PER_GEN, complete=True)
    write_checkpoint(gen_p)
    return gen_p


def case_bad_config(top_level, clone_index, structure):
    """A run record that will not parse."""
    gen_p = gen_dir(top_level, clone_index, 0)
    (gen_p / hv.CONFIG_NAME).write_text('{"gen_index": 0,')
    write_hconfig(gen_p, structure)
    write_traj(gen_p / f'{TRAJ_NAME}{TRAJ_SUFFIX}', FRAMES_PER_GEN)
    return gen_p


def case_no_hconfig(top_level, clone_index, structure):
    """A generation whose clone holds no harvester config to borrow."""
    gen_p = gen_dir(top_level, clone_index, 0)
    write_config(gen_p, top_level, clone_index, 0, structure)
    write_traj(gen_p / f'{TRAJ_NAME}{TRAJ_SUFFIX}', FRAMES_PER_GEN)
    write_status(gen_p, target_step=STEPS_PER_GEN,
                 reached_step=STEPS_PER_GEN, complete=True)
    write_checkpoint(gen_p)
    return gen_p


# label -> (builder, category the classifier must reach)
CASES = (
    ('already harvested', case_harvested, hv.CATEGORY_HARVESTED),
    ('gromacs, finished mid-chain', case_gromacs_complete, hv.CATEGORY_COMPLETE),
    ('openmm, finished mid-chain', case_openmm_complete, hv.CATEGORY_COMPLETE),
    ('harvest died after the swap', case_repairable, hv.CATEGORY_REPAIRABLE),
    ('harvest died before the swap', case_partly_harvested,
     hv.CATEGORY_PARTLY_HARVESTED),
    ('gromacs, still owes steps', case_gromacs_unfinished,
     hv.CATEGORY_UNFINISHED),
    ('openmm, checkpoint behind the trajectory', case_openmm_unfinished,
     hv.CATEGORY_UNFINISHED),
    ('frames fit no plan', case_frame_count, hv.CATEGORY_INCONSISTENT),
    ('complete against the wrong target', case_off_chain_target,
     hv.CATEGORY_INCONSISTENT),
    ('no engine record at all', case_no_witness, hv.CATEGORY_UNPROVEN),
    ('generation length unrecorded', case_no_steps_per_gen,
     hv.CATEGORY_UNPROVEN),
    ('run record will not parse', case_bad_config, hv.CATEGORY_UNREADABLE),
    ('no harvester config anywhere', case_no_hconfig, hv.CATEGORY_UNREADABLE),
)


def build_campaign(top_level, structure, cases=CASES):
    """One clone per case. Returns label -> (gen_dir, expected category)."""
    return {label: (build(top_level, clone_index, structure), expected)
            for clone_index, (label, build, expected) in enumerate(cases)}


def tree(top_level):
    """Every path under a campaign, for proving the report changed nothing."""
    return sorted(str(p) for p in Path(top_level).rglob('*'))


def forced_harvest(source, dest_root):
    """Harvest a copy of a generation as if the classifier had cleared it.

    Returns the exception the harvest raised, or None when it went through --
    which is how a case the harvest refuses on its own is told from one the
    classifier is the only guard against.
    """
    copy = Path(dest_root) / Path(source).name
    shutil.rmtree(copy, ignore_errors=True)
    copy.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(source, copy)
    row = dict(hv.classify_gen_dir(source), gen_dir=str(copy),
               category=hv.CATEGORY_COMPLETE)
    try:
        hv.harvest_recovered([row])
    except Exception as exc:
        return exc
    return None


def row_for(rows, gen_p):
    """The classification row for one generation directory."""
    for row in rows:
        if Path(row['gen_dir']) == Path(gen_p):
            return row
    return None


def main():
    suite = Suite('harvest_recovery')
    work = harness.workdir('harvest_recovery')
    structure = write_structure(work / STRUCTURE_NAME)
    top = work / 'trajectories'
    built = build_campaign(top, structure)

    suite.section('every category is reached, and only by its own generation')
    before = tree(top)
    rows = hv.classify_campaign(top)
    reported = hv.format_report(rows)
    print(reported, flush=True)
    suite.check('the report modifies nothing', tree(top) == before)
    suite.check('the harvested generation is not in the sweep',
                row_for(rows, built['already harvested'][0]) is None)
    for label, (gen_p, expected) in built.items():
        row = row_for(rows, gen_p)
        if row is None:
            direct = hv.classify_gen_dir(gen_p)
            suite.check(f'{label} -> {expected}', direct['category'] == expected,
                        f"-> {direct['category']}: {direct['reason'][:60]}")
            continue
        suite.check(f'{label} -> {expected}', row['category'] == expected,
                    f"-> {row['category']}: {row['reason'][:60]}")

    suite.section('the evidence behind the two safe verdicts')
    gmx_row = row_for(rows, built['gromacs, finished mid-chain'][0])
    suite.check('the seam frame is seen and dropped',
                gmx_row['n_orig'] == FRAMES_PER_GEN + 1 and gmx_row['skip_first'],
                f"-> {gmx_row['n_orig']} frames, skip_first={gmx_row['skip_first']}")
    suite.check('the chain places it at the right global frame and step',
                (gmx_row['first_global_index'], gmx_row['target_step'])
                == (FRAMES_PER_GEN, 2 * STEPS_PER_GEN),
                f"-> {gmx_row['first_global_index']}, {gmx_row['target_step']}")
    omm_row = row_for(rows, built['openmm, finished mid-chain'][0])
    suite.check('the openmm generation is judged on its state.xml',
                omm_row['witness'] == OMM_RESTART_NAME,
                f"-> {omm_row['witness']}")
    suite.check('a generation with no harvester config borrows a sibling one',
                omm_row['hconfig_source'].endswith(hv.HARVESTER_CONFIG_NAME)
                and omm_row['hconfig_source'] != hv.HARVESTER_CONFIG_NAME,
                f"-> {omm_row['hconfig_source']}")

    suite.section('acting touches the safe generations and nothing else')
    unsafe = {label: gen_p for label, (gen_p, expected) in built.items()
              if expected not in hv.SAFE_CATEGORIES
              and expected != hv.CATEGORY_HARVESTED}
    unsafe_before = {label: tree(gen_p) for label, gen_p in unsafe.items()}
    records = hv.harvest_recovered(rows)
    suite.check('one result per safe generation', len(records) == 3,
                f'-> {len(records)}')
    for label in ('gromacs, finished mid-chain', 'openmm, finished mid-chain'):
        gen_p = built[label][0]
        record = [r for r in records if r.get('n_dry') and r['status'] == 'harvested'
                  and Path(gen_p, r['dry']).is_file()]
        suite.check(f'{label} is harvested', (gen_p / hv.SENTINEL_NAME).is_file()
                    and bool(record))
        suite.check(f'{label} keeps its original',
                    (gen_p / f'{TRAJ_NAME}{TRAJ_SUFFIX}').is_file()
                    and not (gen_p / f'{TRAJ_NAME}{TRAJ_SUFFIX}').is_symlink())
        suite.check(f'{label} writes the frames the plan predicted',
                    record[0]['n_dry'] == GEN1_DRY
                    and record[0]['n_down'] == GEN1_DOWN,
                    f"-> {record[0]['n_dry']}, {record[0]['n_down']}")
        suite.check(f'{label} records that it did not unlink',
                    record[0]['unlinked'] is False)
    repaired = built['harvest died after the swap'][0]
    suite.check('the interrupted harvest gets only its sentinel',
                (repaired / hv.SENTINEL_NAME).is_file()
                and json.loads((repaired / hv.SENTINEL_NAME).read_text())
                ['status'] == 'repaired')
    for label, gen_p in unsafe.items():
        suite.check(f'{label} is left exactly as it was',
                    tree(gen_p) == unsafe_before[label])

    suite.section('the guard is what stops an unfinished generation')
    # Its trajectory is full, so nothing downstream of the classifier objects:
    # widen the safe set and the same generation is harvested.
    mutant = work / 'mutant'
    guarded = case_openmm_unfinished(mutant, 0, structure)
    guarded_before = tree(guarded)
    hv.harvest_recovered(hv.classify_campaign(mutant))
    suite.check('the default safe set refuses it', tree(guarded) == guarded_before)
    loosened = work / 'loosened'
    widened = case_openmm_unfinished(loosened, 0, structure)
    hv.harvest_recovered(
        hv.classify_campaign(loosened),
        safe_categories=hv.SAFE_CATEGORIES + (hv.CATEGORY_UNFINISHED,))
    suite.check('widening the safe set does harvest it, so the guard is load '
                'bearing', (widened / hv.SENTINEL_NAME).is_file())

    suite.section('what the harvest catches on its own, and what it cannot')
    short = forced_harvest(built['frames fit no plan'][0], work / 'forced-short')
    suite.check('a frame count fitting no plan is refused by the harvest too',
                isinstance(short, hv.HarvestError), f'-> {short}')
    off_chain = forced_harvest(built['complete against the wrong target'][0],
                               work / 'forced-off-chain')
    suite.check('a step target off the chain is invisible to the harvest, so '
                'the classifier is the only guard', off_chain is None,
                f'-> {off_chain}')

    suite.section('unlink is honoured when it is asked for')
    asked = work / 'asked'
    gen_p = case_gromacs_complete(asked, 0, structure)
    hv.harvest_recovered(hv.classify_campaign(asked), unlink=True)
    traj_p = gen_p / f'{TRAJ_NAME}{TRAJ_SUFFIX}'
    suite.check('the original becomes a symlink to the dry copy',
                traj_p.is_symlink()
                and traj_p.resolve().name.startswith(hv.DRY_PREFIX),
                f'-> {traj_p.is_symlink()}')
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
