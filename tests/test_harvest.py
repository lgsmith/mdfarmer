"""Suite for mdfarmer.harvester: drives real short GROMACS generations of the
water box, harvests them, and checks the guarantees harvest_generation makes.
"""
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import mdtraj as md

import harness

import mdfarmer
from mdfarmer import harvester
from mdfarmer import gmx_simulate

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

N_GENS = 3
STEPS_PER_GEN = 1000
WRITE_INTERVAL = 100
DOWNSAMPLE_FRQ = 5
HARVESTER_SELECTION = 'resid <= 20'          # pretend solute
SEED_INDEX = 0
CLONE_INDEX = 0
DIRNAME_PAD = 2
SEP = '-'
GEN_TEMPERATURE = 300.0
TITLE = 'harvest-test'
CONFIG_NAME = 'config.json'
HCONFIG_NAME = 'hconfig.json'
RESTART_NAME = 'state.cpt'
GROMPP_MAXWARN = 5
GMX_BIN = harness.GMX_BIN
MDRUN_ARGS = ()                              # let mdrun auto-pick hardware
REPAIR_TARGET_GEN = 1                        # which chain gen to break/repair

# Round-trip through a freshly written XTC loses a bit of precision to
# compression; allow a little more than one compression step's worth.
COORD_ATOL_NM = 2e-3


# ---------------------------------------------------------------------------
# Driving real generations
# ---------------------------------------------------------------------------

def run_generations(top_level, structure, topology, mdp, n_gens,
                    seed_index=SEED_INDEX, clone_index=CLONE_INDEX,
                    dirname_pad=DIRNAME_PAD, sep=SEP,
                    steps_per_gen=STEPS_PER_GEN, write_interval=WRITE_INTERVAL,
                    downsample_frq=DOWNSAMPLE_FRQ,
                    selection=HARVESTER_SELECTION,
                    temperature=GEN_TEMPERATURE, gmx_bin=GMX_BIN,
                    mdrun_args=MDRUN_ARGS, grompp_maxwarn=GROMPP_MAXWARN,
                    title=TITLE, restart_name=RESTART_NAME,
                    config_name=CONFIG_NAME, hconfig_name=HCONFIG_NAME):
    """Run `n_gens` chained GROMACS generations of the water box.

    Mirrors what `Clone` does around `gmx_generation`: generation 0 is
    grompp'd from the base .mdp/.gro, and every later generation has the
    previous one's checkpoint copied in under `restart_name` before it runs,
    so convert-tpr + mdrun -cpi continue it exactly. Writes config.json and
    hconfig.json into each generation directory. Returns the list of
    generation directories, in order.
    """
    gen_dirs = []
    for gen_index in range(n_gens):
        gen_dir = mdfarmer.dir_seeds_clones_gens(
            Path(top_level), seed_index, clone_index, gen_index, dirname_pad,
            sep=sep)
        prev_state_cpt = None
        if gen_index > 0:
            shutil.copy(gen_dirs[-1] / restart_name, gen_dir / restart_name)
            prev_state_cpt = gen_dir / restart_name
        config = gmx_simulate.gmx_config_template(
            traj_dir_top_level=str(top_level),
            top_fn=str(topology),
            seed_index=seed_index,
            clone_index=clone_index,
            gen_index=gen_index,
            title=title,
            seed_fn=(str(structure) if gen_index == 0
                     else str(prev_state_cpt)),
            structure_fn=(str(structure) if gen_index == 0 else None),
            mdp_fn=(str(mdp) if gen_index == 0 else None),
            steps=steps_per_gen,
            steps_per_gen=steps_per_gen,
            target_step=harness.target_step(gen_index, steps_per_gen),
            write_interval=write_interval,
            temperature=(temperature if gen_index == 0 else None),
            new_velocities=(gen_index == 0),
            gmx_bin=gmx_bin,
            mdrun_args=mdrun_args,
            grompp_maxwarn=grompp_maxwarn,
            dirname_pad=dirname_pad,
            sep=sep,
            restart_name=restart_name)
        (gen_dir / config_name).write_text(json.dumps(config))
        gmx_simulate.gmx_generation(**config)
        hconfig = dict(harvester_subset=selection,
                      downsample_frq=downsample_frq,
                      harvester_structure=str(structure),
                      steps_per_gen=steps_per_gen)
        (gen_dir / hconfig_name).write_text(json.dumps(hconfig))
        gen_dirs.append(gen_dir)
    return gen_dirs


def copy_gen(src_dir, dest_dir):
    """A fresh, independent copy of a harvested-or-not generation directory."""
    shutil.rmtree(dest_dir, ignore_errors=True)
    shutil.copytree(src_dir, dest_dir)
    return dest_dir


def traj_path(gen_dir, config_name=CONFIG_NAME):
    """The trajectory name harvest_generation would compute for this gen."""
    config = json.loads((Path(gen_dir) / config_name).read_text())
    return Path(gen_dir) / f"{config['traj_name']}{config['traj_suffix']}"


def read_xtc(path):
    """(xyz, time, step, box) for every frame, without holding a topology."""
    with md.open(str(path)) as fh:
        return fh.read()


def catch(fn, *args, **kwargs):
    """The exception `fn` raises, or None if it returns normally."""
    try:
        fn(*args, **kwargs)
    except Exception as exc:
        return exc
    return None


def fake_loos_undercount(traj_fn, structure_fn, subset_spec, dry_out, down_out,
                         first_global_index, downsample_frq, skip_first,
                         timing=None,
                         dry_topology_name=harvester.DRY_TOPOLOGY_NAME):
    """Like `_harvest_loos`, but reports one dry frame fewer than it wrote."""
    n_orig, n_dry, n_down = harvester._harvest_loos(
        traj_fn, structure_fn, subset_spec, dry_out, down_out,
        first_global_index, downsample_frq, skip_first, timing=timing,
        dry_topology_name=dry_topology_name)
    return n_orig, n_dry - 1, n_down


# ---------------------------------------------------------------------------
# Suite
# ---------------------------------------------------------------------------

def main():
    suite = harness.Suite('harvester')
    work = harness.workdir('harvest')
    structure, topology, mdp = harness.build_water_system(work)

    suite.section('single generation (template)')
    single_gen, = run_generations(work / 'single', structure, topology, mdp,
                                  n_gens=1)
    source_xyz, source_time, source_step, source_box = read_xtc(
        single_gen / 'prod.xtc')

    # -- 1. the .top blocker -------------------------------------------------
    blocker_dir = copy_gen(single_gen, work / 'top-blocker')
    hconfig_p = blocker_dir / HCONFIG_NAME
    hconfig = json.loads(hconfig_p.read_text())
    del hconfig['harvester_structure']
    hconfig_p.write_text(json.dumps(hconfig))
    blocker_traj = traj_path(blocker_dir)
    exc = catch(harvester.harvest_generation, blocker_dir / CONFIG_NAME,
               hconfig_p)
    suite.check('harvesting with only a .top structure raises', exc is not None,
               detail=repr(exc))
    suite.check('original trajectory intact after the .top blocker',
               blocker_traj.is_file() and not blocker_traj.is_symlink())

    # -- 2. truncated-harvest guard ------------------------------------------
    trunc_dir = copy_gen(single_gen, work / 'truncated')
    trunc_traj = traj_path(trunc_dir)
    trunc_size_before = trunc_traj.stat().st_size
    trunc_backends = {harvester.BACKEND_LOOS: fake_loos_undercount,
                      harvester.BACKEND_MDTRAJ: harvester._harvest_mdtraj}
    exc = catch(harvester.harvest_generation, trunc_dir / CONFIG_NAME,
               trunc_dir / HCONFIG_NAME, backend=harvester.BACKEND_LOOS,
               backends=trunc_backends)
    suite.check('an undercounting backend raises HarvestError',
               isinstance(exc, harvester.HarvestError), detail=repr(exc))
    suite.check('original trajectory intact after the truncated-harvest guard',
               trunc_traj.is_file() and not trunc_traj.is_symlink()
               and trunc_traj.stat().st_size == trunc_size_before)
    suite.check('no sentinel written after the truncated-harvest guard',
               not (trunc_dir / harvester.SENTINEL_NAME).is_file())

    # -- 3. backend agreement; 6. time axis; 7. box --------------------------
    mdtraj_dir = copy_gen(single_gen, work / 'backend-mdtraj')
    auto_dir = copy_gen(single_gen, work / 'backend-auto')
    record_mdtraj = harvester.harvest_generation(
        mdtraj_dir / CONFIG_NAME, mdtraj_dir / HCONFIG_NAME,
        backend=harvester.BACKEND_MDTRAJ)
    record_auto = harvester.harvest_generation(
        auto_dir / CONFIG_NAME, auto_dir / HCONFIG_NAME)

    suite.check('auto backend resolves to loos for this rectangular box',
               record_auto['backend'] == harvester.BACKEND_LOOS,
               detail=record_auto['backend'])

    xyz_m, time_m, step_m, box_m = read_xtc(mdtraj_dir / record_mdtraj['dry'])
    xyz_a, time_a, step_a, box_a = read_xtc(auto_dir / record_auto['dry'])
    suite.check('backends agree on dry frame count',
               xyz_m.shape[0] == xyz_a.shape[0],
               detail=f'{xyz_m.shape[0]} vs {xyz_a.shape[0]}')
    suite.check('backends agree on dry atom count',
               xyz_m.shape[1] == xyz_a.shape[1],
               detail=f'{xyz_m.shape[1]} vs {xyz_a.shape[1]}')
    suite.check('backends agree on dry frame times',
               time_m.shape == time_a.shape and np.allclose(time_m, time_a))
    suite.check('backends agree on dry coordinates',
               xyz_m.shape == xyz_a.shape
               and np.allclose(xyz_m, xyz_a, atol=COORD_ATOL_NM))

    for label, gen_dir, record in (('mdtraj', mdtraj_dir, record_mdtraj),
                                   ('loos', auto_dir, record_auto)):
        xyz, time, step, box = read_xtc(gen_dir / record['dry'])
        suite.check(f'{label} backend: dry frame times match the source',
                   time.shape == source_time.shape
                   and np.allclose(time, source_time))
        suite.check(f'{label} backend: dry MD steps match the source',
                   step.shape == source_step.shape
                   and np.array_equal(step, source_step))

    suite.check('dry frames carry a non-zero unit cell',
               bool(np.all(np.linalg.det(box_a) > 0)))
    suite.check('dry frames carry the source unit cell',
               box_a.shape == source_box.shape
               and np.allclose(box_a, source_box))

    # -- 8. idempotency -------------------------------------------------------
    idem_dry_p = auto_dir / record_auto['dry']
    idem_size_before = idem_dry_p.stat().st_size
    idem_record = harvester.harvest_generation(auto_dir / CONFIG_NAME,
                                               auto_dir / HCONFIG_NAME)
    suite.check('re-harvesting reports already-harvested',
               idem_record['status'] == 'already-harvested')
    suite.check('re-harvesting leaves the dry file size unchanged',
               idem_dry_p.stat().st_size == idem_size_before)

    # -- generation chain: 4. seams, 5. wet spacing, 9. repair, 10. sweep -----
    suite.section('generation chain')
    chain_top = work / 'chain'
    gen_dirs = run_generations(chain_top, structure, topology, mdp,
                               n_gens=N_GENS)
    chain_records = [harvester.harvest_generation(d / CONFIG_NAME,
                                                   d / HCONFIG_NAME)
                     for d in gen_dirs]

    chain_check = harvester.verify_dry_chain(gen_dirs)
    suite.check('dry chain is contiguous across seams',
               chain_check['contiguous'], detail=str(chain_check))
    suite.check('dry chain time is strictly increasing',
               chain_check['strictly_increasing'])
    suite.check('dry chain has no duplicate times',
               chain_check['n_duplicate_times'] == 0)
    suite.check('dry chain is uniformly spaced',
               chain_check['uniform_spacing'])

    down_times = [read_xtc(d / rec['downsample'])[1]
                 for d, rec in zip(gen_dirs, chain_records)]
    down_time = np.concatenate(down_times)
    down_spacing = np.diff(down_time)
    suite.check('downsampled stream has no duplicate times across seams',
               bool((down_spacing > 0).all()))
    suite.check('downsampled stream is uniformly spaced across seams',
               len(down_spacing) > 0
               and np.allclose(down_spacing, down_spacing[0]),
               detail=str(np.unique(np.round(down_spacing, 6))))

    target_dir = gen_dirs[REPAIR_TARGET_GEN]
    target_record = chain_records[REPAIR_TARGET_GEN]
    (target_dir / harvester.SENTINEL_NAME).unlink()

    stale = [Path(p).resolve() for p in
            harvester.unharvested_gen_dirs(chain_top)]
    suite.check('unharvested_gen_dirs finds the gen with a removed sentinel',
               target_dir.resolve() in stale, detail=str(stale))
    suite.check('unharvested_gen_dirs does not flag still-harvested gens',
               all(d.resolve() not in stale
                   for i, d in enumerate(gen_dirs) if i != REPAIR_TARGET_GEN))

    repaired = harvester.harvest_generation(target_dir / CONFIG_NAME,
                                            target_dir / HCONFIG_NAME)
    suite.check('repair after a missing sentinel reports repaired',
               repaired['status'] == 'repaired', detail=str(repaired))
    suite.check('repair reproduces the original frame counts',
               repaired['n_dry'] == target_record['n_dry']
               and repaired['n_down'] == target_record['n_down'])

    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
