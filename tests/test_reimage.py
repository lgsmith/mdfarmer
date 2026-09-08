"""Reimaging keeps the clock, and reads each format's own length units."""
import shutil
import sys

import numpy as np

import harness
from harness import Suite

from mdfarmer import reimage

STEPS = 400
WRITE_INTERVAL = 100
ANGSTROM_PER_NM = 10.0


def _refuses(call, *args):
    try:
        call(*args)
    except ValueError:
        return True
    return False


def main(steps=STEPS, write_interval=WRITE_INTERVAL,
         angstrom_per_nm=ANGSTROM_PER_NM, gmx_bin=harness.GMX_BIN):
    import mdtraj as md
    from mdfarmer import gmx_simulate as gs

    suite = Suite('reimage')
    work = harness.workdir('reimage')
    structure, topology, mdp = harness.build_water_system(work, gmx_bin=gmx_bin)

    traj = gs.gmx_generation(
        traj_dir_top_level=str(work / 'farm'), top_fn=str(topology),
        seed_index=0, clone_index=0, gen_index=0, title='ri',
        seed_fn=str(structure), structure_fn=str(structure), mdp_fn=str(mdp),
        dirname_pad=2, sep='-', traj_name='prod', traj_suffix='.xtc',
        restart_name='state.cpt', steps=steps, steps_per_gen=steps,
        target_step=harness.target_step(0, steps),
        write_interval=write_interval, temperature=300, gen_seed_base=7,
        gmx_bin=gmx_bin, grompp_maxwarn=3, new_velocities=True, append=False,
        mdrun_args=('-nb', 'cpu', '-pme', 'cpu', '-ntomp', '2'))

    source_time, source_step = (np.asarray(a) for a in
                                md.open(str(traj)).read()[1:3])

    # The force field lives in the GROMACS install, not beside the .top, and
    # reimage_with_loos has no include_dir of its own; resolve them here.
    gmx_top = harness.require_gmx(gmx_bin=gmx_bin)
    ranges = reimage.molecule_ranges(str(topology), include_dir=str(gmx_top))

    suite.section('the reimaged trajectory keeps the source clock')
    out, n_written = reimage.reimage_with_loos(
        traj, structure_fn=str(structure), out_fn=str(work / 'whole.xtc'),
        ranges=ranges, verify=False)
    whole_time, whole_step = (np.asarray(a) for a in md.open(str(out)).read()[1:3])
    print(f'   source time {source_time[:4]} step {source_step[:4]}', flush=True)
    print(f'   whole  time {whole_time[:4]} step {whole_step[:4]}', flush=True)
    suite.check('every frame is written', n_written == len(source_time),
                f'-> {n_written}')
    suite.check('frame times survive reimaging',
                np.allclose(whole_time, source_time, atol=1e-5))
    suite.check('MD step numbers survive reimaging',
                np.array_equal(whole_step, source_step))

    suite.section('length units follow the trajectory format')
    suite.check('an .xtc is already nanometres',
                reimage.length_scale('prod.xtc') == 1.0)
    suite.check('a .dcd is converted from Angstroms',
                np.isclose(reimage.length_scale('prod.dcd'), 1 / angstrom_per_nm))

    # The same coordinates in both formats must give the same bond lengths.
    frames = md.load(str(traj), top=str(structure))
    frames.save_dcd(str(work / 'prod.dcd'))
    pairs = reimage.bond_pairs(str(topology), include_dir=str(gmx_top))
    n_xtc, _ = reimage.check_bond_lengths(traj, pairs=pairs)
    n_dcd, _ = reimage.check_bond_lengths(work / 'prod.dcd', pairs=pairs)
    print(f'   long bonds: {n_xtc} in the .xtc, {n_dcd} in the same frames '
          f'as .dcd', flush=True)
    suite.check('a .dcd is not falsely flagged as one long bond per bond',
                n_dcd == n_xtc, f'-> {n_dcd} vs {n_xtc}')
    suite.section('the topology is reachable wherever its force field lives')
    out_inc, _ = reimage.reimage_with_loos(
        traj, structure_fn=str(structure), out_fn=str(work / 'inc.xtc'),
        top_fn=str(topology), include_dir=str(gmx_top), verify=True)
    suite.check('include_dir reaches molecule_ranges and the bond check',
                out_inc.is_file())
    suite.check('gromacs_topology refuses anything but a .top',
                _refuses(reimage.gromacs_topology, str(structure)))

    suite.section('molecules come back whole, never imaged atom by atom')
    # GROMACS writes whatever the integrator holds, so its own output already
    # has molecules cut across box faces.
    n_raw, worst = reimage.check_bond_lengths(traj, pairs=pairs)
    longest = max((v[3] for v in worst), default=0.0)
    print(f'   raw GROMACS output: {n_raw} bonds over {reimage.MAX_BOND} nm, '
          f'longest {longest:.3f} nm', flush=True)
    suite.check('the raw trajectory really is broken by the boundary', n_raw > 0)
    n_whole, _ = reimage.check_bond_lengths(out, pairs=pairs)
    suite.check('reimaging leaves no bond longer than a bond can be',
                n_whole == 0, f'-> {n_whole}')

    # The worst case: shift by half a box so every molecule straddles a face,
    # then wrap each atom on its own, which is what breaks a trajectory beyond
    # repair if the molecules are ever lost.
    frames = md.load(str(traj), top=str(structure))
    box = frames.unitcell_lengths[0]
    shifted = frames.xyz + box / 2.0
    frames.xyz = shifted - box * np.floor(shifted / box)
    frames.save_xtc(str(work / 'per-atom.xtc'))
    n_split, _ = reimage.check_bond_lengths(work / 'per-atom.xtc', pairs=pairs)
    print(f'   after wrapping every atom on its own: {n_split} long bonds',
          flush=True)
    suite.check('wrapping atom by atom breaks more molecules', n_split >= n_raw)
    healed, _ = reimage.reimage_with_loos(
        work / 'per-atom.xtc', structure_fn=str(structure),
        out_fn=str(work / 'per-atom-whole.xtc'), ranges=ranges, verify=False)
    n_healed, _ = reimage.check_bond_lengths(healed, pairs=pairs)
    suite.check('reimaging puts every one of them back together',
                n_healed == 0, f'-> {n_healed}')

    suite.section('an unsafe anchor margin is refused, bonds or no bonds')
    # Declaring the whole box one molecule puts its furthest atom well past the
    # half-edge mergeImage() assumes, which is the regime LOOS cannot answer in.
    whole_system = [(0, md.load(str(structure)).n_atoms)]
    try:
        reimage.reimage_with_loos(
            traj, structure_fn=str(structure),
            out_fn=str(work / 'unsafe.xtc'), ranges=whole_system, verify=True)
        suite.check('reimaging outside the safe regime raises', False,
                    '-> no exception')
    except RuntimeError as exc:
        suite.check('reimaging outside the safe regime raises',
                    'safe regime' in str(exc), f'-> {str(exc)[:50]}')
        suite.check('the message names the backend that can do it',
                    reimage.BACKEND_MDTRAJ in str(exc))
    # Measured on the reimaged output: the raw trajectory's molecules are still
    # split, so its atoms really are a box length from their anchors.
    margin = reimage.check_anchor_distances(out, ranges)
    suite.check('whole molecules are comfortably inside it',
                margin['loos_safe'], f"-> {margin['max_anchor_offset']:.3f} nm")

    suite.section('the mdtraj backend never guesses what a molecule is')
    try:
        reimage.reimage_with_mdtraj(traj, None, str(work / 'no-top.xtc'))
        suite.check('no topology at all is reported as such', False,
                    '-> no exception')
    except ValueError as exc:
        suite.check('no topology at all is reported as such',
                    'top_fn' in str(exc), f'-> {str(exc)[:50]}')
    except TypeError as exc:
        suite.check('no topology at all is reported as such', False,
                    f'-> TypeError instead: {exc}')
    # One molecule per atom is what mdtraj falls back to on a bondless
    # topology, and it is the grouping that breaks a trajectory for good.
    per_atom = [(i, i + 1) for i in range(ranges[-1][1])]
    suite.check('a grouping that cuts through bonds is refused',
                _refuses(reimage.check_bonds_within_molecules, pairs,
                         per_atom))

    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
