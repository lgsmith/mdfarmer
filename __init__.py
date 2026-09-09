from . import harvester
from .harvester import (Harvester, HarvestError, harvest_generation,
                        check_commensurability, unharvested_gen_dirs,
                        verify_dry_chain, select_backend, resolve_subset,
                        frame_plan, expected_counts,
                        classify_gen_dir, classify_campaign, format_report,
                        harvest_recovered, completion_witness)
from .farmer import *
from . import utilities
from .utilities import *
# The OpenMM runner is imported only where OpenMM is: a GROMACS campaign drives
# gmx as a subprocess, so a GROMACS-only site must not need OpenMM to import
# this package. Absent, the two OpenMM entry points are still exported, as
# stubs that say what to install when they are called.
if utilities.openmm_available():
    from .simulate import *
else:
    omm_generation = utilities.missing_openmm('omm_generation')
    omm_basic_sim_block_json = utilities.missing_openmm(
        'omm_basic_sim_block_json')
from .seeder import *
# GROMACS runner imported by name (not *) so its module-local Preempted does
# not shadow simulate.Preempted in the package namespace.
from . import gmx_simulate
from .gmx_simulate import (gmx_generation, gmx_basic_sim_block_json,
                           gmx_try_recover_gen, gmx_gen_progress,
                           gmx_config_template, default_gmx_run_script)
# Reimaging is imported as a module, not splatted: it defines short names like
# box_vectors and bond_pairs that would collide unhelpfully in the package
# namespace, and callers read better as reimage.reimage_trajectory(...).
from . import reimage
# MPS packing: K replicas sharing one GPU inside a single job.
from . import gmx_pack
from .gmx_pack import (gmx_pack_sim_block_json, replica_mdrun_args,
                       write_pack_manifest, default_gmx_pack_run_script)
from .seeder import ClonePack

__all__ = ['fdir', 'dir_seeds_clones', 'dir_seeds_clones_gens', 'omm_generation',
           'omm_basic_sim_block_json', 'Clone', 'Farmer', 'basic_scheduler_reports',
           'basic_scheduler_fstrings', 'basic_scheduler_fstrings_preempt',
           'basic_scheduler_assoc_reports', 'basic_gpu_lines',
           'BadNodeRegistry', 'default_bad_node_patterns',
           # Used by the README's own example, so they have to be exported.
           'Harvester', 'default_harvest_shellscript',
           'default_harvest_shellscript_slurm',
           'harvester', 'HarvestError', 'harvest_generation',
           'check_commensurability', 'unharvested_gen_dirs',
           'verify_dry_chain', 'select_backend', 'resolve_subset',
           'frame_plan', 'expected_counts',
           # How a finished campaign is inspected and repaired, which the
           # README walks through and unharvested_gen_dirs alone cannot do.
           'classify_gen_dir', 'classify_campaign', 'format_report',
           'harvest_recovered', 'completion_witness',
           'strip_and_downsample',
           'default_straight_sampling_config_template',
           'default_straight_sampling_init_config',
           'merge_args_defaults_dict', 'calx_remaining_steps', 'get_traj_len',
           'frame_timing',
           'gmx_simulate', 'gmx_generation', 'gmx_basic_sim_block_json',
           'gmx_try_recover_gen', 'gmx_gen_progress', 'gmx_config_template',
           'default_gmx_run_script', 'reimage',
           'gmx_pack', 'gmx_pack_sim_block_json', 'replica_mdrun_args',
           'write_pack_manifest', 'default_gmx_pack_run_script',
           'ClonePack', 'basic_scheduler_fstrings_mps',
           'missing_seed_inputs', 'ready_seed_count', 'check_seed_map']


def __getattr__(name):
    """Names that exist only with OpenMM installed, refused by name.

    mdfarmer.simulate is the OpenMM runner module and is not imported without
    OpenMM, so asking for it raises an ImportError naming the missing package
    rather than a bare AttributeError about the package.
    """
    if name == 'simulate':
        utilities.require_openmm('mdfarmer.simulate')
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')
