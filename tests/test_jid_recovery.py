"""A job of ours the queue no longer names in a way we can read.

The tender learns which jobs are its own by reading their names, so a site that
rewrites job names, or a report that mangles one, makes a live job look dead --
and a dead-looking job gets a second one launched into its generation
directory. The jids file each tick writes is the safety net: an id recorded
there, or bound to a clone, is ours whatever the queue now calls it. A name
that parses to somebody else's title is still somebody else's.
"""
import contextlib
import io
import sys

import harness
from harness import Suite

from mdfarmer import utilities as util
from mdfarmer import gmx_simulate as gs
from mdfarmer import farmer as fm

TITLE = 'fm'
SEP = '_'
STEPS_PER_GEN = 1000
WRITE_INTERVAL = 100
# One of ours, named the way we named it.
OUR_JID = 55
# One of ours, whose queue entry no longer carries a name we can split.
MANGLED_JID = 66
# One that parses cleanly to a title that is not this campaign's.
STRANGER_JID = 77
# One still queued at boot under a name that stops one index short.
BOOT_JID = 88


def captured(call):
    """Run call(), returning (result, everything it printed)."""
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        result = call()
    return result, out.getvalue()


def report_cmd(*lines):
    """A scheduler report command printing the given '<jid> <name>' lines."""
    return "printf '{}'".format(''.join(line + r'\n' for line in lines))


class BoundStub:
    """Stands in a queue for a unit the tender has already submitted."""

    def __init__(self, job_number):
        self.job_number = job_number


def make_farmer(work, steps_per_gen=STEPS_PER_GEN,
                write_interval=WRITE_INTERVAL, title=TITLE, sep=SEP):
    template = gs.gmx_config_template(
        traj_dir_top_level=str(work / 'farm'), title=title,
        structure_fn=str(work / 'a.gro'), mdp_fn=str(work / 'base.mdp'),
        dirname_pad=2, sep=sep, traj_name='prod', traj_suffix='.xtc',
        restart_name='state.cpt', steps=steps_per_gen,
        steps_per_gen=steps_per_gen, write_interval=write_interval,
        traj_list=str(work / 'tl.txt'))
    farmer, _ = captured(lambda: fm.Farmer(
        n_seeds=1, n_clones=1, n_gens=1, config_template=template,
        seed_structure_fns=[str(work / 'a.gro')],
        system_fns=[str(work / 'base.mdp')],
        top_fns=[str(work / 'topol.top')],
        scheduler='sbatch',
        scheduler_fstring=util.basic_scheduler_fstrings['slurm'],
        scheduler_kws=dict(gpu_line='', queue_name='gpu', exclude_nodes='',
                           cpus=8, run_script_name='run.py'),
        scheduler_report_cmd='true', scheduler_assoc_rep_cmd='true',
        sep=sep, dirname_pad=2, runner=gs.gmx_generation, dry_run=True,
        overwrite=True, jids_file=work / 'jids.txt'))
    return farmer


def tick(farmer, *lines):
    """Run one update_jids over a queue report, returning (jids, log)."""
    farmer.scheduler_report_cmd = report_cmd(*lines)
    _, log = captured(farmer.update_jids)
    return farmer.current_jids, log


def main(our_jid=OUR_JID, mangled_jid=MANGLED_JID, stranger_jid=STRANGER_JID,
         boot_jid=BOOT_JID, title=TITLE, sep=SEP):
    suite = Suite('jid_recovery')
    work = harness.workdir('jid_recovery')
    for name in ('a.gro', 'topol.top', 'base.mdp'):
        (work / name).write_text('placeholder\n')
    ours = f'{title}{sep}0{sep}0{sep}0'
    theirs = f'other{sep}0{sep}0{sep}0'

    suite.section('what the campaign remembers between ticks')
    farmer = make_farmer(work)
    suite.check('a boot over an empty queue remembers nothing',
                farmer.recorded_jids() == set(),
                f'-> {farmer.recorded_jids()}')
    farmer.jids_file.write_text(f'{our_jid} {mangled_jid}')
    suite.check('the ids the last tick wrote are read back',
                farmer.recorded_jids() == {our_jid, mangled_jid},
                f'-> {farmer.recorded_jids()}')
    farmer.jids_file.write_text(f'nonsense {our_jid}')
    suite.check('a line that is not a job id is passed over',
                farmer.recorded_jids() == {our_jid},
                f'-> {farmer.recorded_jids()}')
    farmer.jids_file.unlink()
    suite.check('no file at all is not an error',
                farmer.recorded_jids() == set(),
                f'-> {farmer.recorded_jids()}')

    suite.section('a recorded id whose queue name no longer parses')
    farmer = make_farmer(work)
    farmer.jids_file.write_text(str(mangled_jid))
    jids, log = tick(farmer, f'{mangled_jid} rewritten.by.the.site')
    suite.check('it is still counted as live', jids == {mangled_jid},
                f'-> {jids}')
    suite.check('and the tick says why it kept it',
                str(mangled_jid) in log and 'NOTE' in log,
                f'-> {log.strip()[:70]}')
    suite.check('so the next tick still remembers it',
                farmer.recorded_jids() == {mangled_jid},
                f'-> {farmer.recorded_jids()}')

    suite.section('the same job with nothing to tie it to this campaign')
    farmer = make_farmer(work)
    jids, log = tick(farmer, f'{mangled_jid} rewritten.by.the.site')
    suite.check('an unrecorded, unbound id is not claimed', jids == set(),
                f'-> {jids}')

    suite.section("another campaign's job, whose id we once recorded")
    farmer = make_farmer(work)
    farmer.jids_file.write_text(str(stranger_jid))
    jids, _ = tick(farmer, f'{stranger_jid} {theirs}')
    suite.check('a name that parses to another title is left alone',
                jids == set(), f'-> {jids}')

    suite.section('a recorded id that has left the queue')
    farmer = make_farmer(work)
    farmer.jids_file.write_text(f'{our_jid} {mangled_jid}')
    jids, _ = tick(farmer, f'{our_jid} {ours}')
    suite.check('it is dropped, so its clone can start its next generation',
                jids == {our_jid}, f'-> {jids}')

    suite.section('a job this tender submitted itself, before any file')
    farmer = make_farmer(work)
    farmer.jids_file.unlink()
    farmer.priority_ordered_clones = [[BoundStub(mangled_jid)]]
    jids, _ = tick(farmer, f'{mangled_jid} rewritten.by.the.site')
    suite.check('the bound job number is enough to claim it',
                jids == {mangled_jid}, f'-> {jids}')

    suite.section('boot, with a recorded id still queued under a bad name')
    farmer = make_farmer(work)
    farmer.jids_file.write_text(str(boot_jid))
    farmer.scheduler_assoc_rep_cmd = report_cmd(
        f'{boot_jid} {title}{sep}0{sep}0')
    rep_dict, log = captured(farmer.reassociate_running_jobs)
    suite.check('no clone can be bound to it', rep_dict == {},
                f'-> {rep_dict}')
    suite.check('it is still counted as live',
                farmer.current_jids == {boot_jid},
                f'-> {farmer.current_jids}')
    suite.check('and boot warns that a second job may land on top of it',
                f'[{boot_jid}]' in log and 'cancel them by hand' in log,
                f'-> {log.strip()[-90:]}')
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
