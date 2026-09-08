"""Only this campaign's jobs may come back from a scheduler report.

Job names are title, seed, clone and gen joined by the separator, so a campaign
called 'sampling' and one called 'sampling-long' produce names that share a
prefix. A report that matches on the title alone binds the foreign job id to
this campaign's (seed, clone, gen); when the foreign job ends, the farmer reads
its own live generation as finished and submits a second job into that
generation's directory.

The reports are unfiltered queries now, and the selection happens in Python:
the trailing indices are split off the name and the title that remains is
compared by string equality. A title holding the separator or a regex
metacharacter is therefore not a special case, and no job name reaches a shell.
"""
import contextlib
import io
import shlex
import sys

import harness
from harness import Suite

from mdfarmer import utilities as util

TITLE = 'sampling'
# Two campaigns whose names share the title as a prefix, plus an unrelated job.
QUEUE_LINES = ('101 sampling-0-0-5',
               '202 sampling-long-0-0-5',
               '303 othercampaign-0-0-1',
               '404 sampling-1-2-3')
OURS = [(101, (0, 0, 5)), (404, (1, 2, 3))]
# What each report has to keep. bjobs already defaults to the invoking user's
# unfinished jobs; squeue does not, hence --me.
REQUIRED = {'slurm': ('--me', "-o '%i %j'"),
            'lsf': ("-o 'JOBID JOB_NAME'", '-noheader')}
# A squeue that cannot reach the controller, to stand in for a broken query.
DEAD_SQUEUE = ('#!/bin/sh\n'
               'echo "squeue: error: Unable to contact slurm controller" >&2\n'
               'exit 1\n')


def said(fn, *args, **kwargs):
    """Run fn with stdout captured, returning (its result, what it printed)."""
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        result = fn(*args, **kwargs)
    return result, out.getvalue()


def check_reports_are_plain(suite, required):
    """The shipped queries must ask for the whole queue and nothing else."""
    suite.section('the reports are plain, unfiltered queries')
    for family, fragments in sorted(required.items()):
        cmd = util.basic_scheduler_reports[family]
        suite.check(f'{family}: nothing is piped through a filter',
                    '|' not in cmd, f'-> {cmd}')
        suite.check(f'{family}: no title is interpolated into the shell',
                    '{' not in cmd and '}' not in cmd, f'-> {cmd}')
        for fragment in fragments:
            suite.check(f'{family}: keeps {fragment}', fragment in cmd,
                        f'-> {cmd}')
        suite.check(f'{family}: the association report asks the same thing',
                    util.basic_scheduler_assoc_reports[family] == cmd)


def check_parsing(suite, text):
    """Report text becomes (job id, job name) pairs."""
    suite.section('reading a report into (job id, job name) pairs')
    pairs, _ = said(util.parse_scheduler_report, text)
    suite.check('every line comes back as an id and a name',
                pairs == [(101, 'sampling-0-0-5'),
                          (202, 'sampling-long-0-0-5'),
                          (303, 'othercampaign-0-0-1'),
                          (404, 'sampling-1-2-3')], f'-> {pairs}')
    padded, _ = said(util.parse_scheduler_report, '  909   sampling-0-0-1  \n')
    suite.check("LSF's column padding is stripped off the name",
                padded == [(909, 'sampling-0-0-1')], f'-> {padded}')
    spaced, _ = said(util.parse_scheduler_report, '55 a name with spaces-0-0-1')
    suite.check('a name holding a space survives whole',
                spaced == [(55, 'a name with spaces-0-0-1')], f'-> {spaced}')
    junk, warned = said(util.parse_scheduler_report,
                        'slurm_load_jobs error: Socket timed out\n')
    suite.check('a line that is not a job is dropped, and said so',
                junk == [] and 'WARNING' in warned, f'-> {junk}')
    blank, warned = said(util.parse_scheduler_report, '\n\n')
    suite.check('blank lines are not worth a warning',
                blank == [] and warned == '', f'-> {warned!r}')


def check_splitting(suite):
    """A name splits into the title and its trailing indices."""
    suite.section('splitting a name into its title and its indices')
    cases = {
        'sampling-0-0-5': ('sampling', (0, 0, 5)),
        # The title leads and may hold the separator itself.
        'sampling-long-0-0-5': ('sampling-long', (0, 0, 5)),
        'a-b-c-1-22-333': ('a-b-c', (1, 22, 333)),
        # Too few fields, so there is no title left over.
        'sampling-0-0': None,
        '0-0-5': None,
        'sampling': None,
        'sampling-0-0-x': None,
        'sampling-0-0-': None,
        # int() would take these; a job name generated from ints never has them.
        'sampling-0-0-+5': None,
        'sampling-0-0- 5': None,
        'sampling-0-0-٥': None,
    }
    for name, want in cases.items():
        got = util.split_job_name(name)
        suite.check(f'{name!r} -> {want}', got == want, f'-> {got}')
    suite.check("a campaign whose sep is '_' splits on '_'",
                util.split_job_name('run_1_2_3', sep='_') == ('run', (1, 2, 3)))
    suite.check("and its '-' are just part of the title",
                util.split_job_name('a-b_1_2_3', sep='_') == ('a-b', (1, 2, 3)))


def check_selection(suite, title, text, ours):
    """Only names whose parsed title equals ours are selected."""
    suite.section('picking this campaign out of the queue')
    got, warned = said(util.campaign_jobs, text, title)
    suite.check("only this campaign's jobs come back", got == ours, f'-> {got}')
    suite.check("a foreign campaign's longer title is left out",
                202 not in [jid for jid, _ in got])
    suite.check('and is not worth a warning, since it is simply not ours',
                warned == '', f'-> {warned!r}')

    suite.section('a title that itself contains the separator')
    got, _ = said(util.campaign_jobs, text, 'sampling-long')
    suite.check("the longer campaign claims its own job and only it",
                got == [(202, (0, 0, 5))], f'-> {got}')

    suite.section('a title full of regex metacharacters')
    meta = 'a.b*c[d]+e'
    lines = f'11 {meta}-0-0-1\n22 aXbXXcd-e-0-0-1\n33 a.b*c[d]+e-x-0-0-1'
    got, _ = said(util.campaign_jobs, lines, meta)
    suite.check('the title is matched literally, never as a pattern',
                got == [(11, (0, 0, 1))], f'-> {got}')

    suite.section('a name with our title but the wrong indices')
    odd = ('505 sampling-0-0\n606 sampling\n707 sampling-0-0-5\n'
           '808 othercampaign-alpha')
    got, warned = said(util.campaign_jobs, odd, title)
    suite.check('it is left out, since no clone could be bound to it',
                got == [(707, (0, 0, 5))], f'-> {got}')
    suite.check('and both of ours are warned about',
                warned.count('WARNING') == 2, f'-> {warned!r}')
    suite.check("a foreign name that will not split is nobody's business",
                'othercampaign' not in warned, f'-> {warned!r}')

    suite.section('a campaign whose separator is not the default')
    got, _ = said(util.campaign_jobs, '88 run_0_1_2\n99 run-0-1-2', 'run',
                  sep='_')
    suite.check('only names joined by that separator are ours',
                got == [(88, (0, 1, 2))], f'-> {got}')


def check_end_to_end(suite, title, text, ours, dead_squeue=DEAD_SQUEUE):
    """A failed query must stay distinguishable from an empty queue."""
    suite.section('end to end, through scheduler_query')
    trusted, out = util.scheduler_query(f'echo {shlex.quote(text)}')
    got, _ = said(util.campaign_jobs, out, title)
    suite.check('a queue read off the scheduler selects the same jobs',
                trusted and got == ours, f'-> {got}')

    empty = util.scheduler_query('true')
    failed = util.scheduler_query('exit 1')
    suite.check('an empty queue is trusted and holds none of ours',
                empty[0] is True and util.campaign_jobs(empty[1], title) == [])
    suite.check('a failed query is not trusted', failed[0] is False)
    suite.check('the two do not read alike', empty[0] != failed[0])

    work = harness.workdir('job_name_matching')
    (work / 'squeue').write_text(dead_squeue)
    (work / 'squeue').chmod(0o755)
    cmd = (f'PATH={shlex.quote(str(work))}:$PATH '
           + util.basic_scheduler_reports['slurm'])
    trusted, _ = util.scheduler_query(cmd)
    suite.check('with no pipe left to swallow it, a dead squeue is untrusted',
                trusted is False)


def main(title=TITLE, queue_lines=QUEUE_LINES, ours=OURS, required=REQUIRED):
    suite = Suite('job_name_matching')
    text = '\n'.join(queue_lines)
    check_reports_are_plain(suite, required)
    check_parsing(suite, text)
    check_splitting(suite)
    check_selection(suite, title, text, ours)
    check_end_to_end(suite, title, text, ours)
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
