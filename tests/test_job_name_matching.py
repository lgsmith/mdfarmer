"""Only this campaign's jobs may come back from the scheduler reports.

Job names are title, seed, clone and gen joined by the separator, so a campaign
called 'sampling' and one called 'sampling-long' produce names that share a
prefix. A report that matches on the title alone binds the foreign job id to
this campaign's (seed, clone, gen); when the foreign job ends, the farmer reads
its own live generation as finished and submits a second job into that
generation's directory.
"""
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
OURS = ('101 sampling-0-0-5', '404 sampling-1-2-3')


def fake_queue(fstring, title, lines):
    """Run a report fstring against canned queue output instead of the real
    scheduler: everything up to the first pipe is replaced by an echo."""
    awk = fstring.format(title=title).split('|', 1)[1]
    text = '\n'.join(lines)
    trusted, output = util.scheduler_query(f'echo {shlex.quote(text)} |{awk}')
    return trusted, output.split('\n') if output else []


def main(title=TITLE, queue_lines=QUEUE_LINES, ours=OURS):
    suite = Suite('job_name_matching')
    our_ids = [row.split()[0] for row in ours]
    # How each scheduler is asked for our jobs rather than the whole queue.
    narrowing = {'slurm': '--me', 'lsf': f"-J '{title}-*'"}

    for family in ('slurm', 'lsf'):
        report = util.basic_scheduler_reports[family]
        assoc = util.basic_scheduler_assoc_reports[family]

        suite.section(f'{family}: the job-id report')
        trusted, ids = fake_queue(report, title, queue_lines)
        suite.check('the report runs', trusted)
        suite.check("only this campaign's job ids come back", ids == our_ids,
                    f'-> {ids}')

        suite.section(f'{family}: the association report')
        trusted, rows = fake_queue(assoc, title, queue_lines)
        suite.check('the report runs', trusted)
        suite.check("a foreign campaign's longer title is left out",
                    rows == list(ours), f'-> {rows}')

        suite.section(f'{family}: an empty queue')
        trusted, rows = fake_queue(assoc, title, ())
        suite.check('a quiet queue is trusted and empty',
                    trusted and rows == [])

        suite.section(f'{family}: what is asked of the scheduler')
        suite.check('the queue is narrowed before awk ever sees it',
                    narrowing[family] in report.format(title=title))
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
