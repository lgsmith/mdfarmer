"""Run every suite, reporting pass, fail or skip for each.

Exit status is the number of suites that failed. A suite that skipped -- its
dependencies are absent, usually GROMACS -- does not count as a failure.
"""
import subprocess as sp
import sys
from pathlib import Path

import harness

TESTS_DIR = Path(__file__).resolve().parent
SUITE_GLOB = 'test_*.py'

# Exit status a suite uses to mean "dependencies absent", per harness.run_suite.
SKIP_STATUS = 77


def main(tests_dir=TESTS_DIR, suite_glob=SUITE_GLOB, skip_status=SKIP_STATUS):
    outcomes = {}
    for suite in sorted(tests_dir.glob(suite_glob)):
        print(f'\n{"=" * 72}\n{suite.name}\n{"=" * 72}', flush=True)
        result = sp.run([sys.executable, '-u', str(suite)], cwd=tests_dir)
        outcomes[suite.name] = result.returncode

    print(f'\n{"=" * 72}\nsummary\n{"=" * 72}')
    failed = 0
    for name, status in outcomes.items():
        if status == 0:
            label = 'pass'
        elif status == skip_status:
            label = 'skip'
        else:
            label = 'FAIL'
            failed += 1
        print(f'  {label}  {name}')
    return failed


if __name__ == '__main__':
    sys.exit(main())
