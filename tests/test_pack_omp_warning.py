"""A packed job says out loud that it is overriding OMP_NUM_THREADS.

_mdrun_env already rewrites the variable per replica, because mdrun refuses to
start when it disagrees with -ntomp and a batch job inherits whatever the
submitting shell exported. Silently correcting it leaves an operator who set it
deliberately with no way to see that the job ignored them, so the pack says so
once at startup rather than once per replica.
"""
import io
import sys
from contextlib import redirect_stdout

import harness
from harness import Suite

from mdfarmer import gmx_pack as gp

CPUS_PER_TASK = 8
N_REPLICAS = 2
UNEVEN_CORES = (6, 2)


def warn_output(layout, environ):
    """(returned value, everything the warning printed) for one environment."""
    buffer = io.StringIO()
    with redirect_stdout(buffer):
        inherited = gp.warn_inherited_omp_threads(layout, environ=environ)
    return inherited, buffer.getvalue()


def main():
    suite = Suite('pack_omp_warning')
    even = gp.member_core_layout(CPUS_PER_TASK, N_REPLICAS)
    uneven = gp.member_core_layout(CPUS_PER_TASK, N_REPLICAS,
                                   member_cores=UNEVEN_CORES)

    suite.section('a hostile value inherited from the submitting shell')
    inherited, out = warn_output(even, {gp.OMP_THREADS_ENV: '16'})
    suite.check('is reported back to the caller', inherited == '16',
                f'-> {inherited!r}')
    suite.check('and named in the warning', '16' in out, f'-> {out.strip()}')
    suite.check('alongside the cores each replica really gets',
                '4' in out, f'-> {out.strip()}')
    suite.check('saying which variable is being overridden',
                gp.OMP_THREADS_ENV in out)
    suite.check('once, not once per replica', out.count('WARNING') == 1,
                f'-> {out.count("WARNING")}')

    suite.section('a value that already agrees is not worth a warning')
    inherited, out = warn_output(even, {gp.OMP_THREADS_ENV: '4'})
    suite.check('nothing is printed', out == '', f'-> {out.strip()}')
    suite.check('and nothing is reported', inherited is None,
                f'-> {inherited!r}')

    suite.section('an unset variable is the ordinary case')
    inherited, out = warn_output(even, {})
    suite.check('nothing is printed', out == '', f'-> {out.strip()}')
    suite.check('and nothing is reported', inherited is None,
                f'-> {inherited!r}')

    suite.section('a pack whose members get different core counts')
    inherited, out = warn_output(uneven, {gp.OMP_THREADS_ENV: '6'})
    suite.check('warns even though one replica happens to match',
                inherited == '6', f'-> {inherited!r}')
    suite.check('and lists every count it will hand out',
                '6' in out and '2' in out, f'-> {out.strip()}')
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
