"""A gmx that has to be launched by something else, and a tender that never is.

Some sites ship only an MPI GROMACS. Such a binary calls MPI_Init even for
grompp, which is pure serial preprocessing, and under a scheduler it aborts
unless it was launched the way its MPI expects -- so every generation dies
before any MD. The cure is a launcher in front of the binary, which gmx_bin
could not express while it was interpolated as a single argv element.

Whether a launcher is needed, and which, is the site's business. This suite
only pins that the library can carry one.

The second half pins the reason a launcher in the job script is sufficient:
the tender never runs gmx. It reads each generation's progress from the JSON
status file the runner wrote, so nothing on the submitting host needs GROMACS
at all -- and a hook that quietly started shelling out would break exactly the
sites whose gmx is reachable only inside a job.
"""
import inspect
import sys

import harness
from harness import Suite

from mdfarmer import gmx_pack, gmx_simulate as gmx

# What a site with an MPI-only GROMACS has to run.
LAUNCHED = ['mpirun', '-n', '1', 'gmx_mpi']
# Every gmx subcommand the runner issues.
SUBCOMMANDS = ('grompp', 'mdrun', 'convert-tpr', 'dump')
# Hooks the tender calls per tick, which must not need a gmx of their own.
TENDER_HOOKS = ('gmx_gen_progress', 'gmx_try_recover_gen')


def main(launched=LAUNCHED, subcommands=SUBCOMMANDS, hooks=TENDER_HOOKS):
    suite = Suite('gmx_launcher')

    suite.section('a bare binary is still a bare binary')
    suite.check('a plain name leads the argv',
                gmx.gmx_argv('gmx', 'grompp') == ['gmx', 'grompp'],
                f"-> {gmx.gmx_argv('gmx', 'grompp')}")
    spaced = '/opt/my gromacs/gmx'
    suite.check('a path with a space stays one word, not two',
                gmx.gmx_argv(spaced, '-version') == [spaced, '-version'],
                f'-> {gmx.gmx_argv(spaced, "-version")}')
    suite.check('the default is a bare string, so nothing assumes a launcher',
                isinstance(gmx.GMX_BIN, str), f'-> {gmx.GMX_BIN!r}')

    suite.section('a launcher vector reaches the front of every call')
    got = gmx.gmx_argv(launched, 'grompp', '-f', 'gen.mdp')
    suite.check('the whole vector leads, and the subcommand follows it',
                got == [*launched, 'grompp', '-f', 'gen.mdp'], f'-> {got}')
    suite.check('the vector is copied, not aliased',
                gmx.gmx_argv(launched, 'mdrun')[:len(launched)] == launched
                and len(launched) == 4, f'-> {launched}')
    for sub in subcommands:
        argv = gmx.gmx_argv(launched, sub)
        suite.check(f'{sub} is launched too',
                    argv[:len(launched)] == launched and argv[-1] == sub)

    suite.section('every gmx the runner issues goes through it')
    source = inspect.getsource(gmx)
    # A bare [gmx_bin, ...] would bypass the launcher and abort on such a site.
    bare = [line.strip() for line in source.splitlines()
            if '[gmx_bin,' in line or '[gmx_bin ,' in line]
    suite.check('no call site still builds argv around gmx_bin directly',
                not bare, f'-> {bare[:2]}')
    pack_source = inspect.getsource(gmx_pack)
    pack_bare = [line.strip() for line in pack_source.splitlines()
                 if '[gmx_bin,' in line]
    suite.check('nor does the pack, which probes the binary for -ntmpi',
                not pack_bare, f'-> {pack_bare[:2]}')

    suite.section('the tender needs no gmx of its own')
    for name in hooks:
        hook_source = inspect.getsource(getattr(gmx, name))
        shells_out = [token for token in
                      ('_run(', '_run_capture(', 'sp.run(', 'checkpoint_step(',
                       'checkpoint_part_step(', 'gmx_argv(')
                      if token in hook_source]
        suite.check(f'{name} reads what the runner wrote, and runs nothing',
                    not shells_out, f'-> calls {shells_out}')
    return suite.report()


if __name__ == '__main__':
    sys.exit(harness.run_suite(main))
