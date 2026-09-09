"""`python -m mdfarmer <command>`: the parts of the package with a command line.

Running a submodule directly -- `python -m mdfarmer.harvester` -- warns, and the
warning is earned: importing the package has already executed that module, and
runpy then runs it a second time under the name __main__, so its module-level
state exists twice. Dispatching from here runs each module once. The old
spelling still works, warning and all, for the job scripts already on disk.
"""
import sys

from . import harvester

COMMANDS = {'harvest': harvester._main}


def main(argv=None, commands=COMMANDS):
    """Hand argv to the named command, or explain what the names are."""
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ('-h', '--help'):
        print(__doc__.splitlines()[0])
        print('\ncommands:')
        for name, command in sorted(commands.items()):
            summary = (command.__doc__ or '').strip().splitlines()
            print(f'  {name:10s} {summary[0] if summary else ""}')
        print('\nEach takes its own --help.')
        return 0
    if argv[0] not in commands:
        print(f'mdfarmer: {argv[0]!r} is not a command; '
              f'try one of {", ".join(sorted(commands))}', file=sys.stderr)
        return 2
    return commands[argv[0]](argv[1:])


if __name__ == '__main__':
    sys.exit(main())
