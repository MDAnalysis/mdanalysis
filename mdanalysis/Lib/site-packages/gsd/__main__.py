# Copyright (c) 2016-2023 The Regents of the University of Michigan
# Part of GSD, released under the BSD 2-Clause License.

"""The GSD command line interface.

To simplify ad hoc usage of :py:mod:`gsd`, this module provides a command line
interface for interacting with GSD files. The primary entry point is a single
command for starting a Python interpreter with a GSD file pre-loaded::

    $ gsd read trajectory.gsd

The following options are available for the ``read`` subcommand:

.. program:: read

.. option:: -s schema, --schema schema

    The schema of the GSD file. Supported values for ``schema`` are "hoomd" and
    "none".

.. option:: -m mode, --mode mode

    The mode in which to open the file. Valid modes are identical to those
    accepted by :func:`gsd.fl.open`.
"""

import sys
import argparse
import code

from .version import __version__
from .hoomd import open as hoomd_open
from . import fl


def _print_err(msg=None, *args):
    print(msg, *args, file=sys.stderr)


SHELL_BANNER = """Python {python_version}
gsd {gsd_version}

File: {fn}
{extras}
The GSD file handle is available via the "handle" variable.
For supported schema, you may access the trajectory using the "traj" variable.
Type "help(handle)" or "help(traj)" for more information.
The gsd and gsd.fl packages are always loaded.
Schema-specific modules (e.g. gsd.hoomd) are loaded if available."""


def main_read(args):
    """Main function to launch a Python interpreter with an open GSD file."""
    # Default to a new line for well-formatted printing.
    local_ns = {
        'gsd': sys.modules['gsd'],
        'gsd.hoomd': sys.modules['gsd.hoomd'],
        'gsd.fl': sys.modules['gsd.fl'],
    }
    attributes = {}

    if args.schema == 'hoomd':
        traj = hoomd_open(args.file, mode=args.mode)
        handle = traj.file
        local_ns.update({
            'handle': handle,
            'traj': traj,
        })
        attributes.update({"Number of frames": len(traj)})
    else:
        if args.mode not in ['rb', 'rb+', 'ab']:
            raise ValueError("Unsupported schema for creating a file.")
        handle = fl.open(args.file, args.mode)
        local_ns.update({
            'handle': handle,
        })

    extras = "\n".join(
        "{}: {}".format(key, val) for key, val in attributes.items())

    code.interact(local=local_ns,
                  banner=SHELL_BANNER.format(python_version=sys.version,
                                             gsd_version=__version__,
                                             fn=args.file,
                                             extras=extras + "\n"))


def main():
    """Entry point to the GSD command-line interface.

    This function handles parsing command-line arguments and launching the
    appropriate subcommand based on the first argument to ``gsd`` on the
    command line. At present the following commands are supported:

        * read
    """
    parser = argparse.ArgumentParser(
        description="The gsd package encodes canonical readers and writers "
        "for the gsd file format.")
    parser.add_argument('--version',
                        action='store_true',
                        help="Display the version number and exit.")
    parser.add_argument('--debug',
                        action='store_true',
                        help="Show traceback on error for debugging.")
    subparsers = parser.add_subparsers()

    parser_read = subparsers.add_parser('read')
    parser_read.add_argument('file',
                             type=str,
                             nargs='?',
                             help="GSD file to read.")
    parser_read.add_argument('-s',
                             '--schema',
                             type=str,
                             default='hoomd',
                             choices=['hoomd', 'none'],
                             help="The file schema.")
    parser_read.add_argument(
        '-m',
        '--mode',
        type=str,
        default='rb',
        choices=['rb', 'rb+', 'wb', 'wb+', 'xb', 'xb+', 'ab'],
        help="The file mode.")
    parser_read.set_defaults(func=main_read)

    # This is a hack, as argparse itself does not
    # allow to parse only --version without any
    # of the other required arguments.
    if '--version' in sys.argv:
        print('gsd', __version__)
        sys.exit(0)

    args = parser.parse_args()

    if not hasattr(args, 'func'):
        parser.print_usage()
        sys.exit(2)
    try:
        args.func(args)
    except KeyboardInterrupt:
        _print_err()
        _print_err("Interrupted.")
        if args.debug:
            raise
        sys.exit(1)
    except RuntimeWarning as warning:
        _print_err("Warning: {}".format(warning))
        if args.debug:
            raise
        sys.exit(1)
    except Exception as error:
        _print_err('Error: {}'.format(error))
        if args.debug:
            raise
        sys.exit(1)
    else:
        sys.exit(0)


if __name__ == '__main__':
    main()
