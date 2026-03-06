from argparse import ArgumentParser
import datetime
import os
try:
    import readline
except ImportError:
    readline = None
import signal
import sys
import warnings
from . import load_file
from .exceptions import ParmedError
from .tools.logos import Logo
from .tools.exceptions import SeriousParmWarning, InterpreterError, ParmError
from .tools.parmed_cmd import ParmedCmd
from .tools.actions import Action
from .tools.parmlist import ParmList
from . import __version__

def clapp():
    # Set up new excepthook to clean up fatal exception printouts
    def interrupted(*args, **kwargs):
        """ Handle interruptions gracefully """
        sys.stdout.write('Interrupted\n')
        sys.exit(1)

    signal.signal(signal.SIGINT, interrupted)

    # Define our own custom warning printer
    def _print_warnings(message, category, filename, lineno, file=None, line=None):
        """ Override the default showwarning method """
        file = sys.stderr if file is None else file
        try:
            file.write(f'{category.__name__}: {message}\n')
        except IOError:
            pass

    warnings.showwarning = _print_warnings

    # Set up parser
    parser = ArgumentParser()
    parser.add_argument('-v', '--version', action='version', version=f"%(prog)s: Version {__version__}")
    group = parser.add_argument_group('Input Files')
    group.add_argument('-i', '--input', dest='script', default=[],
             metavar='FILE', help='''Script with ParmEd commands to execute. Default
             reads from stdin. Can be specified multiple times to process multiple
             input files.''', action='append')
    group.add_argument('-p', '--parm', dest='prmtop', default=[],
             metavar='<prmtop>', action='append', help='''List of topology files to load into
             ParmEd. Can be specified multiple times to process multiple topologies.''')
    group.add_argument('-c', '--inpcrd', dest='inpcrd', default=[],
             metavar='<inpcrd>', action='append', help='''List of inpcrd files to
             load into ParmEd. They are paired with the topology files in the same
             order that each set of files is specified on the command-line.''')
    group = parser.add_argument_group('Output Files')
    group.add_argument('-O', '--overwrite', dest='overwrite', default=False,
             help='Allow ParmEd to overwrite existing files.', action='store_true')
    group.add_argument('-l', '--logfile', dest='logfile', metavar='FILE',
             default='parmed.log', help='''Log file with every command executed
             during an interactive ParmEd session. Default is parmed.log''')
    group = parser.add_argument_group('Interpreter Options', '''These options affect
             how the ParmEd interpreter behaves in certain cases.''')
    group.add_argument('--prompt', dest='prompt', default='>', metavar='PROMPT',
             help='String to use as a command prompt.')
    group.add_argument('-n', '--no-splash', dest='printlogo', action='store_false',
             help='Prevent printing the greeting logo.', default=True)
    group.add_argument('-e', '--enable-interpreter', dest='interpreter',
             action='store_true', default=False,
             help='''Allow arbitrary single Python commands or blocks of Python code
             to be run. By default Python commands will not be run as a safeguard
             for your system. Make sure you trust the source of the ParmEd command
             before turning this option on.''')
    group = parser.add_argument_group('Error Handling', '''These options control how
             ParmEd handles various errors and warnings that appear occur during the
             course of Action execution''')
    group.add_argument('-s', '--strict', dest='strict', action='store_true',
             default=True, help='''Prevent scripts from running past unrecognized input and actions
             that end with an error. In interactive mode, actions with unrecognized inputs and
             failed actions prevent any changes from being made to the topology, but does not quit
             the interpreter. This is the default behavior.''')
    group.add_argument('-r', '--relaxed', dest='strict', action='store_false',
             help='''Scripts ignore unrecognized input and simply skip over failed
             actions, executing the rest of the script. Unrecognized input in the
             interactive interpreter emits a non-fatal warning.''')
    parser.add_argument('prmtop_cl', nargs='?', metavar='<prmtop>', default=None,
             help='Topology file to analyze.')
    parser.add_argument('script_cl', nargs='?', metavar='<script>', default=None,
             help='File with a series of ParmEd commands to execute.')

    opt = parser.parse_args()

    # If the user specified a prmtop and script in the 'old' way, append them to the
    # list constructed via the --parm and -i flags -- they come at the end
    if opt.script_cl is not None:
        opt.script.append(opt.script_cl)

    if opt.prmtop_cl is not None:
        opt.prmtop.append(opt.prmtop_cl)

    # Load the splash screen
    if opt.printlogo:
        splash = Logo()
        print(splash)

    # Set our warning filter
    if not opt.strict:
        warnings.filterwarnings('always', category=SeriousParmWarning)

    # Set our overwrite preferences
    Action.overwrite = opt.overwrite

    amber_prmtop = ParmList()
    for i, parm in enumerate(opt.prmtop):
        if i < len(opt.inpcrd):
            amber_prmtop.add_parm(parm)
            amber_prmtop.parm.load_rst7(opt.inpcrd[i])
            print(f'Loaded Amber topology file {parm} with coordinates from {opt.inpcrd[i]}\n')
        else:
            amber_prmtop.add_parm(parm)
            print(f'Loaded Amber topology file {parm}')

    if len(opt.script) > 0:
        # Read from the list of scripts
        print(opt.script)
        # Make sure that all scripts exist, quitting if we are strict and not all
        # scripts exist. Don't do anything until we know that all scripts exist.
        for script in opt.script:
            if not os.path.exists(script):
                warnings.warn(f'Script file {script} cannot be found.', SeriousParmWarning)

        # We have already pre-screened the scripts.
        for script in opt.script:
            if not os.path.exists(script):
                continue
            print(f'Reading actions from {script}\n')
            parmed_commands = ParmedCmd(amber_prmtop, stdin=open(script, 'r'))
            parmed_commands.use_rawinput = 0
            parmed_commands.interpreter = opt.interpreter
            parmed_commands.prompt = ''
            # Loop through all of the commands
            try:
                parmed_commands.cmdloop()
            except InterpreterError as err:
                sys.exit(f'{err.__class__.__name__}: {err}')
            except ParmedError:
                # This has already been caught and printed. If it was re-raised,
                # then that means we wanted to exit
                sys.exit(1)

    else:
        close_log_file = False
        parmed_commands = ParmedCmd(amber_prmtop)
        parmed_commands.intro = "Reading input from STDIN..."
        parmed_commands.interpreter = opt.interpreter
        parmed_commands.prompt = opt.prompt.strip() + ' '
        # Set the log file and logging, but only if interactive
        if os.isatty(sys.stdin.fileno()):
            if os.path.exists(opt.logfile) and readline is not None:
                # Load the logfile as a history, and get rid of any timestamps
                # from the history
                try:
                    f = open(opt.logfile, 'r')
                except IOError:
                    pass # Do not have permission or something
                else:
                    for line in f:
                        if line.startswith('# Log started on'):
                            continue
                        readline.add_history(line.rstrip())
                    f.close()
            try:
                logfile = open(opt.logfile, 'a')
            except IOError:
                print("Could not open logfile. Skip logging", file=sys.stderr)
                close_log_file = False
            else:
                now = datetime.datetime.now()
                logfile.write(f'# Log started on {now.month:02d}/{now.day:02d}/{now.year:d} '
                              f'[mm/dd/yyyy] at {now.hour:02d}:{now.minute:02d}:{now.second:02d}\n')
                if len(sys.argv) > 1:
                    logfile.write(f"# Command line arguments: {' '.join(sys.argv[1:])}\n")
                if len(amber_prmtop) > 0:
                    logfile.write(f"# Loaded topologies: {', '.join(parm.name for parm in amber_prmtop)}")
                parmed_commands.setlog(logfile)
                close_log_file = True
        # Loop through all of the commands
        try:
            try:
                parmed_commands.cmdloop()
            except InterpreterError as err:
                sys.exit(f"{err.__class__.__name__}: {err}")
            except ParmedError:
                # This has already been caught and printed. If it was re-raised,
                # then that means we wanted to exit
                sys.exit(1)
        finally:
            if close_log_file:
                logfile.close()

    print('Done!')
