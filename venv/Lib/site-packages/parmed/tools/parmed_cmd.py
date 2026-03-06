"""
This sets up the command interpreter for textual ParmEd (parmed).
"""

# Load some system modules that may be useful for various users in shell mode
import cmd
import os
import sys
import traceback
from glob import glob

from ..exceptions import ParmedError, ParmedWarning
from .actions import COMMANDMAP, Usages
from .argumentlist import ArgumentList
from .exceptions import InterpreterError

try:
    import readline
except ImportError:
    readline = None

_COMMANDLOGS = []

def log_commands(func):
    """ Decorator to write the line to a file """
    global _COMMANDLOGS
    def new_func(self, line, *args, **kwargs):
        if self._logfile is not None and line != 'EOF':
            self._logfile.write('%s\n' % line)
            if readline is None:
                _COMMANDLOGS.append(line)
            try:
                self._logfile.flush()
            except AttributeError:
                pass
        return func(self, line, *args, **kwargs)

    return new_func

class ParmedCmd(cmd.Cmd):
    """
    ParmEd command interpreter. The docstrings for each do_* command are simply
    copied from the docstring for that Action's docstring
    """

    prompt = "> "
    _populated = False

    def __init__(self, parm, stdin=None, stdout=None):
        """ Load a topology file into the interpreter """
        self.parm = parm
        cmd.Cmd.__init__(self, stdin=stdin, stdout=stdout)
        self.continued = ''
        self._logfile = None
        try:
            self._exit_on_fatal = not os.isatty(self.stdin.fileno())
        except AttributeError:
            self._exit_on_fatal = False
        if not self._populated: self.populate_actions()

    def setlog(self, f):
        """ Open up a log file to start logging the commands. """
        if f is None:
            return
        if hasattr(f, 'write'):
            self._logfile = f
        else:
            self._logfile = open(f, 'w')

    def emptyline(self):
        """ Override emptyline so that empty lines are simply ignored """
        pass

    @log_commands
    def precmd(self, line):
        """ Override this to strip out comments, but preserve #s in quotes """
        in_quotes = False
        quote_type = ''
        idx = -1
        for i, char in enumerate(line):
            if in_quotes:
                if char == quote_type:
                    in_quotes = False
                    continue
                continue
            if char == '"':
                in_quotes = True
                quote_type = '"'
                continue
            if char == "'":
                in_quotes = True
                quote_type = "'"
                continue
            if char == '#':
                idx = i
                break

        if idx < 0:
            if line.endswith('\\'):
                self.continued += line[:-1] + ' '
                return ''
            line = self.continued + line
            self.continued = ''
            return line
        else:
            line = self.continued + line[:idx]
            self.continued = ''
            return line

    def parseline(self, line):
        """
        Override parseline so that it will set args as ArgumentList
        """
        line = line.strip()
        if not line:
            return None, None, line
        elif line[0] == '?':
            line = 'help ' + line[1:]
        elif line[-1] == '?' and len(line.split()) == 1:
            line = 'help ' + line[:-1]
        elif line[0] == '!':
            if hasattr(self, 'do_shell'):
                line = 'shell ' + line[1:]
            else:
                return None, None, line
        i, n = 0, len(line)
        while i < n and line[i] in self.identchars:
            i = i+1
        cmd, arg = line[:i], ArgumentList(line[i:].strip())
        if len(line) == 4 and line == 'help':
            arg = ''
        if line[:5] == 'help ':
            arg = line[i:].strip()
        return cmd, arg, line

    def do_parmout(self, line):
        # Store this action for later use. This action is unique
        try:
            self.parmout = COMMANDMAP['parmout'](self.parm, line)
        except (ParmedError, ParmedWarning) as err:
            self.stdout.write('Action parmout failed.\n\t')
            self.stdout.write(f'{type(err).__name__}: {err}\n')
            if self._exit_on_fatal:
                raise
        except Exception as err:
            if self._exit_on_fatal:
                raise
            self.stdout.write(f'Unexpected failure:\n{type(err).__name__}: {err}\n\nTraceback is\n\n')
            traceback.print_tb(sys.exc_info()[2], file=self.stdout)

    def _normaldo(self, ActionClass, line):
        """ The standard action command does this stuff """
        actionname = ActionClass.__name__
        try:
            action = ActionClass(self.parm, line)
            if action.valid:
                self.stdout.write('%s\n' % action)
                action.execute()
        except (ParmedError, ParmedWarning) as err:
            self.stdout.write('Action %s failed\n\t' % actionname)
            self.stdout.write('%s: %s\n' % (type(err).__name__, err))
            if self._exit_on_fatal:
                raise
        except Exception as err:
            if self._exit_on_fatal:
                raise
            self.stdout.write('Unexpected failure:\n%s: %s\n\n'
                              'Traceback is\n\n' % (type(err).__name__, err))
            traceback.print_tb(sys.exc_info()[2], file=self.stdout)

    @classmethod
    def populate_actions(cls):
        """
        This will create all of the do_Command methods to trigger command
        auto-complete. This eliminates the need to modify the ParmedCmd class
        when a new command is added
        """
        for _cmd, cmdclass in COMMANDMAP.items():
            if _cmd in ('source', 'go', 'EOF', 'quit', 'help', 'parmout'):
                continue
            cmdname = cmdclass.__name__
            exec(f'cls.do_{cmdname} = lambda self, line: self._normaldo(COMMANDMAP["{_cmd}"], line)')
        cls._populated = True

    def do_source(self, line):
        action = COMMANDMAP['source'](self.parm, line)
        if not action.valid:
            return
        self.stdout.write(f'{str(action)}\n')
        _cmd = ParmedCmd(self.parm, stdin=open(action.filename, 'r'), stdout=self.stdout)
        _cmd.prompt = ''
        _cmd.interpreter = self.interpreter
        _cmd.use_rawinput = 0
        # If we exit on error, call cmdloop without try protection. Otherwise,
        # we need to protect cmdloop to catch any ParmedError's that are passed
        # through to avoid exiting the top-level interpreter
        if self._exit_on_fatal:
            _cmd.cmdloop()
        else:
            # If we are not exiting on error, just ignore these parm errors
            # since they've already been printed. Now just wait for the next
            # command from the user.
            try:
                _cmd.cmdloop()
            except:
                pass

    def do_go(self, line):
        """
        Stops reading commands and executes any 'parmout' command that had
        been issued
        """
        if hasattr(self, 'parmout'):
            self.stdout.write(f'{str(self.parmout)}\n')
            try:
                self.parmout.execute()
            except ParmedError as err:
                self.stdout.write(f'{type(err).__name__}: {err}\n')
                if self._exit_on_fatal:
                    raise
            except Exception as err:
                if self._exit_on_fatal:
                    raise
                self.stdout.write(f'Unexpected failure:\n{type(err).__name__}: {err}\n\n')
                self.stdout.write('Traceback is\n\n')
                traceback.print_tb(sys.exc_info()[2], file=self.stdout)

        return True

    # EOF is treated the same as "go"
    do_EOF = do_go

    def do_quit(self, line):
        """
        Stops reading commands and quits WITHOUT running the last 'parmout' 
        command that had been issued
        """
        return True

    def do_history(self, line):
        """
        Print the command history
        """
        global _COMMANDLOGS
        if readline is None:
            for line in _COMMANDLOGS:
                self.stdout.write(f'{line}\n')
        else:
            for i in range(readline.get_current_history_length()):
                self.stdout.write(f'{readline.get_history_item(i + 1)}\n')

    def default(self, line):
        words = line.split()
        mycmd = words[0].lower()
        if mycmd not in COMMANDMAP:
            self.stdout.write(f"{words[0]} command not recognized\n")
        else:
            args = ArgumentList(words[1:])
            self._normaldo(COMMANDMAP[mycmd], args)

    def do_shell(self, line):
        """ Support the limited interpreter """
        if not self.interpreter:
            raise InterpreterError("Interpreter not enabled! Use '-e' to enable")
        line = str(line)
        if line.strip() == '!':
            self._python_shell()
        else:
            globals_ = dict(amber_prmtop=self.parm)
            globals_.update(globals())
            try:
                exec(line.strip(), globals_)
            except Exception as err:
                self.stdout.write(f"{type(err).__name__}: {err}\n")

    def completedefault(self, text, line, begidx, endidx):
        partial = line[:endidx]
        idx = max(partial.rfind(' '), partial.rfind('\t')) + 1
        beg_token = partial[idx:begidx]
        return [s.replace(beg_token, '', 1) for s in glob(partial[idx:] + '*')]

    def _python_shell(self):
        """ Drop into a limited interpreter and read until we see !! """
        from . import actions
        python_interpreter = PythonCmd(stdin=self.stdin, stdout=self.stdout)
        python_interpreter.use_rawinput = self.use_rawinput
        python_interpreter.setlog(self._logfile)
        if not self.prompt:
            python_interpreter.prompt = ''
        python_interpreter.cmdloop()
        try:
            globals_ = dict(amber_prmtop=self.parm, actions=actions)
            globals_.update(globals())
            exec(python_interpreter.command_string, globals_)
        except Exception as err:
            self.stdout.write(f"{type(err).__name__}: {err}\n")
      
    def do_help(self, arg):
        " Modify the original do_help to pull docstrings from actions "
        if arg:
            # XXX check arg syntax
            try:
                func = getattr(self, 'help_' + arg)
            except (AttributeError, TypeError):
                try:
                    doc=getattr(self, 'do_' + arg).__doc__
                    if doc:
                        self.stdout.write(f"{str(doc)}\n")
                        return
                except (AttributeError, TypeError):
                    pass
                try:
                    _action = COMMANDMAP[arg.lower()]
                    if _action.needs_parm:
                        doc = '{} [parm <idx>|<name>]\n{}'
                    else:
                        doc = '{}\n{}'
                    doc = doc.format(_fmt_usage(Usages[arg.lower()]), _action.__doc__)
                    if doc:
                        self.stdout.write(f"{str(doc)}\n")
                        return
                except (AttributeError, KeyError):
                    pass
                self.stdout.write(f"No help available for {arg}\n")
                return
                func()
        else:
            names = self.get_names()
            cmds_doc = []
            cmds_undoc = []
            help = {}
            for name in names:
                if name[:5] == 'help_':
                    help[name[5:]]=1
            names.sort()
            # There can be duplicates if routines overridden
            prevname = ''
            for name in names:
                if name[:3] == 'do_':
                    if name == prevname:
                        continue
                    prevname = name
                    cmd=name[3:]
                    if cmd in help:
                        cmds_doc.append(cmd)
                        del help[cmd]
                    elif getattr(self, name).__doc__:
                        cmds_doc.append(cmd)
                    elif cmd.lower() in COMMANDMAP and cmd.lower() in Usages:
                        cmds_doc.append(cmd)
                    else:
                        cmds_undoc.append(cmd)
            self.stdout.write(f"{str(self.doc_leader)}\n")
            self.print_topics(self.doc_header, cmds_doc, 15,80)
            self.print_topics(self.misc_header, help.keys(), 15,80)
            self.print_topics(self.undoc_header, cmds_undoc, 15,80)

class PythonCmd(cmd.Cmd):
    """ Command interpreter for limited Python interpreter """
    prompt = "py >>> "

    def __init__(self, stdin=None, stdout=None):
        cmd.Cmd.__init__(self, stdin=stdin, stdout=stdout)
        self._logfile = None
        self.command_string = ""

    def setlog(self, f):
        """ Open up a log file to start logging the commands """
        if f is None: return
        if hasattr(f, 'write'):
            self._logfile = f
        else:
            self._logfile = open(f, 'w')

    precmd = log_commands(cmd.Cmd.precmd)

    def emptyline(self):
        """ Ignore all empty lines """
        self.command_string += "\n"

    def do_shell(self, line):
        """ Break out of the shell """
        return True

    def do_EOF(self, line):
        """ Break out of the shell """
        raise InterpreterError("EOF hit while in interpreter!")

    def parseline(self, line):
        """ Override parseline """
        if not line.strip():
            return False, False, False
        return line, line, line

    def default(self, line):
        """ Add this onto the command string """
        if line.strip() == '!!':
            return True
        self.command_string += line + '\n'

# To pretty-print usage statements. Some are getting long, so they need to be
# shortened and printed in a helpfully formatted way
def _fmt_usage(usage):
    MAX_LINE_LEN = 79

    def split_options(string):
        whitespace = '\r\n\t '
        openbraces = '[<'
        closebraces = '>]'
        # Splits up the options to keep <> and []'s together
        words, word = [], ''
        i = 0
        # Keep track of when we are inside a <> and when we're not
        nopen = 0
        while i < len(string):
            char = string[i]
            if char in whitespace:
                # Ignore consecutive whitespace
                if not word:
                    i += 1
                    continue
                else:
                    # Only split the string here if we are not in open braces
                    if nopen > 0:
                        word = word.rstrip() + ' '
                        i += 1
                        continue
                    else:
                        words.append(word)
                        word = ''
                        i += 1
                        continue
            # This is not whitespace. If it is an open brace, increment nopen.
            # If it is a closing brace, decrement nopen (check that it does not
            # go below zero)
            if char in openbraces:
                nopen += 1
                word += char
                i += 1
                continue
            if char in closebraces:
                nopen = max(nopen - 1, 0)
                word += char
                i += 1
                continue
            # It is none of these, simply add the character to the word and
            # continue on with the next letter
            word += char
            i += 1
        # Now add the final word to the list of words if it's not blank
        if word: words.append(word)
        return words

    uwords = split_options(usage)
    indent_size = len(uwords[0]) + 1
    lines = []
    line = ''
    for i, uword in enumerate(uwords):
        if not line:
            line = uword
        elif len(line) + len(uword) < MAX_LINE_LEN:
            line += ' ' + uword
        else:
            lines.append(line)
            line = ' ' * indent_size + uword
    if line:
        lines.append(line)
    return '\n'.join(lines)
