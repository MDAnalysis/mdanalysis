"""
Provides a class for reading GROMACS-style files. The key component to these
files is that the ; character is a comment character and everything after ; is
ignored.
"""
from ._cpp import CPreProcessor

class GromacsFile:
    """
    A GROMACS file that recognizes the ";" character as a 'comment' token. It
    can be iterated over and generally treated like a file object, but only
    spits out strings that have been truncated at its first comment character.
    
    There is currently no way to recognize a ; as a _non_ comment character,
    since allowing an escape character does not seem to be common practice and
    would likely introduce negative performance implications.

    Parameters
    ----------
    fname : str or file-like
        Name of the file to parse or file-like object to parse
    defines : dict{str : str}, optional
        List of defines for the preprocessed file, if any
    includes : list of str, optional
        List of include files. Default is taken from environment variables
        GMXDATA or GMXBIN if they are set. Otherwise, it is looked for in /usr,
        /usr/local, /opt, or /opt/local. If it is still not found, it is looked
        for relative to whatever ``mdrun`` executable is in your path
    notfound_fatal : bool, optional
        If True, missing include files are fatal. If False, they are a warning.
        Default is True
    """

    def __init__(self, fname, **kwargs):
        self._handle = CPreProcessor(fname, **kwargs)
        self.closed = False
        self.line_number = 0

    def __iter__(self):
        # Iterate over the file, treating an ending \ as a continuation
        parts = []
        for line in self._handle:
            try:
                idx = line.index(';')
                if not parts:
                    yield '%s\n' % line[:idx]
                else:
                    parts.append('%s' % line[:idx])
                    yield '%s\n' % ''.join(parts)
                    parts = []
            except ValueError:
                # There is no comment...
                if line.rstrip('\r\n').endswith('\\'):
                    chars = list(reversed(line.rstrip('\r\n')))
                    del chars[chars.index('\\')]
                    parts.append('%s ' % ''.join(reversed(chars)))
                elif parts:
                    parts.append(line)
                    yield ''.join(parts)
                    parts = []
                else:
                    yield line

    @property
    def included_files(self):
        return self._handle.included_files

    def readline(self):
        parts = []
        self.line_number += 1
        line = True
        while line:
            line = self._handle.readline()
            try:
                idx = line.index(';')
                if not parts:
                    return '%s\n' % line[:idx]
                else:
                    parts.append('%s' % line[:idx])
                    return '%s\n' % ''.join(parts)
            except ValueError:
                # There is no comment...
                if line.rstrip('\r\n').endswith('\\'):
                    chars = list(reversed(line.rstrip('\r\n')))
                    del chars[chars.index('\\')]
                    parts.append(f"{''.join(reversed(chars))} ")
                elif parts:
                    parts.append(line)
                    return ''.join(parts)
                else:
                    return line

    def readlines(self):
        return [line for line in self]

    def read(self):
        return ''.join(self.readlines())

    def close(self):
        self._handle.close()
        self.closed = True

    def __del__(self):
        try:
            self.closed or self._handle.close()
        except AttributeError:
            # It didn't make it out of the constructor
            pass
