"""
Provides a class for reading CHARMM-style files. The key component to these
files is that the ! character is a comment character and everything after ! is
ignored.
"""
from ..utils.io import genopen

class CharmmFile:
    """
    A CHARMM file that recognizes the "!" character as a 'comment' token. It
    can be iterated over and generally treated like a file object, but only
    spits out strings that have been truncated at its first comment character.
    
    There is currently no way to recognize a ! as a _non_ comment character,
    since allowing an escape character does not seem to be common practice and
    would likely introduce negative performance implications.
    """

    def __init__(self, fname, mode='r'):
        if mode not in ('r', 'w'):
            raise ValueError('Cannot open CharmmFile with mode "%s"' % mode)
        if mode == 'r':
            self.status = 'OLD'
        else:
            self.status = 'NEW'
        self._handle = genopen(fname, mode)
        self.closed = False
        self.line_number = 0
        self.comment = ''

    def __enter__(self):
        self._handle.__enter__()
        return self

    def __exit__(self, *args):
        if not self.closed:
            self.close()

    def tell(self):
        return self._handle.tell()

    def seek(self, value):
        return self._handle.seek(value)

    def write(self, *args, **kwargs):
        return self._handle.write(*args, **kwargs)

    def __iter__(self):
        # Iterate over the file
        parts = []
        for line in self._handle:
            try:
                idx = line.index('!')
            except ValueError:
                # There is no comment...
                idx = None
                end = ''
                self.comment = ''
                if line.rstrip('\r\n').endswith('-'):
                    # Continuation
                    parts.append(line.rstrip('\r\n')[:-1]) # Skip the continuation character
                    continue
            else:
                # Lines with no comment cannot continue
                end = '\n'
                self.comment = line[idx:].rstrip()
            parts.append(line[:idx] + end)
            yield ' '.join(parts)
            # Reset parts
            parts = []

    def readline(self):
        self.line_number += 1
        line = self._handle.readline()
        parts = []
        while line:
            if line.rstrip('\r\n').endswith('-'):
                # Continuation
                parts.append(line.rstrip('\r\n')[:-1]) # Skip the continuation character
                line = self._handle.readline()
                self.line_number += 1
            else:
                parts.append(line)
                break # done with this line
        line = ' '.join(parts)
        try:
            idx = line.index('!')
            self.comment = line[idx:].rstrip()
            end = '\n'
        except ValueError:
            idx = None
            end = ''
            self.comment = ''
        return line[:idx] + end

    def readlines(self):
        return [line for line in self]

    def read(self):
        return ''.join(self.readlines())

    def close(self):
        self._handle.close()
        self.closed = True

    def rewind(self):
        """ Return to the beginning of the file """
        self._handle.seek(0)

    def __del__(self):
        try:
            self.closed or self._handle.close()
        except AttributeError:
            # It didn't make it out of the constructor
            pass

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class CharmmStreamFile(object):
    """
    The stream file is broken down into sections of commands delimited by the
    strings:
        read <section> <options>
        ....
        ....
        end
    This object provides iterators over those sections and a file-like API for
    dealing with the text.

    """
    def __init__(self, fname):
        self.lines = []
        self.comments = []
        with CharmmFile(fname, 'r') as f:
            for line in f:
                self.lines.append(line)
                self.comments.append(f.comment)
        self.line_number = 0

    def __iter__(self):
        return iter(self.lines)

    def rewind(self):
        """ Return to the beginning of the file """
        self.line_number = 0

    def next_section(self):
        """
        Fast-forwards the file to the next CHARMM command section

        Returns
        -------
        name, data, comments : str, list of str, list of str
            name is the line defining the section that's being returned, whereas
            data is a list of all lines in the section, and comments is a list
            of all comments (same size as data) for all those lines

        Notes
        -----
        The line pointer will be set to the line defining the section
        """
        lines = []
        comments = []
        while self.line_number < len(self.lines):
            line = self.lines[self.line_number].strip()
            comment = self.comments[self.line_number].strip()
            if line[:4].lower() == 'read':
                title = line.strip()
                self.line_number += 1
                line = self.lines[self.line_number]
                while line and not line.strip().lower().startswith('end'):
                    lines.append(line)
                    comments.append(comment)
                    self.line_number += 1
                    line = self.lines[self.line_number]
                    comment = self.comments[self.line_number]
                if line[:3].upper() == 'END':
                    lines.append(line)
                    comments.append(comment)
                return title, lines, comments
            self.line_number += 1
        # No sections left
        return None, None, None

    def __del__(self): 
        pass
