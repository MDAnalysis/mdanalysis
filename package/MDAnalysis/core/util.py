# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

"""
Helper functions --- :mod:`MDAnalysis.core.util`
====================================================

Small helper functions that don't fit anywhere else.

Files and directories
---------------------

.. autofunction:: filename
.. function:: openany(directory[,mode='r'])

   Context manager to open a compressed (bzip2, gzip) or plain file
   (uses :func:`anyopen`).

.. autofunction:: anyopen

.. autofunction:: greedy_splitext

.. autofunction:: which

.. autofunction:: realpath

Containers and lists
--------------------

.. autofunction:: iterable
.. autofunction:: asiterable


File parsing
------------

.. autoclass:: FORTRANReader
   :members:
.. autodata:: FORTRAN_format_regex


Data manipulation and handling
------------------------------

.. autofunction:: fixedwidth_bins


Strings
-------

.. autofunction:: convert_aa_code
.. autofunction:: parse_residue


Mathematics and Geometry
------------------------

.. autofunction:: normal
.. autofunction:: norm
.. autofunction:: angle
.. autofunction:: dihedral
.. autofunction:: stp

"""
from __future__ import with_statement

__docformat__ = "restructuredtext en"

import os.path
from contextlib import contextmanager
import bz2, gzip
import re

import numpy

def filename(name,ext=None,keep=False):
    """Return a new name that has suffix attached; replaces other extensions.

    :Arguments:
      *name*
           filename; extension is replaced unless keep=True
      *ext*
           extension
      *keep*
           ``False``: replace existing extension; ``True``: keep if exists
    """
    name = str(name)
    if ext is None:
        return name
    if not ext.startswith(os.path.extsep):
        ext = os.path.extsep + ext
    #if name.find(ext) > 0:    # normally >= 0 but if we start with '.' then we keep it
    #    return name
    root, origext = os.path.splitext(name)
    if len(origext) == 0 or not keep:
        return root + ext
    return name

@contextmanager
def openany(datasource, mode='r'):
    """Open the datasource and close it when the context exits."""
    stream, filename = anyopen(datasource, mode=mode)
    try:
        yield stream
    finally:
        stream.close()

def anyopen(datasource, mode='r'):
    """Open datasource (gzipped, bzipped, uncompressed) and return a stream.

    :Arguments:
     *datasource*
        a file or a stream
     *mode*
        'r' or 'w'
    """
    # TODO: - make this act as ContextManager (and not return filename)
    #       - need to add binary 'b' to mode for compressed files?

    handlers = {'bz2': bz2.BZ2File, 'gz': gzip.open, '': file}

    if mode.startswith('r'):
        if hasattr(datasource,'next') or hasattr(datasource,'readline'):
            stream = datasource
            filename = '(%s)' % stream.name  # maybe that does not always work?
        else:
            stream = None
            filename = datasource
            for ext in ('bz2', 'gz', ''):   # file == '' should be last
                openfunc = handlers[ext]
                stream = _get_stream(datasource, openfunc, mode=mode)
                if not stream is None:
                    break
            if stream is None:
                raise IOError("Cannot open %(filename)r in mode=%(mode)r." % vars())
    elif mode.startswith('w'):
        if hasattr(datasource, 'write'):
            stream = datasource
            filename = '(%s)' % stream.name  # maybe that does not always work?
        else:
            stream = None
            filename = datasource
            name, ext = os.path.splitext(filename)
            if ext.startswith('.'):
                ext = ext[1:]
            if not ext in ('bz2', 'gz'):
                ext = ''   # anything else but bz2 or gz is just a normal file
            openfunc = handlers[ext]
            stream = openfunc(datasource, mode=mode)
            if stream is None:
                raise IOError("Cannot open %(filename)r in mode=%(mode)r with type %(ext)r." % vars())
    else:
        raise NotImplementedError("Sorry, mode=%(mode)r is not implemented for %(datasource)r" % vars())

    return stream, filename

def _get_stream(filename, openfunction=file, mode='r'):
    try:
        stream = openfunction(filename, mode=mode)
    except IOError:
        return None

    try:
        stream.readline()
        stream.close()
        stream = openfunction(filename,'r')
    except IOError:
        stream.close()
        stream = None
    return stream

def greedy_splitext(p):
    """Split extension in path *p* at the left-most separator."""
    path, root = os.path.split(p)
    extension = ''
    while True:
        root, ext = os.path.splitext(root)
        extension = ext + extension
        if not ext:
            break
    return root, extension

def which(program):
    """Determine full path of executable *program* on :envvar:`PATH`.

    (Jay at http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python)
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        real_program = realpath(program)
        if is_exe(real_program):
            return real_program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def realpath(*args):
    """Join all args and return the real path, rooted at /.

    Expands '~', '~user', and environment variables such as $HOME.

    Returns ``None`` if any of the args is ``None``.
    """
    if None in args:
        return None
    return os.path.realpath(os.path.expanduser(os.path.expandvars(os.path.join(*args))))


def iterable(obj):
    """Returns ``True`` if *obj* can be iterated over and is *not* a  string."""
    if isinstance(obj, basestring):
        return False    # avoid iterating over characters of a string

    if hasattr(obj, 'next'):
        return True    # any iterator will do
    try:
        len(obj)       # anything else that might work
    except TypeError:
        return False
    return True

def asiterable(obj):
    """Returns obj so that it can be iterated over; a string is *not* treated as iterable"""
    if not iterable(obj):
        obj = [obj]
    return obj

#: Regular expresssion (see :mod:`re`) to parse a simple `FORTRAN edit descriptor`_.
#: ``(?P<repeat>\d?)(?P<format>[IFELAX])(?P<numfmt>(?P<length>\d+)(\.(?P<decimals>\d+))?)?``
#:
#: .. _FORTRAN edit descriptor: http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html
FORTRAN_format_regex = "(?P<repeat>\d+?)(?P<format>[IFEAX])(?P<numfmt>(?P<length>\d+)(\.(?P<decimals>\d+))?)?"
_FORTRAN_format_pattern = re.compile(FORTRAN_format_regex)

def strip(s):
    """Convert *s* to a string and return it white-space stripped."""
    return str(s).strip()

class FixedcolumnEntry(object):
    """Represent an entry at specific fixed columns.

    Reads from line[start:stop] and converts according to
    typespecifier.
    """
    convertors = {'I': int, 'F': float, 'E': float, 'A': strip}
    def __init__(self, start, stop, typespecifier):
        """
        :Arguments:
         *start*
            first column
         *stop*
            last column + 1
         *typespecifier*
            'I': int, 'F': float, 'E': float, 'A': stripped string

        The start/stop arguments follow standard Python convention in that
        they are 0-based and that the *stop* argument is not included.
        """
        self.start = start
        self.stop = stop
        self.typespecifier = typespecifier
        self.convertor = self.convertors[typespecifier]
    def read(self, line):
        """Read the entry from *line* and convert to appropriate type."""
        try:
            return self.convertor(line[self.start:self.stop])
        except ValueError:
            raise ValueError("%r: Failed to read&convert %r" % (self, line[self.start:self.stop]))
    def __len__(self):
        """Length of the field in columns (stop - start)"""
        return self.stop - self.start
    def __repr__(self):
        return "FixedcolumnEntry(%d,%d,%r)" % (self.start, self.stop, self.typespecifier)

class FORTRANReader(object):
    """FORTRANReader provides a method to parse FORTRAN formatted lines in a file.

    Usage::

       atomformat = FORTRANReader('2I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A8,F20.10')
       for line in open('coordinates.crd'):
           serial,TotRes,resName,name,x,y,z,chainID,resSeq,tempFactor = atomformat.read(line)

    Fortran format edit descriptors; see `Fortran Formats`_ for the syntax.

    Only simple one-character specifiers supported here: *I F E A X* (see
    :data:`FORTRAN_format_regex`).

    Strings are stripped of leading and trailing white space.

    .. _`Fortran Formats`: http://www.webcitation.org/5xbaWMV2x
    .. _`Fortran Formats (URL)`:
       http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html
    """

    def __init__(self, fmt):
        """Set up the reader with the FORTRAN format string.

        The string *fmt* should look like '2I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A8,F20.10'.
        """
        self.fmt = fmt.split(',')
        descriptors = [self.parse_FORTRAN_format(descriptor) for descriptor in self.fmt]
        start = 0
        self.entries = []
        for d in descriptors:
            if d['format'] != 'X':
                for x in range(d['repeat']):
                    stop = start + d['length']
                    self.entries.append(FixedcolumnEntry(start,stop,d['format']))
                    start = stop
            else:
                start += d['totallength']

    def read(self, line):
        """Parse *line* according to the format string and return list of values.

        Values are converted to Python types according to the format specifier.

        :Returns: list of entries with appropriate types
        :Raises: :exc:`ValueError` if any of the conversions cannot be made
                 (e.g. space for an int)

        .. SeeAlso:: :meth:`FORTRANReader.number_of_matches`
        """
        return [e.read(line) for e in self.entries]

    def number_of_matches(self, line):
        """Return how many format entries could be populated with legal values."""
        # not optimal, I suppose...
        matches = 0
        for e in self.entries:
            try:
                e.read(line)
                matches += 1
            except ValueError:
                pass
        return matches

    def parse_FORTRAN_format(self, edit_descriptor):
        """Parse the descriptor.

          parse_FORTRAN_format(edit_descriptor) --> dict

        :Returns: dict with totallength (in chars), repeat, length,
                  format, decimals
        :Raises: :exc:`ValueError` if the *edit_descriptor* is not recognized
                 and cannot be parsed

        .. Note::

           Specifiers: *L ES EN T TL TR / r S SP SS BN BZ* are *not*
           supported, and neither are the scientific notation *Ew.dEe*
           forms.
        """

        m = _FORTRAN_format_pattern.match(edit_descriptor.upper())
        if m is None:
            try:
                m = _FORTRAN_format_pattern.match("1"+edit_descriptor.upper())
                if m is None:
                    raise ValueError  # really no idea what the descriptor is supposed to mean
            except:
                raise ValueError("unrecognized FORTRAN format %r" % edit_descriptor)
        d = m.groupdict()
        if d['repeat'] == '':
            d['repeat'] = 1
        if d['format'] == 'X':
            d['length'] = 1
        for k in ('repeat','length','decimals'):
            try:
                d[k] = int(d[k])
            except ValueError:   # catches ''
                d[k] = 0
            except TypeError:    # keep None
                pass
        d['totallength'] = d['repeat'] * d['length']
        return d

    def __len__(self):
        """Returns number of entries."""
        return len(self.entries)

    def __repr__(self):
        return self.__class__.__name__+"("+",".join(self.fmt)+")"

def fixedwidth_bins(delta,xmin,xmax):
    """Return bins of width delta that cover xmin,xmax (or a larger range).

    dict = fixedwidth_bins(delta,xmin,xmax)

    The dict contains 'Nbins', 'delta', 'min', and 'max'.
    """
    if not numpy.all(xmin < xmax):
        raise ValueError('Boundaries are not sane: should be xmin < xmax.')
    _delta = numpy.asarray(delta,dtype=numpy.float_)
    _xmin = numpy.asarray(xmin,dtype=numpy.float_)
    _xmax = numpy.asarray(xmax,dtype=numpy.float_)
    _length = _xmax - _xmin
    N = numpy.ceil(_length/_delta).astype(numpy.int_)      # number of bins
    dx = 0.5 * (N*_delta - _length)   # add half of the excess to each end
    return {'Nbins':N, 'delta':_delta,'min':_xmin-dx, 'max':_xmax+dx}


# geometric functions
def norm(v):
    """Returns the length of a vector, ``sqrt(v.v)``.

    Faster than :func:`numpy.linalg.norm` because no frills.
    """
    return numpy.sqrt(numpy.dot(v,v))

def normal(vec1, vec2):
    """Returns the unit vector normal to two vectors."""
    normal = numpy.cross(vec1, vec2)
    n = norm(normal)
    if n == 0.0:
        return normal  # returns [0,0,0] instead of [nan,nan,nan]
    return normal/n    # ... could also use numpy.nan_to_num(normal/norm(normal))

def angle(a, b):
    """Returns the angle between two vectors in radians"""
    x = numpy.dot(a, b) / (norm(a)*norm(b))
    # catch roundoffs that lead to nan otherwise
    if x > 1.0:
        return 0.0
    elif x < -1.0:
        return -numpy.pi
    return numpy.arccos(x)

def stp(vec1, vec2, vec3):
    """Takes the scalar triple product of three vectors.

    :Returns: v3 . (v1 x v2)
    """
    return numpy.dot(vec3, numpy.cross(vec1, vec2))

def dihedral(ab, bc, cd):
    """Returns the dihedral angle in radians between vectors connecting A,B,C,D.

    The dihedral measures the rotation around bc::

         ab
       A---->B
              \ bc
              _\'
                C---->D
                  cd

    The dihedral angle is -pi <= x <= pi.

    .. versionadded:: 0.7.8
    """
    x = angle(normal(ab, bc), normal(bc, cd))
    return (x if stp(ab, bc, cd) <= 0.0 else -x)


# String functions
# ----------------

# TODO: make it work for non-default charge state amino acids.
#: translation table for 1-letter codes --> 3-letter codes
#: .. Note: This does not work for HISB and non-default charge state aa!
amino_acid_codes = {'A':'ALA', 'C':'CYS', 'D':'ASP', 'E':'GLU',
                    'F':'PHE', 'G':'GLY', 'H':'HIS', 'I':'ILE',
                    'K':'LYS', 'L':'LEU', 'M':'MET', 'N':'ASN',
                    'P':'PRO', 'Q':'GLN', 'R':'ARG', 'S':'SER',
                    'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR'}
inverse_aa_codes = dict([(three, one) for one,three in amino_acid_codes.items()])

# from GromacsWrapper utilities
def convert_aa_code(x):
    """Converts between 3-letter and 1-letter amino acid codes."""
    if len(x) == 1:
        return amino_acid_codes[x.upper()]
    elif len(x) == 3:
        return inverse_aa_codes[x.upper()]
    else:
        raise ValueError("Can only convert 1-letter or 3-letter amino acid codes, "
                         "not %r" % x)


#: Regular expression to match and parse a residue-atom selection; will match
#: "LYS300:HZ1" or "K300:HZ1" or "K300" or "4GB300:H6O" or "4GB300" or "YaA300".
RESIDUE = re.compile("""
                 (?P<aa>([ACDEFGHIKLMNPQRSTVWY])   # 1-letter amino acid
                        |                          #   or
                        ([0-9A-Z][a-zA-Z][A-Z][A-Z]?)    # 3-letter or 4-letter residue name
                 )
                 \s*                               # white space allowed
                 (?P<resid>\d+)                    # resid
                 \s*
                 (:                                # separator ':'
                   \s*
                   (?P<atom>\w+)                   # atom name
                 )?                                # possibly one
            """, re.VERBOSE | re.IGNORECASE)

# from GromacsWrapper cbook.IndexBuilder
def parse_residue(residue):
    """Process residue string.

    Examples:
     - "LYS300:HZ1" --> ("LYS", 300, "HZ1")
     - "K300:HZ1" --> ("LYS", 300, "HZ1")
     - "K300" --> ("LYS", 300, None)
     - "4GB300:H6O" --> ("4GB", 300, "H6O")
     - "4GB300" --> ("4GB", 300, None)

    :Argument: The *residue* must contain a 1-letter or 3-letter or
               4-letter residue string, a number (the resid) and
               optionally an atom identifier, which must be separate
               from the residue with a colon (":"). White space is
               allowed in between.

    :Returns: `(3-letter aa string, resid, atomname)`; known 1-letter
              aa codes are converted to 3-letter codes
    """

    # XXX: use _translate_residue() ....
    m = RESIDUE.match(residue)
    if not m:
        raise ValueError("Selection %(residue)r is not valid (only 1/3/4 letter resnames, resid required)." % vars())
    resid = int(m.group('resid'))
    residue = m.group('aa')
    if len(residue) == 1:
        resname = convert_aa_code(residue) # only works for AA
    else:
        resname = residue                            # use 3-letter for any resname
    atomname = m.group('atom')
    return (resname, resid, atomname)
