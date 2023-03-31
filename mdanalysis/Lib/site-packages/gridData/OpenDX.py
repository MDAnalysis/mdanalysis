# gridData --- python modules to read and write gridded data
# Copyright (c) 2009-2014 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Lesser General Public License, version 3 or later.

r"""
:mod:`~gridData.OpenDX` --- routines to read and write simple OpenDX files
==========================================================================

The OpenDX format for multi-dimensional grid data. OpenDX is a free
visualization software, see http://www.opendx.org.

.. Note:: This module only implements a primitive subset, sufficient
          to represent n-dimensional regular grids.

The OpenDX scalar file format is specified in Appendix `B.2 Data
Explorer Native Files`_ [#OpenDXformat]_.

If you want to build a dx object from your data you can either use the
convenient :class:`~gridData.core.Grid` class from the top level
module (:class:`gridData.Grid`) or see the lower-level methods
described below.


.. _opendx-read-write:

Reading and writing OpenDX files
--------------------------------

If you have OpenDX files from other software and you just want to
**read** it into a Python array then you do not really need to use the
interface in :mod:`gridData.OpenDX`: just use
:class:`~gridData.core.Grid` and load the file::

  from gridData import Grid
  g = Grid("data.dx")

This should work for files produced by common visualization programs
(VMD_, PyMOL_, Chimera_). The documentation for :mod:`gridData` tells
you more about what to do with the :class:`~gridData.core.Grid`
object.

If you want to **write** an OpenDX file then you just use the
:meth:`gridData.core.Grid.export` method with `file_format="dx"` (or
just use a filename with extension ".dx")::

  g.export("data.dx")

However, some visualization programs do not implement full OpenDX
specifications and only read very specific, "OpenDX-like"
files. :mod:`gridData.OpenDX` tries to be compatible with these
formats. However, sometimes additional help is needed to write an
OpenDX file that can be read by a specific software, as described
below:

Known issues for writing OpenDX files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* APBS require the delta to be written to the seventh significant figure.
  The delta is now written to reflect this increase in precision.

  .. versionchanged:: 0.6.0

* PyMOL_ requires OpenDX files with the type specification "double" in
  the `class array` section (see issue `#35`_). By default (since
  release 0.4.0), the type is set to the one that most closely
  approximates the dtype of the numpy array :attr:`Grid.grid`, which
  holds all data. This is often :class:`numpy.float64`, which will
  create an OpenDX type "double", which PyMOL will read.

  However, if you want to *force* a specific OpenDX type (such as
  "float" or "double", see :attr:`gridData.OpenDX.array.dx_types` for
  available values) then you can use the ``type`` keyword argument::

    g.export("for_pymol.dx", type="double")

  If you always want to be able to read OpenDX files with PyMOL, it is
  suggested to always export with ``type="double"``.

  .. versionadded:: 0.4.0



.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _PyMOL: http://www.pymol.org/
.. _Chimera: https://www.cgl.ucsf.edu/chimera/
.. _`#35`: https://github.com/MDAnalysis/GridDataFormats/issues/35




Building a dx object from a numpy array ``A``
---------------------------------------------

If you have a numpy array ``A`` that represents a density in cartesian
space then you can construct a dx object (named a *field* in OpenDX
parlance) if you provide some additional information that fixes the
coordinate system in space and defines the units along the axes.

The following data are required:

grid
    numpy nD array (typically a nD histogram)
grid.shape
    the shape of the array
origin
    the cartesian coordinates of the center of the (0,0,..,0) grid cell
delta
    :math:`n \times n` array with the length of a grid cell along
    each axis; for regular rectangular grids the off-diagonal
    elements are 0 and the diagonal ones correspond to the
    'bin width' of the histogram, eg ``delta[0,0] = 1.0`` (Angstrom)

The DX data type ("type" in the DX file) is determined from the
:class:`numpy.dtype` of the :class:`numpy.ndarray` that is provided as
the *grid* (or with the *type* keyword argument to
:class:`gridData.OpenDX.array`).

For example, to build a :class:`field`::

  dx = OpenDX.field('density')
  dx.add('positions', OpenDX.gridpositions(1, grid.shape, origin, delta))
  dx.add('connections', OpenDX.gridconnections(2, grid.shape))
  dx.add('data', OpenDX.array(3, grid))

or all with the constructor::

  dx = OpenDX.field('density', components=dict(
            positions=OpenDX.gridpositions(1,grid.shape, d.origin, d.delta),
            connections=OpenDX.gridconnections(2, grid.shape),
            data=OpenDX.array(3, grid)))


Building a dx object from a dx file
-----------------------------------

One can also read data from an existing dx file::

 dx = OpenDX.field(0)
 dx.read('file.dx')

Only simple arrays are read and initially stored as a 1-d
:class:`numpy.ndarray` in the `dx.components['data'].array` with the
:class:`numpy.dtype` determined by the DX type in the file.

The dx :class:`field` object has a method
:meth:`~OpenDX.field.histogramdd` that produces output identical to
the :func:`numpy.histogramdd` function by taking the stored dimension
and deltas into account. In this way, one can store nD histograms in a
portable and universal manner::

  histogram, edges = dx.histogramdd()

.. rubric:; Footnotes

.. [#OpenDXformat] The original link to the OpenDX file format specs
   http://opendx.sdsc.edu/docs/html/pages/usrgu068.htm#HDREDF is dead so I am linking
   to an archived copy at the Internet Archive , `B.2 Data Explorer Native Files`_.

.. _`B.2 Data Explorer Native Files`:
   https://web.archive.org/web/20080808140524/http://opendx.sdsc.edu/docs/html/pages/usrgu068.htm
.. http://opendx.sdsc.edu/docs/html/pages/usrgu068.htm#HDREDF

Classes and functions
---------------------

"""
import numpy
import re
import gzip

import warnings

# Python 2/3 compatibility (see issue #99)
# and https://bugs.python.org/issue30012
import sys
if sys.version_info >= (3, ):
    def _gzip_open(filename, mode="rt"):
        return gzip.open(filename, mode)
else:
    def _gzip_open(filename, mode="rt"):
        return gzip.open(filename)
del sys

class DXclass(object):
    """'class' object as defined by OpenDX"""
    def __init__(self,classid):
        """id is the object number"""
        self.id = classid  # serial number of the object
        self.name = None   # name of the DXclass
        self.component = None   # component type
        self.D = None      # dimensions

    def write(self, stream, optstring="", quote=False):
        """write the 'object' line; additional args are packed in string"""
        classid = str(self.id)
        if quote: classid = '"'+classid+'"'
        # Only use a *single* space between tokens; both chimera's and pymol's DX parser
        # does not properly implement the OpenDX specs and produces garbage with multiple
        # spaces. (Chimera 1.4.1, PyMOL 1.3)
        to_write = 'object '+classid+' class '+str(self.name)+' '+optstring+'\n'
        self._write_line(stream, to_write)

    @staticmethod
    def _write_line(stream, line="", quote=False):
        """write a line to the file"""
        if isinstance(stream, gzip.GzipFile):
            line = line.encode()
        stream.write(line)

    def read(self, stream):
        raise NotImplementedError('Reading is currently not supported.')

    def ndformat(self,s):
        """Returns a string with as many repetitions of s as self
        has dimensions (derived from shape)"""
        return s * len(self.shape)

    def __repr__(self):
        return '<OpenDX.'+str(self.name)+' object, id='+str(self.id)+'>'


class gridpositions(DXclass):
    """OpenDX gridpositions class.

    shape     D-tuplet describing size in each dimension
    origin    coordinates of the centre of the grid cell with index 0,0,...,0
    delta     DxD array describing the deltas
    """
    def __init__(self,classid,shape=None,origin=None,delta=None,**kwargs):
        if shape is None or origin is None or delta is None:
            raise ValueError('all keyword arguments are required')
        self.id = classid
        self.name = 'gridpositions'
        self.component = 'positions'
        self.shape = numpy.asarray(shape)      # D dimensional shape
        self.origin = numpy.asarray(origin)    # D vector
        self.rank = len(self.shape)            # D === rank

        self.delta = numpy.asarray(delta)      # DxD array of grid spacings
        # gridDataFormats  actually provides a simple 1D array with the deltas because only
        # regular grids are used but the following is a reminder that OpenDX should be able
        # to handle more complicated volume elements
        if len(self.delta.shape) == 1:
            self.delta = numpy.diag(delta)
        if self.delta.shape != (self.rank, self.rank):
            # check OpenDX specs for irreg spacing if we want to implement
            # anything more complicated
            raise NotImplementedError('Only regularly spaced grids allowed, '
                                      'not delta={}'.format(self.delta))
    def write(self, stream):
        super(gridpositions, self).write(
            stream, ('counts '+self.ndformat(' %d')) % tuple(self.shape))
        self._write_line(stream, 'origin %f %f %f\n' % tuple(self.origin))
        for delta in self.delta:
            self._write_line(
                stream, ('delta ' +
                         self.ndformat(' {:.7g}').format(*delta) +
                         '\n'))

    def edges(self):
        """Edges of the grid cells, origin at centre of 0,0,..,0 grid cell.

        Only works for regular, orthonormal grids.
        """
        return [self.delta[d,d] * numpy.arange(self.shape[d]+1) + self.origin[d]\
                - 0.5*self.delta[d,d]     for d in range(self.rank)]


class gridconnections(DXclass):
    """OpenDX gridconnections class"""
    def __init__(self,classid,shape=None,**kwargs):
        if shape is None:
            raise ValueError('all keyword arguments are required')
        self.id = classid
        self.name = 'gridconnections'
        self.component = 'connections'
        self.shape = numpy.asarray(shape)      # D dimensional shape

    def write(self, stream):
        super(gridconnections, self).write(
            stream, ('counts '+self.ndformat(' %d')) % tuple(self.shape))


class array(DXclass):
    """OpenDX array class.

    See `Array Objects`_ for details.

    .. _Array Objects:
       https://web.archive.org/web/20080808140524/http://opendx.sdsc.edu/docs/html/pages/usrgu068.htm#Header_440
    """
    #: conversion from :attr:`numpy.dtype.name` to closest OpenDX array type
    #: (round-tripping is not guaranteed to produce identical types); not all
    #: types are supported (e.g., strings are missing)
    np_types = {
        "uint8": "byte",         # DX "unsigned byte" equivalent
        "int8": "signed byte",
        "uint16": "unsigned short",
        "int16": "short",         # DX "signed short" equivalent
        "uint32": "unsigned int",
        "int32": "int",           # DX "signed int"   equivalent
        "uint64": "unsigned int", # not explicit in DX, for compatibility
        "int64": "int",           # not explicit in DX, for compatibility
        # "hyper",                # ?
        "float32": "float",       # default
        "float64": "double",
        "float16": "float",       # float16 not available in DX, use float
        # numpy "float128 not available, raise error
        # "string" not automatically supported
    }
    #: conversion from OpenDX type to closest :class:`numpy.dtype`
    #: (round-tripping is not guaranteed to produce identical types); not all
    #: types are supported (e.g., strings and conversion to int64 are missing)
    dx_types = {
        "byte": "uint8",
        "unsigned byte": "uint8",
        "signed byte": "int8",
        "unsigned short": "uint16",
        "short": "int16",
        "signed short": "int16",
        "unsigned int": "uint32",
        "int": "int32",
        "signed int": "int32",
        # "hyper",                # ?
        "float": "float32",       # default
        "double": "float64",
        # "string" not automatically supported
    }

    def __init__(self, classid, array=None, type=None, typequote='"',
                 **kwargs):
        """
        Parameters
        ----------
        classid : int
        array : array_like
        type : str (optional)
             Set the DX type in the output file and cast `array` to
             the closest numpy dtype.  `type` must be one of the
             allowed types in DX files as defined under `Array
             Objects`_.  The default ``None`` tries to set the type
             from the :class:`numpy.dtype` of `array`.

             .. versionadded:: 0.4.0

        Raises
        ------
        ValueError
             if `array` is not provided; or if `type` is not of the correct
             DX type
        """
        if array is None:
            raise ValueError('array keyword argument is required')
        self.id = classid
        self.name = 'array'
        self.component = 'data'
        # detect type https://github.com/MDAnalysis/GridDataFormats/issues/35
        if type is None:
            self.array = numpy.asarray(array)
            try:
                self.type = self.np_types[self.array.dtype.name]
            except KeyError:
                warnings.warn(("array dtype.name = {0} can not be automatically "
                               "converted to a DX array type. Use the 'type' keyword "
                               "to manually specify the correct type.").format(
                                   self.array.dtype.name))
                self.type = self.array.dtype.name  # will raise ValueError on writing
        else:
            try:
                self.array = numpy.asarray(array, dtype=self.dx_types[type])
            except KeyError:
                raise ValueError(("DX type {0} cannot be converted to an "
                                  "appropriate numpy dtype. Available "
                                  "types are: {1}".format(type,
                                                          list(self.dx_types.values()))))
            self.type = type
        self.typequote = typequote

    def write(self, stream):
        """Write the *class array* section.

        Parameters
        ----------
        stream : stream

        Raises
        ------
        ValueError
             If the `dxtype` is not a valid type, :exc:`ValueError` is
             raised.

        """
        if self.type not in self.dx_types:
            raise ValueError(("DX type {} is not supported in the DX format. \n"
                              "Supported valus are: {}\n"
                              "Use the type=<type> keyword argument.").format(
                                  self.type, list(self.dx_types.keys())))
        typelabel = (self.typequote+self.type+self.typequote)
        super(array, self).write(stream, 'type {0} rank 0 items {1} data follows'.format(
            typelabel, self.array.size))

        # grid data, serialized as a C array (z fastest varying)
        # (flat iterator is equivalent to: for x: for y: for z: grid[x,y,z])
        # VMD's DX reader requires exactly 3 values per line
        fmt_string = "{:d}"
        if (self.array.dtype.kind == 'f' or self.array.dtype.kind == 'c'):
            precision = numpy.finfo(self.array.dtype).precision
            fmt_string = "{:."+"{:d}".format(precision)+"f}"
        values_per_line = 3
        values = self.array.flat
        while 1:
            try:
                for i in range(values_per_line):
                    self._write_line(stream, fmt_string.format(next(values)) + "\t")
                self._write_line(stream, '\n')
            except StopIteration:
                self._write_line(stream, '\n')
                break
        self._write_line(stream, 'attribute "dep" string "positions"\n')

class field(DXclass):
    """OpenDX container class

    The *field* is the top-level object and represents the whole
    OpenDX file. It contains a number of other objects.

    Instantiate a DX object from this class and add subclasses with
    :meth:`add`.

    """
    # perhaps this should not derive from DXclass as those are
    # objects in field but a field cannot contain itself
    def __init__(self,classid='0',components=None,comments=None):
        """OpenDX object, which is build from a list of components.

        Parameters
        ----------

        id : str
               arbitrary string
        components : dict
               dictionary of DXclass instances (no sanity check on the
               individual ids!) which correspond to

               * positions
               * connections
               * data

        comments : list
               list of strings; each string becomes a comment line
               prefixed with '#'. Avoid newlines.


        A field must have at least the components 'positions',
        'connections', and 'data'. Those components are associated
        with objects belonging to the field. When writing a dx file
        from the field, only the required objects are dumped to the file.

        (For a more general class that can use field:
        Because there could be more objects than components, we keep a
        separate object list. When dumping the dx file, first all
        objects are written and then the field object describes its
        components. Objects are referenced by their unique id.)

        .. Note:: uniqueness of the *id* is not checked.


        Example
        -------
        Create a new dx object::

           dx = OpenDX.field('density',[gridpoints,gridconnections,array])

        """
        if components is None:
            components = dict(positions=None,connections=None,data=None)
        if comments is None:
            comments = ['OpenDX written by gridData.OpenDX',
                        'from https://github.com/MDAnalysis/GridDataFormats']
        elif type(comments) is not list:
            comments = [str(comments)]
        self.id = classid       # can be an arbitrary string
        self.name = 'field'
        self.component = None   # cannot be a component of a field
        self.components = components
        self.comments= comments

    def _openfile_writing(self, filename):
        """Returns a regular or gz file stream for writing"""
        if filename.endswith('.gz'):
            return gzip.open(filename, 'wb')
        else:
            return open(filename, 'w')

    def write(self, filename):
        """Write the complete dx object to the file.

        This is the simple OpenDX format which includes the data into
        the header via the 'object array ... data follows' statement.

        Only simple regular arrays are supported.

        The format should be compatible with VMD's dx reader plugin.
        """
        # comments (VMD chokes on lines of len > 80, so truncate)
        maxcol = 80
        with self._openfile_writing(str(filename)) as outfile:
            for line in self.comments:
                comment = '# '+str(line)
                self._write_line(outfile, comment[:maxcol]+'\n')
            # each individual object
            for component, object in self.sorted_components():
                object.write(outfile)
            # the field object itself
            super(field, self).write(outfile, quote=True)
            for component, object in self.sorted_components():
                self._write_line(outfile, 'component "%s" value %s\n' % (
                    component, str(object.id)))

    def read(self, stream):
        """Read DX field from file.

            dx = OpenDX.field.read(dxfile)

        The classid is discarded and replaced with the one from the file.
        """
        DXfield = self
        p = DXParser(stream)
        p.parse(DXfield)

    def add(self,component,DXobj):
        """add a component to the field"""
        self[component] = DXobj

    def add_comment(self,comment):
        """add comments"""
        self.comments.append(comment)

    def sorted_components(self):
        """iterator that returns (component,object) in id order"""
        for component, object in \
                sorted(self.components.items(),
                       key=lambda comp_obj: comp_obj[1].id):
            yield component, object

    def histogramdd(self):
        """Return array data as (edges,grid), i.e. a numpy nD histogram."""
        shape = self.components['positions'].shape
        edges = self.components['positions'].edges()
        hist = self.components['data'].array.reshape(shape)
        return (hist,edges)

    def __getitem__(self,key):
        return self.components[key]

    def __setitem__(self,key,value):
        self.components[key] = value

    def __repr__(self):
        return '<OpenDX.field object, id='+str(self.id)+', with '+\
               str(len(self.components))+' components and '+\
               str(len(self.components))+' objects>'


#------------------------------------------------------------
# DX file parsing
#------------------------------------------------------------

class DXParseError(Exception):
    """general exception for parsing errors in DX files"""
    pass
class DXParserNoTokens(DXParseError):
    """raised when the token buffer is exhausted"""
    pass

class Token:
    # token categories (values of dx_regex must match up with these categories)
    category = {'COMMENT': ['COMMENT'],
                'WORD': ['WORD'],
                'STRING': ['QUOTEDSTRING','BARESTRING','STRING'],
                'WHITESPACE': ['WHITESPACE'],
                'INTEGER': ['INTEGER'],
                'REAL': ['REAL'],
                'NUMBER': ['INTEGER','REAL']}
    # cast functions
    cast = {'COMMENT': lambda s:re.sub(r'#\s*','',s),
            'WORD': str,
            'STRING': str, 'QUOTEDSTRING': str, 'BARESTRING': str,
            'WHITESPACE': None,
            'NUMBER': float, 'INTEGER': int, 'REAL': float}

    def __init__(self,code,text):
        self.code = code    # store raw code
        self.text = text
    def equals(self,v):
        return self.text == v
    def iscode(self,code):
        return self.code in self.category[code]  # use many -> 1 mappings
    def value(self,ascode=None):
        """Return text cast to the correct type or the selected type"""
        if ascode is None:
            ascode = self.code
        return self.cast[ascode](self.text)
    def __repr__(self):
        return '<token '+str(self.code)+','+str(self.value())+'>'

class DXInitObject(object):
    """Storage class that holds data to initialize one of the 'real'
    classes such as OpenDX.array, OpenDX.gridconnections, ...

    All variables are stored in args which will be turned into the
    arguments for the DX class.
    """
    DXclasses = {'gridpositions':gridpositions,
                 'gridconnections':gridconnections,
                 'array':array, 'field':field,
                 }

    def __init__(self,classtype,classid):
        self.type = classtype
        self.id = classid
        self.args = dict()
    def initialize(self):
        """Initialize the corresponding DXclass from the data.

        class = DXInitObject.initialize()
        """
        return self.DXclasses[self.type](self.id,**self.args)
    def __getitem__(self,k):
        return self.args[k]
    def __setitem__(self,k,v):
        self.args[k] = v
    def __repr__(self):
        return '<DXInitObject instance type='+str(self.type)+', id='+str(self.id)+'>'

class DXParser(object):
    """Brain-dead baroque implementation to read a simple (VMD) dx file.

    Requires a OpenDX.field instance.

    1) scan for 'object' lines:
       'object' id 'class' class  [data]
       [data ...]
    2) parse data according to class
    3) construct dx field from classes
    """

    # the regexes must match with the categories defined in the Token class
    # REAL regular expression will catch both integers and floats.
    # Taken from
    # https://docs.python.org/3/library/re.html#simulating-scanf
    dx_regex = re.compile(r"""
    (?P<COMMENT>\#.*$)            # comment (until end of line)
    |(?P<WORD>(object|class|counts|origin|delta|type|counts|rank|items|data))
    |"(?P<QUOTEDSTRING>[^\"]*)"   # string in double quotes  (quotes removed)
    |(?P<WHITESPACE>\s+)          # white space
    |(?P<REAL>[-+]?               # true real number (decimal point or
    (\d+(\.\d*)?|\.\d+)           # scientific notation) and integers
    ([eE][-+]?\d+)?)
    |(?P<BARESTRING>[a-zA-Z_][^\s\#\"]+) # unquoted strings, starting with non-numeric
    """, re.VERBOSE)


    def __init__(self, filename):
        """Setup a parser for a simple DX file (from VMD)

        >>> DXfield_object = OpenDX.field(id)
        >>> p = DXparser('bulk.dx')
        >>> p.parse(DXfield_object)

        The field object will be completely rewritten (including the
        id if one is found in the input file. The input files
        component layout is currently ignored.

        Note that quotes are removed from quoted strings.
        """
        self.filename = str(filename)
        self.field = field('grid data',comments=['filename: {0}'.format(self.filename)])
        # other variables are initialised every time parse() is called

        self.parsers = {'general':self.__general,
                        'comment':self.__comment, 'object':self.__object,
                        'gridpositions':self.__gridpositions,
                        'gridconnections':self.__gridconnections,
                        'array':self.__array, 'field':self.__field,
                        }


    def parse(self, DXfield):
        """Parse the dx file and construct a DX field object with component classes.

        A :class:`field` instance *DXfield* must be provided to be
        filled by the parser::

           DXfield_object = OpenDX.field(*args)
           parse(DXfield_object)

        A tokenizer turns the dx file into a stream of tokens. A
        hierarchy of parsers examines the stream. The level-0 parser
        ('general') distinguishes comments and objects (level-1). The
        object parser calls level-3 parsers depending on the object
        found. The basic idea is that of a 'state machine'. There is
        one parser active at any time. The main loop is the general
        parser.

        * Constructing the dx objects with classtype and classid is
          not implemented yet.
        * Unknown tokens raise an exception.
        """

        self.DXfield = DXfield              # OpenDX.field (used by comment parser)
        self.currentobject = None           # containers for data
        self.objects = []                   # |
        self.tokens = []                    # token buffer

        if self.filename.endswith('.gz'):
            with _gzip_open(self.filename, 'rt') as self.dxfile:
                self.use_parser('general')
        else:
            with open(self.filename, 'r') as self.dxfile:
                self.use_parser('general')      # parse the whole file and populate self.objects

        # assemble field from objects
        for o in self.objects:
            if o.type == 'field':
                # Almost ignore the field object; VMD, for instance,
                # does not write components. To make this work
                # seamlessly I have to think harder how to organize
                # and use the data, eg preping the field object
                # properly and the initializing. Probably should also
                # check uniqueness of ids etc.
                DXfield.id = o.id
                continue
            c = o.initialize()
            self.DXfield.add(c.component,c)

        # free space
        del self.currentobject, self.objects



    def __general(self):
        """Level-0 parser and main loop.

        Look for a token that matches a level-1 parser and hand over control."""
        while 1:                            # main loop
            try:
                tok = self.__peek()         # only peek, apply_parser() will consume
            except DXParserNoTokens:
                # save previous DXInitObject
                # (kludge in here as the last level-2 parser usually does not return
                # via the object parser)
                if self.currentobject and self.currentobject not in self.objects:
                    self.objects.append(self.currentobject)
                return                      # stop parsing and finish
            # decision branches for all level-1 parsers:
            # (the only way to get out of the lower level parsers!)
            if tok.iscode('COMMENT'):
                self.set_parser('comment')  # switch the state
            elif tok.iscode('WORD') and tok.equals('object'):
                self.set_parser('object')   # switch the state
            elif self.__parser is self.__general:
                # Either a level-2 parser screwed up or some level-1
                # construct is not implemented.  (Note: this elif can
                # be only reached at the beginning or after comments;
                # later we never formally switch back to __general
                # (would create inifinite loop)
                raise DXParseError('Unknown level-1 construct at '+str(tok))

            self.apply_parser()     # hand over to new parser
                                    # (possibly been set further down the hierarchy!)

    # Level-1 parser
    def __comment(self):
        """Level-1 parser for comments.

        pattern: #.*
        Append comment (with initial '# ' stripped) to all comments.
        """
        tok = self.__consume()
        self.DXfield.add_comment(tok.value())
        self.set_parser('general')   # switch back to general parser

    def __object(self):
        """Level-1 parser for objects.

        pattern: 'object' id 'class' type ...

        id ::=   integer|string|'"'white space string'"'
        type ::= string
        """
        self.__consume()                    # 'object'
        classid = self.__consume().text
        word = self.__consume().text
        if word != "class":
            raise DXParseError("reserved word %s should have been 'class'." % word)
        # save previous DXInitObject
        if self.currentobject:
            self.objects.append(self.currentobject)
        # setup new DXInitObject
        classtype = self.__consume().text
        self.currentobject = DXInitObject(classtype=classtype,classid=classid)

        self.use_parser(classtype)

    # Level-2 parser (object parsers)
    def __gridpositions(self):
        """Level-2 parser for gridpositions.

        pattern:
        object 1 class gridpositions counts 97 93 99
        origin -46.5 -45.5 -48.5
        delta 1 0 0
        delta 0 1 0
        delta 0 0 1
        """
        try:
            tok = self.__consume()
        except DXParserNoTokens:
            return

        if tok.equals('counts'):
            shape = []
            try:
                while True:
                    # raises exception if not an int
                    self.__peek().value('INTEGER')
                    tok = self.__consume()
                    shape.append(tok.value('INTEGER'))
            except (DXParserNoTokens, ValueError):
                pass
            if len(shape) == 0:
                raise DXParseError('gridpositions: no shape parameters')
            self.currentobject['shape'] = shape
        elif tok.equals('origin'):
            origin = []
            try:
                while (self.__peek().iscode('INTEGER') or
                       self.__peek().iscode('REAL')):
                    tok = self.__consume()
                    origin.append(tok.value())
            except DXParserNoTokens:
                pass
            if len(origin) == 0:
                raise DXParseError('gridpositions: no origin parameters')
            self.currentobject['origin'] = origin
        elif tok.equals('delta'):
            d = []
            try:
                while (self.__peek().iscode('INTEGER') or
                       self.__peek().iscode('REAL')):
                    tok = self.__consume()
                    d.append(tok.value())
            except DXParserNoTokens:
                pass
            if len(d) == 0:
                raise DXParseError('gridpositions: missing delta parameters')
            try:
                self.currentobject['delta'].append(d)
            except KeyError:
                self.currentobject['delta'] = [d]
        else:
            raise DXParseError('gridpositions: '+str(tok)+' not recognized.')


    def __gridconnections(self):
        """Level-2 parser for gridconnections.

        pattern:
        object 2 class gridconnections counts 97 93 99
        """
        try:
            tok = self.__consume()
        except DXParserNoTokens:
            return

        if tok.equals('counts'):
            shape = []
            try:
                while True:
                    # raises exception if not an int
                    self.__peek().value('INTEGER')
                    tok = self.__consume()
                    shape.append(tok.value('INTEGER'))
            except (DXParserNoTokens, ValueError):
                pass
            if len(shape) == 0:
                raise DXParseError('gridconnections: no shape parameters')
            self.currentobject['shape'] = shape
        else:
            raise DXParseError('gridconnections: '+str(tok)+' not recognized.')


    def __array(self):
        """Level-2 parser for arrays.

        pattern:
        object 3 class array type double rank 0 items 12 data follows
        0 2 0
        0 0 3.6
        0 -2.0 1e-12
        +4.534e+01 .34534 0.43654
        attribute "dep" string "positions"
        """
        try:
            tok = self.__consume()
        except DXParserNoTokens:
            return

        if tok.equals('type'):
            tok = self.__consume()
            if not tok.iscode('STRING'):
                raise DXParseError('array: type was "%s", not a string.'%\
                                   tok.text)
            self.currentobject['type'] = tok.value()
        elif tok.equals('rank'):
            tok = self.__consume()
            try:
                self.currentobject['rank'] = tok.value('INTEGER')
            except ValueError:
                raise DXParseError('array: rank was "%s", not an integer.'%\
                                   tok.text)
        elif tok.equals('items'):
            tok = self.__consume()
            try:
                self.currentobject['size'] = tok.value('INTEGER')
            except ValueError:
                raise DXParseError('array: items was "%s", not an integer.'%\
                                   tok.text)
        elif tok.equals('data'):
            tok = self.__consume()
            if not tok.iscode('STRING'):
                raise DXParseError('array: data was "%s", not a string.'%\
                                   tok.text)
            if tok.text != 'follows':
                raise NotImplementedError(\
                            'array: Only the "data follows header" format is supported.')
            if not self.currentobject['size']:
                raise DXParseError("array: missing number of items")
            # This is the slow part.  Once we get here, we are just
            # reading in a long list of numbers.  Conversion to floats
            # will be done later when the numpy array is created.

            # Don't assume anything about whitespace or the number of elements per row
            self.currentobject['array'] = []
            while len(self.currentobject['array']) <self.currentobject['size']:
                 self.currentobject['array'].extend(self.dxfile.readline().strip().split())

            # If you assume that there are three elements per row
            # (except the last) the following version works and is a little faster.
            # for i in range(int(numpy.ceil(self.currentobject['size']/3))):
            #     self.currentobject['array'].append(self.dxfile.readline())
            # self.currentobject['array'] = ' '.join(self.currentobject['array']).split()
        elif tok.equals('attribute'):
            # not used at the moment
            attribute = self.__consume().value()
            if not self.__consume().equals('string'):
                raise DXParseError('array: "string" expected.')
            value = self.__consume().value()
        else:
            raise DXParseError('array: '+str(tok)+' not recognized.')

    def __field(self):
        """Level-2 parser for a DX field object.

        pattern:
        object "site map 1" class field
        component "positions" value 1
        component "connections" value 2
        component "data" value 3
        """
        try:
            tok = self.__consume()
        except DXParserNoTokens:
            return

        if tok.equals('component'):
            component = self.__consume().value()
            if not self.__consume().equals('value'):
                raise DXParseError('field: "value" expected')
            classid = self.__consume().value()
            try:
                self.currentobject['components'][component] = classid
            except KeyError:
                self.currentobject['components'] = {component:classid}
        else:
            raise DXParseError('field: '+str(tok)+' not recognized.')

    # parser routines independent of the dx classes
    # (with ideas from MDAnalysis.Selection and
    # http://effbot.org/zone/xml-scanner.htm)

    def use_parser(self,parsername):
        """Set parsername as the current parser and apply it."""
        self.__parser = self.parsers[parsername]
        self.__parser()
    def set_parser(self,parsername):
        """Set parsername as the current parser."""
        self.__parser = self.parsers[parsername]
    def apply_parser(self):
        """Apply the current parser to the token stream."""
        self.__parser()

    def __tokenize(self,string):
        """Split s into tokens and update the token buffer.

        __tokenize(string)

        New tokens are appended to the token buffer, discarding white
        space.  Based on http://effbot.org/zone/xml-scanner.htm
        """
        for m in self.dx_regex.finditer(string.strip()):
            code = m.lastgroup
            text = m.group(m.lastgroup)
            tok = Token(code,text)
            if not tok.iscode('WHITESPACE'):
                 self.tokens.append(tok)
                 # print "DEBUG tokenize: "+str(tok)

    def __refill_tokenbuffer(self):
        """Add a new tokenized line from the file to the token buffer.

        __refill_tokenbuffer()

        Only reads a new line if the buffer is empty. It is safe to
        call it repeatedly.

        At end of file, method returns empty strings and it is up to
        __peek and __consume to flag the end of the stream.
        """
        if len(self.tokens) == 0:
            self.__tokenize(self.dxfile.readline())

    def __peek(self):
        self.__refill_tokenbuffer()
        try:
            return self.tokens[0]
        except IndexError:
            raise DXParserNoTokens

    def __consume(self,):
        """Get the next token from the buffer and remove it/them.

        try:
          while 1:
             token = __consume()
        except DXParserNoTokens:
          pass
        """
        self.__refill_tokenbuffer()
        #print "DEBUG consume: "+str(self.__parser)+' '+str(self.__peek())
        try:
            return self.tokens.pop(0)  # singlet
        except IndexError:
            raise DXParserNoTokens
