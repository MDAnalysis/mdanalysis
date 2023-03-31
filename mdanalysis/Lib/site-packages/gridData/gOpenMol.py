# gridDataFormats --- python modules to read and write gridded data
# Copyright (c) 2009-2014 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Lesser General Public License, version 3 or later.
#
# Part of the documentation and format specification: Copyright CSC, 2005

"""
:mod:`gOpenMol` --- the gOpenMol plt format
===========================================

.. _gOpenMol: http://www.csc.fi/english/pages/g0penMol

The module provides a simple implementation of a reader for gOpenMol_
*plt* files. Plt files are binary files. The :class:`Plt` reader tries
to guess the endianess of the file, but this can fail (with a
:exc:`TypeError`); you are on your own in this case.

Only the reader is implemented. If you want to write gridded data use a format
that is more standard, such as OpenDX (see :mod:`OpenDX`).


Background
----------

gOpenMol http://www.csc.fi/english/pages/g0penMol plt format.

Used to be documented at http://www.csc.fi/gopenmol/developers/plt_format.phtml but currently this is only accessible through the internet archive at
http://web.archive.org/web/20061011125817/http://www.csc.fi/gopenmol/developers/plt_format.phtml



Grid data plt file format
-------------------------

Copyright CSC, 2005. Last modified: September 23, 2003 09:18:50

Plot file (plt) format The plot files are regular 3D grid files for plotting of
molecular orbitals, electron densities or other molecular properties. The plot
files are produced by several programs. It is also possible to format/unformat
plot files using the pltfile program in the utility directory. It is also
possible to produce plot files with external (own) programs. Produce first a
formatted text file and use then the pltfile program to unformat the file for
gOpenMol. The format for the plot files are very simple and a description of
the format can be found elsewhere in this manual. gOpenMol can read binary plot
files from different hardware platforms independent of the system type (little
or big endian machines).

Format of the binary ``*.plt`` file
...................................

The ``*.plt`` file binary and formatted file formats are very simple but please
observe that unformatted files written with a FORTRAN program are not pure
binary files because there are file records between the values while pure
binary files do not have any records between the values. gOpenMol should be
able to figure out if the file is pure binary or FORTRAN unformatted but it is
not very well tested.

Binary ``*.plt`` (grid) file format
...................................

Record number and meaning::

   #1: Integer, rank value must always be = 3
   #2: Integer, possible values are 1 ... 50. This value is not used but
   it can be used to define the type of surface!
   Values used (you can use your own value between 1... 50):

   1:   VSS surface
   2:   Orbital/density surface
   3:   Probe surface
   200: Gaussian 94/98
   201: Jaguar
   202: Gamess
   203: AutoDock
   204: Delphi/Insight
   205: Grid

   Value 100 is reserved for grid data coming from OpenMol!

   #3: Integer, number of points in z direction
   #4: Integer, number of points in y direction
   #5: Integer, number of points in x direction
   #6: Float, zmin value
   #7: Float, zmax value
   #8: Float, ymin value
   #9: Float, ymax value
   #10: Float, xmin value
   #11: Float, xmax value
   #12 ... Float, grid data values running (x is inner loop, then y and last z):

1.      Loop in the z direction
2.      Loop in the y direction
3.      Loop in the x direction

Example::

   nx=2  ny=1  nz=3

   0,0,0   1,0,0     y=0, z=0
   0,0,1   1,0,0     y=0, z=1
   0,0,2   1,0,2     y=0, z=2

The formatted (the first few lines) file can look like::

   3 2
   65 65 65
   -3.300000e+001 3.200000e+001 -3.300000e+001 3.200000e+001 -3.300000e+001 3.200000e+001
   -1.625609e+001 -1.644741e+001 -1.663923e+001 -1.683115e+001 -1.702274e+001 -1.721340e+001
   -1.740280e+001 -1.759018e+001 -1.777478e+001 -1.795639e+001 -1.813387e+001 -1.830635e+001
   ...

Formatted ``*.plt`` (grid) file format
......................................

Line numbers and variables on the line::

   line #1: Integer, Integer. Rank and type of surface (rank is always = 3)
   line #2: Integer, Integer, Integer. Zdim, Ydim, Xdim (number of points in the z,y,x directions)
   line #3: Float, Float, Float, Float, Float, Float. Zmin, Zmax, Ymin, Ymax, Xmin,Xmax (min and max values)
   line #4: ... Float. Grid data values running (x is inner loop, then y and last z) with one or several values per line:

    1. Loop in the z direction
    2. Loop in the y direction
    3. Loop in the x direction

Classes
-------

"""
import warnings
import struct
import numpy

class Record(object):
    def __init__(self, key, bintype, values=None):
        self.key = key
        self.bintype = bintype
        self.values = values  # dict(value='comment', ...)
    def is_legal(self, value):
        if self.values is None:
            return True
        return value in self.values
    def is_legal_dict(self, d):
        return self.is_legal(d[self.key])
    def __repr__(self):
        return "Record(%(key)r,%(bintype)r,...)" % vars(self)

class Plt(object):
    """A class to represent a gOpenMol_ plt file.

    Only reading is implemented; either supply a filename to the constructor
      >>> G = Plt(filename)
    or load the file with the read method
      >>> G = Plt()
      >>> G.read(filename)

    The data is held in :attr:`GOpenMol.array` and all header information is in
    the dict :attr:`GOpenMol.header`.

    :attr:`Plt.shape`
         D-tuplet describing size in each dimension
    :attr:`Plt.origin`
         coordinates of the centre of the grid cell with index 0,0,...,0
    :attr:`Plt.delta`
         DxD array describing the deltas

    """

    _header_struct =  (Record('rank',   'I', {3:'dimension'}),
                       Record('surface','I', {1: 'VSS surface',
                                              2: 'Orbital/density surface',
                                              3: 'Probe surface',
                                              42: 'gridcount',
                                              100: 'OpenMol',
                                              200: 'Gaussian 94/98',
                                              201: 'Jaguar',
                                              202: 'Gamess',
                                              203: 'AutoDock',
                                              204: 'Delphi/Insight',
                                              205: 'Grid',
                                              }),     # update in init with all user defined values
                       Record('nz',     'I'),
                       Record('ny',     'I'),
                       Record('nx',     'I'),
                       Record('zmin',   'f'),
                       Record('zmax',   'f'),
                       Record('ymin',   'f'),
                       Record('ymax',   'f'),
                       Record('xmin',   'f'),
                       Record('xmax',   'f'))
    _data_bintype = 'f'   #  write(&value,sizeof(float),1L,output);

    def __init__(self, filename=None):
        self.filename = str(filename)
        # fix header_struct because I cannot do {...}.update()
        rec_surf = [r for r in self._header_struct if r.key == 'surface'][0]
        rec_surf.values.update(dict((k,'user-defined') for k in range(4,51) if k != 42))
        # assemble format
        self._headerfmt = "".join([r.bintype for r in self._header_struct])

        if not filename is None:
            self.read(filename)

    def read(self, filename):
        """Populate the instance from the plt file *filename*."""
        from struct import calcsize, unpack
        if not filename is None:
            self.filename = str(filename)
        with open(self.filename, 'rb') as plt:
            h = self.header = self._read_header(plt)
            nentries = h['nx'] * h['ny'] * h['nz']
            # quick and dirty... slurp it all in one go
            datafmt = h['bsaflag']+str(nentries)+self._data_bintype
            a = numpy.array(unpack(datafmt, plt.read(calcsize(datafmt))))
        self.header['filename'] = self.filename
        self.array = a.reshape(h['nz'], h['ny'], h['nx']).transpose()  # unpack plt in reverse!!
        self.delta = self._delta()
        self.origin = numpy.array([h['xmin'], h['ymin'], h['zmin']]) + 0.5*numpy.diagonal(self.delta)
        self.rank = h['rank']

    @property
    def shape(self):
        return self.array.shape

    @property
    def edges(self):
        """Edges of the grid cells, origin at centre of 0,0,..,0 grid cell.

        Only works for regular, orthonormal grids.
        """
        return [self.delta[d,d] * numpy.arange(self.shape[d]+1) + self.origin[d]\
                - 0.5*self.delta[d,d]     for d in range(self.rank)]

    def _delta(self):
        h = self.header
        qmin = numpy.array([h['xmin'],h['ymin'],h['zmin']])
        qmax = numpy.array([h['xmax'],h['ymax'],h['zmax']])
        delta = numpy.abs(qmax - qmin) / self.shape
        return numpy.diag(delta)

    def _read_header(self, pltfile):
        """Read header bytes, try all possibilities for byte order/size/alignment."""
        nheader = struct.calcsize(self._headerfmt)
        names = [r.key for r in self._header_struct]
        binheader = pltfile.read(nheader)
        def decode_header(bsaflag='@'):
            h = dict(zip(names, struct.unpack(bsaflag+self._headerfmt, binheader)))
            h['bsaflag'] = bsaflag
            return h
        for flag in '@=<>':
            # try all endinaness and alignment options until we find something that looks sensible
            header = decode_header(flag)
            if header['rank'] == 3:
                break   # only legal value according to spec
            header = None
        if header is None:
            raise TypeError("Cannot decode header --- corrupted or wrong format?")
        for rec in self._header_struct:
            if not rec.is_legal_dict(header):
                warnings.warn("Key %s: Illegal value %r" % (rec.key, header[rec.key]))
        return header

    def histogramdd(self):
        """Return array data as (edges,grid), i.e. a numpy nD histogram."""
        return (self.array, self.edges)
