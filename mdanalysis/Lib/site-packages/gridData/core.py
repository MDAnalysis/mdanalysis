# gridDataFormats --- python modules to read and write gridded data
# Copyright (c) 2009-2014 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Lesser General Public License, version 3 or later.
r"""
Core functionality for storing n-D grids --- :mod:`gridData.core`
=================================================================

The :mod:`core` module contains classes and functions that are
independent of the grid data format. In particular this module
contains the :class:`Grid` class that acts as a universal constructor
for specific formats::

 g = Grid(**kwargs)           # construct
 g.export(filename, format)   # export to the desired format

Some formats can also be read::

 g = Grid()                   # make an empty Grid
 g.load(filename)             # populate with data from filename


Classes and functions
---------------------

"""
import os
import errno
import pickle

import numpy

# For interpolated grids: need scipy.ndimage but we import it only when needed:
# import scipy

from . import OpenDX
from . import gOpenMol
from . import mrc


def _grid(x):
    """Access the underlying ndarray of a Grid object or return the object itself"""
    try:
        return x.grid
    except AttributeError:
        return x


class Grid(object):
    """A multidimensional grid object with origin and grid spacings.

    :class:`Grid` objects can be used in arithmetical calculations
    just like numpy arrays *if* they are compatible, i.e., they have
    the same shapes and lengths. In order to make arrays compatible,
    they an be resampled (:meth:`resample`) on a common grid.

    The attribute :attr:`grid` that holds the data is a standard numpy
    array and so the data can be directly manipulated.

    Data can be read from a number of molecular volume/density formats
    and written out in different formats with :meth:`export`.


    Parameters
    ----------
    grid : numpy.ndarray or str (optional)
      Build the grid either from a histogram or density (a numpy nD
      array) or read data from a filename.

    edges : list (optional)
      List of arrays, the lower and upper bin edges along the axes
      (same as the output by :func:`numpy.histogramdd`)

    origin : :class:`numpy.ndarray` (optional)
      Cartesian coordinates of the center of grid position at index
      ``[0, 0, ..., 0]``.

    delta : :class:`numpy.ndarray` (optional)
      Either ``n x n`` array containing the cell lengths in each dimension,
      or ``n x 1`` array for rectangular arrays.

    metadata : dict (optional)
      A user defined dictionary of arbitrary key/value pairs
      associated with the density; the class does not touch
      :attr:`metadata` but stores it with :meth:`save`

    interpolation_spline_order : int (optional)
      Order of interpolation function for resampling with
      :func:`resample`; cubic splines = 3 and the default is 3

    file_format : str (optional)
      Name of the file format; only necessary when `grid` is a
      filename (see :meth:`load`) and autodetection of the file
      format fails. The default is ``None`` and normally the file
      format is guessed from the file extension.

    Raises
    ------
    TypeError
      If the dimensions of the various input data do not agree with
      each other.
    ValueError
      If some of the required data are not provided in the keyword
      arguments, e.g., if only the `grid` is supplied as an array but
      not the `edges` or only `grid` and one of `origin` and `delta`.
    NotImplementedError
      If triclinic (non-orthorhombic) boxes are supplied in `delta`

      .. Note:: `delta` can only be a 1D array of length :attr:`grid.ndim`


    Attributes
    ----------
    grid : :class:`numpy.ndarray`
      This array can be any number of dimensions supported by NumPy
      in order to represent high-dimensional data. When used with
      data that represents real space densities then the **axis
      convention in GridDataFormats** is that axis 0 corresponds to
      the Cartesian :math:`x` component, axis 1 corresponds to the
      :math:`y` component, and axis 2 to the :math:`z` component.

    delta : :class:`numpy.ndarray`
      Length of a grid cell (spacing or voxelsize) in :math:`x`,
      :math:`y`, :math:`z` dimensions. This is a *1D array* with
      length :attr:`Grid.grid.ndim`.

    origin : :class:`numpy.ndarray`
      Array with the Cartesian coordinates of the coordinate system
      origin, the *center* of cell ``Grid.grid[0, 0, .., 0]``.

    edges : list
      List of arrays, one for each axis in :attr:`grid`.  Each 1D edge
      array describes the *edges* of the grid cells along the
      corresponding axis. The length of an edge array for axis ``i``
      is ``grid.shape[i] + 1`` because it contains the lower boundary
      for the first cell, the boundaries between all grid cells, and
      the upper boundary for the last cell. The edges are assumed to
      be regular with spacing indicated in :attr:`delta`, namely
      ``Grid.delta[i]`` for axis ``i``.

    midpoints : list
      List of arrays, one for each axis in :attr:`grid`.  Each 1D
      midpoints array contains the *midpoints* of the grid cells along
      the corresponding axis.

    metadata : dict
      A user-defined dictionary that can be used to annotate the
      data. The content is not touched by :class:`Grid`. It is saved
      together with the other data with :meth:`save`.


    Example
    -------
    Create a Grid object from data.

    From :func:`numpy.histogramdd`::

      grid, edges = numpy.histogramdd(...)
      g = Grid(grid, edges=edges)

    From an arbitrary grid::

      g = Grid(grid, origin=origin, delta=delta)

    From a saved file::

      g = Grid(filename)

    or ::

      g = Grid()
      g.load(filename)


    Notes
    -----
    In principle, the dimension (number of axes) is arbitrary but in
    practice many formats only support three and almost all
    functionality is only tested for this special case.

    The :meth:`export` method with ``format='dx'`` always exports a 3D
    object. Other methods might work for an array of any dimension (in
    particular the Python pickle output).


    .. versionchanged:: 0.5.0
       New *file_format* keyword argument.

    .. versionchanged:: 0.7.0
       CCP4 files are now read with :class:`gridData.mrc.MRC` and not anymore
       with the deprecated/buggy `ccp4.CCP4`

    """

    #: Default format for exporting with :meth:`export`.
    default_format = 'DX'

    def __init__(self, grid=None, edges=None, origin=None, delta=None,
                 metadata=None, interpolation_spline_order=3,
                 file_format=None):
        # file formats are guessed from extension == lower case key
        self._exporters = {
            'DX': self._export_dx,
            'PKL': self._export_python,
            'PICKLE': self._export_python,  # compatibility
            'PYTHON': self._export_python,  # compatibility
        }
        self._loaders = {
            'CCP4': self._load_mrc,
            'MRC':  self._load_mrc,
            'DX': self._load_dx,
            'PLT': self._load_plt,
            'PKL': self._load_python,
            'PICKLE': self._load_python,  # compatibility
            'PYTHON': self._load_python,  # compatibility
        }

        self.metadata = metadata if metadata is not None else {}
        self.__interpolated = None  # cache for interpolated grid
        self.__interpolation_spline_order = interpolation_spline_order
        self.interpolation_cval = None  # default to using min(grid)

        if grid is not None:
            if isinstance(grid, str):
                # can probably safely try to load() it...
                filename = grid
            else:
                try:
                    # Can we read this as a file?
                    # Use str(x) to work with py.path.LocalPath and pathlib.Path instances
                    # even for Python < 3.6
                    with open(str(grid), 'rb'):
                        pass
                except (OSError, IOError):
                    # no, this is probably an array-like thingy
                    filename = None
                else:
                    # yes, let's use it as a file
                    filename = str(grid)

            if filename is not None:
                self.load(filename, file_format=file_format)
            else:
                self._load(grid, edges, metadata, origin, delta)

    @property
    def interpolation_spline_order(self):
        """Order of the B-spline interpolation of the data.

        3 = cubic; 4 & 5 are also supported

        Only choose values that are acceptable to
        :func:`scipy.ndimage.spline_filter`!

        See Also
        --------
        interpolated
        """

        return self.__interpolation_spline_order

    @interpolation_spline_order.setter
    def interpolation_spline_order(self, x):
        """Setting the ``interpolation_spline_order`` updates :func:`interpolated`

        Because we cache the interpolation function, we need to rebuild the
        cache whenever the interpolation order changes: this is
        handled by :meth:`_update`

        """
        self.__interpolation_spline_order = x
        self._update()

    def resample(self, edges):
        """Resample data to a new grid with edges *edges*.

        This method creates a new grid with the data from the current
        grid resampled to a regular grid specified by `edges`.  The
        order of the interpolation is set by
        :attr:`Grid.interpolation_spline_order`: change the value
        *before* calling :meth:`resample`.

        Parameters
        ----------
        edges : tuple of arrays or Grid
             edges of the new grid or a :class:`Grid` instance that
             provides :attr:`Grid.edges`

        Returns
        -------
        Grid
             a new :class:`Grid` with the data interpolated over the
             new grid cells


        Examples
        --------

        Providing `edges` (a tuple of three arrays, indicating the
        boundaries of each grid cell)::

          g = grid.resample(edges)

        As a convenience, one can also supply another :class:`Grid` as
        the argument for this method ::

          g = grid.resample(othergrid)

        and the edges are taken from :attr:`Grid.edges`.

        """
        try:
            edges = edges.edges  # can also supply another Grid
        except AttributeError:
            pass
        midpoints = self._midpoints(edges)
        coordinates = ndmeshgrid(*midpoints)
        # feed a meshgrid to generate all points
        newgrid = self.interpolated(*coordinates)
        return self.__class__(newgrid, edges)

    def resample_factor(self, factor):
        """Resample to a new regular grid.


        Parameters
        ----------
        factor : float
            The number of grid cells are scaled with `factor` in each
            dimension, i.e., ``factor * N_i`` cells along each
            dimension i. Must be positive, and cannot result in fewer
            than 2 cells along a dimension.


        Returns
        -------
        interpolated grid : Grid
            The resampled data are represented on a :class:`Grid` with the new
            grid cell sizes.

        See Also
        --------
        resample


        .. versionchanged:: 0.6.0
           Previous implementations would not alter the range of the grid edges
           being resampled on. As a result, values at the grid edges would creep
           steadily inward. The new implementation recalculates the extent of
           grid edges for every resampling.

        """
        if float(factor) <= 0:
            raise ValueError("Factor must be positive")
        # Determine current spacing
        spacing = (numpy.array(self._max_edges()) - numpy.array(self._min_edges())) / (
                  -1 + numpy.array(self._len_edges()))
        # First guess at the new spacing is inversely related to the
        # magnification factor.
        newspacing = spacing / float(factor)
        smidpoints = numpy.array(self._midpoints())
        # We require that the new spacing result in an even subdivision of the
        # existing midpoints
        newspacing = (smidpoints[:, -1] - smidpoints[:, 0]) / (numpy.maximum(
            1, numpy.floor((smidpoints[:, -1] - smidpoints[:, 0]) / newspacing)))
        # How many edge points should there be? It is the number of intervals
        # between midpoints + 2
        edgelength = 2 + \
            numpy.round((smidpoints[:, -1] - smidpoints[:, 0]) / newspacing)
        edges = [numpy.linspace(start, stop, num=int(N), endpoint=True) for (start, stop, N) in zip(
            smidpoints[:, 0] - 0.5 * newspacing, smidpoints[:, -1] + 0.5 * newspacing, edgelength)]
        return self.resample(edges)

    def _update(self):
        """compute/update all derived data

        Can be called without harm and is idem-potent.

        Updates these attributes and methods:
           :attr:`origin`
              the center of the cell with index 0,0,0
           :attr:`midpoints`
              centre coordinate of each grid cell
           :meth:`interpolated`
              spline interpolation function that can generated a value for
              coordinate
        """
        self.delta = numpy.array(list(
            map(lambda e: (e[-1] - e[0]) / (len(e) - 1), self.edges)))
        self.midpoints = self._midpoints(self.edges)
        self.origin = numpy.array(list(map(lambda m: m[0], self.midpoints)))
        if self.__interpolated is not None:
            # only update if we are using it
            self.__interpolated = self._interpolationFunctionFactory()

    @property
    def interpolated(self):
        """B-spline function over the data grid(x,y,z).

        The :func:`interpolated` function allows one to obtain data
        values for any values of the coordinates::

           interpolated([x1,x2,...],[y1,y2,...],[z1,z2,...]) -> F[x1,y1,z1],F[x2,y2,z2],...

        The interpolation order is set in
        :attr:`Grid.interpolation_spline_order`.

        The interpolated function is computed once and is cached for better
        performance. Whenever :attr:`~Grid.interpolation_spline_order` is
        modified, :meth:`Grid.interpolated` is recomputed.

        The value for unknown data is set in :attr:`Grid.interpolation_cval`
        (TODO: also recompute when ``interpolation_cval`` value is changed.)

        Example
        -------
        Example usage for resampling::

           XX, YY, ZZ = numpy.mgrid[40:75:0.5, 96:150:0.5, 20:50:0.5]
           FF = interpolated(XX, YY, ZZ)

        Note
        ----
        Values are interpolated with a spline function. It is possible
        that the spline will generate values that would not normally
        appear in the data. For example, a density is non-negative but
        a cubic spline interpolation can generate negative values,
        especially at the boundary between 0 and high values.

        Internally, the function uses :func:`scipy.ndimage.map_coordinates`
        with ``mode="constant"`` whereby interpolated values outside
        the interpolated grid are determined by filling all values beyond
        the edge with the same constant value, defined by the
        :attr:`interpolation_cval` parameter, which when not set defaults
        to the minimum value in the interpolated grid.


        .. versionchanged:: 0.6.0
           Interpolation outside the grid is now performed with
           ``mode="constant"`` rather than ``mode="nearest"``, eliminating
           extruded volumes when interpolating beyond the grid.

        """
        if self.__interpolated is None:
            self.__interpolated = self._interpolationFunctionFactory()
        return self.__interpolated

    def _map_edges(self, func, edges=None):
        if edges is None:
            edges = self.edges
        return [func(e) for e in edges]

    def _midpoints(self, edges=None):
        return self._map_edges(lambda e: 0.5 * (e[:-1] + e[1:]), edges=edges)

    def _len_edges(self, edges=None):
        return self._map_edges(len, edges=edges)

    def _min_edges(self, edges=None):
        return self._map_edges(numpy.min, edges=edges)

    def _max_edges(self, edges=None):
        return self._map_edges(numpy.max, edges=edges)

    def _guess_format(self, filename, file_format=None, export=True):
        if export:
            available = self._exporters
        else:
            available = self._loaders
        if file_format is None:
            splitted = os.path.splitext(filename)
            if splitted[1][1:] in ('gz', ):
                file_format = os.path.splitext(splitted[0])[1][1:]
            else:
                file_format = splitted[1][1:]
        file_format = file_format.upper()
        if not file_format:
            file_format = self.default_format
        if file_format not in available:
            raise ValueError(
                "File format {} not available, choose one of {}".format(
                    file_format, available.keys()))
        return file_format

    def _get_exporter(self, filename, file_format=None):
        return self._exporters[self._guess_format(filename,
                                                  file_format=file_format,
                                                  export=True)]

    def _get_loader(self, filename, file_format=None):
        return self._loaders[self._guess_format(filename,
                                                file_format=file_format,
                                                export=False)]

    def _load(
            self,
            grid=None,
            edges=None,
            metadata=None,
            origin=None,
            delta=None):
        if edges is not None:
            # set up from histogramdd-type data
            self.grid = numpy.asanyarray(grid)
            self.edges = edges
            self._update()
        elif origin is not None and delta is not None:
            # setup from generic data
            origin = numpy.asanyarray(origin)
            delta = numpy.asanyarray(delta)
            if len(origin) != grid.ndim:
                raise TypeError(
                    "Dimension of origin is not the same as grid dimension.")
            if delta.shape == () and numpy.isreal(delta):
                delta = numpy.ones(grid.ndim) * delta
            elif delta.ndim > 1:
                raise NotImplementedError(
                    "Non-rectangular grids are not supported.")
            elif len(delta) != grid.ndim:
                raise TypeError("delta should be scalar or array-like of"
                                "len(grid.ndim)")
            # note that origin is CENTER so edges must be shifted by -0.5*delta
            self.edges = [origin[dim] +
                          (numpy.arange(m + 1) - 0.5) * delta[dim]
                          for dim, m in enumerate(grid.shape)]
            self.grid = numpy.asanyarray(grid)
            self._update()
        else:
            raise ValueError(
                "Wrong/missing data to set up Grid. Use Grid() or "
                "Grid(grid=<array>, edges=<list>) or "
                "Grid(grid=<array>, origin=(x0, y0, z0), delta=(dx, dy, dz)):\n"
                "grid={0} edges={1} origin={2} delta={3}".format(
                    grid, edges, origin, delta))

    def load(self, filename, file_format=None):
        """Load saved grid and edges from `filename`

        The :meth:`load` method calls the class's constructor method and
        completely resets all values, based on the loaded data.
        """
        filename = str(filename)
        if not os.path.exists(filename):
            # check before we try to detect the file type because
            # _guess_fileformat() does not work well with things that
            # are not really a file
            raise IOError(errno.ENOENT, "file not found", filename)
        loader = self._get_loader(filename, file_format=file_format)
        loader(filename)

    def _load_python(self, filename):
        with open(filename, 'rb') as f:
            saved = pickle.load(f)
        self._load(grid=saved['grid'],
                   edges=saved['edges'],
                   metadata=saved['metadata'])

    def _load_mrc(self, filename):
        """Initializes Grid from a MRC/CCP4 file."""
        mrcfile = mrc.MRC(filename)
        grid, edges = mrcfile.histogramdd()
        self._load(grid=grid, edges=edges, metadata=self.metadata)
        # Store header for access from Grid object (undocumented)
        # https://github.com/MDAnalysis/GridDataFormats/pull/100#discussion_r782604833
        self._mrc_header = mrcfile.header.copy()

    def _load_dx(self, filename):
        """Initializes Grid from a OpenDX file."""
        dx = OpenDX.field(0)
        dx.read(filename)
        grid, edges = dx.histogramdd()
        self._load(grid=grid, edges=edges, metadata=self.metadata)

    def _load_plt(self, filename):
        """Initialize Grid from gOpenMol plt file."""
        g = gOpenMol.Plt()
        g.read(filename)
        grid, edges = g.histogramdd()
        self._load(grid=grid, edges=edges, metadata=self.metadata)

    def export(self, filename, file_format=None, type=None, typequote='"'):
        """export density to file using the given format.

        The format can also be deduced from the suffix of the filename
        although the `file_format` keyword takes precedence.

        The default format for :meth:`export` is 'dx'.  Use 'dx' for
        visualization.

        Implemented formats:

        dx
            :mod:`OpenDX`
        pickle
            pickle (use :meth:`Grid.load` to restore); :meth:`Grid.save`
            is simpler than ``export(format='python')``.

        Parameters
        ----------
        filename : str
            name of the output file

        file_format : {'dx', 'pickle', None} (optional)
            output file format, the default is "dx"

        type : str (optional)
            for DX, set the output DX array type, e.g., "double" or "float".
            By default (``None``), the DX type is determined from the numpy
            dtype of the array of the grid (and this will typically result in
            "double").

            .. versionadded:: 0.4.0

        typequote : str (optional)
            For DX, set the character used to quote the type string;
            by default this is a double-quote character, '"'.
            Custom parsers like the one from NAMD-GridForces (backend for MDFF)
            expect no quotes, and typequote='' may be used to appease them.

            .. versionadded:: 0.5.0

        """
        filename = str(filename)
        exporter = self._get_exporter(filename, file_format=file_format)
        exporter(filename, type=type, typequote=typequote)

    # note: the _export_FORMAT() methods all take the filename as a mandatory
    # argument. They can process kwargs but they are not required to do
    # so. However, they must ignore any kwargs that they are not processing.

    def _export_python(self, filename, **kwargs):
        """Pickle the Grid object

        The object is dumped as a dictionary with grid and edges: This
        is sufficient to recreate the grid object with ``__init__()``.
        """
        data = dict(grid=self.grid, edges=self.edges, metadata=self.metadata)
        with open(filename, 'wb') as f:
            pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)

    def _export_dx(self, filename, type=None, typequote='"', **kwargs):
        """Export the density grid to an OpenDX file.

        The file format is the simplest regular grid array and it is
        also understood by VMD's and Chimera's DX reader; PyMOL
        requires the dx `type` to be set to "double".

        For the file format see
        http://opendx.sdsc.edu/docs/html/pages/usrgu068.htm#HDREDF

        """
        root, ext = os.path.splitext(filename)
        filename = root + '.dx'

        comments = [
            'OpenDX density file written by gridDataFormats.Grid.export()',
            'File format: http://opendx.sdsc.edu/docs/html/pages/usrgu068.htm#HDREDF',
            'Data are embedded in the header and tied to the grid positions.',
            'Data is written in C array order: In grid[x,y,z] the axis z is fastest',
            'varying, then y, then finally x, i.e. z is the innermost loop.']

        # write metadata in comments section
        if self.metadata:
            comments.append('Meta data stored with the python Grid object:')
        for k in self.metadata:
            comments.append('   ' + str(k) + ' = ' + str(self.metadata[k]))
        comments.append(
            '(Note: the VMD dx-reader chokes on comments below this line)')

        components = dict(
            positions=OpenDX.gridpositions(1, self.grid.shape, self.origin,
                                           self.delta),
            connections=OpenDX.gridconnections(2, self.grid.shape),
            data=OpenDX.array(3, self.grid, type=type, typequote=typequote),
        )
        dx = OpenDX.field('density', components=components, comments=comments)
        if ext == '.gz':
            filename = root + ext
        dx.write(filename)

    def save(self, filename):
        """Save a grid object to `filename` and add ".pickle" extension.

        Internally, this calls
        ``Grid.export(filename, format="python")``. A grid can be
        regenerated from the saved data with ::

           g = Grid(filename="grid.pickle")

        .. note::
           The pickle format depends on the Python version and
           therefore it is not guaranteed that a grid saved with, say,
           Python 2.7 can also be read with Python 3.5. The OpenDX format
           is a better alternative for portability.

        """
        self.export(filename, file_format="pickle")

    def centers(self):
        """Returns the coordinates of the centers of all grid cells as an iterator.


        See Also
        --------
        :func:`numpy.ndindex`
        """
        for idx in numpy.ndindex(self.grid.shape):
            yield self.delta * numpy.array(idx) + self.origin

    def check_compatible(self, other):
        """Check if `other` can be used in an arithmetic operation.

        `other` is compatible if

        1) `other` is a scalar
        2) `other` is a grid defined on the same edges

        In order to make `other` compatible, resample it on the same
        grid as this one using :meth:`resample`.

        Parameters
        ----------
        other : Grid or float or int
           Another object to be used for standard arithmetic
           operations with this :class:`Grid`

        Raises
        ------
        TypeError
              if not compatible

        See Also
        --------
        :meth:`resample`
        """

        if not (numpy.isreal(other) or self == other):
            raise TypeError(
                "The argument cannot be arithmetically combined with the grid. "
                "It must be a scalar or a grid with identical edges. "
                "Use Grid.resample(other.edges) to make a new grid that is "
                "compatible with other.")
        return True

    def _interpolationFunctionFactory(self, spline_order=None, cval=None):
        """Returns a function F(x,y,z) that interpolates any values on the grid.

        _interpolationFunctionFactory(self,spline_order=3,cval=None) --> F

        *cval* is set to :meth:`Grid.grid.min`. *cval* cannot be chosen too
        large or too small or NaN because otherwise the spline interpolation
        breaks down near that region and produces wild oscillations.

        .. Note:: Only correct for equally spaced values (i.e. regular edges with
                  constant delta).
        .. SeeAlso:: http://www.scipy.org/Cookbook/Interpolation
        """
        # for scipy >=0.9: should use scipy.interpolate.griddata
        # http://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html#scipy.interpolate.griddata
        # (does it work for nD?)
        import scipy.ndimage

        if spline_order is None:
            # must be compatible with whatever
            # :func:`scipy.ndimage.spline_filter` takes.
            spline_order = self.interpolation_spline_order
        if cval is None:
            cval = self.interpolation_cval

        data = self.grid
        if cval is None:
            cval = data.min()
        try:
            # masked arrays, fill with min: should keep spline happy
            _data = data.filled(cval)
        except AttributeError:
            _data = data

        coeffs = scipy.ndimage.spline_filter(_data, order=spline_order)
        x0 = self.origin
        dx = self.delta

        def _transform(cnew, c0, dc):
            return (numpy.atleast_1d(cnew) - c0) / dc

        def interpolatedF(*coordinates):
            """B-spline function over the data grid(x,y,z).

            interpolatedF([x1,x2,...],[y1,y2,...],[z1,z2,...]) -> F[x1,y1,z1],F[x2,y2,z2],...

            Example usage for resampling::
              >>> XX,YY,ZZ = numpy.mgrid[40:75:0.5, 96:150:0.5, 20:50:0.5]
              >>> FF = _interpolationFunction(XX,YY,ZZ)
            """
            _coordinates = numpy.array(
                [_transform(coordinates[i], x0[i], dx[i]) for i in range(len(
                    coordinates))])
            return scipy.ndimage.map_coordinates(coeffs,
                                                 _coordinates,
                                                 prefilter=False,
                                                 mode='constant',
                                                 cval=cval)
        return interpolatedF

    def __eq__(self, other):
        if not isinstance(other, Grid):
            return False
        return numpy.all(
            other.grid == self.grid) and numpy.all(
            other.origin == self.origin) and numpy.all(
            numpy.all(
                other_edge == self_edge) for other_edge,
            self_edge in zip(
                other.edges,
                self.edges))

    def __ne__(self, other):
        return not self.__eq__(other)

    def __add__(self, other):
        self.check_compatible(other)
        return self.__class__(self.grid + _grid(other), edges=self.edges)

    def __sub__(self, other):
        self.check_compatible(other)
        return self.__class__(self.grid - _grid(other), edges=self.edges)

    def __mul__(self, other):
        self.check_compatible(other)
        return self.__class__(self.grid * _grid(other), edges=self.edges)

    def __truediv__(self, other):
        self.check_compatible(other)
        return self.__class__(self.grid / _grid(other), edges=self.edges)

    def __floordiv__(self, other):
        self.check_compatible(other)
        return self.__class__(self.grid // _grid(other), edges=self.edges)

    def __pow__(self, other):
        self.check_compatible(other)
        return self.__class__(
            numpy.power(
                self.grid,
                _grid(other)),
            edges=self.edges)

    def __radd__(self, other):
        self.check_compatible(other)
        return self.__class__(_grid(other) + self.grid, edges=self.edges)

    def __rsub__(self, other):
        self.check_compatible(other)
        return self.__class__(_grid(other) - self.grid, edges=self.edges)

    def __rmul__(self, other):
        self.check_compatible(other)
        return self.__class__(_grid(other) * self.grid, edges=self.edges)

    def __rtruediv__(self, other):
        self.check_compatible(other)
        return self.__class__(_grid(other) / self.grid, edges=self.edges)

    def __rfloordiv__(self, other):
        self.check_compatible(other)
        return self.__class__(_grid(other) // self.grid, edges=self.edges)

    def __rpow__(self, other):
        self.check_compatible(other)
        return self.__class__(
            numpy.power(
                _grid(other),
                self.grid),
            edges=self.edges)

    def __repr__(self):
        try:
            bins = self.grid.shape
        except AttributeError:
            bins = "no"
        return '<{0} with {1!r} bins>'.format(self.__class__, bins)


def ndmeshgrid(*arrs):
    """Return a mesh grid for N dimensions.

    The input are N arrays, each of which contains the values along one axis of
    the coordinate system. The arrays do not have to have the same number of
    entries. The function returns arrays that can be fed into numpy functions
    so that they produce values for *all* points spanned by the axes *arrs*.

    Original from
    http://stackoverflow.com/questions/1827489/numpy-meshgrid-in-3d and fixed.

    .. SeeAlso: :func:`numpy.meshgrid` for the 2D case.
    """
    # arrs = tuple(reversed(arrs)) <-- wrong on stackoverflow.com
    arrs = tuple(arrs)
    lens = list(map(len, arrs))
    dim = len(arrs)

    sz = 1
    for s in lens:
        sz *= s

    ans = []
    for i, arr in enumerate(arrs):
        slc = [1] * dim
        slc[i] = lens[i]
        arr2 = numpy.asanyarray(arr).reshape(slc)
        for j, sz in enumerate(lens):
            if j != i:
                arr2 = arr2.repeat(sz, axis=j)
        ans.append(arr2)

    return tuple(ans)
