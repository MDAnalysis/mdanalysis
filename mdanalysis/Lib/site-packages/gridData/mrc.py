# gridDataFormats --- python modules to read and write gridded data
# Copyright (c) 2009-2021 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Lesser General Public License, version 3 or later.
#

""":mod:`mrc` --- the MRC/CCP4 volumetric data format
===================================================

.. versionadded:: 0.7.0

Reading of MRC/CCP4 volumetric files (`MRC2014 file format`_) using
the mrcfile_ library [Burnley2017]_.

.. _mrcfile: https://mrcfile.readthedocs.io/
.. _`MRC2014 file format`: http://www.ccpem.ac.uk/mrc_format/mrc2014.php


References
----------

.. [Burnley2017] Burnley T, Palmer C and Winn M (2017) Recent
                 developments in the CCP-EM software suite. *Acta
                 Cryst.* D73:469-477. doi: `10.1107/S2059798317007859`_


.. _`10.1107/S2059798317007859`: https://doi.org/10.1107/S2059798317007859

Classes
-------

"""
import numpy as np
import mrcfile


class MRC(object):
    """Represent a MRC/CCP4 file.

    Load `MRC/CCP4 2014 <MRC2014 file format>`_ 3D volumetric data with
    the mrcfile_ library.

    Parameters
    ----------
    filename : str (optional)
       input file (or stream), can be compressed

    Raises
    ------
    ValueError
       If the unit cell is not orthorhombic or if the data
       are not volumetric.


    Attributes
    ----------
    header : numpy.recarray
       Header data from the MRC file as a numpy record array.

    array : numpy.ndarray
       Data as a 3-dimensional array where axis 0 corresponds to X,
       axis 1 to Y, and axis 2 to Z. This order is always enforced,
       regardless of the order in the mrc file.

    delta : numpy.ndarray
       Diagonal matrix with the voxel size in X, Y, and Z direction
       (taken from the :attr:`mrcfile.mrcfile.voxel_size` attribute)

    origin : numpy.ndarray
       numpy array with coordinates of the coordinate system origin
       (computed from :attr:`header.origin`, the offsets
       :attr:`header.origin.nxstart`, :attr:`header.origin.nystart`,
       :attr:`header.origin.nzstart` and the spacing :attr:`delta`)

    rank : int
       The integer 3, denoting that only 3D maps are read.


    Notes
    -----
    * Only volumetric (3D) densities are read.
    * Only orthorhombic unitcells supported (other raise :exc:`ValueError`)
    * Only reading is currently supported.


    .. versionadded:: 0.7.0

    """

    def __init__(self, filename=None):
        self.filename = filename
        if filename is not None:
            self.read(filename)

    def read(self, filename):
        """Populate the instance from the MRC/CCP4 file *filename*."""
        if filename is not None:
            self.filename = filename
        with mrcfile.open(filename) as mrc:
            if not mrc.is_volume():                           #pragma: no cover
                raise ValueError(
                    "MRC file {} is not a volumetric density.".format(filename))
            self.header = h = mrc.header.copy()
            # check for being orthorhombic
            if not np.allclose([h.cellb.alpha, h.cellb.beta, h.cellb.gamma],
                               [90, 90, 90]):
                raise ValueError("Only orthorhombic unitcells are currently "
                                 "supported, not "
                                 "alpha={0}, beta={1}, gamma={2}".format(
                                     h.cellb.alpha, h.cellb.beta, h.cellb.gamma))
            # mrc.data[z, y, x] indexed: convert to x,y,z as used in GridDataFormats
            # together with the axes orientation information in mapc/mapr/maps.
            # mapc, mapr, maps = 1, 2, 3 for Fortran-ordering and 3, 2, 1 for C-ordering.
            # Other combinations are possible. We reorder the data for the general case
            # by sorting mapc, mapr, maps in ascending order, i.e., to obtain x,y,z.
            # mrcfile provides the data in zyx shape (without regard to map*) so we first
            # transpose it to xyz and then reorient with axes_c_order.
            #
            # All other "xyz" quantitities are also reordered.
            axes_order = np.hstack([h.mapc, h.mapr, h.maps])
            axes_c_order = np.argsort(axes_order)
            transpose_order = np.argsort(axes_order[::-1])
            self.array = np.transpose(mrc.data, axes=transpose_order)
            self.delta = np.diag(np.array([mrc.voxel_size.x, mrc.voxel_size.y, mrc.voxel_size.z]))
            # the grid is shifted to the MRC origin by offset
            # (assume orthorhombic)
            offsets = np.hstack([h.nxstart, h.nystart, h.nzstart])[axes_c_order] * np.diag(self.delta)
            # GridData origin is centre of cell at x=col=0, y=row=0 z=seg=0
            self.origin = np.hstack([h.origin.x, h.origin.y, h.origin.z]) + offsets
            self.rank = 3

    @property
    def shape(self):
        """Shape of the :attr:`array`"""
        return self.array.shape

    @property
    def edges(self):
        """Edges of the grid cells, origin at centre of 0,0,0 grid cell.

        Only works for regular, orthonormal grids.
        """
        # TODO: Add triclinic cell support.
        return [self.delta[d, d] * np.arange(self.shape[d] + 1) +
                self.origin[d] - 0.5 * self.delta[d, d]
                for d in range(self.rank)]

    def histogramdd(self):
        """Return array data as (edges,grid), i.e. a numpy nD histogram."""
        return (self.array, self.edges)



