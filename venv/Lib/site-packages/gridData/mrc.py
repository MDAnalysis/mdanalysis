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

.. autoclass:: MRC
   :members:


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
    assume_volumetric : bool (optional)
      If ``False`` (default), check the file header to determine whether
      the data in `grid` is a 3D volume. If ``True``, assume `grid` is volumetric.

      .. versionadded:: 1.1.0
        

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
    * Reading and writing are supported.


    .. versionadded:: 0.7.0

    """

    def __init__(self, filename=None, assume_volumetric=False):
        self.filename = filename
        if filename is not None:
            self.read(filename, assume_volumetric=assume_volumetric)

    def read(self, filename, assume_volumetric=False):
        """Populate the instance from the MRC/CCP4 file *filename*."""
        if filename is not None:
            self.filename = filename
        with mrcfile.open(filename) as mrc:
            if assume_volumetric:
                # non 3D volumes should always fail, regardless of assume_volumetric value
                is_volume = mrc.data is not None and len(mrc.data.shape) == 3
            else:
                is_volume = mrc.is_volume()

            if not is_volume:
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

    def write(self, filename):
        """Write grid data to MRC/CCP4 file format.
        
        Parameters
        ----------
        filename : str
            Output filename for the MRC file
        
        Notes
        -----
        The data array should be in xyz order (axis 0=X, axis 1=Y, axis 2=Z).
        
        If the MRC object was created by reading an existing file, the original
        header information (including mapc, mapr, maps ordering) is preserved.
        Otherwise, standard ordering (mapc=1, mapr=2, maps=3) is used.
        
        
        .. versionadded:: 1.1.0
        """
        if filename is not None:
            self.filename = filename
        
        # Preserve header if it exists, otherwise use defaults
        if hasattr(self, 'header'):
            # File was read - preserve original ordering
            h = self.header
            axes_order = np.hstack([h.mapc, h.mapr, h.maps])
            mapc, mapr, maps = int(h.mapc), int(h.mapr), int(h.maps)
        else:
            # New file - use standard ordering
            axes_order = np.array([1, 2, 3])
            mapc, mapr, maps = 1, 2, 3
            h = None
        
        # Reverse the transformation done in read()
        transpose_order = np.argsort(axes_order[::-1])
        inverse_transpose_order = np.argsort(transpose_order)
        
        # Transform our xyz array back to the file's native ordering
        data_for_file = np.transpose(self.array, axes=inverse_transpose_order)
        
        # Ensure proper data type (float32 is standard for mode 2)
        data_for_file = data_for_file.astype(np.float32)
        
        # Create new MRC file
        with mrcfile.new(filename, overwrite=True) as mrc:
            mrc.set_data(data_for_file)
            
            # Set voxel size from delta (diagonal elements)
            voxel_size = np.diag(self.delta).astype(np.float32)
            mrc.voxel_size = tuple(voxel_size)
            
            # Set map ordering
            mrc.header.mapc = mapc
            mrc.header.mapr = mapr
            mrc.header.maps = maps
            
            # Handle nstart and origin
            if h is not None:
                # Preserve original header values
                nxstart = int(h.nxstart)
                nystart = int(h.nystart)
                nzstart = int(h.nzstart)
                header_origin_xyz = np.array([h.origin.x, h.origin.y, h.origin.z], dtype=np.float32)
                
                mrc.header.mx = int(h.mx)
                mrc.header.my = int(h.my)
                mrc.header.mz = int(h.mz)
                
                # Preserve cell dimensions
                if hasattr(h, 'cella'):
                    mrc.header.cella.x = float(h.cella.x)
                    mrc.header.cella.y = float(h.cella.y)
                    mrc.header.cella.z = float(h.cella.z)
                if hasattr(h, 'cellb'):
                    mrc.header.cellb.alpha = float(h.cellb.alpha)
                    mrc.header.cellb.beta = float(h.cellb.beta)
                    mrc.header.cellb.gamma = float(h.cellb.gamma)
                # Copy space group if available
                if hasattr(h, 'ispg'):
                    mrc.header.ispg = int(h.ispg)
            else:
                # For new files, calculate nstart from origin
                if np.any(voxel_size <= 0):
                    raise ValueError(f"Voxel size must be positive, got {voxel_size}")
                
                # Set header.origin = 0 and encode everything in nstart
                header_origin_xyz = np.zeros(3, dtype=np.float32)
                nxstart = int(np.round(self.origin[0] / voxel_size[0]))
                nystart = int(np.round(self.origin[1] / voxel_size[1]))
                nzstart = int(np.round(self.origin[2] / voxel_size[2]))
            
            # Set the start positions
            mrc.header.nxstart = nxstart
            mrc.header.nystart = nystart
            mrc.header.nzstart = nzstart
            
            # Set explicit origin
            mrc.header.origin.x = float(header_origin_xyz[0])
            mrc.header.origin.y = float(header_origin_xyz[1])
            mrc.header.origin.z = float(header_origin_xyz[2])
            
            # Update statistics only
            mrc.update_header_stats()
