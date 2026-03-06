# Copyright (c) 2016, Science and Technology Facilities Council
# This software is distributed under a BSD licence. See LICENSE.txt.
"""
mrcobject
---------

Module which exports the :class:`MrcObject` class.

Classes:
    :class:`MrcObject`: An object representing image or volume data in the MRC
    format.

"""

# Import Python 3 features for future-proofing
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from datetime import datetime
import warnings

import numpy as np

from . import utils
from .dtypes import (HEADER_DTYPE, VOXEL_SIZE_DTYPE, NSTART_DTYPE,
                     get_ext_header_dtype)
from .constants import (MAP_ID, IMAGE_STACK_SPACEGROUP, VOLUME_SPACEGROUP,
                        VOLUME_STACK_SPACEGROUP)


class MrcObject(object):

    """An object representing image or volume data in the MRC format.

    The header, extended header and data are stored as numpy arrays and
    exposed as read-only attributes. To replace the data or extended header,
    call :meth:`set_data` or :meth:`set_extended_header`. The header cannot be
    replaced but can be modified in place.

    Voxel size is exposed as a writeable attribute, but is calculated
    on-the-fly from the header's ``cella`` and ``mx``/``my``/``mz`` fields.

    Three-dimensional data can represent either a stack of 2D images, or a 3D
    volume. This is indicated by the header's ``ispg`` (space group) field,
    which is set to 0 for image data and >= 1 for volume data. The
    :meth:`is_single_image`, :meth:`is_image_stack`, :meth:`is_volume` and
    :meth:`is_volume_stack` methods can be used to identify the type of
    information stored in the data array. For 3D data, the
    :meth:`set_image_stack` and :meth:`set_volume` methods can be used to
    switch between image stack and volume interpretations of the data.

    If the data contents have been changed, you can use the
    :meth:`update_header_from_data` and :meth:`update_header_stats` methods to
    make the header consistent with the data. These methods are called
    automatically if the data array is replaced by calling :meth:`set_data`.
    :meth:`update_header_from_data` is fast, even with very large data arrays,
    because it only examines the shape and type of the data array.
    :meth:`update_header_stats` calculates statistics from all items in the
    data array and so can be slow for very large arrays. If necessary, the
    :meth:`reset_header_stats` method can be called to set the header fields to
    indicate that the statistics are undetermined.

    Attributes:

    * :attr:`header`
    * :attr:`extended_header`
    * :attr:`indexed_extended_header`
    * :attr:`data`
    * :attr:`voxel_size`
    * :attr:`nstart`

    Methods:

    * :meth:`set_extended_header`
    * :meth:`set_data`
    * :meth:`is_single_image`
    * :meth:`is_image_stack`
    * :meth:`is_volume`
    * :meth:`is_volume_stack`
    * :meth:`set_image_stack`
    * :meth:`set_volume`
    * :meth:`update_header_from_data`
    * :meth:`update_header_stats`
    * :meth:`reset_header_stats`
    * :meth:`print_header`
    * :meth:`get_labels`
    * :meth:`add_label`

    Attributes and methods relevant to subclasses:

    * ``_read_only``
    * :meth:`_check_writeable`
    * :meth:`_create_default_attributes`
    * :meth:`_close_data`
    * :meth:`_set_new_data`

    """

    def __init__(self, **kwargs):
        """Initialise a new :class:`MrcObject`.

        This initialiser deliberately avoids creating any arrays and simply
        sets the header, extended header and data attributes to :data:`None`.
        This allows subclasses to call :meth:`__init__` at the start of their
        initialisers and then set the attributes themselves, probably by
        reading from a file, or by calling :meth:`_create_default_attributes`
        for a new empty object.

        Note that this behaviour might change in future: this initialiser could
        take optional arguments to allow the header and data to be provided
        by the caller, or might create the standard empty defaults rather than
        setting the attributes to :data:`None`.
        """
        super(MrcObject, self).__init__(**kwargs)

        # Set empty default attributes
        self._header = None
        self._extended_header = None
        self._data = None
        self._read_only = False

    def _check_writeable(self):
        """Check that this MRC object is writeable.

        Raises:
            :exc:`ValueError`: If this object is read-only.
        """
        if self._read_only:
            raise ValueError('MRC object is read-only')

    def _create_default_attributes(self):
        """Set valid default values for the header and data attributes."""
        self._create_default_header()
        self._extended_header = np.empty(0, dtype='V1')
        self._set_new_data(np.empty(0, dtype=np.int8))

    def _create_default_header(self):
        """Create a default MRC file header.

        The header is initialised with standard file type and version
        information, default values for some essential fields, and zeros
        elsewhere. The first text label is also set to indicate the file was
        created by this module.
        """
        self._header = np.zeros(shape=(), dtype=HEADER_DTYPE).view(np.recarray)
        header = self._header
        header.map = MAP_ID
        header.nversion = 20141  # current MRC 2014 format version
        header.machst = utils.machine_stamp_from_byte_order(header.mode.dtype.byteorder)

        # Default space group is P1
        header.ispg = VOLUME_SPACEGROUP

        # Standard cell angles all 90.0 degrees
        default_cell_angle = 90.0
        header.cellb.alpha = default_cell_angle
        header.cellb.beta = default_cell_angle
        header.cellb.gamma = default_cell_angle
        # (this can also be achieved by assigning a 3-tuple to header.cellb
        # directly but using the sub-fields individually is easier to read and
        # understand)

        # Standard axes: columns = X, rows = Y, sections = Z
        header.mapc = 1
        header.mapr = 2
        header.maps = 3

        time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        header.label[0] = '{0:40s}{1:>39s} '.format('Created by mrcfile.py',
                                                    time)
        header.nlabl = 1

        self.reset_header_stats()

    @property
    def header(self):
        """Get the header as a :class:`numpy record array <numpy.recarray>`."""
        return self._header

    @property
    def extended_header(self):
        """Get the extended header as a :class:`numpy array <numpy.ndarray>`.

        The dtype will be void (raw data, dtype ``V'``). If the actual data type
        of the extended header is known, the dtype of the array can be changed
        to match. For supported types (e.g. ``'FEI1'`` and ``'FEI2'``), the
        indexed part of the extended header (excluding any zero padding) can be
        accessed using :meth:`indexed_extended_header`.

        The extended header may be modified in place. To replace it completely,
        call :meth:`set_extended_header`.
        """
        return self._extended_header

    @property
    def indexed_extended_header(self):
        """Get the indexed part of the extended header as a
        :class:`numpy array <numpy.ndarray>` with the appropriate dtype set.

        Currently only ``'FEI1'`` and ``'FEI2'`` extended headers are supported.
        Modifications to the indexed extended header will not change the
        extended header data recorded in this :class:`MrcObject`. If the
        extended header type is unrecognised or extended header data is not of
        sufficient length a warning will be produced and the indexed extended
        header will be None.
        """
        # Use the header's byte order for the extended header
        dtype = get_ext_header_dtype(self.header.exttyp,
                                     self.header.mode.dtype.byteorder)

        # Interpret one element
        try:
            if self.extended_header.nbytes < dtype.itemsize:
                raise ValueError
            first = self.extended_header[0:dtype.itemsize]
            first.dtype = dtype
            if first["Metadata size"][0] != dtype.itemsize:
                raise ValueError
        except ValueError:
            warnings.warn("The header has exttyp '{}' but the extended header "
                            "cannot be interpreted as that type"
                            .format(self.header.exttyp), RuntimeWarning)
            return None

        nbytes = int(self.header["nz"]) * dtype.itemsize
        try:
            if self.extended_header.nbytes < nbytes:
                raise ValueError
            full = self.extended_header[0:nbytes]
            full.dtype = dtype
        except ValueError:
            warnings.warn("The header has exttyp '{}' but the extended header "
                            "cannot be interpreted as that type"
                            .format(self.header.exttyp), RuntimeWarning)
            return None

        return full

    def set_extended_header(self, extended_header):
        """Replace the extended header.

        If you set the extended header you should also set the
        ``header.exttyp`` field to indicate the type of extended header.

        Raises:
            :exc:`ValueError`: If the new extended header has more than
                2,147,483,647 bytes (and therefore its size cannot be stored
                in the header).
        """
        self._check_writeable()
        if extended_header.nbytes > np.iinfo(np.int32).max:
            raise ValueError("New extended header is too large! It has {} "
                             "bytes. The maximum allowed is {}."
                             .format(extended_header.nbytes,
                                     np.iinfo(np.int32).max))
        self._extended_header = extended_header
        self.header.nsymbt = extended_header.nbytes

    @property
    def data(self):
        """Get the data as a :class:`numpy array <numpy.ndarray>`."""
        return self._data

    def set_data(self, data):
        """Replace the data array.

        This replaces the current data with the given array (or a copy of it),
        and updates the header to match the new data dimensions. The data
        statistics (min, max, mean and rms) stored in the header will also be
        updated.

        Raises:
            :exc:`ValueError`: if the new data has a dimension larger than
                2,147,483,647 (and therefore its size cannot be stored in the
                header).

        Warns:
            RuntimeWarning: If the data array contains Inf or NaN values.
        """
        self._check_writeable()

        # Check if the new data's dtype is valid without changes
        mode = utils.mode_from_dtype(data.dtype)
        new_dtype = (utils.dtype_from_mode(mode)
                     .newbyteorder(data.dtype.byteorder))
        
        for dim in data.shape:
            if dim > np.iinfo(np.int32).max:
                raise ValueError("New data array is too large! Found a "
                                 "dimension of size {}. The maximum allowed "
                                 "is {}."
                                 .format(dim, np.iinfo(np.int32).max))
        
        # Set new_dtype to None if it is the same type as the original array,
        # to avoid numpy >= 1.24 copying the array unnecessarily
        if new_dtype == data.dtype:
            new_dtype = None

        # Copy the data if necessary to ensure correct dtype and C ordering
        new_data = np.ascontiguousarray(data, dtype=new_dtype)

        # Replace the old data array with the new one, and update the header
        self._close_data()
        self._set_new_data(new_data)
        self.update_header_from_data()
        self.update_header_stats()

    def _close_data(self):
        """Close the data array."""
        self._data = None

    def _set_new_data(self, data):
        """Replace the data array with a new one.

        The new data array is not checked - it must already be valid for use in
        an MRC file.
        """
        self._data = data

    @property
    def voxel_size(self):
        """Get or set the voxel size in angstroms.

        The voxel size is returned as a structured NumPy :class:`record array
        <numpy.recarray>` with three fields (x, y and z). For example:

        >>> mrc.voxel_size
        rec.array((0.44825, 0.3925, 0.45874998),
          dtype=[('x', '<f4'), ('y', '<f4'), ('z', '<f4')])
        >>> mrc.voxel_size.x
        array(0.44825, dtype=float32)

        Note that changing the voxel_size array in-place will *not* change the
        voxel size in the file -- to prevent this being overlooked
        accidentally, the writeable flag is set to :data:`False` on the
        voxel_size array.

        To set the voxel size, assign a new value to the voxel_size attribute.
        You may give a single number, a 3-tuple ``(x, y ,z)`` or a modified
        version of the voxel_size array. The following examples are all
        equivalent:

        >>> mrc.voxel_size = 1.0

        >>> mrc.voxel_size = (1.0, 1.0, 1.0)

        >>> vox_sizes = mrc.voxel_size
        >>> vox_sizes.flags.writeable = True
        >>> vox_sizes.x = 1.0
        >>> vox_sizes.y = 1.0
        >>> vox_sizes.z = 1.0
        >>> mrc.voxel_size = vox_sizes
        """
        x = self.header.cella.x / self.header.mx
        y = self.header.cella.y / self.header.my
        z = self.header.cella.z / self.header.mz
        sizes = np.rec.array((x, y, z), VOXEL_SIZE_DTYPE)
        sizes.flags.writeable = False
        return sizes

    @voxel_size.setter
    def voxel_size(self, voxel_size):
        self._check_writeable()
        try:
            # First, assume we have a single numeric value
            sizes = (float(voxel_size),) * 3
        except TypeError:
            try:
                # Not a single value. Next, if voxel_size is an array (as
                # produced by the voxel_size getter), item() gives a 3-tuple
                sizes = voxel_size.item()
            except AttributeError:
                # If the item() method doesn't exist, assume we have a 3-tuple
                sizes = voxel_size
        self._set_voxel_size(*sizes)

    def _set_voxel_size(self, x_size, y_size, z_size):
        """Set the voxel size.

        Args:
            x_size: The voxel size in the X direction, in angstroms
            y_size: The voxel size in the Y direction, in angstroms
            z_size: The voxel size in the Z direction, in angstroms
        """
        self.header.cella.x = x_size * self.header.mx
        self.header.cella.y = y_size * self.header.my
        self.header.cella.z = z_size * self.header.mz

    @property
    def nstart(self):
        """Get or set the grid start locations.

        This provides a convenient way to get and set the values of the
        header's ``nxstart``, ``nystart`` and ``nzstart`` fields. Note that
        these fields are integers and are measured in voxels, not angstroms.
        The start locations are returned as a structured NumPy :class:`record
        array <numpy.recarray>` with three fields (x, y and z). For example:

        >>> mrc.header.nxstart
        array(0, dtype=int32)
        >>> mrc.header.nystart
        array(-21, dtype=int32)
        >>> mrc.header.nzstart
        array(-12, dtype=int32)
        >>> mrc.nstart
        rec.array((0, -21, -12),
          dtype=[('x', '<i4'), ('y', '<i4'), ('z', '<i4')])
        >>> mrc.nstart.y
        array(-21, dtype=int32)

        Note that changing the nstart array in-place will *not* change the
        values in the file -- to prevent this being overlooked accidentally,
        the writeable flag is set to :data:`False` on the nstart array.

        To set the start locations, assign a new value to the nstart
        attribute. You may give a single number, a 3-tuple ``(x, y ,z)`` or a
        modified version of the nstart array. The following examples are all
        equivalent:

        >>> mrc.nstart = -150

        >>> mrc.nstart = (-150, -150, -150)

        >>> starts = mrc.nstart
        >>> starts.flags.writeable = True
        >>> starts.x = -150
        >>> starts.y = -150
        >>> starts.z = -150
        >>> mrc.nstart = starts
        """
        x = self.header.nxstart
        y = self.header.nystart
        z = self.header.nzstart
        nstart = np.rec.array((x, y, z), NSTART_DTYPE)
        nstart.flags.writeable = False
        return nstart

    @nstart.setter
    def nstart(self, nstart):
        self._check_writeable()
        try:
            # First, assume we have a single numeric value
            starts = (int(nstart),) * 3
        except TypeError:
            try:
                # Not a single value. Next, if nstart is an array (as
                # produced by the nstart getter), item() gives a 3-tuple
                starts = nstart.item()
            except AttributeError:
                # If the item() method doesn't exist, assume we have a 3-tuple
                starts = nstart
        self._set_nstart(*starts)

    def _set_nstart(self, nxstart, nystart, nzstart):
        """Set the grid start locations.

        Args:
            nxstart: The location of the first column in the unit cell
            nystart: The location of the first row in the unit cell
            nzstart: The location of the first section in the unit cell
        """
        self.header.nxstart = nxstart
        self.header.nystart = nystart
        self.header.nzstart = nzstart

    def is_single_image(self):
        """Identify whether the file represents a single image.

        Returns:
            :data:`True` if the data array is two-dimensional.
        """
        return self.data.ndim == 2

    def is_image_stack(self):
        """Identify whether the file represents a stack of images.

        Returns:
            :data:`True` if the data array is three-dimensional and the space group
            is zero.
        """
        return (self.data.ndim == 3
                and self.header.ispg == IMAGE_STACK_SPACEGROUP)

    def is_volume(self):
        """Identify whether the file represents a volume.

        Returns:
            :data:`True` if the data array is three-dimensional and the space
            group is not zero.
        """
        return (self.data.ndim == 3
                and self.header.ispg != IMAGE_STACK_SPACEGROUP)

    def is_volume_stack(self):
        """Identify whether the file represents a stack of volumes.

        Returns:
            :data:`True` if the data array is four-dimensional.
        """
        return self.data.ndim == 4

    def set_image_stack(self):
        """Change three-dimensional data to represent an image stack.

        This method changes the space group number (``header.ispg``) to zero.

        Raises:
            :exc:`ValueError`: If the data array is not three-dimensional.
        """
        self._check_writeable()
        if self.data.ndim != 3:
            raise ValueError('Only 3D data can be changed into an image stack')
        self.header.ispg = IMAGE_STACK_SPACEGROUP
        self.header.mz = 1

    def set_volume(self):
        """Change three-dimensional data to represent a volume.

        If the space group was previously zero (representing an image stack),
        this method sets it to one. Otherwise the space group is not changed.

        Raises:
            :exc:`ValueError`: If the data array is not three-dimensional.
        """
        self._check_writeable()
        if self.data.ndim != 3:
            raise ValueError('Only 3D data can be changed into a volume')
        if self.is_image_stack():
            self.header.ispg = VOLUME_SPACEGROUP
            self.header.mz = self.header.nz

    def update_header_from_data(self):
        """Update the header from the data array.

        This function updates the header byte order and machine stamp to match
        the byte order of the data. It also updates the file mode, space group
        and the dimension fields ``nx``, ``ny``, ``nz``, ``mx``, ``my`` and
        ``mz``.

        If the data is 2D, the space group is set to 0 (image stack). For 3D
        data the space group is not changed, and for 4D data the space group is
        set to 401 (simple P1 volume stack) unless it is already in the volume
        stack range (401--630).

        This means that new 3D data will be treated as an image stack if the
        previous data was a single image or image stack, or as a volume if the
        previous data was a volume or volume stack.

        Note that this function does *not* update the data statistics fields in
        the header (``dmin``, ``dmax``, ``dmean`` and ``rms``). Use the
        :meth:`update_header_stats` function to update the statistics.
        (This is for performance reasons -- updating the statistics can take a
        long time for large data sets, but updating the other header
        information is always fast because only the type and shape of the data
        array need to be inspected.)
        """
        self._check_writeable()

        # Check the dtype is one we can handle and update mode to match
        header = self.header
        header.mode = utils.mode_from_dtype(self.data.dtype)

        # Ensure header byte order and machine stamp match the data's byte order
        data_byte_order = self.data.dtype.byteorder
        header_byte_order = header.mode.dtype.byteorder
        if (data_byte_order != '|'
            and not utils.byte_orders_equal(data_byte_order, header_byte_order)):
            header.byteswap(True)
            header.dtype = header.dtype.newbyteorder(data_byte_order)
        header.machst = utils.machine_stamp_from_byte_order(header.mode.dtype
                                                            .byteorder)

        shape = self.data.shape
        axes = len(shape)
        if axes == 2:
            # Single image. Space group 0, nz = mz = 1
            header.ispg = IMAGE_STACK_SPACEGROUP
            header.nx = header.mx = shape[1]
            header.ny = header.my = shape[0]
            header.nz = header.mz = 1
        elif axes == 3:
            header.nx = header.mx = shape[2]
            header.ny = header.my = shape[1]
            if header.ispg == IMAGE_STACK_SPACEGROUP:
                # Image stack. mz = 1, nz = sections in the volume
                header.mz = 1
                header.nz = shape[0]
            else:
                # Volume. nz = mz = sections in the volume
                header.nz = header.mz = shape[0]
        elif axes == 4:
            # Volume stack. Space group 401, mz = secs per vol, nz = total sections
            if not utils.spacegroup_is_volume_stack(header.ispg):
                header.ispg = VOLUME_STACK_SPACEGROUP
            header.nx = header.mx = shape[3]
            header.ny = header.my = shape[2]
            header.mz = shape[1]
            header.nz = shape[0] * shape[1]
        else:
            raise ValueError('Data must be 2-, 3- or 4-dimensional')

    def update_header_stats(self):
        """Update the header's ``dmin``, ``dmax``, ``dmean`` and ``rms`` fields
        from the data.

        Note that this can take some time with large files, particularly with
        files larger than the currently available memory.

        Warns:
            RuntimeWarning: If the data array contains Inf or NaN values.
        """
        self._check_writeable()

        if self.data.size > 0:
            # Header stats are always in float32. If we have complex data, this doesn't
            # make sense so just set rms and leave min, max and mean at their default
            # un-set values
            if np.iscomplexobj(self.data):
                # Avoid ComplexWarning by explicitly taking the real part
                self.header.rms = np.float32(self.data.std().real)
            else:
                min = self.data.min()
                max = self.data.max()

                if np.isnan(min):
                    warnings.warn("Data array contains NaN values", RuntimeWarning)
                if np.isinf(min) or np.isinf(max):
                    warnings.warn("Data array contains infinite values", RuntimeWarning)

                self.header.dmin = np.float32(min)
                self.header.dmax = np.float32(max)
                self.header.dmean = self.data.mean(dtype=np.float32)
                self.header.rms = self.data.std(dtype=np.float32)
        else:
            self.reset_header_stats()

    def reset_header_stats(self):
        """Set the header statistics to indicate that the values are unknown."""
        self._check_writeable()

        self.header.dmin = 0
        self.header.dmax = -1
        self.header.dmean = -2
        self.header.rms = -1

    def print_header(self, print_file=None):
        """Print the contents of all header fields.

        Args:
            print_file: The output text stream to use for printing the header.
                This is passed directly to the ``file`` argument of Python's
                :func:`print` function. The default is :data:`None`, which
                means output will be printed to :data:`sys.stdout`.
        """
        for item in self.header.dtype.names:
            print('{0:15s} : {1}'.format(item, self.header[item]),
                  file=print_file)

    def get_labels(self):
        """Get the labels from the MRC header.

        Up to ten labels are stored in the header as arrays of 80 bytes. This method
        returns the labels as Python strings, filtered to remove non-printable
        characters. To access the raw bytes (including any non-printable characters) use
        the ``header.label`` attribute (and note that ``header.nlabl`` stores the number
        of labels currently set).

        Returns:
            The labels, as a list of strings. The list will contain between 0 and 10
            items, each containing up to 80 characters.
        """
        return [
            utils.printable_string_from_bytes(label)
            for label in self.header.label[:self.header.nlabl]
        ]

    def add_label(self, label):
        """Add a label to the MRC header.

        The new label will be stored after any labels already in the header. If all ten
        labels are already in use, an exception will be raised.

        Future versions of this method might add checks to ensure that labels containing
        valid text are not overwritten even if the ``nlabl`` value is incorrect.

        Args:
            label: The label value to store, as a string containing only printable
                ASCII characters.

        Raises:
            :exc:`ValueError`: If the label is longer than 80 bytes or contains
                non-printable or non-ASCII characters.
            :exc:`IndexError`: If the file already contains 10 labels and so an
                additional label cannot be stored.
        """
        if not utils.is_printable_ascii(label):
            raise ValueError("Label contains non-printable or non-ASCII characters")
        label_bytes = utils.bytes_from_string(label)
        if len(label_bytes) > 80:
            raise ValueError("Label value has more than 80 bytes")
        self.header.label[self.header.nlabl] = label
        self.header.nlabl += 1

    def validate(self, print_file=None):
        """Validate this MrcObject.

        This method runs a series of tests to check whether this object
        complies strictly with the MRC2014 format specification:

        #. MRC format ID string: The header's ``map`` field must contain
           "MAP ".
        #. Machine stamp: The machine stamp should contain one of
           ``0x44 0x44 0x00 0x00``, ``0x44 0x41 0x00 0x00`` or
           ``0x11 0x11 0x00 0x00``.
        #. MRC mode: the ``mode`` field should be one of the supported mode
           numbers: 0, 1, 2, 4, 6 or 12. (Note that MRC modes 3 and 101 are
           also valid according to the MRC 2014 specification but are not
           supported by mrcfile.)
        #. Map and cell dimensions: The header fields ``nx``, ``ny``, ``nz``,
           ``mx``, ``my``, ``mz``, ``cella.x``, ``cella.y`` and ``cella.z``
           must all be positive numbers.
        #. Axis mapping: Header fields ``mapc``, ``mapr`` and ``maps`` must
           contain the values 1, 2, and 3 (in any order).
        #. Volume stack dimensions: If the spacegroup is in the range 401--630,
           representing a volume stack, the ``nz`` field should be exactly
           divisible by ``mz`` to represent the number of volumes in the stack.
        #. Header labels: The ``nlabl`` field should be set to indicate the
           number of labels in use, and the labels in use should appear first
           in the label array (that is, there should be no blank labels between
           text-filled ones).
        #. MRC format version: The ``nversion`` field should be 20140 or 20141
           for compliance with the MRC2014 standard.
        #. Extended header type: If an extended header is present, the
           ``exttyp`` field should be set to indicate the type of extended
           header.
        #. Data statistics: The statistics in the header should be correct for
           the actual data, or marked as undetermined.

        Args:
            print_file: The output text stream to use for printing messages
                about the validation. This is passed directly to the ``file``
                argument of Python's :func:`print` function. The default is
                :data:`None`, which means output will be printed to
                :data:`sys.stdout`.

        Returns:
            :data:`True` if this MrcObject  is valid, or :data:`False` if it
            does not meet the MRC format specification in any way.
        """
        valid = True

        def log(message):
            print(message, file=print_file)

        # Check map ID string
        if self.header.map != MAP_ID:
            log("Map ID string is incorrect: found {0}, should be {1}"
                .format(self.header.map, MAP_ID))
            valid = False

        # Check machine stamp
        try:
            utils.byte_order_from_machine_stamp(self.header.machst)
        except ValueError:
            pretty_bytes = utils.pretty_machine_stamp(self.header.machst)
            log("Invalid machine stamp: " + pretty_bytes)
            valid = False

        # Check mode is valid
        try:
            utils.dtype_from_mode(self.header.mode)
        except ValueError:
            log("Invalid mode: {0}".format(self.header.mode))
            valid = False

        # Check map dimensions and other fields are non-negative
        for field in ['nx', 'ny', 'nz', 'mx', 'my', 'mz', 'ispg', 'nlabl']:
            if self.header[field] < 0:
                log("Header field '{0}' is negative".format(field))
                valid = False

        # Check cell dimensions are non-negative
        for field in ['x', 'y', 'z']:
            if self.header.cella[field] < 0:
                log("Cell dimension '{0}' is negative".format(field))
                valid = False

        # Check axis mapping is valid
        axes = set()
        for field in ['mapc', 'mapr', 'maps']:
            axes.add(int(self.header[field]))
        if axes != {1, 2, 3}:
            log("Invalid axis mapping: found {0}, should be [1, 2, 3]"
                .format(sorted(list(axes))))
            valid = False

        # Check mz value for volume stacks
        if utils.spacegroup_is_volume_stack(self.header.ispg):
            if self.header.nz % self.header.mz != 0:
                log("Error in dimensions for volume stack: nz should be "
                    "divisible by mz. Found nz = {0}, mz = {1})"
                    .format(self.header.nz, self.header.mz))
                valid = False

        # Check nlabl is correct
        count = 0
        seen_empty_label = False
        for label in self.header.label:
            if len(label.strip()) > 0:
                count += 1
                if seen_empty_label:
                    log("Error in header labels: empty labels appear between "
                        "text-containing labels")
                    valid = False
            else:
                seen_empty_label = True
        if count != self.header.nlabl:
            log("Error in header labels: nlabl is {0} "
                "but {1} labels contain text".format(self.header.nlabl, count))
            valid = False

        # Check MRC format version
        if self.header.nversion not in (20140, 20141):
            log("File does not declare MRC format version 20140 or 20141: nversion ="
                " {0}".format(self.header.nversion))
            valid = False

        # Check extended header type is set to a known value
        valid_exttypes = [b'CCP4', b'MRCO', b'SERI', b'AGAR', b'FEI1', b'FEI2', b'HDF5']
        if self.header.nsymbt > 0 and self.header.exttyp not in valid_exttypes:
            log("Extended header type is undefined or unrecognised: exttyp = "
                "'{0}'".format(self.header.exttyp.item().decode('ascii')))
            valid = False

        # Check data statistics
        real_rms = real_min = real_max = real_mean = 0
        if self.header.rms >= 0:
            if self.data is not None and len(self.data > 0):
                real_rms = self.data.std()
            if not np.isclose(real_rms, self.header.rms, rtol=0.01):
                log("Data statistics appear to be inaccurate: RMS deviation is {0} but"
                    " the value in the header is {1}".format(real_rms, self.header.rms))
                valid = False
        if self.header.dmin < self.header.dmax:
            if self.data is not None and len(self.data > 0):
                real_min = self.data.min()
                real_max = self.data.max()
            if self.header.dmin != real_min:
                log("Data statistics appear to be inaccurate: minimum is {0} but the"
                    " value in the header is {1}".format(real_min, self.header.dmin))
                valid = False
            if self.header.dmax != real_max:
                log("Data statistics appear to be inaccurate: maximum is {0} but the"
                    " value in the header is {1}".format(real_max, self.header.dmax))
                valid = False
        if self.header.dmean > min(self.header.dmin, self.header.dmax):
            if self.data is not None and len(self.data > 0):
                real_mean = self.data.mean(dtype=np.float64)
            if not np.isclose(real_mean, self.header.dmean, rtol=0.01):
                log("Data statistics appear to be inaccurate: mean is {0} but the"
                    " value in the header is {1}".format(real_mean, self.header.dmean))
                valid = False

        return valid
