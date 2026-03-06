"""
This module contains classes for reading and writing Amber NetCDF-style files,
including both restarts and trajectories. The NetCDF engine used here is pulled
from the scipy distribution, and depends *only* on numpy. It is available
through parmed.utils.netcdf.

This module contains objects relevant to Amber NetCDF files. The use() function
is responsible for selecting the API based on a default choice or user-selection
(the latter is really only helpful for development to ensure that all packages
work correctly---there is no difference from a user perspective). ALL
NetCDF-file manipulation that the parmed/amber package does should be
contained in this module.
"""
try:
    # netCDF4 is *much* faster to write NetCDF files, and can make a huge
    # difference when using it as an OpenMM reporter. So do what we can to use
    # the faster library when available
    import netCDF4 as nc
except ImportError:
    nc = None
import numpy as np
from .. import __version__
from ..formats.registry import FileFormatType
from .. import unit as u
from ..utils.netcdf import netcdf_file as NetCDFFile
import warnings

class NetCDFRestart(metaclass=FileFormatType):
    """ Class to read or write NetCDF restart files """

    @staticmethod
    def id_format(filename):
        """ Identifies the file type as an Amber NetCDF restart file

        Parameters
        ----------
        filename : str
            Name of the file to check format for

        Returns
        -------
        is_fmt : bool
            True if it is an Amber NetCDF restart file. False otherwise

        Notes
        -----
        Remote NetCDF files cannot be loaded
        """
        if filename.startswith('http://') or filename.startswith('https://')\
                or filename.startswith('ftp://'):
            return False
        try:
            f = NetCDFFile(filename, 'r', mmap=False)
        except (TypeError, OSError):
            return False
        try:
            try:
                if f.Conventions.decode() != 'AMBERRESTART':
                    return False
            except AttributeError:
                return False
            # Passed all our tests
            return True
        finally:
            f.close()

    def __init__(self, fname, mode='r'):
        """
        Opens a NetCDF File. The main constructor should never be called
        directly.  The alternative "open_old" and "open_new" constructors should
        be called instead for reading existing or writing new NetCDF files,
        respectively.
        """
        self.closed = False
        if mode.startswith('w') and nc is not None:
            self._ncfile = nc.Dataset(fname, mode, format='NETCDF3_64BIT')
        else:
            if mode.startswith('w'):
                warnings.warn('Could not find netCDF4 module. Falling back on '
                              'scipy implementation, which can significantly '
                              'slow down simulations if used as a reporter')
            self._ncfile = NetCDFFile(fname, mode, mmap=False)

    @classmethod
    def open_new(cls, fname, natom, box, vels, title='',
                 remd=None, temp=None, remd_dimtypes=None):
        """
        Opens a new NetCDF file and sets the attributes

        Parameters
        ----------
        fname : str
            Name of the new file to open (overwritten)
        natom : int
            The number of atoms in the system
        box : bool
            Whether unit cell information is written or not
        vels : bool
            Whether velocity information is written or not
        title : str=''
            The title to write to the NetCDF restart file
        remd : str=None
            None -- No REMD information is written
            'T[emperature]' -- target temperature (or pH) will be written
            'M[ulti-D]' -- remd_dimtypes will be written
        remd_dimtypes : iterable of int=None
            Array of exchange types for each group. The length will be the REMD
            dimension (if `remd` above is "M[ulti-D]")

        Notes
        -----
        `remd` is case-insensitive, and done based on first-letter matching
        """
        if remd is not None:
            if remd[0] in 'tT':
                remd_type = 'TEMPERATURE'
                if temp is None:
                    raise ValueError('temp must be specified for T-REMD restarts.')
            elif remd[0] in 'mM':
                remd_type = 'MULTI'
                if remd_dimtypes is None:
                    raise ValueError('remd_dimtypes must be given for multi-D '
                                     'REMD, and must have the same length.')
                for dt in remd_dimtypes:
                    if dt not in (1, 3):
                        raise ValueError(
                            'remd_dimtypes only supports dimension types 1 and 3 currently'
                        )
                remd_dimension = len(remd_dimtypes)
            else:
                raise ValueError('remd must be None, T[emperature] or M[ulti]')
        else:
            remd_type = None
        inst = cls(fname, 'w')
        ncfile = inst._ncfile
        inst.hasbox = bool(box)
        inst.hasvels = bool(vels)
        # Assign the main attributes
        ncfile.Conventions = 'AMBERRESTART'
        ncfile.ConventionVersion = "1.0"
        ncfile.title = str(title) # Cast to avoid ScientificPython segfault
        ncfile.application = "AmberTools"
        ncfile.program = "ParmEd"
        ncfile.programVersion = str(__version__)
        # Make all of the dimensions
        ncfile.createDimension('spatial', 3)
        ncfile.createDimension('atom', natom)
        inst.spatial = 3
        inst.atom = natom
        inst.title = ncfile.title
        if box:
            ncfile.createDimension('cell_spatial', 3)
            ncfile.createDimension('label', 5)
            ncfile.createDimension('cell_angular', 3)
            inst.cell_spatial = 3
            inst.label = 5
            inst.cell_angular = 3
        if remd_type == 'MULTI':
            ncfile.createDimension('remd_dimension', remd_dimension)
            inst.remd_dimension = remd_dimension
        ncfile.createDimension('time', 1)
        # Now make the variables
        v = ncfile.createVariable('time', 'd', ('time',))
        v.units = 'picosecond'
        v = ncfile.createVariable('spatial', 'c', ('spatial',))
        v[:] = np.asarray(list('xyz'))
        if inst.hasbox:
            v = ncfile.createVariable('cell_angular', 'c', ('cell_angular', 'label'))
            v[0] = np.asarray(list('alpha'))
            v[1] = np.asarray(list('beta '))
            v[2] = np.asarray(list('gamma'))
            v = ncfile.createVariable('cell_spatial', 'c', ('cell_spatial',))
            v[0], v[1], v[2] = 'a', 'b', 'c'
            v = ncfile.createVariable('cell_lengths', 'd', ('cell_spatial',))
            v.units = 'angstrom'
            v = ncfile.createVariable('cell_angles', 'd', ('cell_angular',))
            v.units = 'degree'
        v = ncfile.createVariable('coordinates', 'd', ('atom', 'spatial'))
        v.units = 'angstrom'
        if inst.hasvels:
            v = ncfile.createVariable('velocities', 'd', ('atom', 'spatial'))
            v.units = 'angstrom/picosecond'
            v.scale_factor = np.float32(20.455)
            inst.velocity_scale = 20.455
            if nc is not None:
                v.set_auto_maskandscale(False)

        if remd_type == 'TEMPERATURE':
            v = ncfile.createVariable('temp0', 'd', ('time',))
            v.units = 'kelvin'
            v[0] = temp
        elif remd_type == 'MULTI':
            ncfile.createVariable('remd_indices', 'i', ('remd_dimension',))
            v = ncfile.createVariable('remd_dimtype', 'i', ('remd_dimension',))
            v[:] = remd_dimtypes

        return inst

    @classmethod
    def open_old(cls, fname):
        """
        Opens the NetCDF file and sets the global attributes that the file sets

        Parameters
        ----------
        fname : str
            Name of the file to read
        """
        inst = cls(fname, 'r')
        ncfile = inst._ncfile
        inst.Conventions = ncfile.Conventions.decode()
        inst.ConventionVersion = ncfile.ConventionVersion.decode()
        inst.program = ncfile.program.decode()
        inst.programVersion = ncfile.programVersion.decode()
        if hasattr(ncfile, 'application'):
            inst.application = ncfile.application.decode()
        else:
            inst.application = None
        if hasattr(ncfile, 'title'):
            inst.title = ncfile.title.decode()
        else:
            inst.title = None
        # Set up the dimensions as attributes
        for dim in ncfile.dimensions:
            # Exception for ParmEd-created ncrst files
            if dim == 'time': continue
            setattr(inst, dim, ncfile.dimensions[dim])
        inst.hasvels = 'velocities' in ncfile.variables
        inst.hasbox = ('cell_lengths' in ncfile.variables and 'cell_angles' in ncfile.variables)
        if inst.hasvels:
            vels = ncfile.variables['velocities']
            inst.velocity_scale = vels.scale_factor
        return inst

    parse = open_old

    @property
    def coordinates(self):
        coords = self._ncfile.variables['coordinates'][:]
        return coords.reshape((-1, self.atom, 3))

    @coordinates.setter
    def coordinates(self, stuff):
        stuff = np.asarray(stuff).reshape((self.atom, 3))
        self._ncfile.variables['coordinates'][:] = stuff
        self.flush()

    @property
    def velocities(self):
        if 'velocities' in self._ncfile.variables:
            vels = self._ncfile.variables['velocities'][:]
            return (vels.reshape((-1, self.atom, 3)) * self.velocity_scale)

    @velocities.setter
    def velocities(self, stuff):
        self._ncfile.variables['velocities'][:] = \
                    np.reshape(stuff, (self.atom, 3)) / self.velocity_scale
        self.flush()

    @property
    def cell_lengths(self):
        if 'cell_lengths' in self._ncfile.variables:
            return self._ncfile.variables['cell_lengths'][:]

    @cell_lengths.setter
    def cell_lengths(self, stuff):
        self._ncfile.variables['cell_lengths'][:] = np.asarray(stuff)
        self.flush()

    @property
    def cell_angles(self):
        if 'cell_angles' in self._ncfile.variables:
            return self._ncfile.variables['cell_angles'][:]

    @cell_angles.setter
    def cell_angles(self, stuff):
        self._ncfile.variables['cell_angles'][:] = np.asarray(stuff)
        self.flush()

    @property
    def box(self):
        if self.cell_lengths is not None and self.cell_angles is not None:
            leng, ang = self.cell_lengths, self.cell_angles
            return np.concatenate((leng, ang))

    @box.setter
    def box(self, stuff):
        self.cell_lengths = stuff[:3]
        self.cell_angles = stuff[3:]

    @property
    def time(self):
        return self._ncfile.variables['time'].getValue()

    @time.setter
    def time(self, stuff):
        self._ncfile.variables['time'][0] = float(stuff)
        self.flush()

    @property
    def temp0(self):
        return self._ncfile.variables['temp0'].getValue()

    @property
    def remd_indices(self):
        return self._ncfile.variables['remd_indices'][:]

    @remd_indices.setter
    def remd_indices(self, stuff):
        self._ncfile.variables['remd_indices'][:] = np.asarray(stuff, dtype='i')
        self.flush()

    @property
    def remd_dimtype(self):
        return self._ncfile.variables['remd_dimtype'][:]

    def close(self):
        if not self.closed:
            self.closed = True
            self._ncfile.close()

    def __del__(self):
        self.closed or (hasattr(self, '_ncfile') and self._ncfile.close())

    def flush(self):
        try:
            self._ncfile.flush()
        except AttributeError:
            # netCDF4.Dataset's flush method is called sync :-P
            self._ncfile.sync()

class NetCDFTraj(metaclass=FileFormatType):
    """ Class to read or write NetCDF restart files

    Parameters
    ----------
    fname : str
        Name of the file to open
    mode : str
        Mode to open in:
            - 'w' means write-mode
            - 'r' means read-mode

    Notes
    -----
    You should use the open_new and open_old alternative constructors instead of
    the default constructor
    """

    @staticmethod
    def id_format(filename):
        """ Identifies the file type as an Amber NetCDF trajectory file

        Parameters
        ----------
        filename : str
            Name of the file to check format for

        Returns
        -------
        is_fmt : bool
            True if it is an Amber NetCDF trajectory file. False otherwise

        Notes
        -----
        Remote NetCDF files cannot be loaded
        """
        if filename.startswith('http://') or filename.startswith('https://')\
                or filename.startswith('ftp://'):
            return False
        try:
            f = NetCDFFile(filename, 'r', mmap=False)
        except (TypeError, OSError):
            return False
        try:
            try:
                if f.Conventions.decode() != 'AMBER':
                    return False
            except AttributeError:
                return False
            # Passed all our tests
            return True
        finally:
            f.close()

    def __init__(self, fname, mode='r'):
        """ Opens a NetCDF File """
        self.closed = False
        if mode.startswith('w') and nc is not None:
            self._ncfile = nc.Dataset(fname, mode, format='NETCDF3_64BIT')
        else:
            if mode.startswith('w'):
                warnings.warn('Could not find netCDF4 module. Falling back on '
                              'scipy implementation, which can significantly '
                              'slow down simulations if used as a reporter')
            self._ncfile = NetCDFFile(fname, mode, mmap=False)

    @classmethod
    def open_new(cls, fname, natom, box, crds=True, vels=False, frcs=False,
                 remd=None, remd_dimension=None, title=''):
        """
        Opens a new NetCDF file and sets the attributes

        Parameters
        ----------
        fname : str
            Name of the new file to open (overwritten)
        natom : int
            Number of atoms in the restart
        box : bool
            Indicates if cell lengths and angles are written to the NetCDF file
        crds : bool=True
            Indicates if coordinates are written to the NetCDF file
        vels : bool=False
            Indicates if velocities are written to the NetCDF file
        frcs : bool=False
            Indicates if forces are written to the NetCDF file
        remd : str=None
            'T[emperature]' if replica temperature is written
            'M[ulti]' if Multi-D REMD information is written
            None if no REMD information is written
        remd_dimension : int=None
            If remd above is 'M[ulti]', this is how many REMD dimensions exist
        title : str=''
            The title of the NetCDF trajectory file
        """
        inst = cls(fname, 'w')
        ncfile = inst._ncfile
        if remd is not None:
            if remd[0] in 'Tt':
                inst.remd = 'TEMPERATURE'
            elif remd[0] in 'Mm':
                inst.remd = 'MULTI'
                if remd_dimension is None:
                    raise ValueError('remd_dimension must be given for multi-D REMD')
                inst.remd_dimension = int(remd_dimension)
            else:
                raise ValueError('remd must be T[emperature] or M[ultiD]')
        else:
            inst.remd = None
        inst.hasbox = bool(box)
        inst.hasvels = bool(vels)
        inst.hascrds = bool(crds)
        inst.hasfrcs = bool(frcs)
        # Assign the main attributes
        ncfile.Conventions = "AMBER"
        ncfile.ConventionVersion = "1.0"
        ncfile.application = "AmberTools"
        ncfile.program = "ParmEd"
        ncfile.programVersion = __version__
        ncfile.title = "ParmEd-created trajectory"
        inst.Conventions = "AMBER"
        inst.ConventionVersion = "1.0"
        inst.application = "AmberTools"
        inst.program = "ParmEd"
        inst.programVersion = __version__
        inst.title = ncfile.title
        # Create the dimensions
        ncfile.createDimension('frame', None)
        ncfile.createDimension('spatial', 3)
        ncfile.createDimension('atom', natom)
        if inst.remd == 'MULTI':
            ncfile.createDimension('remd_dimension', inst.remd_dimension)
        inst.frame, inst.spatial, inst.atom = None, 3, natom
        if inst.hasbox:
            ncfile.createDimension('cell_spatial', 3)
            ncfile.createDimension('cell_angular', 3)
            ncfile.createDimension('label', 5)
            inst.cell_spatial, inst.cell_angular, inst.label = 3, 3, 5
        # Create the variables and assign units and scaling factors
        v = ncfile.createVariable('spatial', 'c', ('spatial',))
        v[:] = np.asarray(list('xyz'))
        if inst.hasbox:
            v = ncfile.createVariable('cell_spatial', 'c', ('cell_spatial',))
            v[:] = np.asarray(list('abc'))
            v = ncfile.createVariable('cell_angular', 'c', ('cell_angular', 'label',))
            v[:] = np.asarray([list('alpha'), list('beta '), list('gamma')])
        v = ncfile.createVariable('time', 'f', ('frame',))
        v.units = 'picosecond'
        if inst.hascrds:
            v = ncfile.createVariable('coordinates', 'f', ('frame', 'atom', 'spatial'))
            v.units = 'angstrom'
            inst._last_crd_frame = 0
        if inst.hasvels:
            v = ncfile.createVariable('velocities', 'f', ('frame', 'atom', 'spatial'))
            v.units = 'angstrom/picosecond'
            inst.velocity_scale = v.scale_factor = 20.455
            inst._last_vel_frame = 0
            if nc is not None:
                v.set_auto_maskandscale(False)
        if inst.hasfrcs:
            v = ncfile.createVariable('forces', 'f', ('frame', 'atom', 'spatial'))
            v.units = 'kilocalorie/mole/angstrom'
            inst._last_frc_frame = 0
        if inst.hasbox:
            v = ncfile.createVariable('cell_lengths', 'd', ('frame', 'cell_spatial'))
            v.units = 'angstrom'
            v = ncfile.createVariable('cell_angles', 'd', ('frame', 'cell_angular'))
            v.units = 'degree'
            inst._last_box_frame = 0
        if inst.remd == 'TEMPERATURE':
            v = ncfile.createVariable('temp0', 'd', ('frame',))
            v.units = 'kelvin'
            inst._last_remd_frame = 0
        elif inst.remd == 'MULTI':
            ncfile.createVariable('remd_indices', 'i', ('frame', 'remd_dimension'))
            ncfile.createVariable('remd_dimtype', 'i', ('remd_dimension',))
            inst._last_remd_frame = 0

        inst._last_time_frame = 0

        return inst

    @classmethod
    def open_old(cls, fname):
        """
        Opens the NetCDF file and sets the global attributes that the file sets

        Parameters
        ----------
        fname : str
            File name of the trajectory to open. It must exist
        """
        inst = cls(fname, 'r')
        ncfile = inst._ncfile
        inst.Conventions = ncfile.Conventions.decode()
        inst.ConventionVersion = ncfile.ConventionVersion.decode()
        inst.program = ncfile.program.decode()
        inst.programVersion = ncfile.programVersion.decode()
        if hasattr(ncfile, 'application'):
            inst.application = ncfile.application.decode()
        else:
            inst.application = None
        if hasattr(ncfile, 'title'):
            inst.title = ncfile.title.decode()
        else:
            inst.title = None
        # Set up the dimensions as attributes
        for dim in ncfile.dimensions:
            setattr(inst, dim, ncfile.dimensions[dim])
        inst.hascrds = 'coordinates' in ncfile.variables
        inst.hasvels = 'velocities' in ncfile.variables
        inst.hasfrcs = 'forces' in ncfile.variables
        inst.hasbox = ('cell_lengths' in ncfile.variables and 'cell_angles' in ncfile.variables)
        if inst.hascrds:
            inst._coordinates = np.array(ncfile.variables['coordinates'][:])
        if inst.hasvels:
            try:
                scale = ncfile.variables['velocities'].scale_factor
            except AttributeError:
                scale = 1
            inst._velocities = np.array(ncfile.variables['velocities'][:])*scale
            inst.velocity_scale = scale
        if inst.hasfrcs:
            inst._forces = np.array(ncfile.variables['forces'][:])
        if inst.frame is None:
            if 'time' in ncfile.variables:
                inst.frame = len(ncfile.variables['time'][:])
            elif inst.hascrds:
                inst.frame = inst._coordinates.shape[0]
            elif inst.hasvels:
                inst.frame = inst._velocities.shape[0]
            elif inst.hasfrcs:
                inst.frame = inst._forces.shape[0]
        return inst

    @property
    def coordinates(self):
        return self._coordinates

    def add_coordinates(self, stuff):
        """
        Adds a new coordinate frame to the end of a NetCDF trajectory. This
        should only be called on objects created with the "open_new"
        constructor.

        Parameters
        ----------
        stuff : iterable of floats or distance Quantity
            This array of floats is converted into a numpy array of shape
            (natom, 3). It can be passed either in the 2-D format of
            [ [x1, y1, z1], [x2, y2, z2], ... ] or in the 1-D format of
            [x1, y1, z1, x2, y2, z2, ... ].
        """
        if u.is_quantity(stuff):
            stuff = stuff.value_in_unit(u.angstroms)
        stuff = np.asarray(stuff, dtype='f')
        self._ncfile.variables['coordinates'][self._last_crd_frame] = \
                np.reshape(stuff, (self.atom, 3))
        self._last_crd_frame += 1
        self.flush()

    @property
    def velocities(self):
        return self._velocities

    def add_velocities(self, stuff):
        """
        Adds a new velocities frame to the end of a NetCDF trajectory. This
        should only be called on objects created with the "open_new"
        constructor.

        Parameters
        ----------
        stuff : iterable of floats or distance/time Quantity
            This array of floats is converted into a numpy array of shape
            (natom, 3). It can be passed either in the 2-D format of
            [ [x1, y1, z1], [x2, y2, z2], ... ] or in the 1-D format of
            [x1, y1, z1, x2, y2, z2, ... ].
        """
        if u.is_quantity(stuff):
            stuff = stuff.value_in_unit(u.angstrom/u.picosecond)
        stuff = np.asarray(stuff)
        self._ncfile.variables['velocities'][self._last_vel_frame] = \
                np.reshape(stuff, (self.atom, 3)) / self.velocity_scale
        self._last_vel_frame += 1
        self.flush()

    @property
    def forces(self):
        return self._forces

    def add_forces(self, stuff):
        """
        Adds a new coordinate frame to the end of a NetCDF trajectory. This
        should only be called on objects created with the "open_new"
        constructor.

        Parameters
        ----------
        stuff : iterable of floats or energy/distance Quantity
            This array of floats is converted into a numpy array of shape
            (natom, 3). It can be passed either in the 2-D format of
            [ [x1, y1, z1], [x2, y2, z2], ... ] or in the 1-D format of
            [x1, y1, z1, x2, y2, z2, ... ].
        """
        if u.is_quantity(stuff):
            stuff.value_in_unit(u.kilocalories_per_mole/u.angstroms)
        self._ncfile.variables['forces'][self._last_frc_frame] = \
                np.reshape(stuff, (self.atom, 3))
        self._last_frc_frame += 1
        self.flush()

    @property
    def cell_lengths_angles(self):
        try:
            return np.hstack(
                (
                    self._ncfile.variables['cell_lengths'][:],
                    self._ncfile.variables['cell_angles'][:],
                )
            )
        except KeyError:
            return None

    box = cell_lengths_angles

    def add_cell_lengths_angles(self, lengths, angles=None):
        """
        Adds a new cell length and angle frame to the end of a NetCDF
        trajectory.  This should only be called on objects created with the
        "open_new" constructor.

        Parameters
        ----------
        lengths : array of 3 (or 6) floats (or Quantities)
            This should be a 1-D array of 3 or 6 elements.  If 6 elements,
            `angles` should be None (below) and the first 3 elements are the box
            lengths (angstroms) and the last 3 are the box angles (degrees).
        angles : 3-item iterable = None
            These are the box angles (if lengths contains only 3 elements) in
            degrees. Must be a 1-D array of 3 elements or None if lengths
            includes angles as well.
        """
        def strip_units(x, desired_units):
            if u.is_quantity(x): return x.value_in_unit(desired_units)
            return x
        if len(lengths) == 3 and angles is None:
            raise ValueError('Both lengths and angles are required.')
        if len(lengths) == 6 and angles is not None:
            raise ValueError('Angles can be provided only once.')
        if len(lengths) != 6 and (len(lengths) != 3 or len(angles) != 3):
            raise ValueError('6 numbers expected -- 3 lengths and 3 angles.')
        if angles is None:
            angles = [strip_units(x, u.degrees) for x in lengths[3:]]
        lengths = [strip_units(x, u.angstroms) for x in lengths[:3]]
        self._ncfile.variables['cell_lengths'][self._last_box_frame] = np.asarray(lengths)
        self._ncfile.variables['cell_angles'][self._last_box_frame] = np.asarray(angles)
        self._last_box_frame += 1
        self.flush()

    add_box = add_cell_lengths_angles

    @property
    def time(self):
        return self._ncfile.variables['time'][:]

    def add_time(self, stuff):
        """ Adds the time to the current frame of the NetCDF file

        Parameters
        ----------
        stuff : float or time-dimension Quantity
            The time to add to the current frame
        """
        if u.is_quantity(stuff): stuff = stuff.value_in_unit(u.picoseconds)
        self._ncfile.variables['time'][self._last_time_frame] = float(stuff)
        self._last_time_frame += 1
        self.flush()

    @property
    def remd_indices(self):
        return self._ncfile.variables['remd_indices'][:]

    def add_remd_indices(self, stuff):
        """ Add REMD indices to the current frame of the NetCDF file

        Parameters
        ----------
        stuff : iterable of int
            The indices in each REMD dimension
        """
        self._ncfile.variables['remd_indices'][self._last_remd_frame] = \
                np.asarray(stuff, dtype='i')
        self._last_remd_frame += 1
        self.flush()

    @property
    def temp0(self):
        return self._ncfile.variables['temp0'][:]

    def add_temp0(self, stuff):
        """ The temperature to add to the current frame of the NetCDF file

        Parameters
        ----------
        stuff : float or temperature Quantity
            The temperature to add to the current NetCDF file
        """
        if u.is_quantity(stuff): stuff = stuff.value_in_unit(u.kelvin)
        self._ncfile.variables['temp0'][self._last_remd_frame] = float(stuff)
        self._last_remd_frame += 1
        self.flush()

    @property
    def remd_dimtype(self):
        return self._ncfile.variables['remd_dimtype'][:]

    @remd_dimtype.setter
    def remd_dimtype(self, stuff):
        self._ncfile.variables['remd_dimtype'][:] = np.asarray(stuff, dtype='i')
        self.flush()

    def close(self):
        """ Closes the NetCDF file """
        if not self.closed:
            self._ncfile.close()
            self.closed = True

    def __del__(self):
        self.closed or (hasattr(self, '_ncfile') and self._ncfile.close())

    def flush(self):
        if nc is None:
            # netCDF4.Dataset does not have a flush method
            self._ncfile.flush()
