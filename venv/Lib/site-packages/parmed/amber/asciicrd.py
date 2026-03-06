"""
This is a pure-python module for reading and writing ASCII Amber structure
files, like trajectories and restarts.  Because it is pure-python and
formatted-ASCII, parsing and writing these files are expected to be slow. Binary
alternatives (like DCD and NetCDF, provided in netcdffiles.py) are strongly
encouraged, but these are provided for more complete compatibility and for
instances where the prequisites may not be installed.
"""
from math import ceil
import numpy as np
from ..formats.registry import FileFormatType
from ..utils.io import genopen
from ..structure import Structure
from ..topologyobjects import Atom
from .. import unit as u
from ..vec3 import Vec3
import warnings as _warnings

VELSCALE = 20.455
ONEVELSCALE = 1 / VELSCALE

class _FileEOF(Exception):
    """ For control flow """

class _AmberAsciiCoordinateFile(metaclass=FileFormatType):
    """
    Abstract base class for interacting with ASCII coordinate files.
    Opens a new ASCII coordinate file and either parses it (loading
    everything into memory) or sets it up for writing.

    Parameters
    ----------
    fname : str
        File name to open
    natom : int
        Number of atoms in the system
    hasbox : bool
        Does the system have PBCs?
    mode : str={'r', 'w'}
        Whether to open this file for 'r'eading or 'w'riting
    title : str, optional
        Title to write to a new trajectory (when mode='w')
    """

    DEFAULT_TITLE = None
    CRDS_PER_LINE = None

    def __init__(self, fname, natom, hasbox, mode='r', title=None):
        if mode == 'r':
            self._status = 'old'
        elif mode == 'w':
            self._status = 'new'
            # We need to have some way to know whether we need to write the
            # coordinates or the box for this particular frame. Each frame must
            # be written as coordinates first, then box.
            self._writebox = False
        else:
            raise ValueError(f"{self.__class__.__name__} mode must be 'r' or 'w'")
        if isinstance(fname, str):
            self._file = genopen(fname, mode)
            self._own_handle = True
        elif hasattr(fname, 'read'):
            self._file = fname
            self._own_handle = False
        else:
            raise TypeError(f"Unsupported type for {fname}")

        self.natom = natom
        self.hasbox = hasbox
        if self.CRDS_PER_LINE is None:
            raise NotImplementedError("CRDS_PER_LINE must be set by the subclass")
        self._full_lines_per_frame = self.natom * 3 // self.CRDS_PER_LINE
        self._nextras = self.natom * 3 - (self._full_lines_per_frame * self.CRDS_PER_LINE)
        if self._own_handle:
            self.closed = False
        if self._status == 'old':
            self._parse()
        elif self._status == 'new':
            if title is None:
                self._file.write('%s\n' % self.DEFAULT_TITLE)
            else:
                self._file.write(title.rstrip() + '\n')

    def _parse(self):
        """ Handles actual file parsing """
        raise NotImplementedError('virtual method not overwritten')

    @property
    def positions(self):
        return [Vec3(*xyz) for xyz in self._coordinates[0]] * u.angstroms

    def close(self):
        """ Close the open file handler """
        if self._own_handle:
            self.closed or self._file.close()
            self.closed = True

    def __del__(self):
        """ Make sure the open file handler is closed """
        try:
            self.closed or self._file.close()
        except AttributeError:
            pass

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AmberAsciiRestart(_AmberAsciiCoordinateFile):
    """
    Parser for the Amber ASCII inpcrd/restart file format

    Parameters
    ----------
        fname : str
            File name to open
        mode : str={'r', 'w'}
            Whether to open this file for 'r'eading or 'w'riting
        natom : int, optional
            Number of atoms in the system (necessary when mode='w')
        hasbox : bool, optional
            Does the system have PBCs? Necessary when mode='w'
        title : str, optional
            Title to write to a new trajectory (when mode='w')
        time : float, optional
            The time to write to the restart file in ps. Default is 0.
    """
    @staticmethod
    def id_format(filename):
        """ Identifies the file type as an Amber restart/inpcrd file

        Parameters
        ----------
        filename : str
            Name of the file to check format for

        Returns
        -------
        is_fmt : bool
            True if it is an Amber restart/inpcrd file. False otherwise
        """
        if isinstance(filename, str):
            f = genopen(filename, 'r')
            lines = [f.readline() for i in range(5)]
            f.close()
        elif all(hasattr(filename, attr) for attr in ['readline', 'seek', 'tell']):
            cur = filename.tell()
            lines = [filename.readline() for i in range(5)]
            filename.seek(cur)
        # Look for natom
        words = lines[1].split()
        if len(words) > 2 or len(words) < 1:
            return False
        try:
            natom = int(words[0])
            float(words[1])
        except ValueError:
            return False
        except IndexError:
            pass
        # Next 3 lines, make sure we have %12.7f format. This only works if we
        # have at least 6 atoms. Any fewer than that means the restart file is
        # shorter than that.
        try:
            if natom <= 0:
                return False
            i = 0
            for line in lines[2:]:
                i += 1
                if i > natom: break
                for j in range(3):
                    j12 = j * 12
                    if line[j12+4] != '.': return False
                    float(line[j12:j12+12])
                    if line[j12+11] not in '0123456789':
                        return False
                i += 1
                if i > natom: break
                for j in range(3):
                    j12 = j * 12 + 36
                    if line[j12+4] != '.': return False
                    float(line[j12:j12+12])
                    if line[j12+11] not in '0123456789':
                        return False
        except (IndexError, ValueError):
            return False

        # Must be a restart...
        return True

    CRDS_PER_LINE = 6
    DEFAULT_TITLE = 'restart created by ParmEd'

    def __init__(self, fname, mode='r', natom=0, hasbox=None, title=None, time=0.0):
        """
        For restart files, natom and hasbox are determined automatically for
        mode='r', and can be determined at write-time when the coordinates are
        set.
        """
        self._coords_written = False
        self._cell_lengths_written = False
        self._cell_angles_written = False
        self._vels_written = False
        self.time = float(time)
        super().__init__(fname, natom, hasbox, mode, title)

    @classmethod
    def parse(cls, filename, structure=False):
        self = cls(filename)
        if structure:
            obj = Structure()
            for _ in range(self.natom):
                # fake
                obj.add_atom(Atom(), resname='XXX', resnum=0)
            obj.box = self.box
            obj.coordinates = self.coordinates[0]
            obj._coordinates = self.coordinates
            return obj
        else:
            return self

    def _parse(self):
        """
        This method parses the data out of the ASCII restart file and creates
        self._coordinates and self._velocities as np.ndarray(natom*3) arrays

        This method is called automatically for 'old' restart files and should
        not be called by external callers
        """

        lines = self._file.readlines()
        self._file.close()
        self.title = lines[0].strip()
        self.natom = int(lines[1].strip().split()[0])
        try:
            self.time = float(lines[1].strip().split()[1])
        except IndexError:
            self.time = 0.0
        # Get rid of any trailing newlines
        while not lines[-1].strip():
            lines.pop()
        # Determine what information we have based on the number of lines
        # present
        if len(lines) == int(ceil(self.natom / 2.0) + 2):
            self.hasbox = self.hasvels = False
        elif len(lines) == int(ceil(self.natom / 2.0) + 3):
            self.hasbox = True
            self.hasvels = False
        elif len(lines) == int(2 * ceil(self.natom / 2.0) + 2):
            self.hasbox = False
            self.hasvels = True
        elif len(lines) == int(2 * ceil(self.natom / 2.0) + 3):
            self.hasbox = self.hasvels = True
        else:
            raise RuntimeError(f'Badly formatted restart file. Has {len(lines)} lines for {self.natom} atoms.')
        self._coordinates = np.zeros((self.natom, 3))
        if self.hasvels:
            self._velocities = np.zeros((self.natom, 3))
        if self.hasbox:
            self._cell_lengths = np.zeros(3)
            self._cell_angles = np.zeros(3)
        # Now it's time to parse. Coordinates first
        startline = 2
        endline = startline + int(ceil(self.natom / 2.0))
        idx = 0
        for i in range(startline, endline):
            line = lines[i]
            x1 = float(line[ 0:12])
            y1 = float(line[12:24])
            z1 = float(line[24:36])
            self._coordinates[idx] = [x1, y1, z1]
            idx += 1
            try:
                x2 = float(line[36:48])
                y2 = float(line[48:60])
                z2 = float(line[60:72])
            except ValueError:
                pass
            else:
                self._coordinates[idx] = [x2, y2, z2]
                idx += 1
        self._coordinates = self._coordinates.reshape((1, self.natom, 3))
        startline = endline
        # Now it's time to parse the velocities if we have them
        if self.hasvels:
            endline = startline + int(ceil(self.natom / 2.0))
            idx = 0
            for i in range(startline, endline):
                line = lines[i]
                x1 = float(line[ 0:12]) * VELSCALE
                y1 = float(line[12:24]) * VELSCALE
                z1 = float(line[24:36]) * VELSCALE
                self._velocities[idx] = [x1, y1, z1]
                idx += 1
                try:
                    x2 = float(line[36:48]) * VELSCALE
                    y2 = float(line[48:60]) * VELSCALE
                    z2 = float(line[60:72]) * VELSCALE
                except ValueError:
                    pass
                else:
                    self._velocities[idx] = [x2,y2,z2]
                    idx += 1
            startline = endline
            self._velocities = self._velocities.reshape((1, self.natom, 3))
        # Now it's time to parse the box info if we have it
        if self.hasbox:
            line = lines[startline]
            self._cell_lengths[0:3] = [float(line[0:12]), float(line[12:24]), float(line[24:36])]
            self._cell_angles[0:3] = [float(line[36:48]), float(line[48:60]), float(line[60:72])]

    @property
    def coordinates(self):
        if self._status == 'new' and not hasattr(self, '_coordinates'):
            raise RuntimeError('Coordinates not yet set')
        return self._coordinates

    @coordinates.setter
    def coordinates(self, stuff):
        if self._status == 'old':
            raise RuntimeError('Cannot set coordinates on an old restart')
        stuff = np.asarray(stuff).ravel()
        if self.natom > 0 and len(stuff) != 3 * self.natom:
            raise ValueError(f'Got {len(stuff)} coordinates for {self.natom} atoms')
        if self._coords_written:
            raise RuntimeError('Coordinates have already been written.')
        # Error checking done. If we didn't already set our number of atoms,
        # set that now
        self.natom = len(stuff) // 3
        self._coordinates = stuff.reshape((-1, self.natom, 3))
        self._file.write('%5d%15.7e\n' % (self.natom, self.time))
        numwrit = 0
        fmt = '%12.7f%12.7f%12.7f'
        for i in range(self.natom):
            i3 = i * 3
            self._file.write(fmt % (stuff[i3], stuff[i3+1], stuff[i3+2]))
            numwrit += 1
            if numwrit % 2 == 0:
                self._file.write('\n')
        if self.natom % 2 == 1:
            self._file.write('\n')
        self._coords_written = True

    @property
    def velocities(self):
        if self._status == 'new' and not hasattr(self, '_velocities'):
            raise RuntimeError('Velocities not set yet')
        if not self.hasvels:
            return None
        return self._velocities

    @velocities.setter
    def velocities(self, stuff):
        if self._status == 'old':
            raise RuntimeError('Cannot set velocities on an old restart')
        stuff = np.asarray(stuff).ravel()
        if not self._coords_written:
            raise RuntimeError('Coordinates must be set before velocities')
        if self._cell_lengths_written or self._cell_angles_written:
            raise RuntimeError('Velocities must be written before the box info')
        if self._vels_written:
            raise RuntimeError('Can only write velocities once')
        if len(stuff) != 3 * self.natom:
            raise ValueError(f'Got {len(stuff)} velocities for {self.natom} atoms.')
        self._velocities = stuff.reshape((-1, self.natom, 3))
        fmt = '%12.7f%12.7f%12.7f'
        numwrit = 0
        for i in range(self.natom):
            i3 = i * 3
            self._file.write(fmt % tuple(stuff[i3+i]*ONEVELSCALE for i in range(3)))
            numwrit += 1
            if numwrit % 2 == 0:
                self._file.write('\n')
        if self.natom % 2 == 1:
            self._file.write('\n')
        self._vels_written = self.hasvels = True

    @property
    def cell_lengths(self):
        if self._status == 'new' and not hasattr(self, '_cell_lengths'):
            raise RuntimeError('Cell lengths not yet available')
        if not self.hasbox:
            return None
        return self._cell_lengths

    @property
    def cell_angles(self):
        if self._status == 'new' and not hasattr(self, '_cell_angles'):
            raise RuntimeError('Cell angles not yet available')
        if not self.hasbox:
            return None
        return self._cell_angles

    @property
    def box(self):
        """ Combined cell lengths and cell angles """
        if self._status == 'new' and not (hasattr(self, '_cell_lengths') and
                hasattr(self, '_cell_angles')):
            raise RuntimeError('Cell parameters not yet set')
        if not self.hasbox:
            return None
        box = np.zeros(6)
        lengths, angles = self.cell_lengths, self.cell_angles
        box[0:3] = [lengths[0], lengths[1], lengths[2]]
        box[3:6] = [angles[0], angles[1], angles[2]]
        return box

    @cell_lengths.setter
    def cell_lengths(self, stuff):
        if self._status == 'old':
            raise RuntimeError('Cannot set cell lengths on old restart')
        if not self._coords_written:
            raise RuntimeError('Coordinates must be written before box')
        if self._cell_lengths_written:
            raise RuntimeError('Can only write cell lengths once')
        if len(stuff) != 3:
            raise ValueError('Expected 3 numbers for cell lengths')
        self._cell_lengths = np.asarray(stuff)
        self._file.write('%12.7f%12.7f%12.7f' % (stuff[0], stuff[1], stuff[2]))
        self._cell_lengths_written = True

    @cell_angles.setter
    def cell_angles(self, stuff):
        if self._status == 'old':
            raise RuntimeError('Cannot set cell angles on old restart')
        if not self._coords_written:
            raise RuntimeError('Coordinates must be written before box')
        if not self._cell_lengths_written:
            raise RuntimeError('Must write cell lengths before angles')
        if self._cell_angles_written:
            raise RuntimeError('Can only write cell angles once')
        if len(stuff) != 3:
            raise ValueError('Expected 3 numbers for cell angles')
        self._cell_angles = np.asarray(stuff)
        self._file.write('%12.7f%12.7f%12.7f\n' % (stuff[0],stuff[1],stuff[2]))
        self._cell_angles_written = True

    @box.setter
    def box(self, stuff):
        """ Writes both the cell lengths and cell angles """
        self.cell_lengths = stuff[:3]
        self.cell_angles = stuff[3:]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AmberMdcrd(_AmberAsciiCoordinateFile):
    """
    A class to parse Amber ASCII trajectory files. This is *much* slower than
    parsing NetCDF files (or the equivalent parsing done in a compiled language
    like C or C++). For large trajectories, this may be significant.
    """
    extra_args = ('natom', 'hasbox')

    CRDS_PER_LINE = 10
    DEFAULT_TITLE = 'trajectory created by ParmEd'

    @staticmethod
    def id_format(filename):
        """ Identifies the file type as an Amber mdcrd file

        Parameters
        ----------
        filename : str
            Name of the file to check format for

        Returns
        -------
        is_fmt : bool
            True if it is an Amber mdcrd file. False otherwise
        """
        f = genopen(filename, 'r')
        lines = [f.readline() for i in range(5)]
        f.close()
        # Next 4 lines, make sure we have %8.3f format
        try:
            for i in range(4):
                i += 1
                for j in range(10):
                    j8 = j * 8
                    if lines[i][j8+4] != '.': return False
                    float(lines[i][j8:j8+8])
                    if lines[i][j8+7] not in '0123456789':
                        return False
        except (IndexError, ValueError):
            return False

        # Must be a mdcrd
        return True

    def _parse(self):
        """
        This method parses the data out of the mdcrd file and creates
        self._coordinates as a list of np.ndarray(1, natom, 3) and
        self.cell_lengths as a list of np.ndarray(3) for each frame in the
        trajectory. This method is called automatically for 'old' trajectory
        files and should not be called by external callers.
        """
        line = self._file.readline()
        self.title = line.strip()
        self._coordinates = np.ndarray(0)
        self.cell_lengths = np.ndarray(0)
        mainiter = range(0, 8*self.CRDS_PER_LINE, 8)
        extraiter = range(0, 8*self._nextras, 8)
        try:
            while line:
                frame = np.zeros(self.natom*3)
                if not line: raise _FileEOF()
                cell = np.zeros(3)
                idx = 0
                line = self._file.readline()
                if not line: raise StopIteration()
                for i in range(self._full_lines_per_frame):
                    if not line: raise _FileEOF()
                    frame[idx:idx+10] = [float(line[j:j+8]) for j in mainiter]
                    idx += 10
                    line = self._file.readline()

                if self._nextras:
                    frame[idx:idx+self._nextras] = [float(line[j:j+8]) for j in extraiter]

                if self.hasbox:
                    line = self._file.readline()
                    if not line: raise _FileEOF()
                    cell[0] = float(line[:8])
                    cell[1] = float(line[8:16])
                    cell[2] = float(line[16:24])

                self._coordinates = np.concatenate((self._coordinates, frame))
                self.cell_lengths = np.concatenate((self.cell_lengths, cell))

        except _FileEOF:
            _warnings.warn('Unexpected EOF in parsing mdcrd. natom and/or hasbox are likely wrong',
                           RuntimeWarning)
        except StopIteration:
            pass

        self._coordinates = self._coordinates.reshape((-1, self.natom, 3))
        self.cell_lengths = self.cell_lengths.reshape((-1, 3))
        self._file.close()

    @property
    def coordinates(self):
        return self._coordinates

    @property
    def frame(self):
        return self._coordinates.shape[0]

    @property
    def box(self):
        if not self._status == 'old':
            raise RuntimeError('Cannot access box of a new mdcrd')
        return self.cell_lengths

    def add_coordinates(self, stuff):
        """
        Prints 'stuff' (which must be either an iterable of 3*natom or have an
        attribute 'flatten' that converts it into an iterable of 3*natom) to the
        open file handler. Can only be called on a 'new' mdcrd, and adds these
        coordinates to the current end of the file.

        Parameters
        ----------
        stuff : array or iterable
            This must be an iterable of length 3*natom or a numpy array that can
            be flattened to a 3*natom-length array

        Raises
        ------
        If the coordinate file is an old one being parsed or if you are
        currently expected to provide unit cell dimensions, a RuntimeError is
        raised. If the provided coordinate data does not have length 3*natom, or
        cannot be ``flatten()``ed to create a 3*natom array, a ValueError is
        raised.
        """
        # Make sure we can write the coordinates right now
        if not self._status == 'new':
            raise RuntimeError('Cannot print frames to an old mdcrd')
        try:
            stuff = stuff.flatten()
        except AttributeError:
            pass
        if self._writebox:
            raise RuntimeError('Box information not written for last frame')
        if len(stuff) != 3 * self.natom:
            raise ValueError('add_coordinates requires an array of length natom*3')

        # If we can, write the coordinates
        i = 0
        j = -1
        for i in range(self._full_lines_per_frame):
            i10 = i * 10
            for j in range(10):
                self._file.write('%8.3f' % stuff[i10+j])
            self._file.write('\n')
        if self._nextras:
            extra = i*10+j + 1
            while extra < self.natom*3:
                self._file.write('%8.3f' % stuff[extra])
                extra += 1
            self._file.write('\n')
        # Now it's time to write the box info if necessary
        self._writebox = self.hasbox
        self._file.flush()

    def add_box(self, stuff):
        """
        Prints 'stuff' (which must be a 3-element list, array.array, tuple, or
        np.ndarray) as the box lengths for this frame

        Parameters
        ----------
        stuff : array or iterable
            This must be an iterable of length 3 with the box lengths

        Raises
        ------
        If the coordinate file is an old one being parsed or if you are
        currently expected to provide coordinates, a RuntimeError is raised.
        raised. If the provided box lengths are not length 3, a ValueError is
        raised.
        """
        # First make sure we should be writing our box now
        if not self._status == 'new':
            raise RuntimeError('Cannot print box to an old mdcrd')
        if not self._writebox:
            raise RuntimeError('Should not be writing box info right now')
        if len(stuff) != 3:
            raise ValueError('add_box requires an array of length 3')

        self._file.write('%8.3f%8.3f%8.3f\n' % (stuff[0], stuff[1], stuff[2]))
        self._writebox = False
        self._file.flush()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
