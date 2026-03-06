"""
Provides a class for parsing CHARMM-style coordinate files, namely CHARMM .crd
(coordinate) files and CHARMM .rst (restart) file. Uses CharmmFile class in
_charmmfile.py for reading files

Author: Jason Deckman
Contributors: Jason Swails
Date: June 19, 2015
"""
from contextlib import closing
import numpy as np
from ..utils import io
from ..formats.registry import FileFormatType
from ..exceptions import CharmmError
from .. import unit as u
from ..vec3 import Vec3

CHARMLEN = 22
TIMESCALE = 4.888821E-14 * 1e12 # AKMA time units to picoseconds
ONE_TIMESCALE = 1 / TIMESCALE

class CharmmCrdFile(metaclass=FileFormatType):
    """
    Reads and parses a CHARMM coordinate file (.crd) into its components,
    namely the coordinates, CHARMM atom types, resid, resname, etc.

    Parameters
    ----------
    fname : str
        Name of the restart file to parse

    Attributes
    ----------
    natom : int
        Number of atoms in the system
    resname : list of str
        List of all residue names in the system
    coordinates : np.ndarray with shape (1, natom, 3)
        Atomic coordinates in a numpy array
    positions : natom x 3 distance Quantity
        2-D list of all coordinates with the appropriate distance unit attached.
        Has the format [ [x1, y1, z1], [x2, y2, z2], ... ]
    """

    @staticmethod
    def id_format(filename):
        """ Identifies the file type as a CHARMM coordinate file

        Parameters
        ----------
        filename : str
            Name of the file to check format for

        Returns
        -------
        is_fmt : bool
            True if it is a CHARMM coordinate file
        """
        with closing(io.genopen(filename)) as f:
            line = f.readline()
            while line and len(line.strip()) == 0:   # Skip whitespace
                line = f.readline()

            intitle = True
            while intitle:
                line = f.readline()
                if len(line.strip()) == 0:
                    intitle = False
                elif line[0] != '*':
                    intitle = False
                else:
                    intitle = True

            while line and len(line.strip()) == 0:      # Skip whitespace
                line = f.readline()

            try:
                natom = int(line.split()[0])
                for row in range(min(natom, 3)):
                    line = f.readline().split()
                    int(line[0])
                    int(line[1])
                    float(line[4])
                    float(line[5])
                    float(line[6])
                    float(line[9])
            except (IndexError, ValueError):
                return False

            return True

    def __init__(self, fname):
        self.atomno = []                   # Atom number
        self.resno = []                    # Residue number
        self.resname = []                  # Residue name
        self.resid = []                    # Residue ID
        self.atname = []                   # Atom type
        self.coords = []                   # 3N atomic coordinates
        self.title = []                    # .crd file title block
        self.segid = []                    # Segment ID
        self.weighting = []                # Atom weighting

        self.natom = 0                     # Number of atoms specified in file
        self._parse(fname)

    @property
    def positions(self):
        """
        Atomic coordinates with units attached to them with the shape (natom, 3)
        """
        return [Vec3(*xyz) for xyz in self.coordinates[0]] * u.angstroms

    @property
    def coordinates(self):
        return self.coords

    @property
    def box(self):
        return None

    def _parse(self, fname):

        with closing(io.genopen(fname, 'r')) as crdfile:
            line = crdfile.readline().strip()

            while len(line) == 0:   # Skip whitespace, as a precaution
                line = crdfile.readline().strip()

            intitle = True
            while intitle:
                self.title.append(line)
                line = crdfile.readline().strip()
                if len(line) == 0:
                    intitle = False
                elif line[0] != '*':
                    intitle = False
                else:
                    intitle = True

            while len(line) == 0:      # Skip whitespace
                line = crdfile.readline().strip()

            try:
                self.natom = int(line.split()[0])
                for row in range(self.natom):
                    line = crdfile.readline().split()
                    self.atomno.append(int(line[0]))
                    self.resno.append(int(line[1]))
                    self.resname.append(line[2])
                    self.atname.append(line[3])
                    self.coords.append(float(line[4]))
                    self.coords.append(float(line[5]))
                    self.coords.append(float(line[6]))
                    self.segid.append(line[7])
                    self.resid.append(line[8])
                    self.weighting.append(float(line[9]))

                if len(self.coords) != 3 * self.natom:
                    raise RuntimeError(f"{len(self.coords)} coordinates are incompatible with {self.natom} atoms")
            except (ValueError, IndexError):
                raise CharmmError('Error parsing CHARMM coordinate file')
            self.coords = np.array(self.coords).reshape((-1, self.natom, 3))

    @staticmethod
    def write(struct, dest):
        """ Writes a CHARMM coordinate file from a structure

        Parameters
        ----------
        struct : :class:`parmed.structure.Structure`
            The input structure to write the CHARMM coordinate file from
        dest : str or file-like object
            The file name or file object to write the coordinate file to
        """
        if isinstance(dest, str):
            dest = io.genopen(dest, 'w')
            own_handle = True
        else:
            own_handle = False

        dest.write('* GENERATED BY PARMED (HTTPS://GITHUB.COM/PARMED/PARMED)\n')
        dest.write('*\n')
        dest.write(f"{len(struct.atoms):10d}  EXT\n")
        add = 0 if struct.residues[0].number > 0 else 1-struct.residues[0].number
        for i, atom in enumerate(struct.atoms):
            res = atom.residue
            segid = res.segid.strip() or res.chain.strip() or 'SYS'
            dest.write(
                f'{i + 1:10d}{atom.residue.number + add:10d}  {atom.residue.name:<8s}  {atom.name:<8s}'
                f'{atom.xx:20.10f}{atom.xy:20.10f}{atom.xz:20.10f}  {segid:<8s}  '
                f'{str(atom.residue.number):<8s}{0:20.10f}\n'
            )
        if own_handle:
            dest.close()

class CharmmRstFile(metaclass=FileFormatType):
    """
    Reads and parses data, velocities and coordinates from a CHARMM restart
    file (.rst) of file name 'fname' into class attributes

    Parameters
    ----------
    fname : str
        Name of the restart file to parse

    Attributes
    ----------
    natom : int
        Number of atoms in the system
    resname : list of str
        Names of all residues in the system
    coordinates : np.ndarray shape(1, natom, 3)
        List of all coordinates in the format [x1, y1, z1, x2, y2, z2, ...]
    coordinatesold : np.ndarray shape(1, natom, 3)
        List of all old coordinates in the format [x1, y1, z1, x2, y2, z2, ...]
    velocities : np.ndarray shape(1, natom, 3)
        List of all velocities in the format [x1, y1, z1, x2, y2, z2, ...]
    positions : natom x 3 distance Quantity
        2-D list of all coordinates with the appropriate distance unit attached.
        Has the format [ [x1, y1, z1], [x2, y2, z2], ... ]
    positionsold : natom x 3 distance Quantity
        2-D list of all old coordinates with the appropriate distance unit
        attached.  Has the format [ [x1, y1, z1], [x2, y2, z2], ... ]
    """

    @staticmethod
    def id_format(filename):
        """ Identifies the file type as a CHARMM restart file

        Parameters
        ----------
        filename : str
            Name of the file to check format for

        Returns
        -------
        is_fmt : bool
            True if it is a CHARMM restart file
        """
        with closing(io.genopen(filename)) as f:
            line = f.readline()
        return line.startswith('REST')

    def __init__(self, fname):
        self.header = []
        self.title = []
        self.enrgstat = []
        self.coordsold = []
        self.coords = []
        self.vels = []

        self.ff_version = 0
        self.natom = 0
        self.npriv = 0
        self.nstep = 0
        self.nsavc = 0
        self.nsavv = 0
        self.jhstrt = 0

        self._parse(fname)

    @property
    def coordinates(self):
        return self.coords

    @property
    def coordinatesold(self):
        return self.coordsold

    @property
    def positions(self):
        """ Atomic positions with units """
        return [Vec3(*xyz) for xyz in self.coords[0]] * u.angstroms

    @property
    def positionsold(self):
        """ Old atomic positions with units """
        return [Vec3(*xyz) for xyz in self.coordsold[0]] * u.angstroms

    @property
    def velocities(self):
        """ Atomic velocities in Angstroms/picoseconds """
        return self.vels

    @property
    def box(self):
        return None

    def _parse(self, fname):

        with closing(io.genopen(fname, 'r')) as crdfile:
            readingHeader = True
            while readingHeader:
                line = crdfile.readline()
                if not len(line):
                    raise CharmmError('Premature end of file')
                line = line.strip()
                words = line.split()
                if len(line) != 0:
                    if words[0] == 'ENERGIES' or words[0] == '!ENERGIES':
                        readingHeader = False
                    else:
                        self.header.append(line.strip())
                else:
                    self.header.append(line.strip())

            for row in range(len(self.header)):
                if len(self.header[row].strip()) != 0:
                    line = self.header[row].strip().split()
                    if line[0][0:5] == 'NATOM' or line[0][0:6] == '!NATOM':
                        try:
                            line = self.header[row+1].strip().split()
                            self.natom = int(line[0])
                            self.npriv = int(line[1])     # num. previous steps
                            self.nstep = int(line[2])     # num. steps in file
                            self.nsavc = int(line[3])     # coord save frequency
                            self.nsavv = int(line[4])     # velocities "
                            self.jhstrt = int(line[5])    # Num total steps?
                            break
                        except (ValueError, IndexError):
                            raise CharmmError('Problem parsing CHARMM restart')

            self.scan(crdfile, '!XOLD')
            self._get_formatted_crds(crdfile, self.coordsold)
            self.coordsold = np.array(self.coordsold).reshape((-1,self.natom,3))

            self.scan(crdfile, '!VX')
            self._get_formatted_crds(crdfile, self.vels)
            self.vels = np.array(self.vels).reshape((-1, self.natom, 3))
            # Convert velocities to angstroms/ps
            self.vels *= ONE_TIMESCALE

            self.scan(crdfile, '!X')
            self._get_formatted_crds(crdfile, self.coords)
            self.coords = np.array(self.coords).reshape((-1, self.natom, 3))

    def scan(self, handle, str, r=0): # read lines in file till 'str' is found
        scanning = True

        if r:
            handle.seek(0)

        while scanning:
            line = handle.readline()
            if not line:
                raise CharmmError('Premature end of file')

            if len(line.strip()) != 0:
                if line.strip().split()[0][0:len(str)] == str:
                    scanning = False

    def _get_formatted_crds(self, crdfile, crds):
        for row in range(self.natom):
            line = crdfile.readline()

            if not line:
                raise CharmmError('Premature end of file')

            if len(line) < 3 * CHARMLEN:
                raise CharmmError(
                    "Less than 3 coordinates present in coordinate row or coords may be truncated."
                )

            line = line.replace('D','E')     # CHARMM uses 'D' for exponentials

            # CHARMM uses fixed format (len = CHARMLEN = 22) for crds in .rst's
            crds.append(float(line[0:CHARMLEN]))
            crds.append(float(line[CHARMLEN:2*CHARMLEN]))
            crds.append(float(line[2*CHARMLEN:3*CHARMLEN]))
