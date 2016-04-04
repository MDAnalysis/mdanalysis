# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""
PDBQT structure files in MDAnalysis --- :mod:`MDAnalysis.coordinates.PDBQT`
===========================================================================

MDAnalysis reads coordinates from PDBQT_ files and additional optional
data such as B-factors, partial charge and AutoDock_ atom types.  It
is also possible to substitute a PDBQT file for a PSF file in order to
define the list of atoms (but no connectivity information will be
available in this case).

.. _PDBQT:
   http://autodock.scripps.edu/faqs-help/faq/what-is-the-format-of-a-pdbqt-file
.. _AutoDock:
   http://autodock.scripps.edu/
"""

import os
import errno
import numpy as np
import warnings

from ..lib import util
from . import base


class PDBQTReader(base.SingleFrameReader):
    """PDBQTReader that reads a PDBQT-formatted file, no frills.

    Records read:
     - CRYST1 for unitcell A,B,C, alpha,beta,gamm
     - ATOM. HETATM for x,y,z

    Original `PDB format documentation`_    with  `AutoDOCK extensions`_


    .. _PDB format documentation:
       http://www.wwpdb.org/documentation/format32/sect9.html
    .. _AutoDOCK extensions:
       http://autodock.scripps.edu/faqs-help/faq/what-is-the-format-of-a-pdbqt-file

    =============  ============  ===========  =============================================
    COLUMNS        DATA  TYPE    FIELD        DEFINITION
    =============  ============  ===========  =============================================
    1 -  6         Record name   "CRYST1"
    7 - 15         Real(9.3)     a              a (Angstroms).
    16 - 24        Real(9.3)     b              b (Angstroms).
    25 - 33        Real(9.3)     c              c (Angstroms).
    34 - 40        Real(7.2)     alpha          alpha (degrees).
    41 - 47        Real(7.2)     beta           beta (degrees).
    48 - 54        Real(7.2)     gamma          gamma (degrees).

    1 -  6         Record name   "ATOM  "
    7 - 11         Integer       serial       Atom  serial number.
    13 - 16        Atom          name         Atom name.
    17             Character     altLoc       Alternate location indicator. IGNORED
    18 - 21        Residue name  resName      Residue name.
    22             Character     chainID      Chain identifier.
    23 - 26        Integer       resSeq       Residue sequence number.
    27             AChar         iCode        Code for insertion of residues. IGNORED
    31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    55 - 60        Real(6.2)     occupancy    Occupancy.
    61 - 66        Real(6.2)     tempFactor   Temperature  factor.
    67 - 76        Real(10.4)    partialChrg  Gasteiger PEOE partial charge *q*.
    79 - 80        LString(2)    atomType     AutoDOCK atom type *t*.
    =============  ============  ===========  =============================================

    We ignore torsion notation and just pull the partial charge and atom type columns::

         COMPND    NSC7810
         REMARK  3 active torsions:
         REMARK  status: ('A' for Active; 'I' for Inactive)
         REMARK    1  A    between atoms: A7_7  and  C22_23
         REMARK    2  A    between atoms: A9_9  and  A11_11
         REMARK    3  A    between atoms: A17_17  and  C21_21
         ROOT
         123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789. (column reference)
         ATOM      1  A1  INH I           1.054   3.021   1.101  0.00  0.00     0.002 A
         ATOM      2  A2  INH I           1.150   1.704   0.764  0.00  0.00     0.012 A
         ATOM      3  A3  INH I          -0.006   0.975   0.431  0.00  0.00    -0.024 A
         ATOM      4  A4  INH I           0.070  -0.385   0.081  0.00  0.00     0.012 A
         ATOM      5  A5  INH I          -1.062  -1.073  -0.238  0.00  0.00     0.002 A
         ATOM      6  A6  INH I          -2.306  -0.456  -0.226  0.00  0.00     0.019 A
         ATOM      7  A7  INH I          -2.426   0.885   0.114  0.00  0.00     0.052 A
         ATOM      8  A8  INH I          -1.265   1.621   0.449  0.00  0.00     0.002 A
         ATOM      9  A9  INH I          -1.339   2.986   0.801  0.00  0.00    -0.013 A
         ATOM     10  A10 INH I          -0.176   3.667   1.128  0.00  0.00     0.013 A
         ENDROOT
         BRANCH   9  11
         ATOM     11  A11 INH I          -2.644   3.682   0.827  0.00  0.00    -0.013 A
         ATOM     12  A16 INH I          -3.007   4.557  -0.220  0.00  0.00     0.002 A
         ATOM     13  A12 INH I          -3.522   3.485   1.882  0.00  0.00     0.013 A
         ATOM     14  A15 INH I          -4.262   5.209  -0.177  0.00  0.00    -0.024 A
         ATOM     15  A17 INH I          -2.144   4.784  -1.319  0.00  0.00     0.052 A
         ATOM     16  A14 INH I          -5.122   4.981   0.910  0.00  0.00     0.012 A
         ATOM     17  A20 INH I          -4.627   6.077  -1.222  0.00  0.00     0.012 A
         ATOM     18  A13 INH I          -4.749   4.135   1.912  0.00  0.00     0.002 A
         ATOM     19  A19 INH I          -3.777   6.285  -2.267  0.00  0.00     0.002 A
         ATOM     20  A18 INH I          -2.543   5.650  -2.328  0.00  0.00     0.019 A
         BRANCH  15  21
         ATOM     21  C21 INH I          -0.834   4.113  -1.388  0.00  0.00     0.210 C
         ATOM     22  O1  INH I          -0.774   2.915  -1.581  0.00  0.00    -0.644 OA
         ATOM     23  O3  INH I           0.298   4.828  -1.237  0.00  0.00    -0.644 OA
         ENDBRANCH  15  21
         ENDBRANCH   9  11
         BRANCH   7  24
         ATOM     24  C22 INH I          -3.749   1.535   0.125  0.00  0.00     0.210 C
         ATOM     25  O2  INH I          -4.019   2.378  -0.708  0.00  0.00    -0.644 OA
         ATOM     26  O4  INH I          -4.659   1.196   1.059  0.00  0.00    -0.644 OA
         ENDBRANCH   7  24
         TORSDOF 3
         123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789. (column reference)

    .. versionchanged:: 0.11.0
       Frames now 0-based instead of 1-based
    """
    format = 'PDBQT'
    units = {'time': None, 'length': 'Angstrom'}

    def _read_first_frame(self):
        coords = []
        unitcell = np.zeros(6, dtype=np.float32)
        with util.openany(self.filename, 'r') as pdbfile:
            for line in pdbfile:
                # Should only break at the 'END' of a model definition
                # and prevent premature exit for a torsion termination
                # , eg, ENDBRANCH
                if line.startswith('END\n'):
                    break
                if line.startswith('CRYST1'):
                    # lengths
                    x, y, z = map(
                        float, (line[6:15], line[15:24], line[24:33]))
                    # angles
                    A, B, G = map(
                        float, (line[33:40], line[40:47], line[47:54]))
                    unitcell[:] = x, y, z, A, B, G
                if line.startswith(('ATOM', 'HETATM')):
                    x, y, z = map(
                        float, (line[30:38], line[38:46], line[46:54]))
                    coords.append((x, y, z))
        self.n_atoms = len(coords)
        self.ts = self._Timestep.from_coordinates(
            np.array(coords, dtype=np.float32),
            **self._ts_kwargs)
        self.ts._unitcell[:] = unitcell
        self.ts.frame = 0  # 0-based frame number
        if self.convert_units:
            # in-place !
            self.convert_pos_from_native(self.ts._pos)
            self.convert_pos_from_native(self.ts._unitcell[:3])

    def Writer(self, filename, **kwargs):
        """Returns a permissive (simple) PDBQTWriter for *filename*.

        :Arguments:
          *filename*
              filename of the output PDBQT file

        :Returns: :class:`PDBQTWriter`

        """
        return PDBQTWriter(filename, **kwargs)


class PDBQTWriter(base.Writer):
    """PDBQT writer that implements a subset of the PDB_ 3.2 standard and the PDBQT_ spec.

    .. _PDB: http://www.wwpdb.org/documentation/format32/v3.2.html
    .. _PDBQT: http://autodock.scripps.edu/faqs-help/faq/what-is-the-format-of-a-pdbqt-file
    """
    #          1         2         3         4         5         6         7         8
    # 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.
    # ATOM__seria nameAres CressI   xxxxxxxxyyyyyyyyzzzzzzzzOCCUPAtempft          elCH
    # ATOM  %5d   %-4s %-3s %4d %1s %8.3f   %8.3f   %8.3f   %6.2f %6.2f           %2s
    #                 %1s  %1s                                                      %2d
    #            =        =      ===                                    ==========
    # ATOM  %5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2d
    # ATOM  %(serial)5d %(name)-4s%(altLoc)1s%(resName)-3s %(chainID)1s%(resSeq)4d%(iCode)1s   %(x)8.3f%(y)8.3f%(
    # z)8.3f%(occupancy)6.2f%(tempFactor)6.2f          %(element)2s%(charge)2d

    # Strict PDB format:
    #fmt = {'ATOM':   "ATOM  %(serial)5d %(name)-4s%(altLoc)1s%(resName)-3s %(chainID)1s%(resSeq)4d%(iCode)1s   %(
    # x)8.3f%(y)8.3f%(z)8.3f%(occupancy)6.2f%(tempFactor)6.2f          %(element)2s%(charge)2d\n",
    # PDB format as used by NAMD/CHARMM: 4-letter resnames and segID, altLoc ignored
    fmt = {
        #'ATOM':   "ATOM  %(serial)5d %(name)-4s %(resName)-4s%(chainID)1s%(resSeq)4d%(iCode)1s   %(x)8.3f%(y)8.3f%(
        # z)8.3f%(occupancy)6.2f%(tempFactor)6.2f      %(segID)-4s%(element)2s%(charge)2d\n",
        'ATOM': "ATOM  %(serial)5d %(name)-4s %(resName)-4s%(chainID)1s%(resSeq)4d%(iCode)1s   %(x)8.3f%(y)8.3f%("
                "z)8.3f%(occupancy)6.2f%(tempFactor)6.2f      %(partialCharge)-1.4f %(atomtype)-2s\n",
        'REMARK': "REMARK     %s\n",
        'TITLE': "TITLE    %s\n",
        'CRYST1': "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
    }
    format = 'PDBQT'
    units = {'time': None, 'length': 'Angstrom'}
    pdb_coor_limits = {"min": -999.9995, "max": 9999.9995}

    def __init__(self, filename, **kwargs):
        self.filename = util.filename(filename, ext='pdbqt')
        self.pdb = util.anyopen(self.filename, 'w')

    def close(self):
        self.pdb.close()

    def write(self, selection, frame=None):
        """Write selection at current trajectory frame to file.

        write(selection,frame=FRAME)

        :Arguments:
          *selection*
            a :class:`~MDAnalysis.core.groups.AtomGroup`
          *frame*
            optionally move to frame *FRAME*

        .. Note::

           The first letter of the
           :attr:`~MDAnalysis.core.groups.Atom.segid` is used as the PDB
           chainID.

        .. versionchanged:: 0.11.0
           Frames now 0-based instead of 1-based
        """
        u = selection.universe
        if frame is not None:
            u.trajectory[frame]  # advance to frame
        else:
            try:
                frame = u.trajectory.ts.frame
            except AttributeError:
                frame = 0  # should catch cases when we are analyzing a single PDB (?)

        self.TITLE("FRAME " + str(frame) + " FROM " + str(u.trajectory.filename))
        self.CRYST1(self.convert_dimensions_to_unitcell(u.trajectory.ts))
        atoms = selection.atoms  # make sure to use atoms (Issue 46)
        coor = atoms.coordinates()  # can write from selection == Universe (Issue 49)

        # check if any coordinates are illegal (coordinates are already in Angstroem per package default)
        if not self.has_valid_coordinates(self.pdb_coor_limits, coor):
            self.close()
            try:
                os.remove(self.filename)
            except OSError as err:
                if err.errno == errno.ENOENT:
                    pass
            raise ValueError(
                "PDB files must have coordinate values between {0:.3f} and {1:.3f} Angstroem: No file was written.".format(self.pdb_coor_limits["min"], self.pdb_coor_limits["max"]))

        for i, atom in enumerate(atoms):
            self.ATOM(serial=i + 1, name=atom.name.strip(), resName=atom.resname.strip(), resSeq=atom.resid,
                      chainID=atom.segid.strip(), partialCharge=atom.charge,
                      x=coor[i, 0], y=coor[i, 1], z=coor[i, 2], atomtype=atom.type)
            # get bfactor, too, and add to output?
            # 'element' is auto-guessed from atom.name in ATOM()
        self.close()

    def TITLE(self, *title):
        """Write TITLE record.
        http://www.wwpdb.org/documentation/format32/sect2.html
        """
        line = " ".join(title)  # should do continuation automatically
        self.pdb.write(self.fmt['TITLE'] % line)

    def REMARK(self, *remark):
        """Write generic REMARK record (without number).
        http://www.wwpdb.org/documentation/format32/remarks1.html
        http://www.wwpdb.org/documentation/format32/remarks2.html
        """
        line = " ".join(remark)
        self.pdb.write(self.fmt['REMARK'] % line)

    def CRYST1(self, dimensions, spacegroup='P 1', zvalue=1):
        """Write CRYST1 record.
        http://www.wwpdb.org/documentation/format32/sect8.html
        """
        self.pdb.write(self.fmt['CRYST1'] % (tuple(dimensions) + (spacegroup, zvalue)))

    def ATOM(self, serial=None, name=None, altLoc=None, resName=None, chainID=None,
             resSeq=None, iCode=None, x=None, y=None, z=None, occupancy=1.0, tempFactor=0.0,
             segID=None, partialCharge=None, atomtype=None):
        """Write ATOM record.
        http://www.wwpdb.org/documentation/format32/sect9.html
        Only some keyword args are optional (altLoc, iCode, chainID), for some defaults are set.

        All inputs are cut to the maximum allowed length. For integer
        numbers the highest-value digits are chopped (so that the
        serial and reSeq wrap); for strings the trailing characters
        are chopped.

        Note: Floats are not checked and can potentially screw up the format.
        """
        for arg in (
                'serial', 'name', 'resName', 'resSeq', 'x', 'y', 'z',
                'occupancy', 'tempFactor', 'partialCharge', 'atomtype'):
            if locals()[arg] is None:
                raise ValueError('parameter ' + arg + ' must be defined for PDBQT.')
        serial = int(str(serial)[-5:])  # check for overflow here?
        name = name[:4]
        if len(name) < 4:
            name = " " + name  # customary to start in column 14
        altLoc = altLoc or " "
        altLoc = altLoc[:1]
        resName = resName[:4]
        chainID = chainID or ""  # or should we provide a chainID such as 'A'?
        chainID = chainID.strip()[-1:]  # take the last character
        resSeq = int(str(resSeq)[-4:])  # check for overflow here?
        iCode = iCode or ""
        iCode = iCode[:1]
        atomtype = str(atomtype)[:2]  # make sure that is a string for user input
        partialCharge = float(partialCharge)  # Already pre-formatted when passed to this method
        self.pdb.write(self.fmt['ATOM'] % vars())

