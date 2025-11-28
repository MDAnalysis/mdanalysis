# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
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
   https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/AutoDock4.2.6_UserGuide.pdf
.. _AutoDock:
   http://autodock.scripps.edu/
"""

import os
import errno
import itertools
import numpy as np
import warnings

from ..lib import util
from . import base


class PDBQTReader(base.SingleFrameReaderBase):
    """PDBQTReader that reads a PDBQT-formatted file, no frills.

    Records read:
     - CRYST1 for unitcell A,B,C, alpha,beta,gamm
     - ATOM. HETATM for x,y,z

    Original `PDB format documentation`_ with `AutoDOCK extensions`_


    .. _PDB format documentation:
       http://www.wwpdb.org/documentation/file-format-content/format32/v3.2.html
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
    67 - 70        LString(4)    footnote     Usually blank. IGNORED.
    71 - 76        Real(6.4)     partialChrg  Gasteiger PEOE partial charge *q*.
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

    format = "PDBQT"
    units = {"time": None, "length": "Angstrom"}

    def _read_first_frame(self):
        coords = []
        unitcell = np.zeros(6, dtype=np.float32)
        with util.openany(self.filename) as pdbfile:
            for line in pdbfile:
                # Should only break at the 'END' of a model definition
                # and prevent premature exit for a torsion termination
                # , eg, ENDBRANCH
                if line.startswith("END\n"):
                    break
                if line.startswith("CRYST1"):
                    # lengths
                    x, y, z = np.float32(
                        (line[6:15], line[15:24], line[24:33])
                    )
                    # angles
                    A, B, G = np.float32(
                        (line[33:40], line[40:47], line[47:54])
                    )
                    unitcell[:] = x, y, z, A, B, G
                if line.startswith(("ATOM", "HETATM")):
                    # convert all entries at the end once for optimal speed
                    coords.append([line[30:38], line[38:46], line[46:54]])
        self.n_atoms = len(coords)
        self.ts = self._Timestep.from_coordinates(coords, **self._ts_kwargs)
        self.ts.dimensions = unitcell
        self.ts.frame = 0  # 0-based frame number
        if self.convert_units:
            # in-place !
            self.convert_pos_from_native(self.ts._pos)
            if self.ts.dimensions is not None:
                self.convert_pos_from_native(self.ts.dimensions[:3])

    def Writer(self, filename, **kwargs):
        """Returns a permissive (simple) PDBQTWriter for *filename*.

        Parameters
        ----------
        filename : str
            filename of the output PDBQT file

        Returns
        -------
        :class:`PDBQTWriter`

        """
        return PDBQTWriter(filename, **kwargs)


class PDBQTWriter(base.WriterBase):
    """PDBQT writer that implements a subset of the PDB_ 3.2 standard and the PDBQT_ spec.

    .. _PDB: http://www.wwpdb.org/documentation/file-format-content/format32/v3.2.html
    .. _PDBQT: https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/AutoDock4.2.6_UserGuide.pdf

    .. versionchanged:: 2.6.0
       Files are now written in `wt` mode, and keep extensions, allowing
       for files to be written under compressed formats
    """

    fmt = {
        "ATOM": (
            "ATOM  {serial:5d} {name:<4.4s} {resName:<4.4s}"
            "{chainID:1.1s}{resSeq:4d}{iCode:1.1s}"
            "   {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}{occupancy:6.2f}"
            "{tempFactor:6.2f}    {charge:< 1.3f} {element:<2.2s}\n"
        ),
        "REMARK": "REMARK     {0}\n",
        "TITLE": "TITLE     {0}\n",
        "CRYST1": (
            "CRYST1{box[0]:9.3f}{box[1]:9.3f}{box[2]:9.3f}"
            "{ang[0]:7.2f}{ang[1]:7.2f}{ang[2]:7.2f} "
            "{spacegroup:<11s}{zvalue:4d}\n"
        ),
    }
    format = "PDBQT"
    units = {"time": None, "length": "Angstrom"}
    pdb_coor_limits = {"min": -999.9995, "max": 9999.9995}

    def __init__(self, filename, **kwargs):
        self.filename = util.filename(filename, ext="pdbqt", keep=True)
        self.pdb = util.anyopen(self.filename, "wt")

    def close(self):
        self.pdb.close()

    def write(self, selection, frame=None):
        """Write selection at current trajectory frame to file.

        Parameters
        ----------
        selection : AtomGroup
            The selection to be written
        frame : int (optional)
            optionally move to frame index `frame` before writing; the default
            is to write the current trajectory frame

        Note
        ----
        The last letter of the
        :attr:`~MDAnalysis.core.groups.Atom.segid` is used as the PDB
        chainID.


        .. versionchanged:: 0.11.0
           Frames now 0-based instead of 1-based

        """
        try:
            u = selection.universe
        except AttributeError:
            errmsg = "Input obj is neither an AtomGroup or Universe"
            raise TypeError(errmsg) from None

        if frame is not None:
            u.trajectory[frame]  # advance to frame
        else:
            try:
                frame = u.trajectory.ts.frame
            except AttributeError:
                frame = 0  # should catch cases when we are analyzing a single PDB (?)

        atoms = selection.atoms  # make sure to use atoms (Issue 46)
        coor = (
            atoms.positions
        )  # can write from selection == Universe (Issue 49)

        # Check attributes
        attrs = {}
        missing_topology = []
        for attr, dflt in (
            ("altLocs", " "),
            ("charges", 0.0),
            ("icodes", " "),
            ("names", "X"),
            ("occupancies", 1.0),
            ("resids", 1),
            ("resnames", "UNK"),
            ("tempfactors", 0.0),
            ("types", "  "),
        ):
            try:
                attrs[attr] = getattr(atoms, attr)
            except AttributeError:
                attrs[attr] = itertools.cycle((dflt,))
                missing_topology.append(attr)
        # Order of preference: chainids -> segids -> blank string
        try:
            attrs["chainids"] = atoms.chainids
        except AttributeError:
            try:
                attrs["chainids"] = atoms.segids
            except AttributeError:
                attrs["chainids"] = itertools.cycle((" ",))
                missing_topology.append("chainids")
        if missing_topology:
            warnings.warn(
                "Supplied AtomGroup was missing the following attributes: "
                "{miss}. These will be written with default values. "
                "".format(miss=", ".join(missing_topology))
            )

        # check if any coordinates are illegal (coordinates are already
        # in Angstroem per package default)
        if not self.has_valid_coordinates(self.pdb_coor_limits, coor):
            self.close()
            try:
                os.remove(self.filename)
            except OSError as err:
                if err.errno == errno.ENOENT:
                    pass
            raise ValueError(
                "PDB files must have coordinate values between {0:.3f}"
                " and {1:.3f} Angstroem: No file was written."
                "".format(
                    self.pdb_coor_limits["min"], self.pdb_coor_limits["max"]
                )
            )

        # Write title record
        # http://www.wwpdb.org/documentation/file-format-content/format32/sect2.html
        line = "FRAME " + str(frame) + " FROM " + str(u.trajectory.filename)
        self.pdb.write(self.fmt["TITLE"].format(line))

        # Write CRYST1 record
        # http://www.wwpdb.org/documentation/file-format-content/format32/sect8.html
        box = self.convert_dimensions_to_unitcell(u.trajectory.ts)
        self.pdb.write(
            self.fmt["CRYST1"].format(
                box=box[:3], ang=box[3:], spacegroup="P 1", zvalue=1
            )
        )

        # Write atom records
        # http://www.wwpdb.org/documentation/file-format-content/format32/sect9.html
        for serial, (
            pos,
            name,
            resname,
            chainid,
            resid,
            icode,
            occupancy,
            tempfactor,
            charge,
            element,
        ) in enumerate(
            zip(
                coor,
                attrs["names"],
                attrs["resnames"],
                attrs["chainids"],
                attrs["resids"],
                attrs["icodes"],
                attrs["occupancies"],
                attrs["tempfactors"],
                attrs["charges"],
                attrs["types"],
            ),
            start=1,
        ):
            serial = util.ltruncate_int(serial, 5)  # check for overflow here?
            resid = util.ltruncate_int(resid, 4)
            name = name[:4]
            if len(name) < 4:
                name = " " + name  # customary to start in column 14
            chainid = chainid.strip()[-1:]  # take the last character

            self.pdb.write(
                self.fmt["ATOM"].format(
                    serial=serial,
                    name=name,
                    resName=resname,
                    chainID=chainid,
                    resSeq=resid,
                    iCode=icode,
                    pos=pos,
                    occupancy=occupancy,
                    tempFactor=tempfactor,
                    charge=charge,
                    element=element,
                )
            )

        self.close()
