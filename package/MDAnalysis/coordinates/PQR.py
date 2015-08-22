# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
PQR file format --- :mod:`MDAnalysis.coordinates.PQR`
=====================================================

Read atoms with charges from a PQR_ file (as written by PDB2PQR_). The
following is adopted from the description of the PQR_ format as used by APBS_:

*MDAnalysis* reads very loosely-formatted PQR files: all fields are
**whitespace-delimited** rather than the strict column formatting mandated
by the PDB_ format. This more liberal formatting allows coordinates
which are larger/smaller than ±999 Å.

MDAnalysis reads data on a per-line basis from PQR files using the following format::

   recordName serial atomName residueName chainID residueNumber X Y Z charge radius

If this fails it is assumed that the *chainID* was omitted and the shorter
format is read::

   recordName serial atomName residueName residueNumber X Y Z charge radius

Anything else will raise a :exc:`ValueError`.

The whitespace is the most important feature of this format: fields
*must* be separated by at least one space or tab character. The fields
are:

*recordName*
    A string which specifies the type of PQR entry and should either be ATOM or
    HETATM.
*serial*
    An integer which provides the atom index (but note that MDAnalysis renumbers
    atoms so one cannot rely on the *serial*)
*atomName*
    A string which provides the atom name.
*residueName*
    A string which provides the residue name.
*chainID*
    An optional string which provides the chain ID of the atom.
*residueNumber*
    An integer which provides the residue index.
*X Y Z*
    Three floats which provide the atomic coordiantes.
*charge*
    A float which provides the atomic charge (in electrons).
*radius*
    A float which provides the atomic radius (in Å).

Clearly, this format can deviate wildly from PDB_ due to the use of whitespaces
rather than specific column widths and alignments. This deviation can be
particularly significant when large coordinate values are used.

.. Warning:: Fields *must be white-space separated* or data are read
             incorrectly. PDB formatted files are *not* guaranteed to be
             white-space separated so extra care should be taken when quickly
             converting PDB files into PQR files using simple scripts.

For example, PQR files created with PDB2PQR_ and the `--whitespace`
option are guaranteed to conform to the above format::

   pdb2pqr --ff=charmm --whitespace 4ake.pdb 4ake.pqr

.. _PQR:     http://www.poissonboltzmann.org/file-formats/biomolecular-structurw/pqr
.. _APBS:    http://www.poissonboltzmann.org/apbs
.. _PDB2PQR: http://www.poissonboltzmann.org/pdb2pqr
.. _PDB:     http://www.rcsb.org/pdb/info.html#File_Formats_and_Standards
"""

import numpy as np

from ..core import flags
from ..lib import util
from . import base


class PQRReader(base.SingleFrameReader):
    """Read a PQR_ file into MDAnalysis.

    The :mod:`~MDAnalysis.topology.PQRParser` takes charges from the
    PQR file in order to populate the
    :attr:`MDAnalysis.core.AtomGroup.Atom.charge` attribute. Radii are
    accessible through the :meth:`get_radii` method of the reader, the
    :meth:`MDAnalysis.core.AtomGroup.AtomGroup.radii` method and the
    :attr:`MDAnalysis.core.AtomGroup.Atom.radius` attribute.

    .. _PQR:
        http://www.poissonboltzmann.org/file-formats/biomolecular-structurw/pqr

    .. versionchanged:: 0.11.0
       Frames now 0-based instead of 1-based
    """
    format = 'PQR'
    units = {'time': None, 'length': 'Angstrom'}

    def _read_first_frame(self):
        coords = []
        atoms = []
        unitcell = np.zeros(6, dtype=np.float32)
        segID = ''  # use empty string (not in PQR), PQRParsers sets it to SYSTEM
        with util.openany(self.filename, 'r') as pqrfile:
            for line in pqrfile:
                if line[:6] in ('ATOM  ', 'HETATM'):
                    fields = line.split()
                    try:
                        recordName, serial, name, resName, chainID, resSeq, x, y, z, charge, radius = fields
                    except ValueError:
                        # files without the chainID
                        recordName, serial, name, resName, resSeq, x, y, z, charge, radius = fields
                        chainID = ''
                    coords.append((float(x), float(y), float(z)))
                    atoms.append(
                        (int(serial), name, resName, chainID, int(resSeq), float(charge), float(radius), segID))
        self.n_atoms = len(coords)
        self.ts = self._Timestep.from_coordinates(np.array(coords, dtype=np.float32),
                                                  **self._ts_kwargs)
        self.ts._unitcell[:] = unitcell
        self.ts.frame = 0  # 0-based frame number
        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)  # in-place !
            self.convert_pos_from_native(self.ts._unitcell[:3])  # in-place ! (only lengths)

        # hack for PQRParser:
        self._atoms = np.rec.fromrecords(atoms, names="serial,name,resName,chainID,resSeq,charge,radius,segID")

    def get_radii(self):
        """Return an array of atom radii in atom order."""
        return self._atoms.radius

    def get_charges(self):
        """Return an array of charges in atom order."""
        return self._atoms.charge

    def Writer(self, filename, **kwargs):
        """Returns a PQRWriter for *filename*.

        :Arguments:
           *filename*
              filename of the output PQR file

        :Returns: :class:`PQRWriter`
        """
        return PQRWriter(filename, **kwargs)


class PQRWriter(base.Writer):
    """Write a single coordinate frame in whitespace-separated PQR format.

    Charges ("Q") are taken from the
    :attr:`MDAnalysis.core.AtomGroup.Atom.charge` attribute while
    radii are obtaine from the
    :attr:`MDAnalysis.core.AtomGroup.Atom.radius` attribute.

    * If the segid is 'SYSTEM' then it will be set to the empty
      string. Otherwise the first letter will be used as the chain ID.
    * The serial number always starts at 1 and increments sequentially
      for the atoms.

    The output format is similar to ``pdb2pqr --whitespace``.

    .. versionadded:: 0.9.0
    """
    format = 'PQR'
    units = {'time': None, 'length': 'Angstrom'}

    fmt = {
        'ATOM_nochain':
        "ATOM {0:6d} {1:<4}  {2:<3} {4:4d}   {5[0]:-8.3f} {5[1]:-8.3f} {5[2]:-8.3f} {6:-7.4f} {7:6.4f}\n",
        # serial, atomName, residueName, (chainID), residueNumber, XYZ, charge, radius
        'ATOM_chain':
        "ATOM {0:6d} {1:<4}  {2:<3} {3:1.1} {4:4d}   {5[0]:-8.3f} {5[1]:-8.3f} {5[2]:-8.3f} {6:-7.4f} {7:6.4f}\n",
        # serial, atomName, residueName, chainID, residueNumber, XYZ, charge, radius
    }

    def __init__(self, filename, convert_units=None, **kwargs):
        """Set up a PQRWriter with full whitespace separation.

        :Arguments:
          *filename*
             output filename
          *remarks*
             remark lines (list of strings) or single string to be added to the
             top of the PQR file
        """
        self.filename = util.filename(filename, ext='pqr')

        if convert_units is None:
            convert_units = flags['convert_lengths']
        self.convert_units = convert_units  # convert length and time to base units

        self.remarks = kwargs.pop('remarks', "PQR file written by MDAnalysis")

    def write(self, selection, frame=None):
        """Write selection at current trajectory frame to file.

        :Arguments:
            *selection*
                MDAnalysis AtomGroup (selection or Universe.atoms)
                or also Universe
         :Keywords:
             *frame*
                optionally move to frame number *frame*

        .. versionchanged:: 0.11.0
           Frames now 0-based instead of 1-based
        """
        # write() method that complies with the Trajectory API
        u = selection.universe
        if frame is not None:
            u.trajectory[frame]  # advance to frame
        else:
            try:
                frame = u.trajectory.ts.frame
            except AttributeError:
                frame = 0  # should catch cases when we are analyzing a single frame(?)

        atoms = selection.atoms  # make sure to use atoms (Issue 46)
        coordinates = atoms.coordinates()  # can write from selection == Universe (Issue 49)
        if self.convert_units:
            self.convert_pos_to_native(coordinates)  # inplace because coordinates is already a copy

        with util.openany(self.filename, 'w') as pqrfile:
            # Header
            self._write_REMARK(pqrfile, self.remarks)
            self._write_REMARK(pqrfile, "Input: frame {0} of {1}".format(frame, u.trajectory.filename), 5)
            self._write_REMARK(pqrfile, "total charge: {0:+8.4f} e".format(atoms.total_charge()), 6)
            # Atom descriptions and coords
            for atom_index, atom in enumerate(atoms):
                XYZ = coordinates[atom_index]
                self._write_ATOM(pqrfile, atom_index + 1, atom.name, atom.resname, atom.segid, atom.resid, XYZ,
                                 atom.charge, atom.radius)

    def _write_REMARK(self, fh, remarks, remarknumber=1):
        """Write REMARK record.

        The *remarknumber* is typically 1 but :program:`pdb2pgr`
        also uses 6 for the total charge and 5 for warnings.
        """
        for line in util.asiterable(remarks):  # either one line or multiple lines
            fh.write("REMARK   {0} {1}\n".format(remarknumber, line))

    def _write_ATOM(self, fh, serial, atomName, residueName, chainID, residueNumber, XYZ, charge, radius):
        """Write ATOM record.

        Output should look like this (although the only real
        requirement is *whitespace* separation between *all*
        entries). The chainID is optional and can be omitted::

              ATOM      1  N    MET     1     -11.921   26.307   10.410 -0.3000 1.8500
              ATOM     36  NH1  ARG     2      -6.545   25.499    3.854 -0.8000 1.8500
              ATOM     37 HH11  ARG     2      -6.042   25.480    4.723  0.4600 0.2245
        """
        ATOM = self.fmt['ATOM_nochain'] if (chainID == "SYSTEM" or not chainID) else self.fmt['ATOM_chain']
        atomName = (" " + atomName) if len(
            atomName) < 4 else atomName  # pad so that only 4-letter atoms are left-aligned
        fh.write(ATOM.format(serial, atomName, residueName, chainID, residueNumber, XYZ, charge, radius))
