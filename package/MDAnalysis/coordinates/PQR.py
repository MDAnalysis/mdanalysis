# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
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
    A float which provides the atomic radius (in Ã…).

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

import numpy

import MDAnalysis.core
import MDAnalysis.core.util as util
import base
from base import Timestep
import pdb.extensions

class PQRReader(base.Reader):
    """Read a PQR_ file into MDAnalysis.

    The :mod:`~MDAnalysis.topology.PQRParser` takes charges from the
    PQR file in order to populate the
    :attr:`MDAnalysis.core.AtomGroup.Atom.charge` attribute. Radii are
    accessible through the :meth:`get_radii` method of the reader, the
    :meth:`MDAnalysis.core.AtomGroup.AtomGroup.radii` method and the
    :attr:`MDAnalysis.core.AtomGroup.Atom.radius` attribute.
    """
    format = 'PQR'
    units = {'time': None, 'length': 'Angstrom'}
    _Timestep = Timestep

    def __init__(self, filename, convert_units=None, **kwargs):
        """Read coordinates from *filename*.

        *filename* can be a gzipped or bzip2ed compressed PQR_ file.

        .. _PQR:
           http://www.poissonboltzmann.org/file-formats/biomolecular-structurw/pqr
        """
        self.filename = filename
        if convert_units is None:
            convert_units = MDAnalysis.core.flags['convert_gromacs_lengths']
        self.convert_units = convert_units  # convert length and time to base units

        coords = []
        atoms = []
        unitcell = numpy.zeros(6, dtype=numpy.float32)
        segID = ''  # use empty string (not in PQR), PQRParsers sets it to SYSTEM
        with util.openany(filename, 'r') as pqrfile:
            for line in pqrfile:
                if line[:6] in ('ATOM  ', 'HETATM'):
                    fields = line.split()
                    try:
                        recordName,serial,name,resName,chainID,resSeq,x,y,z,charge,radius = fields
                    except ValueError:
                        # files without the chainID
                        recordName,serial,name,resName,resSeq,x,y,z,charge,radius = fields
                        chainID = 'A'
                    coords.append((float(x),float(y),float(z)))
                    atoms.append((int(serial), name, resName, chainID, int(resSeq), float(charge), float(radius),segID))
        self.numatoms = len(coords)
        self.ts = self._Timestep(numpy.array(coords, dtype=numpy.float32))
        self.ts._unitcell[:] = unitcell
        self.ts.frame = 1  # 1-based frame number
        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)             # in-place !
            self.convert_pos_from_native(self.ts._unitcell[:3])    # in-place ! (only lengths)
        self.numframes = 1
        self.fixed = 0
        self.skip = 1
        self.periodic = False
        self.delta = 0
        self.skip_timestep = 1
        # hack for PQRParser:
        self._atoms = numpy.rec.fromrecords(atoms, names="serial,name,resName,chainID,resSeq,charge,radius,segID")

    def get_radii(self):
        """Return an array of atom radii in atom order."""
        return self._atoms.radius

    def get_charges(self):
        """Return an array of charges in atom order."""
        return self._atoms.charge

    def __iter__(self):
        yield self.ts  # Just a single frame
        raise StopIteration

    def _read_frame(self, frame):
        if frame != 0:
            raise IndexError("PQR only contains a single frame at frame index 0")
        return self.ts

    def _read_next_timestep(self):
        # PQR file only contains a single frame
        raise IOError
