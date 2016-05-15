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


"""CRD structure files in MDAnalysis --- :mod:`MDAnalysis.coordinates.CRD`
===========================================================================

Read and write coordinates in CHARMM CARD coordinate format (suffix
"crd"). The CHARMM "extended format" is handled automatically.

"""

import numpy as np

from ..lib import util
from . import base


class CRDReader(base.SingleFrameReader):
    """CRD reader that implements the standard and extended CRD coordinate formats

    .. versionchanged:: 0.11.0
       Now returns a ValueError instead of FormatError
       Frames now 0-based instead of 1-based
    """
    format = 'CRD'
    units = {'time': None, 'length': 'Angstrom'}

    def _read_first_frame(self):
        # EXT:
        #      (i10,2x,a)  natoms,'EXT'
        #      (2I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A8,F20.10)
        #      iatom,ires,resn,typr,x,y,z,segid,rid,wmain
        # standard:
        #      (i5) natoms
        #      (2I5,1X,A4,1X,A4,3F10.5,1X,A4,1X,A4,F10.5)
        #      iatom,ires,resn,typr,x,y,z,segid,orig_resid,wmain
        coords_list = []
        with util.openany(self.filename, 'r') as crdfile:
            extended = False
            natoms = 0
            for linenum, line in enumerate(crdfile):
                if line.strip().startswith('*') or line.strip() == "":
                    continue  # ignore TITLE and empty lines
                fields = line.split()
                if len(fields) <= 2:
                    # should be the natoms line
                    natoms = int(fields[0])
                    extended = (fields[-1] == 'EXT')
                    continue
                # process coordinates
                try:
                    if extended:
                        coords_list.append(np.array(line[45:100].split()[0:3], dtype=float))
                    else:
                        coords_list.append(np.array(line[20:50].split()[0:3], dtype=float))
                except:
                    raise ValueError("Check CRD format at line {0}: {1}"
                                     "".format(linenum, line.rstrip()))

        self.n_atoms = len(coords_list)

        self.ts = self._Timestep.from_coordinates(np.array(coords_list),
                                                  **self._ts_kwargs)
        self.ts.frame = 0  # 0-based frame number
        # if self.convert_units:
        #    self.convert_pos_from_native(self.ts._pos)             # in-place !

        # sanity check
        if self.n_atoms != natoms:
            raise ValueError("Found %d coordinates in %r but the header claims that there "
                              "should be %d coordinates." % (self.n_atoms, self.filename, natoms))

    def Writer(self, filename, **kwargs):
        """Returns a CRDWriter for *filename*.

        :Arguments:
          *filename*
              filename of the output CRD file

        :Returns: :class:`CRDWriter`

        """
        return CRDWriter(filename, **kwargs)


class CRDWriter(base.Writer):
    """CRD writer that implements the CHARMM CRD coordinate format.

    It automatically writes the CHARMM EXT extended format if there
    are more than 99,999 atoms.

    .. versionchanged:: 0.11.0
       Frames now 0-based instead of 1-based
    """
    format = 'CRD'
    units = {'time': None, 'length': 'Angstrom'}

    fmt = {
        #crdtype = 'extended'
        #fortran_format = '(2I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A8,F20.10)'
        'ATOM_EXT': "%(serial)10d%(TotRes)10d  %(resName)-8s  %(name)-8s%(x)20.10f%(y)20.10f%(z)20.10f  %(chainID)-8s "
                    " %(resSeq)-8d%(tempFactor)20.10f\n",
        'NUMATOMS_EXT': "%10d  EXT\n",
        #crdtype = 'standard'
        #fortran_format = '(2I5,1X,A4,1X,A4,3F10.5,1X,A4,1X,A4,F10.5)'
        'ATOM': "%(serial)5d%(TotRes)5d %(resName)-4s %(name)-4s%(x)10.5f%(y)10.5f%(z)10.5f %(chainID)-4s %("
                "resSeq)-4d%(tempFactor)10.5f\n",
        'TITLE': "*%s\n",
        'NUMATOMS': "%5d\n",
    }

    def __init__(self, filename, **kwargs):
        self.filename = util.filename(filename, ext='crd')
        self.crd = None

    def write(self, selection, frame=None):
        """Write selection at current trajectory frame to file.

        write(selection,frame=FRAME)

        selection         MDAnalysis AtomGroup
        frame             optionally move to frame FRAME
        """
        u = selection.universe
        if frame is not None:
            u.trajectory[frame]  # advance to frame
        else:
            try:
                frame = u.trajectory.ts.frame
            except AttributeError:
                frame = 0  # should catch cases when we are analyzing a single PDB (?)

        atoms = selection.atoms  # make sure to use atoms (Issue 46)
        coor = atoms.positions  # can write from selection == Universe (Issue 49)
        with util.openany(self.filename, 'w') as self.crd:
            self._TITLE("FRAME " + str(frame) + " FROM " + str(u.trajectory.filename))
            self._TITLE("")
            self._NUMATOMS(len(atoms))
            current_resid = 0
            for i, atom in enumerate(atoms):
                if atoms[i].resid != atoms[i - 1].resid:
                    # note that this compares first and LAST atom on first iteration... but it works
                    current_resid += 1
                self._ATOM(serial=i + 1, resSeq=atom.resid, resName=atom.resname, name=atom.name,
                           x=coor[i, 0], y=coor[i, 1], z=coor[i, 2], chainID=atom.segid,
                           tempFactor=atom.bfactor, TotRes=current_resid, n_atoms=len(atoms))

    def _TITLE(self, *title):
        """Write TITLE record.
        """
        line = " ".join(title)  # should do continuation automatically
        line = line.strip()
        if len(line) > 0:
            line = " " + line
        self.crd.write(self.fmt['TITLE'] % line)

    def _NUMATOMS(self, n_atoms):
        """Write generic total number of atoms in system)
        """
        if n_atoms > 99999:
            self.crd.write(self.fmt['NUMATOMS_EXT'] % n_atoms)
        else:
            self.crd.write(self.fmt['NUMATOMS'] % n_atoms)

    def _ATOM(self, serial=None, resSeq=None, resName=None, name=None, x=None, y=None, z=None, chainID=None,
              tempFactor=0.0, TotRes=None, n_atoms=None):
        """Write ATOM record.

        All inputs are cut to the maximum allowed length. For integer
        numbers the highest-value digits are chopped (so that the
        serial and reSeq wrap); for strings the trailing characters
        are chopped.

        .. Warning:: Floats are not checked and can potentially screw up the format.
        """
        if tempFactor is None:
            tempFactor = 0.0  # atom.bfactor is None by default
        for arg in ('serial', 'name', 'resName', 'resSeq', 'x', 'y', 'z', 'tempFactor'):
            if locals()[arg] is None:
                raise ValueError('parameter ' + arg + ' must be defined.')

        chainID = chainID or ""  # or should we provide a chainID such as 'A'?
        if n_atoms > 99999:
            serial = int(str(serial)[-10:])  # check for overflow here?
            name = name[:8]
            resName = resName[:8]
            chainID = chainID[:8]
            resSeq = int(str(resSeq)[-8:])  # check for overflow here?
            TotRes = int(str(TotRes)[-10:])
            self.crd.write(self.fmt['ATOM_EXT'] % vars())
        else:
            serial = int(str(serial)[-5:])  # check for overflow here?
            name = name[:4]
            resName = resName[:4]
            chainID = chainID[:4]
            resSeq = int(str(resSeq)[-4:])  # check for overflow here?
            TotRes = int(str(TotRes)[-5:])
            self.crd.write(self.fmt['ATOM'] % vars())

