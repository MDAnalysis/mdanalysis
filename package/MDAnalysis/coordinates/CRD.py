# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://mdanalysis.googlecode.com
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

import MDAnalysis
import MDAnalysis.core.util as util
import base
from MDAnalysis import FormatError
import numpy
import base


class Timestep(base.Timestep):
    @property
    def dimensions(self):
        """unitcell dimensions (*A*, *B*, *C*, *alpha*, *beta*, *gamma*)

        CRD files do not contain unitcell information but in order to
        allow interoperability (and give the use a chance to set the
        simulation box themselves for e.g. writing out to different
        formats) we add an empty unit cell, i.e. when reading a CRD
        file this will only contain zeros.

        lengths *a*, *b*, *c* are in the MDAnalysis length unit (Å), and
        angles are in degrees.
        """
        return self._unitcell

    @dimensions.setter
    def dimensions(self, box):
        self._unitcell[:] = box

class CRDReader(base.Reader):
    """CRD reader that implements the standard and extended CRD coordinate formats
    """
    format = 'CRD'
    units = {'time': None, 'length': 'Angstrom'}
    _Timestep = Timestep

    def __init__(self, crdfilename, convert_units=None, **kwargs):
        # EXT:
        #      (i10,2x,a)  natoms,'EXT'
        #      (2I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A8,F20.10)
        #      iatom,ires,resn,typr,x,y,z,segid,rid,wmain
        # standard:
        #      (i5) natoms
        #      (2I5,1X,A4,1X,A4,3F10.5,1X,A4,1X,A4,F10.5)
        #      iatom,ires,resn,typr,x,y,z,segid,orig_resid,wmain

        self.crdfilename = crdfilename
        self.filename = self.crdfilename
        if convert_units is None:
            # Note: not used at the moment in CRDReader/Writer
            convert_units = MDAnalysis.core.flags['convert_lengths']
        self.convert_units = convert_units  # convert length and time to base units
        coords_list = []
        with util.openany(crdfilename, 'r') as crdfile:
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
                        coords_list.append(numpy.array(map(float, line[45:100].split()[0:3])))
                    else:
                        coords_list.append(numpy.array(map(float, line[20:50].split()[0:3])))
                except:
                    raise FormatError("Check CRD format at line %d: %s" % (linenum, line.rstrip()))

        self.numatoms = len(coords_list)
        self.numframes = 1
        self.fixed = 0  # parse wmain field for fixed atoms?
        self.skip = 1
        self.periodic = False
        self.delta = 0
        self.skip_timestep = 1
        self.ts = self._Timestep(numpy.array(coords_list))
        self.ts.frame = 1  # 1-based frame number
        # if self.convert_units:
        #    self.convert_pos_from_native(self.ts._pos)             # in-place !

        # sanity check
        if self.numatoms != natoms:
            raise FormatError("Found %d coordinates in %r but the header claims that there "
                              "should be %d coordinates." % (self.numatoms, self.filename, natoms))

    def Writer(self, filename, **kwargs):
        """Returns a CRDWriter for *filename*.

        :Arguments:
          *filename*
              filename of the output CRD file

        :Returns: :class:`CRDWriter`

        """
        return CRDWriter(filename, **kwargs)

    def __iter__(self):
        yield self.ts  # Just a single frame
        raise StopIteration

    def _read_frame(self, frame):
        if frame != 0:
            raise IndexError("CRD only contains a single frame at frame index 0")
        return self.ts

    def _read_next_timestep(self):
        # CRD files only contain a single frame
        raise IOError


class CRDWriter(base.Writer):
    """CRD writer that implements the CHARMM CRD coordinate format.

    It automatically writes the CHARMM EXT extended format if there
    are more than 99,999 atoms.
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
                frame = 1  # should catch cases when we are analyzing a single PDB (?)

        atoms = selection.atoms  # make sure to use atoms (Issue 46)
        coor = atoms.coordinates()  # can write from selection == Universe (Issue 49)
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
                           tempFactor=atom.bfactor, TotRes=current_resid, numatoms=len(atoms))

    def _TITLE(self, *title):
        """Write TITLE record.
        """
        line = " ".join(title)  # should do continuation automatically
        line = line.strip()
        if len(line) > 0:
            line = " " + line
        self.crd.write(self.fmt['TITLE'] % line)

    def _NUMATOMS(self, numatoms):
        """Write generic total number of atoms in system)
        """
        if numatoms > 99999:
            self.crd.write(self.fmt['NUMATOMS_EXT'] % numatoms)
        else:
            self.crd.write(self.fmt['NUMATOMS'] % numatoms)

    def _ATOM(self, serial=None, resSeq=None, resName=None, name=None, x=None, y=None, z=None, chainID=None,
              tempFactor=0.0, TotRes=None, numatoms=None):
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
        if numatoms > 99999:
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

