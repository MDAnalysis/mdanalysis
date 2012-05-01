# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
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
PDB structure files in MDAnalysis --- :mod:`MDAnalysis.coordinates.PDB`
========================================================================

MDAnalysis reads coordinates from PDB files and additional optional
data such as B-factors. It is also possible to substitute a PDB file
instead of PSF file in order to define the list of atoms (but no
connectivity information will be  available in this case).

The :mod:`PDB` module makes heavy use of Biopython's :mod:`Bio.PDB`:

  Hamelryck, T., Manderick, B. (2003) PDB parser and structure class
  implemented in Python. Bioinformatics, 19, 2308-2310.

  http://biopython.org

but replaces the standard PDB file parser with one that uses the
:class:`MDAnalysis.coordinates.pdb.extensions.SloppyStructureBuilder`
to cope with very large pdb files as commonly encountered in MD
simulations.
"""
from __future__ import with_statement

try:
    # BioPython is overkill but potentially extensible (altLoc etc)
    import Bio.PDB
except ImportError:
    # TODO: fall back to PrimitivePDBReader
    raise ImportError("No full-feature PDB I/O functionality. Install biopython.")

import os, errno
import numpy

import MDAnalysis.core
import MDAnalysis.core.util as util
import base
import pdb.extensions

from MDAnalysis.topology.core import guess_atom_element
from MDAnalysis.core.AtomGroup import Universe, AtomGroup

import warnings

class Timestep(base.Timestep):
        @property
        def dimensions(self):
                """unitcell dimensions (`A, B, C, alpha, beta, gamma`)

                - `A, B, C` are the lengths of the primitive cell vectors `e1, e2, e3`
                - `alpha` = angle(`e1, e2`)
                - `beta` = angle(`e1, e3`)
                - `gamma` = angle(`e2, e3`)
                """
                # Layout of unitcell is [A,B,C,90,90,90] with the primitive cell vectors
                return self._unitcell

class PDBReader(base.Reader):
    """Read a pdb file into a BioPython pdb structure.

    The coordinates are also supplied as one numpy array and wrapped
    into a Timestep object; attributes are set so that the PDBReader
    object superficially resembles the DCDReader object.
    """
    format = 'PDB'
    units = {'time': None, 'length': 'Angstrom'}
    _Timestep = Timestep

    def __init__(self, pdbfilename, convert_units=None, **kwargs):
        self.pdbfilename = pdbfilename
        self.filename = self.pdbfilename
        if convert_units is None:
            convert_units = MDAnalysis.core.flags['convert_gromacs_lengths']
        self.convert_units = convert_units  # convert length and time to base units

        pdb_id = "0UNK"
        self.pdb = pdb.extensions.get_structure(pdbfilename, pdb_id)
        pos = numpy.array([atom.coord for atom in self.pdb.get_atoms()])
        self.numatoms = pos.shape[0]
        self.numframes = 1
        self.fixed = 0          # parse B field for fixed atoms?
        self.skip = 1
        self.periodic = False
        self.delta = 0
        self.skip_timestep = 1
        #self.ts._unitcell[:] = ??? , from CRYST1?
        self.ts = self._Timestep(pos)
        del pos
        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)             # in-place !

    def get_bfactors(self):
        """Return an array of bfactors (tempFactor) in atom order."""
        warnings.warn("get_bfactors() will be removed in MDAnalysis 0.8; "
                      "use AtomGroup.bfactors [which will become AtomGroup.bfactors()]",
                      DeprecationWarning)
        return numpy.array([a.get_bfactor() for a in self.pdb.get_atoms()])

    def Writer(self, filename, **kwargs):
        """Returns a strict PDBWriter for *filename*.

        :Arguments:
          *filename*
              filename of the output PDB file

        :Returns: :class:`PDBWriter`

        .. Note::

           This :class:`PDBWriter` 's :meth:`~PDBWriter.write` method
           always requires a :class:`Timestep` as an argument (it is
           not optional anymore when the Writer is obtained through
           this method of :class:`PDBReader`.)
        """
        # This is messy; we cannot get a universe from the Reader, which would be
        # also needed to be fed to the PDBWriter (which is a total mess...).
        # Hence we ignore the problem and document it in the doc string... --- the
        # limitation is simply that PDBWriter.write() must always be called with an argument.
        kwargs['BioPDBstructure'] = self.pdb   # make sure that this Writer is
        kwargs.pop('universe', None)           # always linked to this reader, don't bother with Universe
        return PDBWriter(filename, **kwargs)

    def __len__(self):
        return self.numframes
    def __iter__(self):
        def iterPDB():
            yield self.ts   # just a single frame available
            raise StopIteration
        return iterPDB()
    def __getitem__(self, frame):
        if frame != 0:
            raise IndexError('PDBReader only contains a single frame at index 0')
        return self.ts


class PDBWriter(base.Writer):
    """Write out the current time step as a pdb file.

    This is not cleanly implemented at the moment. One must supply a
    universe, even though this is nominally an optional argument. The
    class behaves slightly differently depending on if the structure
    was loaded from a PDB (then the full-fledged :mod:`Bio.PDB` writer is
    used) or if this is really only an atom selection (then a less
    sophistiocated writer is employed).

    .. Note::

      The standard PDBWriter can only write the *whole system*.  In
      order to write a selection, use the :class:`PrimitivePDBWriter`,
      which happens automatically when the
      :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.write` method of a
      :class:`~MDAnalysis.core.AtomGroup.AtomGroup` instance is used.
    """
    format = 'PDB'
    units = {'time': None, 'length': 'Angstrom'}

    # PDBWriter is a bit more complicated than the DCDWriter in the
    # sense that a DCD frame only contains coordinate information. The
    # PDB contains atom data as well and hence it MUST access the
    # universe. In order to present a unified (and backwards
    # compatible) interface we must keep the universe argument an
    # optional keyword argument even though it really is required.

    def __init__(self, pdbfilename, universe=None, multi=False, **kwargs):
        """pdbwriter = PDBWriter(<pdbfilename>,universe=universe,**kwargs)
        :Arguments:
        pdbfilename     filename; if multi=True, embed a %%d formatstring
                        so that write_next_timestep() can insert the frame number
        universe        supply a universe [really REQUIRED; optional only for compatibility]
        multi           False: write a single structure to a single pdb
                        True: write all frames to multiple pdb files
        """
        import Bio.PDB.Structure
        self.universe = universe
        self.PDBstructure = kwargs.pop('BioPDBstructure', None)  # hack for PDBReader.Writer()
        if not self.PDBstructure:
            try:
                self.PDBstructure = universe.trajectory.pdb
            except AttributeError:
                pass
        self.filename = pdbfilename
        self.multi = multi
        if self.multi:
            raise NotImplementedError('Sorry, multi=True does not work yet.')
        if self.PDBstructure is not None and not isinstance(self.PDBstructure, Bio.PDB.Structure.Structure):
            raise TypeError('If defined, PDBstructure must be a Bio.PDB.Structure.Structure, eg '
                            'Universe.trajectory.pdb.')
    def write_next_timestep(self, ts=None):
        self.write(ts)
    def write(self, ts=None):
        """Write timestep as a pdb file.

        If ts=None then we try to get the current one from the universe.
        """
        if self.PDBstructure is None:
            if self.universe is None:
                warnings.warn("PDBWriter: Not writing frame as neither Timestep nor Universe supplied.")
                return
            # primitive PDB writing (ignores timestep argument)
            ppw = PrimitivePDBWriter(self.filename)
            ppw.write(self.universe.selectAtoms('all'))
            ppw.close()
        else:
            # full fledged PDB writer
            # Let's cheat and use universe.pdb.pdb: modify coordinates
            # and save...
            if ts is None:
                try:
                    ts = self.universe.trajectory.ts
                except AttributeError:
                    warnings.warn("PDBWriter: Not writing frame as neither universe nor timestep supplied.")
                    return
            if not hasattr(ts, '_pos'):
                raise TypeError("The PDBWriter can only process a Timestep as optional argument, not "
                                "e.g. a selection. Use the PrimitivePDBWriter instead and see the docs.")
            for a, pos in zip(self.PDBstructure.get_atoms(), ts._pos):
                a.set_coord(pos)
            io = pdb.extensions.SloppyPDBIO()
            io.set_structure(self.PDBstructure)
            io.save(self.filename)


class PrimitivePDBReader(base.Reader):
    """PDBReader that reads a PDB-formatted file, no frills.

    Records read:
     - CRYST1 for unitcell A,B,C, alpha,beta,gamm
     - ATOM. HETATM for x,y,z

    http://www.wwpdb.org/documentation/format32/sect9.html

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
    67 - 76        String        segID        (unofficial CHARMM extension ?)
    77 - 78        LString(2)    element      Element symbol, right-justified. IGNORED
    79 - 80        LString(2)    charge       Charge  on the atom. IGNORED
    =============  ============  ===========  =============================================
    """
    format = 'PDB'
    units = {'time': None, 'length': 'Angstrom'}
    _Timestep = Timestep

    def __init__(self, filename, convert_units=None, **kwargs):
        """Read coordinates from *filename*.

        *filename* can be a gzipped or bzip2ed compressed PDB file.

        If the pdb file contains multiple MODEL records then it is
        read as a trajectory where the MODEL numbers correspond to
        frame numbers. Therefore, the MODEL numbers must be a sequence
        of integeres (typically starting at 1 or 0).
        """
        self.filename = filename
        if convert_units is None:
            convert_units = MDAnalysis.core.flags['convert_gromacs_lengths']
        self.convert_units = convert_units  # convert length and time to base units

        # add this number to the MODEL number to get frame: typically
        # PDB starts at 1 but Gromacs trjconv exports starting with
        # MODEL 0. (MODEL 0 is autodetected in the code, though!)
        self.model_offset = kwargs.pop("model_offset", -1)

        header = ""
        compound = []
        remark = []

        frames = {}
        coords = []
        atoms = []

        unitcell = numpy.zeros(6, dtype=numpy.float32)

        with util.openany(filename, 'r') as pdbfile:
            for i, line in enumerate(pdbfile):
                def _c(start, stop, typeclass=float):
                    return self._col(line, start, stop, typeclass=typeclass)
                if line.split()[0] == 'END':
                    break
                elif line[:6] == 'CRYST1':
                    A, B, C = _c(7, 15), _c(16, 24), _c(25, 33)
                    alpha, beta, gamma = _c(34, 40), _c(41, 47), _c(48, 54)
                    unitcell[:] = A, B, C, alpha, beta, gamma
                    continue
                elif line[:6] == 'HEADER':
                    header = line[6:-1]
                    continue
                elif line[:6] == 'COMPND':
                    l = line[6:-1]
                    compound.append(l)
                    continue
                elif line[:6] == 'REMARK':
                    content = line[6:-1]
                    remark.append(content)
                elif line[:5] == 'MODEL':
                    frameno = int(line.split()[1])
                    if frameno == 0:
                        # detect if MODEL starts at 0 or at 1; switch model_offset appropriately
                        # Will break/loose frames if MODELs are renumbered in the PDB file...
                        self.model_offset = 0
                    frames[frameno + self.model_offset] = i # 0-based indexing
                elif line[:6] in ('ATOM  ', 'HETATM'):
                    # skip atom/hetatm for frames other than the first - they will be read in when next() is called on the trajectory reader
                    if len(frames) > 1:
                        continue
                    # directly use COLUMNS from PDB spec
                    serial = _c(7, 11, int)
                    name = _c(13, 16, str).strip()
                    resName = _c(18, 21, str).strip()
                    chainID = _c(22, 22, str)  # empty chainID is a single space ' '!
                    resSeq = _c(23, 26, int)
                    x, y, z = _c(31, 38), _c(39, 46), _c(47, 54)
                    occupancy = _c(55, 60)
                    tempFactor = _c(61, 66)
                    segID = _c(67, 76, str).strip()
                    element = _c(77, 78, str).strip()
                    coords.append((x, y, z))
                    atoms.append((serial, name, resName, chainID, resSeq, occupancy, tempFactor, segID, element))
                    continue
        self.header = header
        self.compound = compound
        self.remark = remark
        self.numatoms = len(coords)
        self.ts = self._Timestep(numpy.array(coords, dtype=numpy.float32))
        self.ts._unitcell[:] = unitcell
        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)             # in-place !
            self.convert_pos_from_native(self.ts._unitcell[:3])    # in-place ! (only lengths)

        # No 'MODEL' entries
        if len(frames) == 0:
          frames[0] = None

        self.frames = frames
        self.numframes = len(frames) if frames else 1
        self.ts.frame = 0
        self.fixed = 0
        self.skip = 1
        self.periodic = False
        self.delta = 0
        self.skip_timestep = 1
        # hack for PrimitivePDBParser:
        self._atoms = numpy.rec.fromrecords(atoms, names="serial,name,resName,chainID,resSeq,occupancy,tempFactor,segID,element")

    def _col(self, line, start, stop, typeclass=float):
        """Pick out and convert the columns start-stop.

        Numbering starts at column 1 with *start* and includes *stop*;
        this is the convention used in FORTRAN (and also in the PDB format).

        :Returns: ``typeclass(line[start-1:stop])`` or
                  ``typeclass(0)`` if conversion fails
        """
        x = line[start - 1:stop]
        try:
            return typeclass(x)
        except ValueError:
            return typeclass(0)

    def get_bfactors(self):
        """Return an array of bfactors (tempFactor) in atom order."""
        warnings.warn("get_bfactors() will be removed in MDAnalysis 0.8; "
                      "use AtomGroup.bfactors [which will become AtomGroup.bfactors()]",
                      DeprecationWarning)
        return self._atoms.tempFactor

    def get_occupancy(self):
        """Return an array of occupancies in atom order."""
        return self._atoms.occupancy

    def Writer(self, filename, **kwargs):
        """Returns a permissive (simple) PDBWriter for *filename*.

        :Arguments:
          *filename*
              filename of the output PDB file

        :Returns: :class:`PrimitivePDBWriter`

        """
        return PrimitivePDBWriter(filename, **kwargs)

    def __del__(self):
        pass

    def close(self):
        pass

    def rewind(self):
        #self._currentframe = -1
        self.ts.frame = -1
        self._read_next_timestep()

    def next(self):
        self._read_next_timestep()

    def __len__(self):
        return self.numframes
    def __iter__(self):
        def iterPDB():
            for i in xrange(0, self.numframes, self.skip):  # FIXME: skip is not working!!!
                try: yield self._read_next_timestep()
                except IOError: raise StopIteration
        return iterPDB()

    def _read_next_timestep(self):
        #frame = self._currentframe + 1
        frame = self.frame + 1
        if not self.frames.has_key(frame):
            return False
        self.__getitem__(frame)
        #self._currentframe = frame
        self.ts.frame = frame
        return True

    def __getitem__(self, frame):
        if numpy.dtype(type(frame)) != numpy.dtype(int):
            raise TypeError
        if not self.frames.has_key(frame):
            IndexError('PrimitivePDBReader found only %d frames in this trajectory' % len(self.frames))
        line = self.frames[frame]
        if line is None:
          return
        coords = []
        atoms = []
        unitcell = numpy.zeros(6, dtype=numpy.float32)

        def _c(start, stop, typeclass=float):
            return self._col(line, start, stop, typeclass=typeclass)

        with util.openany(self.filename, 'r') as f:
            for i in xrange(line):
                f.next()        # forward to frame
            for line in f:
              if line.split()[0] == 'ENDMDL':
                  break
              elif line[:6] == 'CRYST1':
                  A, B, C = _c(7, 15), _c(16, 24), _c(25, 33)
                  alpha, beta, gamma = _c(34, 40), _c(41, 47), _c(48, 54)
                  unitcell[:] = A, B, C, alpha, beta, gamma
                  continue
              elif line[:6] in ('ATOM  ', 'HETATM'):
                  # we only care about coordinates
                  x, y, z = _c(31, 38), _c(39, 46), _c(47, 54)
                  # TODO import occupancy, bfactors - might these change?
                  # OB: possibly, at least that would be a useful feature of a PDB trajectory
                  coords.append((x, y, z))
                  continue

        # check if atom number chaneged
        if len(coords) != len(self.ts._pos):
            raise Exception("PrimitivePDBReader assumes that the number of atoms remains unchanged between frames; the current frame has %d, the next frame has %d atoms" % (len(self.ts._pos), len(coords)))

        self.ts = self._Timestep(numpy.array(coords, dtype=numpy.float32))
        self.ts._unitcell[:] = unitcell
        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)             # in-place !
            self.convert_pos_from_native(self.ts._unitcell[:3])    # in-place ! (only lengths)

        self.ts.frame = frame


class PrimitivePDBWriter(base.Writer):
    """PDB writer that implements a subset of the `PDB 3.2 standard`_.

    PDB format as used by NAMD/CHARMM: 4-letter resnames and segID, altLoc
    is written.

    .. _`PDB 3.2 standard`:
       http://www.wwpdb.org/documentation/format32/v3.2.html
    """
    #          1         2         3         4         5         6         7         8
    # 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.
    # ATOM__seria nameAres CressI   xxxxxxxxyyyyyyyyzzzzzzzzOCCUPAtempft          elCH
    # ATOM  %5d   %-4s %-3s %4d %1s %8.3f   %8.3f   %8.3f   %6.2f %6.2f           %2s
    #                 %1s  %1s                                                      %2d
    #            =        =      ===                                    ==========
    # ATOM  %5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2d
    # ATOM  %(serial)5d %(name)-4s%(altLoc)1s%(resName)-3s %(chainID)1s%(resSeq)4d%(iCode)1s   %(x)8.3f%(y)8.3f%(z)8.3f%(occupancy)6.2f%(tempFactor)6.2f          %(element)2s%(charge)2d

    # Strict PDB format:
    #fmt = {'ATOM':   "ATOM  %(serial)5d %(name)-4s%(altLoc)1s%(resName)-3s %(chainID)1s%(resSeq)4d%(iCode)1s   %(x)8.3f%(y)8.3f%(z)8.3f%(occupancy)6.2f%(tempFactor)6.2f          %(element)2s%(charge)2d\n",
    # PDB format as used by NAMD/CHARMM: 4-letter resnames and segID, altLoc ignored
    fmt = {'ATOM':   "ATOM  %(serial)5d %(name)-4s %(resName)-4s%(chainID)1s%(resSeq)4d%(iCode)1s   %(x)8.3f%(y)8.3f%(z)8.3f%(occupancy)6.2f%(tempFactor)6.2f      %(segID)-4s%(element)2s%(charge)2d\n",
           'REMARK': "REMARK%s\n",
           'COMPND': "COMPND%s\n",
           'HEADER': "HEADER%s\n",
           'TITLE':  "TITLE    %s\n",
           'MODEL' : "MODEL %8d\n",
           'ENDMDL' : "ENDMDL\n",
           'END'    : "END",
           'CRYST1': "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
           'CONECT': "CONECT%s\n"
           }
    format = 'PDB'
    units = {'time': None, 'length': 'Angstrom'}
    pdb_coor_limits = {"min":-999.9995, "max":9999.9995}

    def __init__(self, filename, numatoms, start=0, step=1, remarks="Created by PrimitvePDBWriter", convert_units=None):
        """Create a new PDBWriter

        :Arguments:
         *filename*
           name of output file
         *start*
           starting timestep
         *step*
           skip between subsequent timesteps
         *remarks*
           comments to annotate pdb file
         *convert_units*
           units are converted to the MDAnalysis base format; ``None`` selects
           the value of :data:`MDAnalysis.core.flags`['convert_gromacs_lengths']
        """

        self.filename = filename
        if convert_units is None:
            convert_units = MDAnalysis.core.flags['convert_gromacs_lengths']
        self.convert_units = convert_units  # convert length and time to base units

        self.frames_written = 0
        if start < 0:
          raise ValueError, "'Start' must be a positive value"

        self.start = start
        self.step = step
        self.pdbfile = file(self.filename, 'w')  # open file on init
        self.remarks = remarks
        self._write_pdb_title(start, step, remarks)

    def __del__(self):
        if hasattr(self, 'pdbfile') and self.pdbfile is not None:
            self.close()

    def close_trajectory(self):
        # Do i really need this?
        self.pdbfile.write("END")
        self.pdbfile.close()
        self.pdbfile = None

    def _write_pdb_title(self, start, step, remarks):
        self.TITLE("FRAME(S) FORM %d, SKIP %d; %s" % (start, step, remarks))

    def _write_pdb_header(self):
        if not self.obj or not hasattr(self.obj, 'universe'):
          return
        u = self.obj.universe
        self.HEADER(u.trajectory)
        self.COMPND(u.trajectory)
        self.REMARK(u.trajectory)
        self.CRYST1(self.convert_dimensions_to_unitcell(u.trajectory.ts))

    def _check_pdb_coordinates(self):
        """
        Check if the coordinate values fall within the range allowed for PDB files.
        Raises :exc:`ValueError` if not and attempts to delete the output file.
        """
        atoms = self.obj.atoms    # make sure to use atoms (Issue 46)
        coor = atoms.coordinates() # can write from selection == Universe (Issue 49)

        # check if any coordinates are illegal (coordinates are already in Angstroem per package default)
        if self.has_valid_coordinates(self.pdb_coor_limits, coor):
          return True
        # note the precarious close() here: we know that the file is open and we now
        # prepare to remove what we have already written (header and such)
        self.close()
        try:
            os.remove(self.filename)
        except OSError, err:
            if err.errno == errno.ENOENT:
                pass
        raise ValueError("PDB files must have coordinate values between %.3f and %.3f Angstroem: No file was written." %
                         (self.pdb_coor_limits["min"], self.pdb_coor_limits["max"]))

    def _write_pdb_bonds(self):
        """
        Writes out all the bond records; works only for Universe objects.
        """
        """
        TODO all bonds are written out, using the old atom numbers - this is
        incorrect. Once a selection is made, the atom numbers have to be updated
        (currently they are unmodified) and bonds have to be seleceted for,
        only if all the atoms for a bond are still present.

        The bonds should not be a list of ints, as are now, but a list of
        Atom objects.
        """

        if not self.obj or not hasattr(self.obj, 'universe'):
          return

        if not hasattr(self.obj.universe, 'bonds'):
          return

        bonds = self.obj.universe.bonds
        for i, conect in enumerate(bonds):
            self.CONECT(conect)

    def write(self, obj):
        """Write selection at current trajectory frame to file.

        write(selection,frame=FRAME)

        :Arguments:
          *selection*
            a :class:`~MDAnalysis.core.AtomGroup.AtomGroup`
          *frame*
            optionally move to frame *FRAME*

        .. Note::

           The first letter of the :attr:`~MDAnalysis.core.AtomGroup.Atom.segid`
           is used as the PDB chainID.
        """

        if isinstance(obj, Timestep):
            raise TypeError("PrimitivePDBWriter cannot write Timestep objects directly, since they lack topology information (atom names and types) required in PDB files")
        self.obj = obj  # remember obj for some of other methods  --- NOTE: this is an evil/lazy hack...
        ts, traj = None, None
        if hasattr(obj, 'universe') and not isinstance(obj, Universe):
            # For AtomGroup and children (Residue, ResidueGroup, Segment)
            ts = obj.ts
            traj = obj.universe.trajectory
        else:
            # For Universe only
            ts = obj.trajectory.ts
            traj = obj.trajectory

        if not (ts and traj):
            raise AssertionError("PrimitivePDBWriter couldn't extract trajectory and timestep information from an object; inheritance problem.")

        self._check_pdb_coordinates()
        self.trajectory, self.timestep = traj, ts
        self._write_pdb_header()
        self.write_all_timesteps()
        self._write_pdb_bonds()
        self.END()

    def write_all_timesteps(self):
        start, step = self.start, self.step

        traj = self.trajectory

        # Start from trajectory[0]/frame 1, if there are more than 1 frame.
        # If there is onyl 1 frame, the traj.frames is not like a python list:
        # accessing trajectory[-1] raises key error.
        if not start and traj.numframes > 1:
          start = traj.frame - 1

        for framenumber in xrange(start, len(traj), step):
            traj[framenumber]
            ts = self.timestep
            self.write_next_timestep(ts)
        # Set the trajectory to the starting position
        traj[start]

    def write_next_timestep(self, ts=None):
        ''' write a new timestep to the PDB file

            ts - timestep object containing coordinates to be written to dcd file
        '''
        if ts is None:
            if not hasattr(self, "ts"):
                raise Exception("PBDWriter: no coordinate data to write to trajectory file")
            else:
                ts = self.ts

        self._write_timestep(ts)

    def _write_timestep(self, ts):
        """Write a new timestep to file

        - does unit conversion if necessary
        - writes MODEL_ ... ENDMDL_ records to represent trajectory frames

        .. Note:: At the moment we do *not* write the NUMMDL_ record.

        .. _MODEL: http://www.wwpdb.org/documentation/format32/sect9.html#MODEL
        .. _ENDMDL: http://www.wwpdb.org/documentation/format32/sect9.html#ENDMDL
        .. _NUMMDL: http://www.wwpdb.org/documentation/format32/sect2.html#NUMMDL
        """

        traj = self.trajectory
        atoms = self.obj.atoms
        if self.convert_units:
            coor = self.convert_pos_to_native(ts._pos, inplace=False)
        else:
            coor = ts._pos

        if len(atoms) != len(coor):
            raise ValueError("Length of the atoms array is %d, this is different form the Timestep coordinate array %d" % (len(atoms), len(ts._pos)))

        # dont write model records, if only one frame
        if len(traj) > 1:
            self.MODEL(self.frames_written + 1)

        for i, atom in enumerate(atoms):
            # TODO Jan: see description in ATOM for why this check has to be made
            if not atom.bfactor:
                atom.bfactor = 0.
            self.ATOM(serial=i + 1, name=atom.name.strip(), resName=atom.resname.strip(),
                      resSeq=atom.resid, chainID=atom.segid.strip(), segID=atom.segid.strip(),
                      tempFactor=atom.bfactor,
                      x=coor[i, 0], y=coor[i, 1], z=coor[i, 2])
            # get bfactor, too, and add to output?
            # 'element' is auto-guessed from atom.name in ATOM()
        # dont write model records, if only one frame
        if len(traj) > 1:
            self.ENDMDL()
        self.frames_written += 1

    def HEADER(self, trajectory):
        """Write HEADER_ record.

        .. _HEADER: http://www.wwpdb.org/documentation/format32/sect2.html#HEADER
        """
        if not hasattr(trajectory, 'header'):
          return
        header = trajectory.header
        self.pdbfile.write(self.fmt['HEADER'] % header)

    def TITLE(self, *title):
        """Write TITLE_ record.

        .. _TITLE: http://www.wwpdb.org/documentation/format32/sect2.html
        """
        line = " ".join(title)    # should do continuation automatically
        self.pdbfile.write(self.fmt['TITLE'] % line)

    def REMARK(self, trajectory):
        """Write generic REMARK_ record (without number).

        See also `REMARK (update)`_.

        .. _REMARK: http://www.wwpdb.org/documentation/format32/remarks1.html
        .. _REMARK (update): http://www.wwpdb.org/documentation/format32/remarks2.html
        """
        if not hasattr(trajectory, 'remark'):
          return
        remarks = trajectory.remark

        for remark in remarks:
          self.pdbfile.write(self.fmt['REMARK'] % (remark))

    def COMPND(self, trajectory):
        if not hasattr(trajectory, 'compound'):
          return
        compound = trajectory.compound
        for c in compound:
          self.pdbfile.write(self.fmt['COMPND'] % c)

    def CRYST1(self, dimensions, spacegroup='P 1', zvalue=1):
        """Write CRYST1_ record.

        .. _CRYST1: http://www.wwpdb.org/documentation/format32/sect8.html
        """
        self.pdbfile.write(self.fmt['CRYST1'] % (tuple(dimensions) + (spacegroup, zvalue)))

    def MODEL(self, modelnumber):
        """Write the MODEL_ record.

        .. _MODEL: http://www.wwpdb.org/documentation/format32/sect9.html#MODEL
        """
        self.pdbfile.write(self.fmt['MODEL'] % modelnumber)

    def END(self):
        """Write END_ record.

        .. _END: http://www.wwpdb.org/documentation/format32/sect11.html#END
        """
        self.pdbfile.write(self.fmt['END'])

    def ENDMDL(self):
        """Write the ENDMDL_ record.

        .. _ENDMDL: http://www.wwpdb.org/documentation/format32/sect9.html#ENDMDL
        """
        self.pdbfile.write(self.fmt['ENDMDL'])

    def ATOM(self, serial=None, name=None, altLoc=None, resName=None, chainID=None,
             resSeq=None, iCode=None, x=None, y=None, z=None, occupancy=1.0, tempFactor=0.0,
             segID=None, element=None, charge=0):
        """Write ATOM_ record.

        Only some keword args are optional (altLoc, iCode, chainID), for some defaults are set.

        All inputs are cut to the maximum allowed length. For integer numbers
        the highest-value digits are chopped (so that the serial and reSeq
        wrap); for strings the trailing characters are chopped. The *last*
        character of *chainID* becomes the PDB *chainID* (unless it has the
        value "SYSTEM" (assigned by MDAnalysis if neither *segID* nor *chainID*
        were available), in which case the PDB will have an empty *chainID*).

        .. Warning: Floats are not checked and can potentially screw up the format.

        .. _ATOM: http://www.wwpdb.org/documentation/format32/sect9.html

        .. versionchanged:: 0.7.6
           If the *chainID* has the special value "SYSTEM" (case insensitive)
           then the chain is set to the empty string "".

        """

        # TODO Jan: PDBReader sets the bfactor value corretly to 0.0 if not
        # defined, DCD error does not. Thus when saving a Universe(DCDfile), the
        # bfactor is None rather than 0.0.
        #
        # Now in the `for` loop below we could remove the offending variable but
        # this gets us when the fromat `ATOM` expects float and not NoneType.
        #
        # Provisional solution for now: custom check if tempfactor is None, in
        # :meth:`_write_timestep`
        #
        # OB: The "provisional" solution is correct. ATOM requires the calling code
        #     to provide a sane value for tempFactor. This is because a tempFactor of
        #     0 might actually mean to some programs "ignore this atom". Hence we
        #     don't want to make this decision in the general I/O routine.

        for arg in ('serial', 'name', 'resName', 'resSeq', 'x', 'y', 'z',
                    'occupancy', 'tempFactor', 'charge'):
            if locals()[arg] is None:
                raise ValueError('parameter ' + arg + ' must be defined.')
        serial = int(str(serial)[-5:])  # check for overflow here?
        name = name[:4]
        if len(name) < 4:
            name = " " + name   # customary to start in column 14
        altLoc = altLoc or " "
        altLoc = altLoc[:1]
        resName = resName[:4]
        chainID = chainID or ""   # or should we provide a chainID such as 'A'?
        chainID = chainID if not chainID.upper() == "SYSTEM" else ""  # special case, new in 0.7.6
        chainID = chainID.strip()[-1:] # take the last character
        resSeq = int(str(resSeq)[-4:]) # check for overflow here?
        iCode = iCode or ""
        iCode = iCode[:1]
        element = element or guess_atom_element(name.strip())  # element == 0|False|None will be guessed
        element = str(element).strip()[:2]                     # make sure that is a string for user input
        segID = segID or chainID
        segID = segID[:4]
        self.pdbfile.write(self.fmt['ATOM'] % vars())

    def CONECT(self, conect):
        """Write CONECT_ record.

        .. _CONECT: http://www.wwpdb.org/documentation/format32/sect10.html#CONECT
        """
        conect = ["%5d" % entry for entry in conect]
        conect = "".join(conect)
        self.pdbfile.write(self.fmt['CONECT'] % conect)

