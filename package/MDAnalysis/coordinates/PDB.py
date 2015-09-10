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
PDB structure files in MDAnalysis --- :mod:`MDAnalysis.coordinates.PDB`
========================================================================

MDAnalysis reads coordinates from PDB files and additional optional
data such as B-factors. It is also possible to substitute a PDB file
instead of PSF file in order to define the list of atoms (but no
connectivity information will be available in this case).

PDB files contain both coordinate and atom information. It is also possible to
write trajectories as multi-frame (or multi-model) PDB files. This is not very
space efficient but is sometimes the lowest common denominator for exchanging
trajectories. Single frame and multi-frame PDB files are automatically
recognized; the only difference is thath the single-frame PDB is represented as
a trajectory with only one frame.

In order to write a selection to a PDB file one typically uses the
:meth:`MDAnalysis.core.AtomGroup.AtomGroup.write` method of the selection::

  calphas = universe.select_atoms("name CA")
  calphas.write("calpha_only.pdb")

This uses the coordinates from the current timestep of the trajectory.

In order to write a PDB trajectory one needs to first obtain a multi-frame
writer (keyword *multiframe* = ``True``) and then iterate through the
trajectory, while writing each frame::

  calphas = universe.select_atoms("name CA")
  W = MDAnalysis.Writer("calpha_traj.pdb", multiframe=True)
  for ts in u.trajectory:
      W.write(calphas)
  W.close()

It is important to *always close the trajectory* when done because only at this
step is the final END_ record written, which is required by the `PDB
standard`_.


Implementations
---------------

Two different implementations of PDB I/O are available: the ":ref:`permissive<permissive>`"
and the ":ref:`strict<strict>`" Reader/Writers. The default are the "permissive" ones
but this can be changed by setting the flag "permissive_pdb_reader" in
:data:`MDAnalysis.core.flags` (see :ref:`flags-label`) to ``False``::

   MDAnalysis.core.flags["permissive_pdb_reader"] = False

The *default for MDAnalysis* is to use the
":ref:`permissive<permissive>`" :class:`PrimitivePDBReader` and
:class:`PrimitivePDBWriter`, corresponding to ::

   MDAnalysis.core.flags["permissive_pdb_reader"] = True

On a case-by-case basis one kind of reader can be selected with the
*permissive* keyword to :class:`~MDAnalysis.core.AtomGroup.Universe`, e.g. ::

  u = MDAnalysis.Universe(PDB, permissive=False)

would select :class:`PDBReader` instead of the default
:class:`PrimitivePDBReader`.

.. _permissive:

Simple (permissive) PDB Reader and Writer
-----------------------------------------

A pure-Python implementation for PDB files commonly encountered in MD
simulations comes under the names :class:`PrimitivePDBReader` and
:class:`PrimitivePDBWriter`. It only implements a subset of the `PDB standard`_
(for instance, it does not deal with insertion codes) and also allows some
typical enhancements such as 4-letter resids (introduced by CHARMM/NAMD). The
"primitive PDB Reader/Writer" are the *default* in MDAnalysis (equivalent to
supplying the *permissive* = ``True`` keyword to
:class:`~MDAnalysis.core.AtomGroup.Universe`).

The :class:`PrimitivePDBReader` can read multi-frame PDB files and represents
them as a trajectory. The :class:`PrimitivePDBWriter` can write single and
multi-frame PDB files as specified by the *multiframe* keyword. By default, it
writes single frames. On the other hand, the :class:`MultiPDBWriter` is set up
to write a PDB trajectory by default (equivalent to using *multiframe* =
``True``).

Examples
~~~~~~~~

In order to write a **multi-frame PDB trajectory** from a universe *u* one can
do the following::

  pdb = MDAnalysis.Writer("all.pdb", multiframe=True)
  for ts in u.trajectory:
      pdb.write(u)
  pdb.close()

Similarly, writing only a protein::

  pdb = MDAnalysis.Writer("protein.pdb", multiframe=True)
  protein = u.select_atoms("protein")
  for ts in u.trajectory:
      pdb.write(protein)
  pdb.close()

A single frame can be written with the
:meth:`~MDAnalysis.core.AtomGroup.AtomGroup.write` method of any
:class:`~MDAnalysis.core.AtomGroup.AtomGroup`::

   protein = u.select_atoms("protein")
   protein.write("protein.pdb")

Alternatively, get the single frame writer and supply the
:class:`~MDAnalysis.core.AtomGroup.AtomGroup`::

  pdb = MDAnalysis.Writer("protein.pdb")
  protein = u.select_atoms("protein")
  pdb.write(protein)
  pdb.close()


Classes
~~~~~~~

.. autoclass:: PrimitivePDBReader
   :members:

.. autoclass:: PrimitivePDBWriter
   :members:

   .. automethod:: _check_pdb_coordinates
   .. automethod:: _write_pdb_bonds
   .. automethod:: _update_frame
   .. automethod:: _write_timestep

.. autoclass:: MultiPDBWriter
   :members:

.. _strict:

Biopython (strict) PDB Reader and Writer
----------------------------------------

The :mod:`PDB` module can make use of Biopython's :mod:`Bio.PDB`
[Hamelryck2003]_ but replaces the standard PDB file parser with one that uses
the :class:`MDAnalysis.coordinates.pdb.extensions.SloppyStructureBuilder` to
cope with very large pdb files as commonly encountered in MD simulations. The
Biopython-based :class:`PDBReader` has the advantage that it implements the
`PDB standard`_ rigorously but this comes at the cost of flexibility and
performance. It is also difficult to write out selections using this
implementation (:class:`PDBWriter`) and multi frame PDB files are not
implemented. The Biopython Reader/Writer can be selected when loading data into
a :class:`~MDAnalysis.core.AtomGroup.Universe` by providing the keyword
*permissive* = ``False``.

The Biopython PDB parser :class:`Bio.PDB.PDBParser` is fairly strict and even
in its own permissive mode (which MDAnalysis employs) typically warns about
missing element names with a
:exc:`Bio.PDB.PDBExceptions.PDBConstructionWarning` . Such warnings, however,
are generally harmless and therefore are filtered (and ignored) by MDAnalysis
with the help of :func:`warnings.filterwarnings`.


Classes
~~~~~~~

.. autoclass:: PDBReader
   :members:

.. autoclass:: PDBWriter
   :members:

References
----------

.. [Hamelryck2003]  Hamelryck, T., Manderick, B. (2003) PDB parser and structure class
                    implemented in Python. Bioinformatics, 19, 2308-2310. http://biopython.org

.. _PDB standard: http://www.wwpdb.org/documentation/format32/v3.2.html
.. _END: http://www.wwpdb.org/documentation/format32/sect11.html#END

"""

try:
    # BioPython is overkill but potentially extensible (altLoc etc)
    import Bio.PDB
    import pdb.extensions
    # disable PDBConstructionWarning from picky builder
    import warnings

    warnings.filterwarnings("ignore", category=Bio.PDB.PDBExceptions.PDBConstructionWarning,
                            message="Could not assign element|Used element .* for Atom")
except ImportError:
    # TODO: fall back to PrimitivePDBReader
    raise ImportError("No full-feature PDB I/O functionality. Install biopython.")

import os
import errno
import textwrap
import warnings
import logging
import numpy as np

from ..core import flags
from ..lib import util
from . import base
from ..topology.core import guess_atom_element
from ..core.AtomGroup import Universe
from ..exceptions import NoDataError


logger = logging.getLogger("MDAnalysis.coordinates.PBD")


class PDBReader(base.SingleFrameReader):
    """Read a pdb file into a :mod:`BioPython.PDB` structure.

    The coordinates are also supplied as one numpy array and wrapped
    into a Timestep object.

    .. Note:: The Biopython.PDB reader does not parse the ``CRYST1``
              record and hence the unitcell dimensions are not set.
              Use the :class:`PrimitivePDBReader` instead (i.e.  use
              the ``primitive=True`` keyword for :class:`Universe`).

    .. versionchanged:: 0.11.0
       * Frames now 0-based instead of 1-based.
       * All PDB header metadata parsed by the reader is available in
         the dict :attr:`metadata`.

    """
    format = 'PDB'
    units = {'time': None, 'length': 'Angstrom'}

    def _read_first_frame(self):
        pdb_id = "0UNK"
        self.pdb = pdb.extensions.get_structure(self.filename, pdb_id)
        pos = np.array([atom.coord for atom in self.pdb.get_atoms()])
        self.n_atoms = pos.shape[0]
        self.fixed = 0  # parse B field for fixed atoms?
        #self.ts._unitcell[:] = ??? , from CRYST1? --- not implemented in Biopython.PDB
        self.ts = self._Timestep.from_coordinates(pos, **self._ts_kwargs)
        self.ts.frame = 0
        del pos
        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)  # in-place !
        # metadata
        self.metadata = self.pdb.header

    def get_bfactors(self):
        """Return an array of bfactors (tempFactor) in atom order."""
        warnings.warn("get_bfactors() will be removed in MDAnalysis 0.8; "
                      "use AtomGroup.bfactors [which will become AtomGroup.bfactors()]",
                      DeprecationWarning)
        return np.array([a.get_bfactor() for a in self.pdb.get_atoms()])

    def Writer(self, filename, **kwargs):
        """Returns a strict PDBWriter for *filename*.

        :Arguments:
          *filename*
              filename of the output PDB file

        :Returns: :class:`PDBWriter`

        .. Note::

           This :class:`PDBWriter` 's :meth:`~PDBWriter.write` method
           always requires a :class:`base.Timestep` as an argument (it is
           not optional anymore when the Writer is obtained through
           this method of :class:`PDBReader` .)
        """
        # This is messy; we cannot get a universe from the Reader, which would be
        # also needed to be fed to the PDBWriter (which is a total mess...).
        # Hence we ignore the problem and document it in the doc string... --- the
        # limitation is simply that PDBWriter.write() must always be called with an argument.
        kwargs['BioPDBstructure'] = self.pdb  # make sure that this Writer is
        kwargs.pop('universe', None)  # always linked to this reader, don't bother with Universe
        return PDBWriter(filename, **kwargs)


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
       order to write a selection, use the :class:`PrimitivePDBWriter` ,
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
        if self.PDBstructure is not None \
           and not isinstance(self.PDBstructure, Bio.PDB.Structure.Structure):
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
            ppw.write(self.universe.select_atoms('all'))
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
    """PDBReader that reads a `PDB-formatted`_ file, no frills.

    The following *PDB records* are parsed (see `PDB coordinate section`_ for details):

     - *CRYST1* for unitcell A,B,C, alpha,beta,gamma
     - *ATOM* or *HETATM* for serial,name,resName,chainID,resSeq,x,y,z,occupancy,tempFactor
     - *HEADER* (:attr:`header`), *TITLE* (:attr:`title`), *COMPND*
       (:attr:`compound`), *REMARK* (:attr:`remarks`)
     - all other lines are ignored

    Reads multi-`MODEL`_ PDB files as trajectories.

    .. _PDB-formatted:
       http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html
    .. _PDB coordinate section:
       http://www.wwpdb.org/documentation/format32/sect9.html
    .. _MODEL:
       http://www.wwpdb.org/documentation/format32/sect9.html#MODEL

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


    .. SeeAlso:: :class:`PrimitivePDBWriter`; :class:`PDBReader`
                 implements a larger subset of the header records,
                 which are accessible as :attr:`PDBReader.metadata`.

    .. versionchanged:: 0.11.0
       * Frames now 0-based instead of 1-based
       * New :attr:`title` (list with all TITLE lines).

    """
    format = 'PDB'
    units = {'time': None, 'length': 'Angstrom'}

    def __init__(self, filename, **kwargs):
        """Read coordinates from *filename*.

        *filename* can be a gzipped or bzip2ed compressed PDB file.

        If the pdb file contains multiple MODEL records then it is
        read as a trajectory where the MODEL numbers correspond to
        frame numbers. Therefore, the MODEL numbers must be a sequence
        of integers (typically starting at 1 or 0).
        """
        super(PrimitivePDBReader, self).__init__(filename, **kwargs)

        try:
            self._n_atoms = kwargs['n_atoms']
        except KeyError:
            raise ValueError("PrimitivePDBReader requires the n_atoms keyword")

        self.model_offset = kwargs.pop("model_offset", 0)

        header = ""
        title = []
        compound = []
        remarks = []

        frames = {}

        self.ts = self._Timestep(self._n_atoms, **self._ts_kwargs)

        pos = 0  # atom position for filling coordinates array
        occupancy = np.ones(self._n_atoms)
        with util.openany(filename, 'r') as pdbfile:
            for i, line in enumerate(pdbfile):
                line = line.strip()  # Remove extra spaces
                if len(line) == 0:  # Skip line if empty
                    continue
                record = line[:6].strip()

                if record == 'END':
                    break
                elif record == 'CRYST1':
                    A, B, C = map(float, [line[6:15], line[15:24], line[24:33]])
                    alpha, beta, gamma = map(float, [line[33:40], line[40:47], line[47:54]])
                    self.ts._unitcell[:] = A, B, C, alpha, beta, gamma
                    continue
                elif record == 'HEADER':
                    # classification = line[10:50]
                    # date = line[50:59]
                    # idCode = line[62:66]
                    header = line[10:66]
                    continue
                elif record == 'TITLE':
                    l = line[8:80].strip()
                    title.append(l)
                    continue
                elif record == 'COMPND':
                    l = line[7:80].strip()
                    compound.append(l)
                    continue
                elif record == 'REMARK':
                    content = line[6:].strip()
                    remarks.append(content)
                elif record == 'MODEL':
                    frames[len(frames)] = i  # 0-based indexing
                elif line[:6] in ('ATOM  ', 'HETATM'):
                    # skip atom/hetatm for frames other than the first
                    # they will be read in when next() is called
                    # on the trajectory reader
                    if len(frames) > 1:
                        continue
                    self.ts._pos[pos] = map(float, [line[30:38], line[38:46], line[46:54]])
                    try:
                        occupancy[pos] = float(line[54:60])
                    except ValueError:
                        pass
                    pos += 1

        self.header = header
        self.title = title
        self.compound = compound
        self.remarks = remarks

        if pos != self._n_atoms:
            raise ValueError("Read an incorrect number of atoms\n"
                             "Expected {expected} got {actual}"
                             "".format(expected=self._n_atoms, actual=pos))
        self.n_atoms = pos

        self.ts.frame = 0  # 0-based frame number as starting frame
        self.ts.data['occupancy'] = occupancy

        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)  # in-place !
            self.convert_pos_from_native(self.ts._unitcell[:3])  # in-place ! (only lengths)

        # No 'MODEL' entries
        if len(frames) == 0:
            frames[0] = 0

        self.frames = frames
        self.n_frames = len(frames) if frames else 1

    def get_occupancy(self):
        """Return an array of occupancies in atom order."""
        return np.array(self._occupancy)

    def Writer(self, filename, **kwargs):
        """Returns a permissive (simple) PDBWriter for *filename*.

        :Arguments:
          *filename*
              filename of the output PDB file

        :Returns: :class:`PrimitivePDBWriter`

        """
        kwargs.setdefault('multiframe', self.n_frames > 1)
        return PrimitivePDBWriter(filename, **kwargs)

    def rewind(self):
        self._read_frame(0)

    def _reopen(self):
        # Pretend the current TS is -1 (in 0 based) so "next" is the
        # 0th frame
        self.ts.frame = -1

    def _read_next_timestep(self, ts=None):
        if ts is None:
            ts = self.ts
        else:
            # TODO: cleanup _read_frame() to use a "free" Timestep
            raise NotImplementedError("PrimitivePDBReader cannot assign to a timestep")
        # frame is 1-based. Normally would add 1 to frame before calling
        # self._read_frame to retrieve the subsequent ts. But self._read_frame
        # assumes it is being passed a 0-based frame, and adjusts.
        frame = self.frame + 1
        return self._read_frame(frame)

    def _read_frame(self, frame):
        try:
            line = self.frames[frame]
        except KeyError:
            raise IOError
        if line is None:
            # single frame file, we already have the timestep
            return self.ts

        # TODO: only open file once and leave the file open; then seek back and
        #       forth; should improve performance substantially
        pos = 0
        occupancy = np.ones(self._n_atoms)
        with util.openany(self.filename, 'r') as f:
            for i in xrange(line):
                f.next()  # forward to frame
            for line in f:
                if line[:6] == 'ENDMDL':
                    break
                # NOTE - CRYST1 line won't be found if it comes before the MODEL
                # line, which is sometimes the case, e.g. output from gromacs
                # trjconv
                elif line[:6] == 'CRYST1':
                    A, B, C = map(float, [line[6:15], line[15:24], line[24:33]])
                    alpha, beta, gamma = map(float, [line[33:40], line[40:47],
                                                     line[47:54]])
                    self.ts._unitcell[:] = A, B, C, alpha, beta, gamma
                    continue
                elif line[:6] in ('ATOM  ', 'HETATM'):
                    # we only care about coordinates
                    self.ts._pos[pos] = map(float, [line[30:38], line[38:46],
                                                    line[46:54]])
                    # TODO import bfactors - might these change?
                    try:
                        occupancy[pos] = float(line[54:60])
                    except:
                        pass
                    pos += 1
                    continue

        # check if atom number changed
        if pos != self._n_atoms:
            raise ValueError("Read an incorrect number of atoms\n"
                             "Expected {expected} got {actual}"
                             "".format(expected=self._n_atoms, actual=pos+1))

        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)  # in-place !
            self.convert_pos_from_native(self.ts._unitcell[:3])  # in-place ! (only lengths)
        self.ts.frame = frame
        self.ts.data['occupancy'] = occupancy
        return self.ts


class PrimitivePDBWriter(base.Writer):
    """PDB writer that implements a subset of the `PDB 3.2 standard`_ .

    PDB format as used by NAMD/CHARMM: 4-letter resnames and segID are allowed,
    altLoc is written.

    The :class:`PrimitivePDBWriter` can be used to either a dump a coordinate
    set to a PDB file (operating as a "single frame writer", selected with the
    constructor keyword *multiframe* = ``False``, the default) or by writing a
    PDB "movie" (multi frame mode, *multiframe* = ``True``), consisting of
    multiple models (using the MODEL_ and ENDMDL_ records).

    .. _`PDB 3.2 standard`:
       http://www.wwpdb.org/documentation/format32/v3.2.html
    .. _MODEL: http://www.wwpdb.org/documentation/format32/sect9.html#MODEL
    .. _ENDMDL: http://www.wwpdb.org/documentation/format32/sect9.html#ENDMDL
    .. _CONECT: http://www.wwpdb.org/documentation/format32/sect10.html#CONECT

    .. SeeAlso::
       This class is identical to :class:`MultiPDBWriter` with the one
       exception that it defaults to writing single-frame PDB files as if
       *multiframe* = ``False`` was selected.

    .. versionchanged:: 0.7.5
       Initial support for multi-frame PDB files.

    .. versionchanged:: 0.7.6
       The *multiframe* keyword was added to select the writing mode. The writing
       of bond information (CONECT_ records) is now disabled by default but can be
       enabled with the *bonds* keyword.

    .. versionchanged:: 0.11.0
       Frames now 0-based instead of 1-based
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
        'ATOM': "ATOM  %(serial)5d %(name)-4s%(altLoc)-1s%(resName)-4s%(chainID)1s%(resSeq)4d%(iCode)1s   %(x)8.3f%("
                "y)8.3f%(z)8.3f%(occupancy)6.2f%(tempFactor)6.2f      %(segID)-4s%(element)2s%(charge)2d\n",
        'REMARK': "REMARK     %s\n",
        'COMPND': "COMPND    %s\n",
        'HEADER': "HEADER    %s\n",
        'TITLE': "TITLE     %s\n",
        'MODEL': "MODEL     %5d\n",
        'NUMMDL': "NUMMDL    %5d\n",
        'ENDMDL': "ENDMDL\n",
        'END': "END\n",
        'CRYST1': "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
        'CONECT': "CONECT%s\n"
    }
    format = 'PDB'
    units = {'time': None, 'length': 'Angstrom'}
    pdb_coor_limits = {"min": -999.9995, "max": 9999.9995}
    #: wrap comments into REMARK records that are not longer than
    # :attr:`remark_max_length` characters.
    remark_max_length = 66
    _multiframe = False

    def __init__(self, filename, bonds="conect", n_atoms=None, start=0, step=1,
                 remarks="Created by PrimitivePDBWriter",
                 convert_units=None, multiframe=None):
        """Create a new PDBWriter

        :Arguments:
         *filename*
           name of output file
         *start*
           starting timestep
         *step*
           skip between subsequent timesteps
         *remarks*
           comments to annotate pdb file (added to the TITLE record); note that
           any remarks from the trajectory that serves as input are
           written to REMARK records with lines longer than :attr:`remark_max_length` (66
           characters) being wrapped.
         *convert_units*
           units are converted to the MDAnalysis base format; ``None`` selects
           the value of :data:`MDAnalysis.core.flags` ['convert_lengths']
         *bonds*
           write bonds to the PDB file as CONECT_ records [``False``]

           .. Note::

              Currently only works when writing a whole :class:`Universe` and
              if bond information is available in the topology. (For selections
              smaller than the whole :class:`Universe`, the atom numbering in
              the CONECT_ records would not match the numbering of the atoms in
              the new PDB file and therefore a :exc:`NotImplementedError` is
              raised.)

         *multiframe*
           ``False``: write a single frame to the file; ``True``: create a
           multi frame PDB file in which frames are written as MODEL_ ... ENDMDL_
           records. If ``None``, then the class default is chosen.    [``None``]

        .. _CONECT: http://www.wwpdb.org/documentation/format32/sect10.html#CONECT
        .. _MODEL: http://www.wwpdb.org/documentation/format32/sect9.html#MODEL
        .. _ENDMDL: http://www.wwpdb.org/documentation/format32/sect9.html#ENDMDL

        """
        # n_atoms = None : dummy keyword argument
        # (not used, but Writer() always provides n_atoms as the second argument)

        # TODO: - remarks should be a list of lines and written to REMARK
        #       - additional title keyword could contain line for TITLE

        self.filename = filename
        if convert_units is None:
            convert_units = flags['convert_lengths']
        self.convert_units = convert_units  # convert length and time to base units
        self.multiframe = self._multiframe if multiframe is None else multiframe
        self.bonds = bonds

        self.frames_written = 0
        if start < 0:
            raise ValueError("'Start' must be a positive value")

        self.start = start
        self.step = step
        self.remarks = remarks

        self.pdbfile = util.anyopen(self.filename, 'w')  # open file on init
        self.has_END = False

    def close(self):
        """Close PDB file and write END record"""
        if hasattr(self, 'pdbfile') and self.pdbfile is not None:
            if not self.has_END:
                self.END()
            else:
                logger.warn("END record has already been written before the final closing of the file")
            self.pdbfile.close()
        self.pdbfile = None

    def _write_pdb_title(self):
        if self.multiframe:
            self.TITLE("MDANALYSIS FRAMES FROM %d, SKIP %d: %s" % (self.start, self.step, self.remarks))
        else:
            self.TITLE("MDANALYSIS FRAME %d: %s" % (self.start, self.remarks))

    def _write_pdb_header(self):
        if not self.obj or not hasattr(self.obj, 'universe'):
            self._write_pdb_title(self)
            return

        u = self.obj.universe
        self.HEADER(u.trajectory)
        self._write_pdb_title()
        self.COMPND(u.trajectory)
        try:
            # currently inconsistent: DCDReader gives a string, PDB*Reader a list, so always get a list
            # split long single lines into chunks of 66 chars
            remarks = []
            for line in util.asiterable(u.trajectory.remarks):
                remarks.extend(textwrap.wrap(line, self.remark_max_length))
            self.REMARK(*remarks)
        except AttributeError:
            pass
        self.CRYST1(self.convert_dimensions_to_unitcell(u.trajectory.ts))

    def _check_pdb_coordinates(self):
        """Check if the coordinate values fall within the range allowed for PDB files.

        Deletes the output file if this is the first frame or if frames have
        already been written (in multi-frame mode) adds a REMARK instead of the
        coordinates and closes the file.

        Raises :exc:`ValueError` if the coordinates fail the check.
        """
        atoms = self.obj.atoms  # make sure to use atoms (Issue 46)
        coor = atoms.coordinates()  # can write from selection == Universe (Issue 49)

        # check if any coordinates are illegal (coordinates are already in Angstroem per package default)
        if self.has_valid_coordinates(self.pdb_coor_limits, coor):
            return True
        # note the precarious close() here: we know that the file is open and we now
        # prepare to remove what we have already written (header and such) or add a REMARK
        # (which allows the user to look at the previously written frames)
        if self.frames_written > 1:
            self.REMARK("Incomplete multi-frame trajectory.",
                        "Coordinates for the current frame cannot be represented in the PDB format.")
            self.close()
        else:
            self.close()
            try:
                os.remove(self.filename)
            except OSError as err:
                if err.errno == errno.ENOENT:
                    pass
        raise ValueError(
            "PDB files must have coordinate values between %.3f and %.3f Angstroem: file writing was aborted." %
            (self.pdb_coor_limits["min"], self.pdb_coor_limits["max"]))

    def _write_pdb_bonds(self):
        """Writes out all the bond records; works only for Universe objects.

        .. Warning::

           All bonds are written out, using the old atom numbers - this is
           incorrect. Once a selection is made, the atom numbers have to be
           updated (currently they are unmodified) and bonds have to be
           selected for, only if all the atoms for a bond are still present.

           Therefore, this method raises a :exc:`NotImplementedError` if CONECT
           records for anything smaller than the :class:`Universe` are written.

        .. versionchanged:: 0.7.6
           Only write CONECT records if :attr:`PrimitivePDBWriter.bonds` ``== True``.
           Raises :exc:`NotImplementedError` if it would produce wrong output.

        """
        if not self.bonds:
            return

        if not self.obj or not hasattr(self.obj, 'universe') or not hasattr(self.obj.universe, 'bonds'):
            return

        if self.obj.atoms.n_atoms != self.obj.universe.atoms.n_atoms:
            pass
            #logger.error("PDB CONECT records not written because this only works correctly for a whole Universe.")
            #raise NotImplementedError("PDB CONECT records not written because this only works correctly for a whole
            # Universe.")

        bonds = set()

        [[bonds.add(b) for b in a.bonds] for a in self.obj.atoms]

        atoms = set([a.index for a in self.obj.atoms])

        mapping = dict([(atom.index, i) for i, atom in enumerate(self.obj.atoms)])

        # Write out only the bonds that were defined in CONECT records
        if self.bonds == "conect":
            bonds = [(bond[0].index, bond[1].index) for bond in bonds if not bond.is_guessed]
        elif self.bonds == "all":
            bonds = [(bond[0].index, bond[1].index) for bond in bonds]
        else:
            raise ValueError("bonds has to be either None, 'conect' or 'all'")
        con = {}

        for a1, a2 in bonds:
            if not (a1 in atoms and a2 in atoms):
                continue
            if not a1 in con:
                con[a1] = []
            if not a2 in con:
                con[a2] = []
            con[a2].append(a1)
            con[a1].append(a2)

        #print con
        atoms = sorted([a.index for a in self.obj.atoms])

        conect = [([a, ] + sorted(con[a])) for a in atoms if a in con]

        conect = [[mapping[e] for e in row] for row in conect]

        for c in conect:
            self.CONECT(c)

    def _update_frame(self, obj):
        """Method to initialize important attributes in writer from a AtomGroup or Universe *obj*.

        Attributes initialized/updated:

        * :attr:`PrimitivePDBWriter.obj` (the entity that provides topology information *and*
          coordinates, either a :class:`~MDAnalysis.core.AtomGroup.AtomGroup` or a whole
          :class:`~MDAnalysis.core.AtomGroup.Universe`)
        * :attr:`PrimitivePDBWriter.trajectory` (the underlying trajectory
          :class:`~MDAnalysis.coordinates.base.Reader`)
        * :attr:`PrimitivePDBWriter.timestep` (the underlying trajectory
          :class:`~MDAnalysis.coordinates.base.Timestep`)

        Before calling :meth:`write_next_timestep` this method **must** be
        called at least once to enable extracting topology information from the
        current frame.
        """

        if isinstance(obj, base.Timestep):
            raise TypeError(
                "PrimitivePDBWriter cannot write Timestep objects directly, since they lack topology information ("
                "atom names and types) required in PDB files")

        self.obj = obj  # remember obj for some of other methods  --- NOTE: this is an evil/lazy hack...
        ts, traj = None, None
        if hasattr(obj, 'universe') and not isinstance(obj, Universe):
            # For AtomGroup and children (Residue, ResidueGroup, Segment)
            ts = obj.universe.trajectory.ts
            traj = obj.universe.trajectory
        else:
            # For Universe only
            ts = obj.trajectory.ts
            traj = obj.trajectory

        if not (ts and traj):
            raise AssertionError(
                "PrimitivePDBWriter couldn't extract trajectory and timestep information from an object; inheritance "
                "problem.")

        self.trajectory = traj  # update trajectory (used by other methods)
        self.ts = ts  # update timestep (used by other methods)
        return traj, ts

    def write(self, obj):
        """Write object *obj* at current trajectory frame to file.

        *obj* can be a selection (i.e. a
        :class:`~MDAnalysis.core.AtomGroup.AtomGroup`) or a whole
        :class:`~MDAnalysis.core.AtomGroup.Universe`.

        The last letter of the :attr:`~MDAnalysis.core.AtomGroup.Atom.segid` is
        used as the PDB chainID (but see :meth:`~PrimitivePDBWriter.ATOM` for
        details).

        :Arguments:
          *obj*
            :class:`~MDAnalysis.core.AtomGroup.AtomGroup` or
            :class:`~MDAnalysis.core.AtomGroup.Universe`
        """

        self._update_frame(obj)
        self._write_pdb_header()
        # Issue 105: with write() ONLY write a single frame; use write_all_timesteps() to dump
        # everything in one go, or do the traditional loop over frames
        self.write_next_timestep(self.ts, multiframe=self.multiframe)
        self._write_pdb_bonds()
        # END record is written when file is being close()d

    def write_all_timesteps(self, obj):
        """Write all timesteps associated with *obj* to the PDB file.

        *obj* can be a :class:`~MDAnalysis.core.AtomGroup.AtomGroup`
        or a :class:`~MDAnalysis.core.AtomGroup.Universe`.

        The method writes the frames from the one specified as *start* until
        the end, using a step of *skip* (*start* and *skip* are set in the
        constructor). Thus, if *u* is a Universe then ::

           u.trajectory[-2]
           pdb = PrimitivePDBWriter("out.pdb", u.atoms.n_atoms)
           pdb.write_all_timesteps(u)

        will write a PDB trajectory containing the last 2 frames and ::

           pdb = PrimitivePDBWriter("out.pdb", u.atoms.n_atoms, start=12, skip=2)
           pdb.write_all_timesteps(u)

        will be writing frames 12, 14, 16, ...

        .. versionchanged:: 0.11.0
           Frames now 0-based instead of 1-based
        """

        self._update_frame(obj)
        self._write_pdb_header()

        start, step = self.start, self.step
        traj = self.trajectory

        # Start from trajectory[0]/frame 0, if there are more than 1 frame.
        # If there is only 1 frame, the traj.frames is not like a python list:
        # accessing trajectory[-1] raises key error.
        if not start and traj.n_frames > 1:
            start = traj.frame

        for framenumber in xrange(start, len(traj), step):
            traj[framenumber]
            self.write_next_timestep(self.ts, multiframe=True)

        self._write_pdb_bonds()
        self.close()

        # Set the trajectory to the starting position
        traj[start]

    def write_next_timestep(self, ts=None, **kwargs):
        '''write a new timestep to the PDB file

        :Keywords:
          *ts*
             :class:`base.Timestep` object containing coordinates to be written to trajectory file;
             if ``None`` then :attr:`PrimitivePDBWriter.ts`` is tried.
          *multiframe*
             ``False``: write a single frame (default); ``True`` behave as a trajectory writer

        .. Note::

           Before using this method with another :class:`base.Timestep` in the *ts*
           argument, :meth:`PrimitivePDBWriter._update_frame` *must* be called
           with the :class:`~MDAnalysis.core.AtomGroup.AtomGroup.Universe` as
           its argument so that topology information can be gathered.
        '''
        if ts is None:
            if not hasattr(self, "ts"):
                raise NoDataError("PBDWriter: no coordinate data to write to trajectory file")
            else:
                ts = self.ts
        self._check_pdb_coordinates()
        self._write_timestep(ts, **kwargs)

    def _write_timestep(self, ts, multiframe=False):
        """Write a new timestep *ts* to file

        Does unit conversion if necessary.

        By setting *multiframe* = ``True``, MODEL_ ... ENDMDL_ records are
        written to represent trajectory frames in a multi-model PDB file. (At
        the moment we do *not* write the NUMMDL_ record.)

        The *multiframe* = ``False`` keyword signals that the
        :class:`PrimitivePDBWriter` is in single frame mode and no MODEL_
        records are written.

        .. _MODEL: http://www.wwpdb.org/documentation/format32/sect9.html#MODEL
        .. _ENDMDL: http://www.wwpdb.org/documentation/format32/sect9.html#ENDMDL
        .. _NUMMDL: http://www.wwpdb.org/documentation/format32/sect2.html#NUMMDL

        .. versionchanged:: 0.7.6
           The *multiframe* keyword was added, which completely determines if
           MODEL records are written. (Previously, this was decided based on the
           underlying trajectory and only if ``len(traj) > 1`` would MODEL records
           have been written.)

        """

        traj = self.trajectory
        atoms = self.obj.atoms
        if self.convert_units:
            coor = self.convert_pos_to_native(ts._pos, inplace=False)
        else:
            coor = ts._pos

        if hasattr(self.obj, "indices"):
            coor = coor[self.obj.indices]

        if len(atoms) != len(coor):
            raise ValueError(
                "Length of the atoms array is %d, this is different form the Timestep coordinate array %d" % (
                    len(atoms), len(ts._pos)))

        if multiframe:
            self.MODEL(self.frames_written + 1)

        for i, atom in enumerate(atoms):
            # TODO Jan: see description in ATOM for why this check has to be made
            if not atom.bfactor:
                atom.bfactor = 0.
            self.ATOM(serial=i + 1, name=atom.name.strip(), resName=atom.resname.strip(),
                      resSeq=atom.resid, chainID=atom.segid.strip(), segID=atom.segid.strip(),
                      tempFactor=atom.bfactor, altLoc=atom.altLoc,
                      x=coor[i, 0], y=coor[i, 1], z=coor[i, 2])
            # get bfactor, too, and add to output?
            # 'element' is auto-guessed from atom.name in ATOM()
        if multiframe:
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
        line = " ".join(title)  # TODO: should do continuation automatically
        self.pdbfile.write(self.fmt['TITLE'] % line)

    def REMARK(self, *remarks):
        """Write generic REMARK_ record (without number).

        Each string provided in *remarks* is written as a separate REMARK_
        record.

        See also `REMARK (update)`_.

        .. _REMARK: http://www.wwpdb.org/documentation/format32/remarks1.html
        .. _REMARK (update): http://www.wwpdb.org/documentation/format32/remarks2.html

        """
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

        Only a single END record is written. Calling END multiple times has no
        effect. Because :meth:`~PrimitivePDBWriter.close` also calls this
        method right before closing the file it is recommended to *not* call
        :meth:`~PrimitivePDBWriter.END` explicitly.

        .. _END: http://www.wwpdb.org/documentation/format32/sect11.html#END

        """
        if not self.has_END:
            self.pdbfile.write(self.fmt['END'])  # only write a single END record
        self.has_END = True

    def ENDMDL(self):
        """Write the ENDMDL_ record.

        .. _ENDMDL: http://www.wwpdb.org/documentation/format32/sect9.html#ENDMDL

        """
        self.pdbfile.write(self.fmt['ENDMDL'])

    def ATOM(self, serial=None, name=None, altLoc=None, resName=None, chainID=None,
             resSeq=None, iCode=None, x=None, y=None, z=None, occupancy=1.0, tempFactor=0.0,
             segID=None, element=None, charge=0):
        """Write ATOM_ record.

        Only some keword args are optional (*altLoc*, *iCode*, *chainID*), for
        some defaults are set.

        All inputs are cut to the maximum allowed length. For integer numbers
        the highest-value digits are chopped (so that the *serial* and *reSeq*
        wrap); for strings the trailing characters are chopped. The *last*
        character of *chainID* becomes the PDB *chainID* (unless it has the
        special value "SYSTEM" (assigned by MDAnalysis if neither *segID* nor
        *chainID* were available), in which case the PDB will have an empty
        *chainID*).

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

        for arg in (
                'serial', 'name', 'resName', 'resSeq', 'x', 'y', 'z',
                'occupancy', 'tempFactor', 'charge'):
            if locals()[arg] is None:
                raise ValueError('parameter ' + arg + ' must be defined.')
        serial = int(str(serial)[-5:])  # check for overflow here?
        name = name[:4]
        if len(name) < 4:
            name = " " + name  # customary to start in column 14
        altLoc = altLoc or " "
        altLoc = altLoc[:1]
        resName = resName[:4]
        chainID = chainID or ""  # or should we provide a chainID such as 'A'?
        chainID = chainID if not chainID.upper() == "SYSTEM" else ""  # special case, new in 0.7.6
        chainID = chainID.strip()[-1:]  # take the last character
        resSeq = int(str(resSeq)[-4:])  # check for overflow here?
        iCode = iCode or ""
        iCode = iCode[:1]
        element = element or guess_atom_element(name.strip())  # element == 0|False|None will be guessed
        element = str(element).strip()[:2]  # make sure that is a string for user input
        segID = segID or chainID
        segID = segID[:4]
        self.pdbfile.write(self.fmt['ATOM'] % vars())

    def CONECT(self, conect):
        """Write CONECT_ record.

        .. _CONECT: http://www.wwpdb.org/documentation/format32/sect10.html#CONECT

        """
        conect = ["%5d" % (entry + 1) for entry in conect]
        conect = "".join(conect)
        self.pdbfile.write(self.fmt['CONECT'] % conect)


class ExtendedPDBReader(PrimitivePDBReader):
    """PDBReader that reads a PDB-formatted file with five-digit residue numbers.

    This reader does not conform to the `PDB standard`_ because it allows
    five-digit residue numbers that may take up columns 23 to 27 (inclusive)
    instead of being confined to 23-26 (with column 27 being reserved for the
    insertion code in the PDB standard). PDB files in this format are written
    by popular programs such as VMD_.

    .. SeeAlso:: :class:`PrimitivePDBReader`

    .. _PDB standard: http://www.wwpdb.org/documentation/format32/sect9.html
    .. _VMD: http://www.ks.uiuc.edu/Research/vmd/

    .. versionadded:: 0.8
    """
    format = "XPDB"


class MultiPDBWriter(PrimitivePDBWriter):
    """PDB writer that implements a subset of the `PDB 3.2 standard`_ .

    PDB format as used by NAMD/CHARMM: 4-letter resnames and segID, altLoc
    is written.

    By default, :class:`MultiPDBWriter` writes a PDB "movie" (multi frame mode,
    *multiframe* = ``True``), consisting of multiple models (using the MODEL_
    and ENDMDL_ records).

    .. _`PDB 3.2 standard`:
       http://www.wwpdb.org/documentation/format32/v3.2.html

    .. _MODEL: http://www.wwpdb.org/documentation/format32/sect9.html#MODEL
    .. _ENDMDL: http://www.wwpdb.org/documentation/format32/sect9.html#ENDMDL
    .. _CONECT: http://www.wwpdb.org/documentation/format32/sect10.html#CONECT


    .. SeeAlso::
       This class is identical to :class:`PrimitivePDBWriter` with the one
       exception that it defaults to writing multi-frame PDB files instead of
       single frames.

    .. versionadded:: 0.7.6

    """
    _multiframe = True
