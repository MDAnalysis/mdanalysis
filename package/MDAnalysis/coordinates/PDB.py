# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""PDB structure files in MDAnalysis --- :mod:`MDAnalysis.coordinates.PDB`
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

PDB I/O is available in the form of the Simple PDB Reader/Writers.

..deprecated:: 0.15.0
Readers and writers solely available in the form of
Simple Readers and Writers, see below.

Simple PDB Reader and Writer
-----------------------------------------
A pure-Python implementation for PDB files commonly encountered in MD
simulations comes under the names :class:`PDBReader` and
:class:`PDBWriter`. It only implements a subset of the `PDB standard`_
(for instance, it does not deal with insertion codes) and also allows some
typical enhancements such as 4-letter resids (introduced by CHARMM/NAMD).

The :class:`PDBReader` can read multi-frame PDB files and represents
them as a trajectory. The :class:`PDBWriter` can write single and
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

.. autoclass:: PDBReader
   :members:

.. autoclass:: PDBWriter
   :members:

   .. automethod:: _check_pdb_coordinates
   .. automethod:: _write_pdb_bonds
   .. automethod:: _update_frame
   .. automethod:: _write_timestep

.. autoclass:: MultiPDBWriter
   :members:


..deprecated:: 0.15.0
    The "permissive" flag is not used anymore (and effectively defaults to True);
    it will be completely removed in 0.16.0.

"""

from six.moves import range, zip

import os
import errno
import textwrap
import warnings
import logging
import collections
import numpy as np

from ..core import flags
from ..lib import util
from . import base
from ..topology.core import guess_atom_element
from ..core.AtomGroup import Universe
from ..exceptions import NoDataError


logger = logging.getLogger("MDAnalysis.coordinates.PBD")

# Pairs of residue name / atom name in use to deduce PDB formatted atom names
Pair = collections.namedtuple('Atom', 'resname name')

class PDBReader(base.Reader):
    """PDBReader that reads a `PDB-formatted`_ file, no frills.

    The following *PDB records* are parsed (see `PDB coordinate section`_ for
    details):

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


    .. SeeAlso:: :class:`PDBWriter`; :class:`PDBReader`
                 implements a larger subset of the header records,
                 which are accessible as :attr:`PDBReader.metadata`.

    .. versionchanged:: 0.11.0
       * Frames now 0-based instead of 1-based
       * New :attr:`title` (list with all TITLE lines).

    """
    format = ['PDB', 'ENT']
    units = {'time': None, 'length': 'Angstrom'}

    def __init__(self, filename, **kwargs):
        """Read coordinates from *filename*.

        *filename* can be a gzipped or bzip2ed compressed PDB file.

        If the pdb file contains multiple MODEL records then it is
        read as a trajectory where the MODEL numbers correspond to
        frame numbers.
        """
        super(PDBReader, self).__init__(filename, **kwargs)

        try:
            self.n_atoms = kwargs['n_atoms']
        except KeyError:
            # hackish, but should work and keeps things DRY
            # regular MDA usage via Universe doesn't follow this route
            from MDAnalysis.topology import PDBParser

            with PDBParser.PDBParser(self.filename) as p:
                top = p.parse()
            self.n_atoms = len(top['atoms'])

        self.model_offset = kwargs.pop("model_offset", 0)

        self.header = header = ""
        self.title = title = []
        self.compound = compound = []
        self.remarks = remarks = []

        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)

        # Record positions in file of CRYST and MODEL headers
        # then build frame offsets to start at the minimum of these
        # This allows CRYST to come either before or after MODEL
        # This assumes that **either**
        # - pdbfile has a single CRYST (NVT)
        # - pdbfile has a CRYST for every MODEL (NPT)
        models = []
        crysts = []

        pdbfile = self._pdbfile = util.anyopen(filename, 'rt')

        line = "magical"
        while line:
            # need to use readline so tell gives end of line
            # (rather than end of current chunk)
            line = pdbfile.readline()

            if line.startswith('MODEL'):
                models.append(pdbfile.tell())
            elif line.startswith('CRYST1'):
                # remove size of line to get **start** of CRYST line
                crysts.append(pdbfile.tell() - len(line))
            elif line.startswith('HEADER'):
                # classification = line[10:50]
                # date = line[50:59]
                # idCode = line[62:66]
                header = line[10:66]
            elif line.startswith('TITLE'):
                title.append(line[8:80].strip())
            elif line.startswith('COMPND'):
                compound.append(line[7:80].strip())
            elif line.startswith('REMARK'):
                remarks.append(line[6:].strip())

        end = pdbfile.tell()  # where the file ends

        if not models:
            # No model entries
            # so read from start of file to read first frame
            models.append(0)
        if len(crysts) == len(models):
            offsets = [min(a, b) for a, b in zip(models, crysts)]
        else:
            offsets = models
        # Position of the start of each frame
        self._start_offsets = offsets
        # Position of the end of each frame
        self._stop_offsets = offsets[1:] + [end]
        self.n_frames = len(offsets)

        self._read_frame(0)

    def Writer(self, filename, **kwargs):
        """Returns a PDBWriter for *filename*.

        :Arguments:
          *filename*
              filename of the output PDB file

        :Returns: :class:`PDBWriter`

        """
        kwargs.setdefault('multiframe', self.n_frames > 1)
        return PDBWriter(filename, **kwargs)

    def rewind(self):
        self._read_frame(0)

    def _reopen(self):
        # Pretend the current TS is -1 (in 0 based) so "next" is the
        # 0th frame
        self.close()
        self._pdbfile = util.anyopen(self.filename, 'rt')
        self.ts.frame = -1

    def _read_next_timestep(self, ts=None):
        if ts is None:
            ts = self.ts
        else:
            # TODO: cleanup _read_frame() to use a "free" Timestep
            raise NotImplementedError("PDBReader cannot assign to a timestep")
        # frame is 1-based. Normally would add 1 to frame before calling
        # self._read_frame to retrieve the subsequent ts. But self._read_frame
        # assumes it is being passed a 0-based frame, and adjusts.
        frame = self.frame + 1
        return self._read_frame(frame)

    def _read_frame(self, frame):
        try:
            start = self._start_offsets[frame]
            stop = self._stop_offsets[frame]
        except IndexError:  # out of range of known frames
            raise IOError

        pos = 0
        occupancy = np.ones(self.n_atoms)

        # Seek to start and read until start of next frame
        self._pdbfile.seek(start)
        chunk = self._pdbfile.read(stop - start)

        for line in chunk.splitlines():
            if line[:6] in ('ATOM  ', 'HETATM'):
                # we only care about coordinates
                self.ts._pos[pos] = [line[30:38],
                                     line[38:46],
                                     line[46:54]]
                # TODO import bfactors - might these change?
                try:
                    occupancy[pos] = line[54:60]
                except ValueError:
                    # Be tolerant for ill-formated or empty occupancies
                    pass
                pos += 1
            elif line[:6] == 'CRYST1':
                self.ts._unitcell[:] = [line[6:15], line[15:24],
                                        line[24:33], line[33:40],
                                        line[40:47], line[47:54]]

        # check if atom number changed
        if pos != self.n_atoms:
            raise ValueError("Read an incorrect number of atoms\n"
                             "Expected {expected} got {actual}"
                             "".format(expected=self.n_atoms, actual=pos+1))

        if self.convert_units:
            # both happen inplace
            self.convert_pos_from_native(self.ts._pos)
            self.convert_pos_from_native(self.ts._unitcell[:3])
        self.ts.frame = frame
        self.ts.data['occupancy'] = occupancy
        return self.ts

    def close(self):
        self._pdbfile.close()


class PDBWriter(base.Writer):
    """PDB writer that implements a subset of the `PDB 3.2 standard`_ .

    PDB format as used by NAMD/CHARMM: 4-letter resnames and segID are allowed,
    altLoc is written.

    The :class:`PDBWriter` can be used to either dump a coordinate
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
       The *multiframe* keyword was added to select the writing mode. The
       writing of bond information (CONECT_ records) is now disabled by default
       but can be enabled with the *bonds* keyword.

    .. versionchanged:: 0.11.0
       Frames now 0-based instead of 1-based

    .. versionchanged:: 0.14.0
       PDB doesn't save charge information

    """
    fmt = {
        'ATOM': (
            "ATOM  {serial:5d} {name:<4s}{altLoc:<1s}{resName:<4s}"
            "{chainID:1s}{resSeq:4d}{iCode:1s}"
            "   {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}{occupancy:6.2f}"
            "{tempFactor:6.2f}      {segID:<4s}{element:>2s}\n"),
        'REMARK': "REMARK     {0}\n",
        'COMPND': "COMPND    {0}\n",
        'HEADER': "HEADER    {0}\n",
        'TITLE': "TITLE     {0}\n",
        'MODEL': "MODEL     {0:5d}\n",
        'NUMMDL': "NUMMDL    {0:5d}\n",
        'ENDMDL': "ENDMDL\n",
        'END': "END\n",
        'CRYST1': ("CRYST1{box[0]:9.3f}{box[1]:9.3f}{box[2]:9.3f}"
                   "{ang[0]:7.2f}{ang[1]:7.2f}{ang[2]:7.2f} "
                   "{spacegroup:<11s}{zvalue:4d}\n"),
        'CONECT': "CONECT{0}\n"
    }
    format = ['PDB', 'ENT']
    units = {'time': None, 'length': 'Angstrom'}
    pdb_coor_limits = {"min": -999.9995, "max": 9999.9995}
    #: wrap comments into REMARK records that are not longer than
    # :attr:`remark_max_length` characters.
    remark_max_length = 66
    multiframe = False

    # These attributes are used to deduce how to format the atom name.
    ions = ('FE', 'AS', 'ZN', 'MG', 'MN', 'CO', 'BR',
            'CU', 'TA', 'MO', 'AL', 'BE', 'SE', 'PT',
            'EU', 'NI', 'IR', 'RH', 'AU', 'GD', 'RU')
    # Mercurial can be confused for hydrogen gamma. Yet, mercurial is
    # rather rare in the PDB. Here are all the residues that contain
    # mercurial.
    special_hg = ('CMH', 'EMC', 'MBO', 'MMC', 'HGB', 'BE7', 'PMB')
    # Chloride can be confused for a carbon. Here are the residues that
    # contain chloride.
    special_cl = ('0QE', 'CPT', 'DCE', 'EAA', 'IMN', 'OCZ', 'OMY', 'OMZ',
                  'UN9', '1N1', '2T8', '393', '3MY', 'BMU', 'CLM', 'CP6',
                  'DB8', 'DIF', 'EFZ', 'LUR', 'RDC', 'UCL', 'XMM', 'HLT',
                  'IRE', 'LCP', 'PCI', 'VGH')
    # In these pairs, the atom name is aligned on the first column
    # (column 13).
    include_pairs = (Pair('OEC', 'CA1'),
                     Pair('PLL', 'PD'),
                     Pair('OEX', 'CA1'))
    # In these pairs, the atom name is aligned on the second column
    # (column 14), but other rules would align them on the first column.
    exclude_pairs = (Pair('C14', 'C14'), Pair('C15', 'C15'),
                     Pair('F9F', 'F9F'), Pair('OAN', 'OAN'),
                     Pair('BLM', 'NI'), Pair('BZG', 'CO'),
                     Pair('BZG', 'NI'), Pair('VNL', 'CO1'),
                     Pair('VNL', 'CO2'), Pair('PF5', 'FE1'),
                     Pair('PF5', 'FE2'), Pair('UNL', 'UNL'))

    def __init__(self, filename, bonds="conect", n_atoms=None, start=0, step=1,
                 remarks="Created by PDBWriter",
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
        # convert length and time to base units
        self.convert_units = convert_units
        self._multiframe = self.multiframe if multiframe is None else multiframe
        self.bonds = bonds

        self.frames_written = 0
        if start < 0:
            raise ValueError("'Start' must be a positive value")

        self.start = start
        self.step = step
        self.remarks = remarks

        self.pdbfile = util.anyopen(self.filename, 'wt')  # open file on init
        self.has_END = False
        self.first_frame_done = False

    def close(self):
        """Close PDB file and write END record"""
        if hasattr(self, 'pdbfile') and self.pdbfile is not None:
            if not self.has_END:
                self.END()
            else:
                logger.warn("END record has already been written"
                            " before the final closing of the file")
            self.pdbfile.close()
        self.pdbfile = None

    def _write_pdb_title(self):
        if self._multiframe:
            self.TITLE("MDANALYSIS FRAMES FROM {0:d}, SKIP {1:d}: {2!s}"
                       "".format(self.start, self.step, self.remarks))
        else:
            self.TITLE("MDANALYSIS FRAME {0:d}: {1!s}"
                       "".format(self.start, self.remarks))

    def _write_pdb_header(self):
        if not self.obj or not hasattr(self.obj, 'universe'):
            self._write_pdb_title(self)
            return
        if self.first_frame_done == True:
            return

        self.first_frame_done = True
        u = self.obj.universe
        self.HEADER(u.trajectory)

        self._write_pdb_title()

        self.COMPND(u.trajectory)
        try:
            # currently inconsistent: DCDReader gives a string,
            # PDB*Reader a list, so always get a list
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
        # can write from selection == Universe (Issue 49)
        coor = atoms.positions

        # check if any coordinates are illegal (coordinates are already in
        # Angstroem per package default)
        if self.has_valid_coordinates(self.pdb_coor_limits, coor):
            return True
        # note the precarious close() here: we know that the file is open and
        # we now prepare to remove what we have already written (header and
        # such) or add a REMARK (which allows the user to look at the
        # previously written frames)
        if self.frames_written > 1:
            self.REMARK("Incomplete multi-frame trajectory.",
                        "Coordinates for the current frame cannot be "
                        "represented in the PDB format.")
            self.close()
        else:
            self.close()
            try:
                os.remove(self.filename)
            except OSError as err:
                if err.errno == errno.ENOENT:
                    pass
        raise ValueError("PDB files must have coordinate values between "
                         "{0:.3f} and {1:.3f} Angstroem: file writing was "
                         "aborted.".format(self.pdb_coor_limits["min"],
                                           self.pdb_coor_limits["max"]))

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
           Only write CONECT records if :attr:`PDBWriter.bonds` ``== True``.
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

        atoms = {a.index for a in self.obj.atoms}

        mapping = {atom.index: i for i, atom in enumerate(self.obj.atoms)}

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

        atoms = sorted([a.index for a in self.obj.atoms])

        conect = [([a, ] + sorted(con[a])) for a in atoms if a in con]

        conect = [[mapping[e] for e in row] for row in conect]

        for c in conect:
            self.CONECT(c)

    def _update_frame(self, obj):
        """Method to initialize important attributes in writer from a AtomGroup or Universe *obj*.

        Attributes initialized/updated:

        * :attr:`PDBWriter.obj` (the entity that provides topology information *and*
          coordinates, either a :class:`~MDAnalysis.core.AtomGroup.AtomGroup` or a whole
          :class:`~MDAnalysis.core.AtomGroup.Universe`)
        * :attr:`PDBWriter.trajectory` (the underlying trajectory
          :class:`~MDAnalysis.coordinates.base.Reader`)
        * :attr:`PDBWriter.timestep` (the underlying trajectory
          :class:`~MDAnalysis.coordinates.base.Timestep`)

        Before calling :meth:`write_next_timestep` this method **must** be
        called at least once to enable extracting topology information from the
        current frame.
        """

        if isinstance(obj, base.Timestep):
            raise TypeError("PDBWriter cannot write Timestep objects "
                            "directly, since they lack topology information ("
                            "atom names and types) required in PDB files")
        # remember obj for some of other methods --- NOTE: this is an evil/lazy
        # hack...
        self.obj = obj
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
            raise AssertionError("PDBWriter couldn't extract "
                                 "trajectory and timestep information "
                                 "from an object; inheritance problem.")

        self.trajectory = traj  # update trajectory (used by other methods)
        self.ts = ts  # update timestep (used by other methods)
        return traj, ts

    def write(self, obj):
        """Write object *obj* at current trajectory frame to file.

        *obj* can be a selection (i.e. a
        :class:`~MDAnalysis.core.AtomGroup.AtomGroup`) or a whole
        :class:`~MDAnalysis.core.AtomGroup.Universe`.

        The last letter of the :attr:`~MDAnalysis.core.AtomGroup.Atom.segid` is
        used as the PDB chainID (but see :meth:`~PDBWriter.ATOM` for
        details).

        :Arguments:
          *obj*
            :class:`~MDAnalysis.core.AtomGroup.AtomGroup` or
            :class:`~MDAnalysis.core.AtomGroup.Universe`
        """

        self._update_frame(obj)
        self._write_pdb_header()
        # Issue 105: with write() ONLY write a single frame; use
        # write_all_timesteps() to dump everything in one go, or do the
        # traditional loop over frames
        self.write_next_timestep(self.ts, multiframe=self._multiframe)
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
           pdb = PDBWriter("out.pdb", u.atoms.n_atoms)
           pdb.write_all_timesteps(u)

        will write a PDB trajectory containing the last 2 frames and ::

           pdb = PDBWriter("out.pdb", u.atoms.n_atoms, start=12, skip=2)
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

        for framenumber in range(start, len(traj), step):
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
             if ``None`` then :attr:`PDBWriter.ts`` is tried.
          *multiframe*
             ``False``: write a single frame (default); ``True`` behave as a trajectory writer

        .. Note::

           Before using this method with another :class:`base.Timestep` in the *ts*
           argument, :meth:`PDBWriter._update_frame` *must* be called
           with the :class:`~MDAnalysis.core.AtomGroup.AtomGroup.Universe` as
           its argument so that topology information can be gathered.
        '''
        if ts is None:
            if not hasattr(self, "ts"):
                raise NoDataError("PBDWriter: no coordinate data to write to "
                                  "trajectory file")
            else:
                ts = self.ts
        self._check_pdb_coordinates()
        self._write_timestep(ts, **kwargs)

    def _deduce_PDB_atom_name(self, atom):
        """Deduce how the atom name should be aligned.

        Atom name format can be deduced from the atom type, yet atom type is
        not always available. This function uses the atom name and residue name
        to deduce how the atom name should be formatted. The rules in use got
        inferred from an analysis of the PDB. See gist at
        <https://gist.github.com/jbarnoud/37a524330f29b5b7b096> for more
        details.
        """
        if len(atom.name) >= 4:
            return atom.name[:4]
        elif len(atom.name) == 1:
            return ' {}  '.format(atom.name)
        elif ((atom.resname == atom.name
               or atom.name[:2] in self.ions
               or atom.name == 'UNK'
               or (atom.resname in self.special_hg and atom.name[:2] == 'HG')
               or (atom.resname in self.special_cl and atom.name[:2] == 'CL')
               or Pair(atom.resname, atom.name) in self.include_pairs)
              and Pair(atom.resname, atom.name) not in self.exclude_pairs):
            return '{:<4}'.format(atom.name)
        return ' {:<3}'.format(atom.name)

    def _write_timestep(self, ts, multiframe=False):
        """Write a new timestep *ts* to file

        Does unit conversion if necessary.

        By setting *multiframe* = ``True``, MODEL_ ... ENDMDL_ records are
        written to represent trajectory frames in a multi-model PDB file. (At
        the moment we do *not* write the NUMMDL_ record.)

        The *multiframe* = ``False`` keyword signals that the
        :class:`PDBWriter` is in single frame mode and no MODEL_
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
        atoms = self.obj.atoms
        pos = atoms.ts.positions
        if self.convert_units:
            pos = self.convert_pos_to_native(pos, inplace=False)

        if multiframe:
            self.MODEL(self.frames_written + 1)

        for i, atom in enumerate(atoms):
            segid = atom.segid if atom.segid is not "SYSTEM" else " "

            vals = {}
            vals['serial'] = int(str(i + 1)[-5:])  # check for overflow here?
            vals['name'] = self._deduce_PDB_atom_name(atom)
            vals['altLoc'] = atom.altLoc[:1] if atom.altLoc is not None else " "
            vals['resName'] = atom.resname[:4]
            vals['chainID'] = segid[:1]
            vals['resSeq'] = int(str(atom.resid)[-4:])
            vals['iCode'] = " "
            vals['pos'] = pos[i]  # don't take off atom so conversion works
            try:
                occ = atom.occupancy
            except NoDataError:
                occ = 1.0
            vals['occupancy'] = occ if occ is not None else 1.0
            temp = atom.bfactor
            vals['tempFactor'] = temp if temp is not None else 0.0
            vals['segID'] = segid[:4]
            vals['element'] = guess_atom_element(atom.name.strip())[:2]

            # .. _ATOM: http://www.wwpdb.org/documentation/format32/sect9.html
            self.pdbfile.write(self.fmt['ATOM'].format(**vals))
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
        self.pdbfile.write(self.fmt['HEADER'].format(header))

    def TITLE(self, *title):
        """Write TITLE_ record.

        .. _TITLE: http://www.wwpdb.org/documentation/format32/sect2.html

        """
        line = " ".join(title)  # TODO: should do continuation automatically
        self.pdbfile.write(self.fmt['TITLE'].format(line))

    def REMARK(self, *remarks):
        """Write generic REMARK_ record (without number).

        Each string provided in *remarks* is written as a separate REMARK_
        record.

        See also `REMARK (update)`_.

        .. _REMARK: http://www.wwpdb.org/documentation/format32/remarks1.html
        .. _REMARK (update): http://www.wwpdb.org/documentation/format32/remarks2.html

        """
        for remark in remarks:
            self.pdbfile.write(self.fmt['REMARK'].format(remark))

    def COMPND(self, trajectory):
        if not hasattr(trajectory, 'compound'):
            return
        compound = trajectory.compound
        for c in compound:
            self.pdbfile.write(self.fmt['COMPND'].format(c))

    def CRYST1(self, dimensions, spacegroup='P 1', zvalue=1):
        """Write CRYST1_ record.

        .. _CRYST1: http://www.wwpdb.org/documentation/format32/sect8.html

        """
        self.pdbfile.write(self.fmt['CRYST1'].format(
            box=dimensions[:3],
            ang=dimensions[3:],
            spacegroup=spacegroup,
            zvalue=zvalue))

    def MODEL(self, modelnumber):
        """Write the MODEL_ record.

        .. _MODEL: http://www.wwpdb.org/documentation/format32/sect9.html#MODEL

        """
        self.pdbfile.write(self.fmt['MODEL'].format(modelnumber))

    def END(self):
        """Write END_ record.

        Only a single END record is written. Calling END multiple times has no
        effect. Because :meth:`~PDBWriter.close` also calls this
        method right before closing the file it is recommended to *not* call
        :meth:`~PDBWriter.END` explicitly.

        .. _END: http://www.wwpdb.org/documentation/format32/sect11.html#END

        """
        if not self.has_END:
            # only write a single END record
            self.pdbfile.write(self.fmt['END'])
        self.has_END = True

    def ENDMDL(self):
        """Write the ENDMDL_ record.

        .. _ENDMDL: http://www.wwpdb.org/documentation/format32/sect9.html#ENDMDL

        """
        self.pdbfile.write(self.fmt['ENDMDL'])

    def CONECT(self, conect):
        """Write CONECT_ record.

        .. _CONECT: http://www.wwpdb.org/documentation/format32/sect10.html#CONECT

        """
        conect = ["{0:5d}".format(entry + 1) for entry in conect]
        conect = "".join(conect)
        self.pdbfile.write(self.fmt['CONECT'].format(conect))


class PrimitivePDBReader(PDBReader):
    def __init__(self, filename, *args, **kwargs):
        warnings.warn('PrimitivePDBReader is identical to the PDBReader,'
                  ' it is deprecated in favor of the shorter name'
                  ' removal targeted for version 0.16.0',
                  category=DeprecationWarning)
        super(PrimitivePDBReader, self).__init__(filename, *args, **kwargs)


class PrimitivePDBWriter(PDBWriter):
    def __init__(self, filename, *args, **kwargs):
        warnings.warn('PrimitivePDBWriter is identical to the Writer,'
                      'it is deprecated in favor of the shorter name'
                      ' removal targeted for version 0.16.0',
                      category=DeprecationWarning)
        super(PrimitivePDBWriter, self).__init__(filename, *args, **kwargs)

class ExtendedPDBReader(PDBReader):
    """PDBReader that reads a PDB-formatted file with five-digit residue numbers.

    This reader does not conform to the `PDB standard`_ because it allows
    five-digit residue numbers that may take up columns 23 to 27 (inclusive)
    instead of being confined to 23-26 (with column 27 being reserved for the
    insertion code in the PDB standard). PDB files in this format are written
    by popular programs such as VMD_.

    .. SeeAlso:: :class:`PDBReader`

    .. _PDB standard: http://www.wwpdb.org/documentation/format32/sect9.html
    .. _VMD: http://www.ks.uiuc.edu/Research/vmd/

    .. versionadded:: 0.8
    """
    format = "XPDB"


class MultiPDBWriter(PDBWriter):
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
       This class is identical to :class:`PDBWriter` with the one
       exception that it defaults to writing multi-frame PDB files instead of
       single frames.

    .. versionadded:: 0.7.6

    """
    format = 'PDB'
    multiframe = True  # For Writer registration
