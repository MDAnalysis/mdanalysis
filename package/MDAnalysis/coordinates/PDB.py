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
:meth:`MDAnalysis.core.groups.AtomGroup.write` method of the selection::

  calphas = universe.select_atoms("name CA")
  calphas.write("calpha_only.pdb")

This uses the coordinates from the current timestep of the trajectory.

In order to write a PDB trajectory one needs to first obtain a multi-frame
writer (keyword *multiframe* = ``True``) and then iterate through the
trajectory, while writing each frame::

  calphas = universe.select_atoms("name CA")
  with MDAnalysis.Writer("calpha_traj.pdb", multiframe=True) as W:
      for ts in u.trajectory:
          W.write(calphas)

It is important to *always close the trajectory* when done because only at this
step is the final END_ record written, which is required by the `PDB 3.3
standard`_. Using the writer as a context manager ensures that this always
happens.


Capabilities
------------

A pure-Python implementation for PDB files commonly encountered in MD
simulations comes under the names :class:`PDBReader` and :class:`PDBWriter`. It
only implements a subset of the `PDB 3.3 standard`_ and also allows some
typical enhancements such as 4-letter resids (introduced by CHARMM/NAMD).

The :class:`PDBReader` can read multi-frame PDB files and represents
them as a trajectory. The :class:`PDBWriter` can write single and
multi-frame PDB files as specified by the *multiframe* keyword. By default, it
writes single frames. On the other hand, the :class:`MultiPDBWriter` is set up
to write a PDB trajectory by default (equivalent to using *multiframe* =
``True``).


Examples for working with PDB files
-----------------------------------

A **single frame PDB** can be written with the
:meth:`~MDAnalysis.core.groups.AtomGroup.write` method of any
:class:`~MDAnalysis.core.groups.AtomGroup`::

   protein = u.select_atoms("protein")
   protein.write("protein.pdb")

Alternatively, get the single frame writer and supply the
:class:`~MDAnalysis.core.groups.AtomGroup`::

  protein = u.select_atoms("protein")
  with MDAnalysis.Writer("protein.pdb") as pdb:
      pdb.write(protein)

In order to write a **multi-frame PDB trajectory** from a universe *u* one can
do the following::

  with MDAnalysis.Writer("all.pdb", multiframe=True) as pdb:
      for ts in u.trajectory:
          pdb.write(u)

Similarly, writing only a protein::

  protein = u.select_atoms("protein")
  with MDAnalysis.Writer("protein.pdb", multiframe=True) as pdb:
      for ts in u.trajectory:
          pdb.write(protein)



Classes
-------

.. versionchanged:: 0.16.0
   PDB readers and writers based on :class:`Bio.PDB.PDBParser` were retired and
   removed.


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

.. autoclass:: ExtendedPDBReader
   :members:
   :inherited-members:


.. _`PDB 3.3 standard`:
    http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html

"""
from io import StringIO, BytesIO
import os
import errno
import itertools
import textwrap
import warnings
import logging
import collections
import numpy as np
import functools

from ..lib import util
from ..lib.util import store_init_arguments
from . import base
from .timestep import Timestep
from ..exceptions import NoDataError


logger = logging.getLogger("MDAnalysis.coordinates.PBD")

# Pairs of residue name / atom name in use to deduce PDB formatted atom names
Pair = collections.namedtuple("Atom", "resname name")


class PDBReader(base.ReaderBase):
    """PDBReader that reads a `PDB-formatted`_ file, no frills.

    The following *PDB records* are parsed (see `PDB coordinate section`_ for
    details):

    - *CRYST1* for unitcell A,B,C, alpha,beta,gamma
    - *ATOM* or *HETATM* for serial,name,resName,chainID,resSeq,x,y,z,occupancy,tempFactor
    - *HEADER* (:attr:`header`), *TITLE* (:attr:`title`), *COMPND*
        (:attr:`compound`), *REMARK* (:attr:`remarks`)
    - all other lines are ignored


    Reads multi-`MODEL`_ PDB files as trajectories.  The `Timestep.data` dictionary
    holds the occupancy and tempfactor (bfactor) values for each atom for a given frame.
    These attributes are commonly appropriated to store other time varying properties
    and so they are exposed here. Note this does not update the `AtomGroup` attributes,
    as the topology does not change with trajectory iteration.

    .. _PDB-formatted:
       http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html
    .. _PDB coordinate section:
       http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html
    .. _MODEL:
       http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#MODEL

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
    17             Character     altLoc       Alternate location indicator.
    18 - 21        Residue name  resName      Residue name.
    22             Character     chainID      Chain identifier.
    23 - 26        Integer       resSeq       Residue sequence number.
    27             AChar         iCode        Code for insertion of residues.
    31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    55 - 60        Real(6.2)     occupancy    Occupancy.
    61 - 66        Real(6.2)     tempFactor   Temperature  factor.
    67 - 72                                   (not used in the official PDB format)
    73 - 76        String        segID        (unofficial PDB format*)
    77 - 78        LString(2)    element      Element symbol, right-justified.
    79 - 80        LString(2)    charge       Charge  on the atom.
    =============  ============  ===========  =============================================

    Notes
    -----
    If a system does not have unit cell parameters (such as in electron
    microscopy structures), the PDB file format requires the CRYST1_ field to
    be provided with unitary values (cubic box with sides of 1 Å) and an
    appropriate REMARK. If unitary values are found within the CRYST1_ field,
    :code:`PDBReader` will not set unit cell dimensions (which will take the
    default value :code:`np.zeros(6)`, see Issue #2698)
    and it will warn the user.

    .. _CRYST1: http://www.wwpdb.org/documentation/file-format-content/format33/sect8.html#CRYST1

    *The columns 73-76 are not part of the official PDB format but are used by
    some programs to store/operate the segment ID. For instance, Chimera_ assigns
    it as the attribute `pdbSegment` to allow command-line specification.

    .. _Chimera:
        https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html#note6

    See Also
    --------
    :class:`PDBWriter`
    :class:`PDBReader`

    If you would like to force the use of chainID as the segID when parsing PDB,
    please use the keyword argument `force_chainids_to_segids=True`
    (:class:`MDAnalysis.topology.PDBParser.PDBParser`).
    This will prioritize the chain ID to the segment ID.


    .. versionchanged:: 0.11.0
       * Frames now 0-based instead of 1-based
       * New :attr:`title` (list with all TITLE lines).
    .. versionchanged:: 0.19.1
       Can now read PDB files with DOS line endings
    .. versionchanged:: 0.20.0
       Strip trajectory header of trailing spaces and newlines
    .. versionchanged:: 1.0.0
       Raise user warning for CRYST1_ record with unitary values
       (cubic box with sides of 1 Å) and do not set cell dimensions.
    .. versionchanged:: 2.5.0
       Tempfactors (aka bfactors) are now read into the ts.data dictionary each
       frame.  Occupancies are also read into this dictionary.
    """

    format = ["PDB", "ENT"]
    units = {"time": None, "length": "Angstrom"}

    @store_init_arguments
    def __init__(self, filename, **kwargs):
        """Read coordinates from *filename*.

        *filename* can be a gzipped or bzip2ed compressed PDB file.

        If the pdb file contains multiple MODEL records then it is
        read as a trajectory where the MODEL numbers correspond to
        frame numbers.
        """
        super(PDBReader, self).__init__(filename, **kwargs)

        try:
            self.n_atoms = kwargs["n_atoms"]
        except KeyError:
            # hackish, but should work and keeps things DRY
            # regular MDA usage via Universe doesn't follow this route
            from MDAnalysis.topology import PDBParser

            with PDBParser.PDBParser(self.filename) as p:
                top = p.parse()
            self.n_atoms = top.n_atoms

        self.model_offset = kwargs.pop("model_offset", 0)
        # dummy/default variables as these are read
        header = ""
        title = []
        compound = []
        remarks = []

        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)

        # Record positions in file of CRYST and MODEL headers
        # then build frame offsets to start at the minimum of these
        # This allows CRYST to come either before or after MODEL
        # This assumes that **either**
        # - pdbfile has a single CRYST (NVT)
        # - pdbfile has a CRYST for every MODEL (NPT)
        models = []
        crysts = []

        # hack for streamIO
        if isinstance(filename, util.NamedStream) and isinstance(
            filename.stream, StringIO
        ):
            filename.stream = BytesIO(filename.stream.getvalue().encode())

        pdbfile = self._pdbfile = util.anyopen(filename, "rb")

        line = "magical"
        while line:
            # need to use readline so tell gives end of line
            # (rather than end of current chunk)
            line = pdbfile.readline()

            if line[:5] == b"MODEL":
                models.append(pdbfile.tell())
            elif line[:5] == b"CRYST":
                # remove size of line to get **start** of CRYST line
                crysts.append(pdbfile.tell() - len(line))
            elif line[:6] == b"HEADER":
                # classification = line[10:50]
                # date = line[50:59]
                # idCode = line[62:66]
                header = line[10:66].strip().decode()
            elif line[:5] == b"TITLE":
                title.append(line[8:80].strip().decode())
            elif line[:6] == b"COMPND":
                compound.append(line[7:80].strip().decode())
            elif line[:6] == b"REMARK":
                remarks.append(line[6:].strip().decode())

        end = pdbfile.tell()  # where the file ends

        self.header = header
        self.title = title
        self.compound = compound
        self.remarks = remarks

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

        Parameters
        ----------
        filename : str
            filename of the output PDB file

        Returns
        -------
        :class:`PDBWriter`

        """
        kwargs.setdefault("multiframe", self.n_frames > 1)
        return PDBWriter(filename, **kwargs)

    def _reopen(self):
        # Pretend the current TS is -1 (in 0 based) so "next" is the
        # 0th frame
        self.close()
        self._pdbfile = util.anyopen(self.filename, "rb")
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
        """
        Read frame from PDB file.

        Notes
        -----
        When the CRYST1_ record has unitary values (cubic box with sides of
        1 Å), cell dimensions are considered fictitious. An user warning is
        raised and cell dimensions are set to
        :code:`np.zeros(6)` (see Issue #2698)


        .. versionchanged:: 1.0.0
           Raise user warning for CRYST1_ record with unitary valuse
           (cubic box with sides of 1 Å) and do not set cell dimensions.
        .. versionchanged:: 2.5.0
           When seen, `ts.data` is populated with tempfactor information at
           each frame read. If any atom has a non-parsable (i.e. non float)
           value in the tempfactor field, the entry is left as `1.0`.
        """
        try:
            start = self._start_offsets[frame]
            stop = self._stop_offsets[frame]
        except IndexError:  # out of range of known frames
            raise IOError from None

        pos = 0
        occupancy = np.zeros(self.n_atoms)
        tempfactor = np.ones(self.n_atoms)
        saw_tempfactor = False

        # Seek to start and read until start of next frame
        self._pdbfile.seek(start)
        chunk = self._pdbfile.read(stop - start).decode()

        tmp_buf = []
        for line in chunk.splitlines():
            if line[:6] in ("ATOM  ", "HETATM"):
                # we only care about coordinates
                tmp_buf.append([line[30:38], line[38:46], line[46:54]])
                try:
                    # does an implicit str -> float conversion
                    occupancy[pos] = line[54:60]
                except ValueError:
                    # Be tolerant for ill-formated or empty occupancies
                    pass
                try:
                    tempfactor[pos] = line[60:66]
                except ValueError:
                    pass
                else:
                    saw_tempfactor = True
                pos += 1
            elif line[:6] == "CRYST1":
                # does an implicit str -> float conversion
                try:
                    cell_dims = np.array(
                        [
                            line[6:15],
                            line[15:24],
                            line[24:33],
                            line[33:40],
                            line[40:47],
                            line[47:54],
                        ],
                        dtype=np.float32,
                    )
                except ValueError:
                    warnings.warn(
                        "Failed to read CRYST1 record, "
                        "possibly invalid PDB file, got:\n{}"
                        "".format(line)
                    )
                    self.ts.dimensions = None
                else:
                    if np.allclose(
                        cell_dims, np.array([1.0, 1.0, 1.0, 90.0, 90.0, 90.0])
                    ):
                        warnings.warn(
                            "1 A^3 CRYST1 record,"
                            " this is usually a placeholder."
                            " Unit cell dimensions will be set to None."
                        )
                        self.ts.dimensions = None
                    else:
                        self.ts.dimensions = cell_dims

        # check if atom number changed
        if pos != self.n_atoms:
            raise ValueError(
                "Inconsistency in file '{}': The number of atoms "
                "({}) in trajectory frame {} differs from the "
                "number of atoms ({}) in the corresponding "
                "topology.\nTrajectories with varying numbers of "
                "atoms are currently not supported."
                "".format(self.filename, pos, frame, self.n_atoms)
            )

        # doing the conversion from list to array at the end is faster
        self.ts.positions = tmp_buf

        if self.convert_units:
            # both happen inplace
            self.convert_pos_from_native(self.ts._pos)
            if not self.ts.dimensions is None:
                self.convert_pos_from_native(self.ts.dimensions[:3])
        self.ts.frame = frame
        self.ts.data["occupancy"] = occupancy
        if saw_tempfactor:
            self.ts.data["tempfactor"] = tempfactor

        return self.ts

    def close(self):
        self._pdbfile.close()


class PDBWriter(base.WriterBase):
    """PDB writer that implements a subset of the `PDB 3.3 standard`_ .

    PDB format as used by NAMD/CHARMM: 4-letter resnames and segID are allowed,
    altLoc is written.

    The :class:`PDBWriter` can be used to either dump a coordinate
    set to a PDB file (operating as a "single frame writer", selected with the
    constructor keyword *multiframe* = ``False``, the default) or by writing a
    PDB "movie" (multi frame mode, *multiframe* = ``True``), consisting of
    multiple models (using the MODEL_ and ENDMDL_ records).

    .. _`PDB 3.3 standard`:
       http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html
    .. _ATOM: http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
    .. _COMPND: http://www.wwpdb.org/documentation/file-format-content/format33/sect2.html#COMPND
    .. _CONECT: http://www.wwpdb.org/documentation/file-format-content/format33/sect10.html#CONECT
    .. _END: http://www.wwpdb.org/documentation/file-format-content/format33/sect11.html#END
    .. _ENDMDL: http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ENDMDL
    .. _HEADER: http://www.wwpdb.org/documentation/file-format-content/format33/sect2.html#HEADER
    .. _HETATM: http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#HETATM
    .. _MODEL: http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#MODEL
    .. _NUMMDL: http://www.wwpdb.org/documentation/file-format-content/format33/sect2.html#NUMMDL
    .. _REMARKS: http://www.wwpdb.org/documentation/file-format-content/format33/remarks.html
    .. _TITLE: http://www.wwpdb.org/documentation/file-format-content/format33/sect2.html#TITLE

    Note
    ----
    Writing bonds currently only works when writing a whole :class:`Universe`
    and if bond information is available in the topology.  (For selections
    smaller than the whole :class:`Universe`, the atom numbering in the CONECT_
    records would not match the numbering of the atoms in the new PDB file and
    therefore a :exc:`NotImplementedError` is raised.)

    The maximum frame number that can be stored in a PDB file is 9999 and it
    will wrap around (see :meth:`MODEL` for further details).

    The CRYST1_ record specifies the unit cell. This record is set to
    unitary values (cubic box with sides of 1 Å) if unit cell dimensions
    are not set (:code:`None` or :code:`np.zeros(6)`,
    see Issue #2698).

    When the :attr:`record_types` attribute is present (e.g. Universe object
    was created by loading a PDB file), ATOM_ and HETATM_ record type
    keywords are written out accordingly. Otherwise, the ATOM_ record type
    is the default output.

    The CONECT_ record is written out, if required, when the output stream
    is closed.

    See Also
    --------
    This class is identical to :class:`MultiPDBWriter` with the one
    exception that it defaults to writing single-frame PDB files as if
    `multiframe` = ``False`` was selected.


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

    .. versionchanged:: 0.20.0
       Strip trajectory header of trailing spaces and newlines

    .. versionchanged:: 1.0.0
       ChainID now comes from the last character of segid, as stated in the documentation.
       An indexing issue meant it previously used the first charater (Issue #2224)

    .. versionchanged:: 2.0.0
       Add the `redindex` argument. Setting this keyword to ``True``
       (the default) preserves the behavior in earlier versions of MDAnalysis.
       The PDB writer checks for a valid chainID entry instead of using the
       last character of segid. Should a chainID not be present, or not
       conform to the PDB standard, the default value of  'X' is used.

    .. versionchanged:: 2.3.0
       Do not write unusable conect records when ag index
       is larger than 100000.
    """

    fmt = {
        "ATOM": (
            "ATOM  {serial:5d} {name:<4s}{altLoc:<1s}{resName:<4s}"
            "{chainID:1s}{resSeq:4d}{iCode:1s}"
            "   {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}{occupancy:6.2f}"
            "{tempFactor:6.2f}      {segID:<4s}{element:>2s}{charge:2s}\n"
        ),
        "HETATM": (
            "HETATM{serial:5d} {name:<4s}{altLoc:<1s}{resName:<4s}"
            "{chainID:1s}{resSeq:4d}{iCode:1s}"
            "   {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}{occupancy:6.2f}"
            "{tempFactor:6.2f}      {segID:<4s}{element:>2s}{charge:2s}\n"
        ),
        "REMARK": "REMARK     {0}\n",
        "COMPND": "COMPND    {0}\n",
        "HEADER": "HEADER    {0}\n",
        "TITLE": "TITLE     {0}\n",
        "MODEL": "MODEL     {0:>4d}\n",
        "NUMMDL": "NUMMDL    {0:5d}\n",
        "ENDMDL": "ENDMDL\n",
        "END": "END\n",
        "CRYST1": (
            "CRYST1{box[0]:9.3f}{box[1]:9.3f}{box[2]:9.3f}"
            "{ang[0]:7.2f}{ang[1]:7.2f}{ang[2]:7.2f} "
            "{spacegroup:<11s}{zvalue:4d}\n"
        ),
        "CONECT": "CONECT{0}\n",
    }
    format = ["PDB", "ENT"]
    units = {"time": None, "length": "Angstrom"}
    pdb_coor_limits = {"min": -999.9995, "max": 9999.9995}
    #: wrap comments into REMARK records that are not longer than
    # :attr:`remark_max_length` characters.
    remark_max_length = 66
    multiframe = False

    # These attributes are used to deduce how to format the atom name.
    ions = (
        "FE",
        "AS",
        "ZN",
        "MG",
        "MN",
        "CO",
        "BR",
        "CU",
        "TA",
        "MO",
        "AL",
        "BE",
        "SE",
        "PT",
        "EU",
        "NI",
        "IR",
        "RH",
        "AU",
        "GD",
        "RU",
    )
    # Mercurial can be confused for hydrogen gamma. Yet, mercurial is
    # rather rare in the PDB. Here are all the residues that contain
    # mercurial.
    special_hg = ("CMH", "EMC", "MBO", "MMC", "HGB", "BE7", "PMB")
    # Chloride can be confused for a carbon. Here are the residues that
    # contain chloride.
    special_cl = (
        "0QE",
        "CPT",
        "DCE",
        "EAA",
        "IMN",
        "OCZ",
        "OMY",
        "OMZ",
        "UN9",
        "1N1",
        "2T8",
        "393",
        "3MY",
        "BMU",
        "CLM",
        "CP6",
        "DB8",
        "DIF",
        "EFZ",
        "LUR",
        "RDC",
        "UCL",
        "XMM",
        "HLT",
        "IRE",
        "LCP",
        "PCI",
        "VGH",
    )
    # In these pairs, the atom name is aligned on the first column
    # (column 13).
    include_pairs = (Pair("OEC", "CA1"), Pair("PLL", "PD"), Pair("OEX", "CA1"))
    # In these pairs, the atom name is aligned on the second column
    # (column 14), but other rules would align them on the first column.
    exclude_pairs = (
        Pair("C14", "C14"),
        Pair("C15", "C15"),
        Pair("F9F", "F9F"),
        Pair("OAN", "OAN"),
        Pair("BLM", "NI"),
        Pair("BZG", "CO"),
        Pair("BZG", "NI"),
        Pair("VNL", "CO1"),
        Pair("VNL", "CO2"),
        Pair("PF5", "FE1"),
        Pair("PF5", "FE2"),
        Pair("UNL", "UNL"),
    )

    def __init__(
        self,
        filename,
        bonds="conect",
        n_atoms=None,
        start=0,
        step=1,
        remarks="Created by PDBWriter",
        convert_units=True,
        multiframe=None,
        reindex=True,
    ):
        """Create a new PDBWriter

        Parameters
        ----------
        filename: str
           name of output file
        start: int (optional)
           starting timestep (the first frame will have MODEL number `start` + 1
           because the PDB standard prescribes MODEL numbers starting at 1)
        step: int (optional)
           skip between subsequent timesteps
        remarks: str (optional)
           comments to annotate pdb file (added to the *TITLE* record); note that
           any remarks from the trajectory that serves as input are
           written to REMARK records with lines longer than :attr:`remark_max_length` (66
           characters) being wrapped.
        convert_units: bool (optional)
           units are converted to the MDAnalysis base format; [``True``]
        bonds : {"conect", "all", None} (optional)
           If set to "conect", then only write those bonds that were already
           defined in an input PDB file as PDB CONECT_ record. If set to "all",
           write all bonds (including guessed ones) to the file. ``None`` does
           not write any bonds. The default is "conect".
        multiframe: bool (optional)
           ``False``: write a single frame to the file; ``True``: create a
           multi frame PDB file in which frames are written as MODEL_ ... ENDMDL_
           records. If ``None``, then the class default is chosen.    [``None``]
        reindex: bool (optional)
            If ``True`` (default), the atom serial is set to be consecutive
            numbers starting at 1. Else, use the atom id.

        """
        # n_atoms = None : dummy keyword argument
        # (not used, but Writer() always provides n_atoms as the second argument)

        # TODO: - remarks should be a list of lines and written to REMARK
        #       - additional title keyword could contain line for TITLE

        self.filename = filename
        # convert length and time to base units
        self.convert_units = convert_units
        self._multiframe = (
            self.multiframe if multiframe is None else multiframe
        )
        self.bonds = bonds
        self._reindex = reindex

        if start < 0:
            raise ValueError("'Start' must be a positive value")

        self.start = self.frames_written = start
        self.step = step
        self.remarks = remarks

        self.pdbfile = util.anyopen(self.filename, "wt")  # open file on init
        self.has_END = False
        self.first_frame_done = False

    def close(self):
        """
        Close PDB file and write CONECT and END record


        .. versionchanged:: 2.0.0
           CONECT_ record written just before END_ record
        """
        if hasattr(self, "pdbfile") and self.pdbfile is not None:
            if not self.has_END:
                self._write_pdb_bonds()
                self.END()
            else:
                logger.warning(
                    "END record has already been written"
                    " before the final closing of the file"
                )
            self.pdbfile.close()
        self.pdbfile = None

    def _write_pdb_title(self):
        if self._multiframe:
            self.TITLE(
                "MDANALYSIS FRAMES FROM {0:d}, STEP {1:d}: {2!s}"
                "".format(self.start, self.step, self.remarks)
            )
        else:
            self.TITLE(
                "MDANALYSIS FRAME {0:d}: {1!s}"
                "".format(self.start, self.remarks)
            )

    def _write_pdb_header(self):
        """
        Write PDB header.

        The HEADER_ record is set to :code: `trajectory.header`.
        The TITLE_ record explicitly mentions MDAnalysis and contains
        information about trajectory frame(s).
        The COMPND_ record is set to :code:`trajectory.compound`.
        The REMARKS_ records are set to :code:`u.trajectory.remarks`
        The CRYST1_ record specifies the unit cell. This record is set to
        unitary values (cubic box with sides of 1 Å) if unit cell dimensions
        are not set.

        .. versionchanged: 1.0.0
           Fix writing of PDB file without unit cell dimensions (Issue #2679).
           If cell dimensions are not found, unitary values (cubic box with
           sides of 1 Å) are used (PDB standard for CRYST1_).
        """

        if self.first_frame_done is True:
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

        if u.trajectory.ts.dimensions is None:
            # Unitary unit cell by default. See PDB standard:
            # http://www.wwpdb.org/documentation/file-format-content/format33/sect8.html#CRYST1
            self.CRYST1(np.array([1.0, 1.0, 1.0, 90.0, 90.0, 90.0]))

            # Add CRYST1 REMARK (285)
            # The SCALE record is not included
            # (We are only implementing a subset of the PDB standard)
            self.REMARK(
                "285 UNITARY VALUES FOR THE UNIT CELL AUTOMATICALLY SET"
            )
            self.REMARK(
                "285 BY MDANALYSIS PDBWRITER BECAUSE UNIT CELL INFORMATION"
            )
            self.REMARK("285 WAS MISSING.")
            self.REMARK("285 PROTEIN DATA BANK CONVENTIONS REQUIRE THAT")
            self.REMARK("285 CRYST1 RECORD IS INCLUDED, BUT THE VALUES ON")
            self.REMARK("285 THIS RECORD ARE MEANINGLESS.")

            warnings.warn(
                "Unit cell dimensions not found. "
                "CRYST1 record set to unitary values."
            )
        else:
            self.CRYST1(self.convert_dimensions_to_unitcell(u.trajectory.ts))

    def _check_pdb_coordinates(self):
        """Check if the coordinate values fall within the range allowed for PDB files.

        Deletes the output file if this is the first frame or if frames have
        already been written (in multi-frame mode) adds a REMARK instead of the
        coordinates and closes the file.

        Raises
        ------
        ValueError
            if the coordinates fail the check.

        .. versionchanged: 1.0.0
            Check if :attr:`filename` is `StringIO` when attempting to remove
            a PDB file with invalid coordinates (Issue #2512)
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
            self.REMARK(
                "Incomplete multi-frame trajectory.",
                "Coordinates for the current frame cannot be "
                "represented in the PDB format.",
            )
            self.close()
        else:
            self.close()
            try:
                os.remove(self.filename)
            except OSError as err:
                if err.errno == errno.ENOENT:
                    pass
                else:
                    raise
            except TypeError:
                if isinstance(self.filename, StringIO):
                    pass
                else:
                    raise

        raise ValueError(
            "PDB files must have coordinate values between "
            "{0:.3f} and {1:.3f} Angstroem: file writing was "
            "aborted.".format(
                self.pdb_coor_limits["min"], self.pdb_coor_limits["max"]
            )
        )

    def _write_pdb_bonds(self):
        """Writes out all the bond records"""
        if self.bonds is None:
            return

        if (
            not hasattr(self, "obj")
            or not self.obj
            or not hasattr(self.obj.universe, "bonds")
        ):
            return

        bondset = set(itertools.chain(*(a.bonds for a in self.obj.atoms)))
        if self._reindex:
            index_attribute = "index"
            mapping = {
                index: i
                for i, index in enumerate(self.obj.atoms.indices, start=1)
            }
            atoms = np.sort(self.obj.atoms.indices)
        else:
            index_attribute = "id"
            mapping = {id_: id_ for id_ in self.obj.atoms.ids}
            atoms = np.sort(self.obj.atoms.ids)
        if self.bonds == "conect":
            # Write out only the bonds that were defined in CONECT records
            bonds = (
                (
                    getattr(bond[0], index_attribute),
                    getattr(bond[1], index_attribute),
                )
                for bond in bondset
                if not bond.is_guessed
            )
        elif self.bonds == "all":
            bonds = (
                (
                    getattr(bond[0], index_attribute),
                    getattr(bond[1], index_attribute),
                )
                for bond in bondset
            )
        else:
            raise ValueError("bonds has to be either None, 'conect' or 'all'")

        con = collections.defaultdict(list)
        for a1, a2 in bonds:
            if not (a1 in mapping and a2 in mapping):
                continue
            if mapping[a1] >= 100000 or mapping[a2] >= 100000:
                warnings.warn(
                    "Atom with index >=100000 cannot write "
                    "bonds to PDB CONECT records."
                )
                return
            con[a2].append(a1)
            con[a1].append(a2)

        conect = (
            [mapping[a]] + sorted([mapping[at] for at in con[a]])
            for a in atoms
            if a in con
        )

        for c in conect:
            self.CONECT(c)

    def _update_frame(self, obj):
        """Method to initialize important attributes in writer from a AtomGroup or Universe *obj*.

        Attributes initialized/updated:

        * :attr:`PDBWriter.obj` (the entity that provides topology information *and*
          coordinates, either a :class:`~MDAnalysis.core.groups.AtomGroup` or a whole
          :class:`~MDAnalysis.core.universe.Universe`)
        * :attr:`PDBWriter.trajectory` (the underlying trajectory
          :class:`~MDAnalysis.coordinates.base.Reader`)
        * :attr:`PDBWriter.timestep` (the underlying trajectory
          :class:`~MDAnalysis.coordinates.timestep.Timestep`)

        Before calling :meth:`_write_next_frame` this method **must** be
        called at least once to enable extracting topology information from the
        current frame.
        """
        if isinstance(obj, Timestep):
            raise TypeError(
                "PDBWriter cannot write Timestep objects "
                "directly, since they lack topology information ("
                "atom names and types) required in PDB files"
            )
        if len(obj.atoms) == 0:
            raise IndexError("Cannot write empty AtomGroup")

        # remember obj for some of other methods --- NOTE: this is an evil/lazy
        # hack...
        self.obj = obj
        self.ts = obj.universe.trajectory.ts

    def write(self, obj):
        """Write object *obj* at current trajectory frame to file.

        *obj* can be a selection (i.e. a
        :class:`~MDAnalysis.core.groups.AtomGroup`) or a whole
        :class:`~MDAnalysis.core.universe.Universe`.

        The last letter of the :attr:`~MDAnalysis.core.groups.Atom.segid` is
        used as the PDB chainID (but see :meth:`~PDBWriter.ATOM` for
        details).

        Parameters
        ----------
        obj
            The :class:`~MDAnalysis.core.groups.AtomGroup` or
            :class:`~MDAnalysis.core.universe.Universe` to write.
        """

        self._update_frame(obj)
        self._write_pdb_header()
        # Issue 105: with write() ONLY write a single frame; use
        # write_all_timesteps() to dump everything in one go, or do the
        # traditional loop over frames
        self._write_next_frame(self.ts, multiframe=self._multiframe)
        # END and CONECT records are written when file is being close()d

    def write_all_timesteps(self, obj):
        """Write all timesteps associated with *obj* to the PDB file.

        *obj* can be a :class:`~MDAnalysis.core.groups.AtomGroup`
        or a :class:`~MDAnalysis.core.universe.Universe`.

        The method writes the frames from the one specified as *start* until
        the end, using a step of *step* (*start* and *step* are set in the
        constructor). Thus, if *u* is a Universe then ::

           u.trajectory[-2]
           pdb = PDBWriter("out.pdb", u.atoms.n_atoms)
           pdb.write_all_timesteps(u)

        will write a PDB trajectory containing the last 2 frames and ::

           pdb = PDBWriter("out.pdb", u.atoms.n_atoms, start=12, step=2)
           pdb.write_all_timesteps(u)

        will be writing frames 12, 14, 16, ...


        .. versionchanged:: 0.11.0
           Frames now 0-based instead of 1-based

        .. versionchanged:: 2.0.0
           CONECT_ record moved to :meth:`close`
        """

        self._update_frame(obj)
        self._write_pdb_header()

        start, step = self.start, self.step
        traj = obj.universe.trajectory

        # Start from trajectory[0]/frame 0, if there are more than 1 frame.
        # If there is only 1 frame, the traj.frames is not like a python list:
        # accessing trajectory[-1] raises key error.
        if not start and traj.n_frames > 1:
            start = traj.frame

        for framenumber in range(start, len(traj), step):
            traj[framenumber]
            self._write_next_frame(self.ts, multiframe=True)

        # CONECT record is written when the file is being close()d
        self.close()

        # Set the trajectory to the starting position
        traj[start]

    def _write_next_frame(self, ts=None, **kwargs):
        """write a new timestep to the PDB file

        :Keywords:
          *ts*
             :class:`timestep.Timestep` object containing coordinates to be written to trajectory file;
             if ``None`` then :attr:`PDBWriter.ts`` is tried.
          *multiframe*
             ``False``: write a single frame (default); ``True`` behave as a trajectory writer

        .. Note::

           Before using this method with another :class:`timestep.Timestep` in the *ts*
           argument, :meth:`PDBWriter._update_frame` *must* be called
           with the :class:`~MDAnalysis.core.groups.AtomGroup.Universe` as
           its argument so that topology information can be gathered.


        .. versionchanged:: 1.0.0
           Renamed from `write_next_timestep` to `_write_next_frame`.
        """
        if ts is None:
            try:
                ts = self.ts
            except AttributeError:
                errmsg = (
                    "PBDWriter: no coordinate data to write to "
                    "trajectory file"
                )
                raise NoDataError(errmsg) from None
        self._check_pdb_coordinates()
        self._write_timestep(ts, **kwargs)

    @functools.cache
    def _deduce_PDB_atom_name(self, atomname, resname):
        """Deduce how the atom name should be aligned.

        Atom name format can be deduced from the atom type, yet atom type is
        not always available. This function uses the atom name and residue name
        to deduce how the atom name should be formatted. The rules in use got
        inferred from an analysis of the PDB. See gist at
        <https://gist.github.com/jbarnoud/37a524330f29b5b7b096> for more
        details.
        """
        if atomname == "":
            return ""
        if len(atomname) >= 4:
            return atomname[:4]
        elif len(atomname) == 1:
            return " {}  ".format(atomname)
        elif (
            resname == atomname
            or atomname[:2] in self.ions
            or atomname == "UNK"
            or (resname in self.special_hg and atomname[:2] == "HG")
            or (resname in self.special_cl and atomname[:2] == "CL")
            or Pair(resname, atomname) in self.include_pairs
        ) and Pair(resname, atomname) not in self.exclude_pairs:
            return "{:<4}".format(atomname)
        return " {:<3}".format(atomname)

    @staticmethod
    def _format_PDB_charges(charges: np.ndarray) -> np.ndarray:
        """Format formal charges to match PDB requirements.

        Formal charge entry is set to empty if charge is 0, otherwise the
        charge is set to a two character ```<charge value><charge sign>``
        entry, e.g. ``1+`` or ``2-``.

        This method also verifies that formal charges can adhere to the PDB
        format (i.e. charge cannot be > 9 or < -9).

        Parameters
        ----------
        charges: np.ndarray
            NumPy array of integers representing the formal charges of
            the atoms being written.

        Returns
        -------
        np.ndarray
            NumPy array of dtype object with strings representing the
            formal charges of the atoms being written.
        """
        if not np.issubdtype(charges.dtype, np.integer):
            raise ValueError("formal charges array should be of `int` type")

        outcharges = charges.astype(object)
        outcharges[outcharges == 0] = ""  # empty strings for no charge case
        # using np.where is more efficient than looping in sparse cases
        for i in np.where(charges < 0)[0]:
            if charges[i] < -9:
                errmsg = "formal charge < -9 is not supported by PDB standard"
                raise ValueError(errmsg)
            outcharges[i] = f"{abs(charges[i])}-"
        for i in np.where(charges > 0)[0]:
            if charges[i] > 9:
                errmsg = "formal charge > 9 is not supported by PDB standard"
                raise ValueError(errmsg)
            outcharges[i] = f"{charges[i]}+"
        return outcharges

    def _write_timestep(self, ts, multiframe=False):
        """Write a new timestep *ts* to file

        Does unit conversion if necessary.

        By setting *multiframe* = ``True``, MODEL_ ... ENDMDL_ records are
        written to represent trajectory frames in a multi-model PDB file. (At
        the moment we do *not* write the NUMMDL_ record.)

        The *multiframe* = ``False`` keyword signals that the
        :class:`PDBWriter` is in single frame mode and no MODEL_
        records are written.

        .. versionchanged:: 0.7.6
           The *multiframe* keyword was added, which completely determines if
           MODEL_ records are written. (Previously, this was decided based on
           the underlying trajectory and only if ``len(traj) > 1`` would
           MODEL records have been written.)

        .. versionchanged:: 1.0.0
           ChainID now comes from the last character of segid, as stated in
           the documentation. An indexing issue meant it previously used the
           first charater (Issue #2224)

        .. versionchanged:: 2.0.0
           When only :attr:`record_types` attribute is present, instead of
           using ATOM_ for both ATOM_ and HETATM_, HETATM_ record
           types are properly written out (Issue #1753).
           Writing now only uses the contents of the elements attribute
           instead of guessing by default. If the elements are missing,
           empty records are written out (Issue #2423).
           Atoms are now checked for a valid chainID instead of being
           overwritten by the last letter of the `segid` (Issue #3144).

        """
        atoms = self.obj.atoms
        pos = atoms.positions
        if self.convert_units:
            pos = self.convert_pos_to_native(pos, inplace=False)

        if multiframe:
            self.MODEL(self.frames_written + 1)

        # Make zero assumptions on what information the AtomGroup has!
        # theoretically we could get passed only indices!
        def get_attr(attrname, default):
            """Try and pull info off atoms, else fake it

            attrname - the field to pull of AtomGroup (plural!)
            default - default value in case attrname not found
            """
            try:
                return getattr(atoms, attrname)
            except AttributeError:
                if self.frames_written == 0:
                    warnings.warn(
                        "Found no information for attr: '{}'"
                        " Using default value of '{}'"
                        "".format(attrname, default)
                    )
                return np.array([default] * len(atoms))

        altlocs = get_attr("altLocs", " ")
        resnames = get_attr("resnames", "UNK")
        icodes = get_attr("icodes", " ")
        segids = get_attr("segids", " ")
        chainids = get_attr("chainIDs", "")
        resids = get_attr("resids", 1)
        occupancies = get_attr("occupancies", 1.0)
        tempfactors = get_attr("tempfactors", 0.0)
        atomnames = get_attr("names", "X")
        elements = get_attr("elements", " ")
        record_types = get_attr("record_types", "ATOM")
        formal_charges = self._format_PDB_charges(get_attr("formalcharges", 0))

        def validate_chainids(chainids, default):
            """Validate each atom's chainID

            chainids - np array of chainIDs
            default - default value in case chainID is considered invalid
            """
            invalid_length_ids = False
            invalid_char_ids = False
            missing_ids = False

            for i, chainid in enumerate(chainids):
                if chainid == "":
                    missing_ids = True
                    chainids[i] = default
                elif len(chainid) > 1:
                    invalid_length_ids = True
                    chainids[i] = default
                elif not chainid.isalnum():
                    invalid_char_ids = True
                    chainids[i] = default

            if invalid_length_ids:
                warnings.warn(
                    "Found chainIDs with invalid length."
                    " Corresponding atoms will use value of '{}'"
                    "".format(default)
                )
            if invalid_char_ids:
                warnings.warn(
                    "Found chainIDs using unnaccepted character."
                    " Corresponding atoms will use value of '{}'"
                    "".format(default)
                )
            if missing_ids:
                warnings.warn(
                    "Found missing chainIDs."
                    " Corresponding atoms will use value of '{}'"
                    "".format(default)
                )
            return chainids

        chainids = validate_chainids(chainids, "X")

        # If reindex == False, we use the atom ids for the serial. We do not
        # want to use a fallback here.
        if not self._reindex:
            try:
                atom_ids = atoms.ids
            except AttributeError:
                raise NoDataError(
                    'The "id" topology attribute is not set. '
                    "Either set the attribute or use reindex=True."
                )
        else:
            atom_ids = np.arange(len(atoms)) + 1

        for i in range(len(atoms)):
            vals = {}
            vals["serial"] = util.ltruncate_int(
                atom_ids[i], 5
            )  # check for overflow here?
            vals["name"] = self._deduce_PDB_atom_name(
                atomnames[i], resnames[i]
            )
            vals["altLoc"] = altlocs[i][:1]
            vals["resName"] = resnames[i][:4]
            vals["resSeq"] = util.ltruncate_int(resids[i], 4)
            vals["iCode"] = icodes[i][:1]
            vals["pos"] = pos[i]  # don't take off atom so conversion works
            vals["occupancy"] = occupancies[i]
            vals["tempFactor"] = tempfactors[i]
            vals["segID"] = segids[i][:4]
            vals["chainID"] = chainids[i]
            vals["element"] = elements[i][:2].upper()
            vals["charge"] = formal_charges[i]

            # record_type attribute, if exists, can be ATOM or HETATM
            try:
                self.pdbfile.write(self.fmt[record_types[i]].format(**vals))
            except KeyError:
                errmsg = (
                    f"Found {record_types[i]} for the record type, but "
                    f"only allowed types are ATOM or HETATM"
                )
                raise ValueError(errmsg) from None

        if multiframe:
            self.ENDMDL()
        self.frames_written += 1

    def HEADER(self, trajectory):
        """Write HEADER_ record.

        .. versionchanged:: 0.20.0
            Strip `trajectory.header` since it can be modified by the user and should be
            sanitized (Issue #2324)
        """
        if not hasattr(trajectory, "header"):
            return
        header = trajectory.header.strip()
        self.pdbfile.write(self.fmt["HEADER"].format(header))

    def TITLE(self, *title):
        """Write TITLE_ record."""
        line = " ".join(title)  # TODO: should do continuation automatically
        self.pdbfile.write(self.fmt["TITLE"].format(line))

    def REMARK(self, *remarks):
        """Write generic REMARKS_ record (without number).

        Each string provided in *remarks* is written as a separate REMARKS_
        record.

        """
        for remark in remarks:
            self.pdbfile.write(self.fmt["REMARK"].format(remark))

    def COMPND(self, trajectory):
        if not hasattr(trajectory, "compound"):
            return
        compound = trajectory.compound
        for c in compound:
            self.pdbfile.write(self.fmt["COMPND"].format(c))

    def CRYST1(self, dimensions, spacegroup="P 1", zvalue=1):
        """Write CRYST1_ record."""
        self.pdbfile.write(
            self.fmt["CRYST1"].format(
                box=dimensions[:3],
                ang=dimensions[3:],
                spacegroup=spacegroup,
                zvalue=zvalue,
            )
        )

    def MODEL(self, modelnumber):
        """Write the MODEL_ record.

        .. note::

           The maximum MODEL number is limited to 9999 in the PDB
           standard (i.e., 4 digits). If frame numbers are larger than
           9999, they will wrap around, i.e., 9998, 9999, 0, 1, 2, ...

        .. versionchanged:: 0.19.0
           Maximum model number is enforced.

        """
        self.pdbfile.write(
            self.fmt["MODEL"].format(int(str(modelnumber)[-4:]))
        )

    def END(self):
        """Write END_ record.

        Only a single END record is written. Calling END multiple times has no
        effect. Because :meth:`~PDBWriter.close` also calls this
        method right before closing the file it is recommended to *not* call
        :meth:`~PDBWriter.END` explicitly.

        """
        if not self.has_END:
            # only write a single END record
            self.pdbfile.write(self.fmt["END"])
        self.has_END = True

    def ENDMDL(self):
        """Write the ENDMDL_ record."""
        self.pdbfile.write(self.fmt["ENDMDL"])

    def CONECT(self, conect):
        """Write CONECT_ record."""
        conect = ["{0:5d}".format(entry) for entry in conect]
        conect = "".join(conect)
        self.pdbfile.write(self.fmt["CONECT"].format(conect))


class ExtendedPDBReader(PDBReader):
    """PDBReader that reads a PDB-formatted file with five-digit residue numbers.

    This reader does not conform to the `PDB 3.3 standard`_ because it allows
    five-digit residue numbers that may take up columns 23 to 27 (inclusive)
    instead of being confined to 23-26 (with column 27 being reserved for the
    insertion code in the PDB standard). PDB files in this format are written
    by popular programs such as VMD_.

    See Also
    --------
    :class:`PDBReader`


    .. _VMD: http://www.ks.uiuc.edu/Research/vmd/

    .. versionadded:: 0.8
    """

    format = "XPDB"


class MultiPDBWriter(PDBWriter):
    """PDB writer that implements a subset of the `PDB 3.3 standard`_ .

    PDB format as used by NAMD/CHARMM: 4-letter resnames and segID, altLoc
    is written.

    By default, :class:`MultiPDBWriter` writes a PDB "movie" (multi frame mode,
    *multiframe* = ``True``), consisting of multiple models (using the MODEL_
    and ENDMDL_ records).


    .. _MODEL: http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#MODEL
    .. _ENDMDL: http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ENDMDL
    .. _CONECT: http://www.wwpdb.org/documentation/file-format-content/format33/sect10.html#CONECT


    See Also
    --------
    This class is identical to :class:`PDBWriter` with the one
    exception that it defaults to writing multi-frame PDB files instead of
    single frames.


    .. versionadded:: 0.7.6

    """

    format = ["PDB", "ENT"]
    multiframe = True  # For Writer registration
    singleframe = False
