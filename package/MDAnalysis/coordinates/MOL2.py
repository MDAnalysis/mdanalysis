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

"""
MOL2 file format --- :mod:`MDAnalysis.coordinates.MOL2`
========================================================

Classes to read Tripos_ molecule structure format (MOL2_)
coordinate and topology files. Used by the DOCK_ docking
code.

Examples
~~~~~~~~

To open a mol2, remove all hydrogens and save as a new file, use the following::

  u = Universe("MDAnalysis/testsuite/MDAnalysisTests/data/mol2/Molecule.mol2")
  gr = u.selectAtoms("not name H*")
  print len(u.atoms), len(gr)
  gr.write("Molecule_noh.mol2")

.. _MOL2: http://www.tripos.com/data/support/mol2.pdf
.. _Tripos: http://www.tripos.com/
.. _DOCK: http://dock.compbio.ucsf.edu/
"""

import numpy as np

import base
from MDAnalysis.core.AtomGroup import Universe, AtomGroup
from MDAnalysis.topology.core import guess_atom_element, guess_bonds
from .. import core
import MDAnalysis.core.util as util


class Timestep(base.Timestep):
    @property
    def dimensions(self):
        """unitcell dimensions (`A, B, C, alpha, beta, gamma`)

        MOL2 does not contain unitcell information but for
        compatibility, an empty unitcell is provided.

        - `A, B, C` are the lengths of the primitive cell vectors `e1, e2, e3`
        - `alpha` = angle(`e1, e2`)
        - `beta` = angle(`e1, e3`)
        - `gamma` = angle(`e2, e3`)

        """
        # Layout of unitcell is [A,B,C,90,90,90] with the primitive cell vectors
        return self._unitcell

    @dimensions.setter
    def dimensions(self, box):
        self._unitcell = box


class MOL2Reader(base.Reader):
    format = 'MOL2'
    units = {'time': None, 'length': 'Angstrom'}
    _Timestep = Timestep

    def __init__(self, filename, convert_units=None, **kwargs):
        """Read coordinates from *filename*.

        *filename* can be a gzipped or bzip2ed compressed PDB file.

        If the pdb file contains multiple MODEL records then it is
        read as a trajectory where the MODEL numbers correspond to
        frame numbers. Therefore, the MODEL numbers must be a sequence
        of integers (typically starting at 1 or 0).
        """
        self.filename = filename
        if convert_units is None:
            convert_units = core.flags['convert_lengths']
        self.convert_units = convert_units  # convert length and time to base units

        # = NOTE to clear up confusion over 0-based vs 1-based frame numbering =
        # self.frame is 1-based for this Reader, which matches the behavior of
        # the MODEL record in a typical multi-model PDB file.  If the MODEL
        # record is 0-based, this is accounted for by __init__.
        # self._read_frame assumes that it is passed a 0-based frame number, so
        # that it functions as expected when slicing is used.

        blocks = []

        with util.openany(filename) as f:
            for i, line in enumerate(f):
                # found new molecules
                if "@<TRIPOS>MOLECULE" in line:
                    blocks.append({"start_line": i, "lines": []})
                blocks[-1]["lines"].append(line)

        block = blocks[0]

        sections, coords = self.parse_block(block)

        self.numatoms = len(coords)
        self.ts = self._Timestep(np.array(coords, dtype=np.float32))
        self.ts.frame = 1  # 1-based frame number as starting frame

        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)  # in-place !
            self.convert_pos_from_native(self.ts._unitcell[:3])  # in-place ! (only lengths)

        self.molecule = {}
        self.substructure = {}
        self.frames = blocks
        self.numframes = len(blocks)
        self.fixed = 0
        self.skip = 1
        self.periodic = False
        self.delta = 0
        self.skip_timestep = 1

    def next(self):
        """Read the next time step."""
        return self._read_next_timestep()

    def __iter__(self):
        for i in xrange(0, self.numframes):
            try:
                yield self._read_frame(i)
            except IOError:
                raise StopIteration

    def parse_block(self, block):
        sections = {}
        cursor = None
        for line in block["lines"]:
            if "@<TRIPOS>" in line:
                cursor = line.split("@<TRIPOS>")[1].strip().lower()
                sections[cursor] = []
                continue
            elif line.startswith("#") or line == "\n":
                continue
            sections[cursor].append(line)

        atom_lines, bond_lines = sections["atom"], sections["bond"]
        if not len(atom_lines):
            raise Exception("The mol2 (starting at line {}) block has no atoms".format(block["start_line"]))
        if not len(bond_lines):
            raise Exception("The mol2 (starting at line {}) block has no bonds".format(block["start_line"]))

        coords = []
        for a in atom_lines:
            aid, name, x, y, z, atom_type, resid, resname, charge = a.split()
            x, y, z = float(x), float(y), float(z)
            coords.append((x, y, z))
        coords = np.array(coords)
        return sections, coords

    def _read_next_timestep(self, ts=None):
        if ts is None:
            ts = self.ts
        else:
            # TODO: cleanup _read_frame() to use a "free" Timestep
            raise NotImplementedError("PrimitivePDBReader cannot assign to a timestep")
        # frame is 1-based. Normally would add 1 to frame before calling
        # self._read_frame to retrieve the subsequent ts. But self._read_frame
        # assumes it is being passed a 0-based frame, and adjusts.
        frame = self.frame
        return self._read_frame(frame)

    def _read_frame(self, frame):
        if np.dtype(type(frame)) != np.dtype(int):
            raise TypeError("frame must be a integer")

        unitcell = np.zeros(6, dtype=np.float32)
        block = self.frames[frame]

        sections, coords = self.parse_block(block)

        self.molecule[frame] = sections["molecule"]
        self.substructure[frame] = sections["substructure"]

        # check if atom number changed
        if len(coords) != len(self.ts._pos):
            raise ValueError(
                "PrimitivePDBReader assumes that the number of atoms remains unchanged between frames; the current "
                "frame has %d, the next frame has %d atoms" % (
                len(self.ts._pos), len(coords)))

        self.ts = self._Timestep(np.array(coords, dtype=np.float32))
        self.ts._unitcell[:] = unitcell
        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)  # in-place !
            self.convert_pos_from_native(self.ts._unitcell[:3])  # in-place ! (only lengths)
        self.ts.frame = frame
        return self.ts


class MOL2Writer(base.Writer):
    """MOL2 writer

    * MOL2 Format Specification:  http://www.tripos.com/data/support/mol2.pdf
    * Example file::

        @<TRIPOS>MOLECULE
        generated by VMD
           94  232    1    0    0
        SMALL
        NO_CHARGES
        ****
        Energy = 0

        @<TRIPOS>ATOM
              1  C1'       -7.7930   17.1890   22.4770  C1'  400  RUN         0.000000
              2  O1'       -7.8030   18.3960   22.1770  O1'  400  RUN         0.000000
              3  C2D       -7.6940   16.1280   21.3860  C2D  400  RUN         0.000000
              4  C3'       -6.5780   16.4880   20.4110  C3'  400  RUN         0.000000
              5  C4D       -9.0400   16.1410   20.6690  C4D  400  RUN         0.000000
              6  N1'       -7.9670   17.8950   24.6830  N1'  400  RUN         0.000000
              7  C2'       -8.0280   17.2320   25.8340  C2'  400  RUN         0.000000
              ...
        @<TRIPOS>BOND
            1     1     2  1
            2     1     3  1
            3     1    10  1
            4     1    67  1
            5     1    68  1
            ...
        @<TRIPOS>SUBSTRUCTURE
        1 ****        1 TEMP                        0 ****  **** 0 ROOT
    """
    format = 'MOL2'
    units = {'time': None, 'length': 'Angstrom'}

    def __init__(self, filename, numatoms=None, start=0, step=1,
                 convert_units=None, multiframe=None):
        """Create a new PDBWriter

        :Arguments:
         *filename*
           name of output file
         *start*
           starting timestep
         *step*
           skip between subsequent timesteps
         *convert_units*
           units are converted to the MDAnalysis base format; ``None`` selects
           the value of :data:`MDAnalysis.core.flags` ['convert_lengths']

        """
        self.filename = filename
        if convert_units is None:
            convert_units = core.flags['convert_lengths']
        self.convert_units = convert_units  # convert length and time to base units

        self.frames_written = 0
        if start < 0:
            raise ValueError("'Start' must be a positive value")

        self.start = start
        self.step = step

        self.file = util.anyopen(self.filename, 'w')  # open file on init

    def close(self):
        self.file.close()

    __del__ = close

    def encode_block(self, obj):
        traj = obj.universe.trajectory

        # Make sure Universe has loaded bonds
        obj.universe.bonds

        bonds = set()
        for a in obj.atoms:
            bonds.update(a.bonds)
        bonds = sorted([(bond[0].id, bond[1].id, bond.order) for bond in bonds])
        mapping = dict([(a.id, i) for i, a in enumerate(obj.atoms)])

        atom_lines = ["{:>4} {:>4} {:>13.4f} {:>9.4f} {:>9.4f} {:>4} {} {} "
                      "{:>7.4f}".format(mapping[a.id] + 1, a.name,
                                        float(a.pos[0]),
                                        float(a.pos[1]),
                                        float(a.pos[2]), a.type,
                                        a.resid, a.resname,
                                        a.charge) for a in obj.atoms]
        atom_lines = ["@<TRIPOS>ATOM"] + atom_lines + ["\n"]
        atom_lines = "\n".join(atom_lines)

        bonds = [(atom1, atom2, order) for atom1, atom2, order in bonds if atom1 in mapping and atom2 in mapping]
        bond_lines = ["{:>5} {:>5} {:>5} {:>2}".format(bid + 1,
                                                       mapping[atom1] + 1,
                                                       mapping[atom2] + 1,
                                                       order) for bid, (atom1, atom2, order) in enumerate(bonds)]
        bond_lines = ["@<TRIPOS>BOND"] + bond_lines + ["\n"]
        bond_lines = "\n".join(bond_lines)

        substructure = traj.substructure[traj.frame]
        substructure = ["@<TRIPOS>SUBSTRUCTURE\n"] + substructure

        molecule = traj.molecule[traj.frame]
        check_sums = molecule[1].split()
        check_sums[0], check_sums[1] = str(len(obj.atoms)), str(len(bonds))
        molecule[1] = "{}\n".format(" ".join(check_sums))
        molecule = ["@<TRIPOS>MOLECULE\n"] + molecule

        return "".join(molecule) + atom_lines + bond_lines + "".join(substructure)

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

        start, step = self.start, self.step
        traj = obj.universe.trajectory

        # Start from trajectory[0]/frame 1, if there are more than 1 frame.
        # If there is onyl 1 frame, the traj.frames is not like a python list:
        # accessing trajectory[-1] raises key error.
        if not start and traj.numframes > 1:
            start = traj.frame - 1

        for framenumber in xrange(start, len(traj), step):
            traj[framenumber]
            self.write_next_timestep(obj)

        self.close()

        # Set the trajectory to the starting position
        traj[start]

    def write_next_timestep(self, obj):
        """Write a new frame to the MOL2 file.

        *obj* can be a :class:`~MDAnalysis.core.AtomGroup.AtomGroup`
        or a :class:`~MDAnalysis.core.AtomGroup.Universe`.
        """
        block = self.encode_block(obj)
        self.file.writelines(block)
