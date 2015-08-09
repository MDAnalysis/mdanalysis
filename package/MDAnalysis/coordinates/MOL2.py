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
MOL2 file format --- :mod:`MDAnalysis.coordinates.MOL2`
========================================================

Classes to read Tripos_ molecule structure format (MOL2_)
coordinate and topology files. Used by the DOCK_ docking
code.

Examples
~~~~~~~~

To open a mol2, remove all hydrogens and save as a new file, use the following::

  u = Universe("MDAnalysis/testsuite/MDAnalysisTests/data/mol2/Molecule.mol2")
  gr = u.select_atoms("not name H*")
  print len(u.atoms), len(gr)
  gr.write("Molecule_noh.mol2")

.. _MOL2: http://www.tripos.com/data/support/mol2.pdf
.. _Tripos: http://www.tripos.com/
.. _DOCK: http://dock.compbio.ucsf.edu/
"""

import numpy as np

from . import base
from ..core import flags
from ..lib import util


class MOL2Reader(base.Reader):
    """Reader for MOL2 structure format.

    .. versionchanged:: 0.11.0
       Frames now 0-based instead of 1-based
       MOL2 now resuses the same Timestep object for every frame
       previously created a new instance of Timestep each frame
    """
    format = 'MOL2'
    units = {'time': None, 'length': 'Angstrom'}

    def __init__(self, filename, **kwargs):
        """Read coordinates from *filename*."""
        super(MOL2Reader, self).__init__(filename, **kwargs)

        blocks = []

        with util.openany(filename) as f:
            for i, line in enumerate(f):
                # found new molecules
                if "@<TRIPOS>MOLECULE" in line:
                    blocks.append({"start_line": i, "lines": []})
                blocks[-1]["lines"].append(line)

        block = blocks[0]

        sections, coords = self.parse_block(block)

        self.n_atoms = len(coords)
        self.ts = self._Timestep.from_coordinates(np.array(coords, dtype=np.float32),
                                                  **self._ts_kwargs)
        self.ts.frame = 0  # 0-based frame number as starting frame

        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)  # in-place !
            self.convert_pos_from_native(self.ts._unitcell[:3])  # in-place ! (only lengths)

        self.molecule = {}
        self.substructure = {}
        self.frames = blocks
        self.n_frames = len(blocks)

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
            raise Exception("The mol2 (starting at line {0}) block has no atoms"
                            "".format(block["start_line"]))
        if not len(bond_lines):
            raise Exception("The mol2 (starting at line {0}) block has no bonds"
                            "".format(block["start_line"]))

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
            raise NotImplementedError("PrimitiveMOL2Reader cannot assign to a timestep")
        frame = self.frame + 1
        return self._read_frame(frame)

    def _read_frame(self, frame):
        unitcell = np.zeros(6, dtype=np.float32)
        block = self.frames[frame]

        sections, coords = self.parse_block(block)

        self.molecule[frame] = sections["molecule"]
        self.substructure[frame] = sections["substructure"]

        # check if atom number changed
        if len(coords) != len(self.ts._pos):
            raise ValueError(
                "PrimitiveMOL2Reader assumes that the number of atoms remains unchanged between frames; the current "
                "frame has %d, the next frame has %d atoms" % (
                len(self.ts._pos), len(coords)))

        self.ts.positions = np.array(coords, dtype=np.float32)
        self.ts._unitcell[:] = unitcell
        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)  # in-place !
            self.convert_pos_from_native(self.ts._unitcell[:3])  # in-place ! (only lengths)
        self.ts.frame = frame
        return self.ts

    def _reopen(self):
        # Make frame think it's before start, so calling next
        # reads first frame
        self.ts.frame = -1


class MOL2Writer(base.Writer):
    """MOL2Writer Limitations ----------- MOL2Writer can only be used to write out previously loaded MOL2 files.  If you're trying to convert, for example, a PDB file to MOL you should use other tools, like rdkit (http://www.rdkit.org/docs/GettingStartedInPython.html).  Here is an example how to use rdkit to convert a PDB to MOL:: from rdkit import Chem mol = Chem.MolFromPDBFile("molecule.pdb", removeHs=False) Chem.MolToMolFile(mol, "molecule.mol" ) MOL2 writer is currently not available for rdkit master. It requires 
    SYBYL atomtype generation.
    See page 7 for list of SYBYL atomtypes (http://tripos.com/tripos_resources/fileroot/pdfs/mol2_format2.pdf).

    * MOL2 Format Specification:  (http://www.tripos.com/data/support/mol2.pdf)
    * Example file (http://www.tripos.com/mol2/mol2_format3.html)::

              #    Name: benzene 
              #    Creating user name: tom 
              #    Creation time: Wed Dec 28 00:18:30 1988 
            
              #    Modifying user name: tom 
              #    Modification time: Wed Dec 28 00:18:30 1988 
            
              @<TRIPOS>MOLECULE 
              benzene 
              12 12 1  0   0 
              SMALL 
              NO_CHARGES 
             
             
              @<TRIPOS>ATOM 
              1   C1  1.207   2.091   0.000   C.ar    1   BENZENE 0.000 
              2   C2  2.414   1.394   0.000   C.ar    1   BENZENE 0.000 
              3   C3  2.414   0.000   0.000   C.ar    1   BENZENE 0.000 
              4   C4  1.207   -0.697  0.000   C.ar    1   BENZENE 0.000 
              5   C5  0.000   0.000   0.000   C.ar    1   BENZENE 0.000 
              6   C6  0.000   1.394   0.000   C.ar    1   BENZENE 0.000 
              7   H1  1.207   3.175   0.000   H   1   BENZENE 0.000 
              8   H2  3.353   1.936   0.000   H   1   BENZENE 0.000 
              9   H3  3.353   -0.542  0.000   H   1   BENZENE 0.000 
              10  H4  1.207   -1.781  0.000   H   1   BENZENE 0.000 
              11  H5  -0.939  -0.542  0.000   H   1   BENZENE 0.000 
              12  H6  -0.939  1.936   0.000   H   1   BENZENE 0.000 
              @<TRIPOS>BOND 
              1   1   2   ar 
              2   1   6   ar 
              3   2   3   ar 
              4   3   4   ar 
              5   4   5   ar 
              6   5   6   ar 
              7   1   7   1 
              8   2   8   1 
              9   3   9   1 
              10  4   10  1 
              11  5   11  1 
              12  6   12  1 
             @<TRIPOS>SUBSTRUCTURE 
              1   BENZENE 1   PERM    0   ****    ****    0   ROOT

    .. versionchanged:: 0.11.0
       Frames now 0-based instead of 1-based

    """
    format = 'MOL2'
    units = {'time': None, 'length': 'Angstrom'}

    def __init__(self, filename, n_atoms=None, start=0, step=1,
                 convert_units=None, multiframe=None):
        """Create a new MOL2Writer

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
            convert_units = flags['convert_lengths']
        self.convert_units = convert_units  # convert length and time to base units

        self.frames_written = 0
        if start < 0:
            raise ValueError("'Start' must be a positive value")

        self.start = start
        self.step = step

        self.file = util.anyopen(self.filename, 'w')  # open file on init

    def close(self):
        self.file.close()

    def encode_block(self, obj):
        traj = obj.universe.trajectory

        # Make sure Universe has loaded bonds
        obj.universe.bonds

        bonds = set()
        for a in obj.atoms:
            bonds.update(a.bonds)
        bonds = sorted([(bond[0].id, bond[1].id, bond.order) for bond in bonds])
        mapping = dict([(a.id, i) for i, a in enumerate(obj.atoms)])

        atom_lines = ["{0:>4} {1:>4} {2:>13.4f} {3:>9.4f} {4:>9.4f} {5:>4} {6} {7} "
                      "{8:>7.4f}".format(mapping[a.id] + 1, a.name,
                                        float(a.pos[0]),
                                        float(a.pos[1]),
                                        float(a.pos[2]), a.type,
                                        a.resid, a.resname,
                                        a.charge) for a in obj.atoms]
        atom_lines = ["@<TRIPOS>ATOM"] + atom_lines + ["\n"]
        atom_lines = "\n".join(atom_lines)

        bonds = [(atom1, atom2, order) for atom1, atom2, order in bonds if atom1 in mapping and atom2 in mapping]
        bond_lines = ["{0:>5} {1:>5} {2:>5} {3:>2}"
                      "".format(bid + 1,
                                mapping[atom1] + 1,
                                mapping[atom2] + 1,
                                order) for bid, (atom1, atom2, order) in enumerate(bonds)]
        bond_lines = ["@<TRIPOS>BOND"] + bond_lines + ["\n"]
        bond_lines = "\n".join(bond_lines)

        try:
            substructure = traj.substructure[traj.frame]
        except AttributeError:
            raise NotImplementedError("No MOL2 substructure type found in traj")

        substructure = ["@<TRIPOS>SUBSTRUCTURE\n"] + substructure

        molecule = traj.molecule[traj.frame]
        check_sums = molecule[1].split()
        check_sums[0], check_sums[1] = str(len(obj.atoms)), str(len(bonds))
        molecule[1] = "{0}\n".format(" ".join(check_sums))
        molecule = ["@<TRIPOS>MOLECULE\n"] + molecule

        return "".join(molecule) + atom_lines + bond_lines + "".join(substructure)

    def write(self, obj):
        """Write object *obj* at current trajectory frame to file.

        *obj* can be a selection (i.e. a
        :class:`~MDAnalysis.core.AtomGroup.AtomGroup`) or a whole
        :class:`~MDAnalysis.core.AtomGroup.Universe`.

        :Arguments:
          *obj*
            :class:`~MDAnalysis.core.AtomGroup.AtomGroup` or
            :class:`~MDAnalysis.core.AtomGroup.Universe`
        """

        start, step = self.start, self.step
        traj = obj.universe.trajectory

        # Start from trajectory[0]/frame 1, if there are more than 1 frame.
        # If there is only 1 frame, the traj.frames is not like a python list:
        # accessing trajectory[-1] raises key error.
        if not start and traj.n_frames > 1:
            start = traj.frame

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
