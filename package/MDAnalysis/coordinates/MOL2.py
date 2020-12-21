# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
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

"""MOL2 file format --- :mod:`MDAnalysis.coordinates.MOL2`
========================================================

Classes to work with Tripos_ molecule structure format (MOL2_) coordinate and
topology files. Used, for instance, by the DOCK_ docking code.


Example for working with mol2 files
-----------------------------------

To open a mol2, remove all hydrogens and save as a new file, use the following::

  u = Universe("Molecule.mol2")
  gr = u.select_atoms("not name H*")
  print(len(u.atoms), len(gr))
  gr.write("Molecule_noh.mol2")

.. _MOL2: http://www.tripos.com/data/support/mol2.pdf
.. _Tripos: http://www.tripos.com/
.. _DOCK: http://dock.compbio.ucsf.edu/


See Also
--------
rdkit: rdkit_ is an open source cheminformatics toolkit with more mol2
       functionality

.. _rdkit: http://www.rdkit.org/docs/GettingStartedInPython.html


Classes
-------

.. autoclass:: MOL2Reader
   :members:

.. autoclass:: MOL2Writer
   :members:


MOL2 format notes
-----------------

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

"""
import numpy as np

from . import base
from ..lib import util


class MOL2Reader(base.ReaderBase):
    """Reader for MOL2 structure format.

    .. versionchanged:: 0.11.0
       Frames now 0-based instead of 1-based.
       MOL2 now reuses the same Timestep object for every frame,
       previously created a new instance of Timestep each frame.
    .. versionchanged:: 0.20.0
       Allows for comments at top of file.
       Ignores status bit strings
    .. versionchanged:: 2.0.0
       Bonds attribute is not added if no bonds are present in MOL2 file
    """
    format = 'MOL2'
    units = {'time': None, 'length': 'Angstrom'}

    def __init__(self, filename, **kwargs):
        """Read coordinates from `filename`.


        Parameters
        ----------
        filename : str or NamedStream
             name of the mol2 file or stream
        """
        super(MOL2Reader, self).__init__(filename, **kwargs)

        self.n_atoms = None

        blocks = []

        with util.openany(filename) as f:
            for i, line in enumerate(f):
                # found new molecules
                if "@<TRIPOS>MOLECULE" in line:
                    blocks.append({"start_line": i, "lines": []})
                if len(blocks):
                    blocks[-1]["lines"].append(line)
        self.n_frames = len(blocks)
        self.frames = blocks

        sections, coords = self.parse_block(blocks[0])
        self.n_atoms = len(coords)

        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)

        self.ts = self._read_frame(0)

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

        atom_lines = sections["atom"]
        if not len(atom_lines):
            raise Exception("The mol2 (starting at line {0}) block has no atoms"
                            "".format(block["start_line"]))
        elif self.n_atoms is None:
            # First time round, remember the number of atoms
            self.n_atoms = len(atom_lines)
        elif len(atom_lines) != self.n_atoms:
            raise ValueError(
                "MOL2Reader assumes that the number of atoms remains unchanged"
                " between frames; the current "
                "frame has {0}, the next frame has {1} atoms"
                "".format(self.n_atoms, len(atom_lines)))

        coords = np.zeros((self.n_atoms, 3), dtype=np.float32)
        for i, a in enumerate(atom_lines):
            aid, name, x, y, z, atom_type, resid, resname, charge = a.split()[:9]

            #x, y, z = float(x), float(y), float(z)
            coords[i, :] = x, y, z

        return sections, coords

    def _read_next_timestep(self, ts=None):
        if ts is None:
            ts = self.ts
        else:
            # TODO: cleanup _read_frame() to use a "free" Timestep
            raise NotImplementedError("MOL2Reader cannot assign to a timestep")
        frame = self.frame + 1
        return self._read_frame(frame)

    def _read_frame(self, frame):
        unitcell = np.zeros(6, dtype=np.float32)
        try:
            block = self.frames[frame]
        except IndexError:
            errmsg = (f"Invalid frame {frame} for trajectory with length "
                      f"{len(self)}")
            raise IOError(errmsg) from None

        sections, coords = self.parse_block(block)

        for sect in ['molecule', 'substructure']:
            try:
                self.ts.data[sect] = sections[sect]
            except KeyError:
                pass

        self.ts.positions = np.array(coords, dtype=np.float32)
        self.ts.unitcell = unitcell
        if self.convert_units:
            # in-place !
            self.convert_pos_from_native(self.ts._pos)
            self.convert_pos_from_native(self.ts._unitcell[:3])
        self.ts.frame = frame

        return self.ts

    def _reopen(self):
        # Make frame think it's before start, so calling next
        # reads first frame
        self.ts.frame = -1


class MOL2Writer(base.WriterBase):
    """mol2 format writer

    Write a file in the Tripos_ molecule structure format (MOL2_).

    Note
    ----
    :class:`MOL2Writer` can only be used to write out previously loaded MOL2 files.
    If you're trying to convert, for example, a PDB file to MOL you should
    use other tools, like rdkit_.

    Here is an example how to use rdkit_ to convert a PDB to MOL::

      from rdkit import Chem
      mol = Chem.MolFromPDBFile("molecule.pdb", removeHs=False)
      Chem.MolToMolFile(mol, "molecule.mol" )

    MOL2 writer is currently not available for rdkit master. It requires SYBYL
    atomtype generation.
    See page 7 for list of SYBYL atomtypes
    (http://tripos.com/tripos_resources/fileroot/pdfs/mol2_format2.pdf).


    .. _rdkit: http://www.rdkit.org/docs/GettingStartedInPython.html

    .. versionchanged:: 0.11.0
       Frames now 0-based instead of 1-based

    """
    format = 'MOL2'
    multiframe = True
    units = {'time': None, 'length': 'Angstrom'}

    def __init__(self, filename, n_atoms=None, convert_units=True):
        """Create a new MOL2Writer

        Parameters
        ----------
        filename: str
            name of output file
        convert_units: bool (optional)
            units are converted to the MDAnalysis base format; [``True``]
        """
        self.filename = filename
        self.convert_units = convert_units  # convert length and time to base units

        self.frames_written = 0

        self.file = util.anyopen(self.filename, 'w')  # open file on init

    def close(self):
        self.file.close()

    def encode_block(self, obj):
        """
        Parameters
        ----------
        obj : AtomGroup or Universe
        """
        # Issue 2717
        try:
            obj = obj.atoms
        except AttributeError:
            errmsg = "Input obj is neither an AtomGroup or Universe"
            raise TypeError(errmsg) from None

        traj = obj.universe.trajectory
        ts = traj.ts

        try:
            molecule = ts.data['molecule']
        except KeyError:
            errmsg = "MOL2Writer cannot currently write non MOL2 data"
            raise NotImplementedError(errmsg) from None

        # Need to remap atom indices to 1 based in this selection
        mapping = {a: i for i, a in enumerate(obj.atoms, start=1)}

        # only write bonds if the Bonds attribute exists (Issue #3057)
        if hasattr(obj, "bonds"):
            # Grab only bonds between atoms in the obj
            # ie none that extend out of it
            bondgroup = obj.bonds.atomgroup_intersection(obj, strict=True)
            bonds = sorted((b[0], b[1], b.order) for b in bondgroup)
            bond_lines = ["@<TRIPOS>BOND"]
            bls = ["{0:>5} {1:>5} {2:>5} {3:>2}".format(bid,
                                                        mapping[atom1],
                                                        mapping[atom2],
                                                        order)
                   for bid, (atom1, atom2, order) in enumerate(bonds, start=1)]
            bond_lines.extend(bls)
            bond_lines.append("\n")
            bond_lines = "\n".join(bond_lines)
        else:
            bondgroup = []
            bond_lines = ""

        atom_lines = ["@<TRIPOS>ATOM"]
        atom_lines.extend("{0:>4} {1:>4} {2:>13.4f} {3:>9.4f} {4:>9.4f}"
                          " {5:>4} {6} {7} {8:>7.4f}"
                          "".format(mapping[a],
                                    a.name,
                                    a.position[0],
                                    a.position[1],
                                    a.position[2],
                                    a.type,
                                    a.resid,
                                    a.resname,
                                    a.charge)
                          for a in obj.atoms)
        atom_lines.append("\n")
        atom_lines = "\n".join(atom_lines)

        try:
            substructure = ["@<TRIPOS>SUBSTRUCTURE\n"] + ts.data['substructure']
        except KeyError:
            substructure = ""

        check_sums = molecule[1].split()
        check_sums[0], check_sums[1] = str(len(obj.atoms)), str(len(bondgroup))

        # prevent behavior change between repeated calls
        # see gh-2678
        molecule_0_store = molecule[0]
        molecule_1_store = molecule[1]

        molecule[1] = "{0}\n".format(" ".join(check_sums))
        molecule.insert(0, "@<TRIPOS>MOLECULE\n")

        return_val = ("".join(molecule) + atom_lines +
                      bond_lines + "".join(substructure))

        molecule[0] = molecule_0_store
        molecule[1] = molecule_1_store
        return return_val

    def _write_next_frame(self, obj):
        """Write a new frame to the MOL2 file.

        Parameters
        ----------
        obj : AtomGroup or Universe


        .. versionchanged:: 1.0.0
            Renamed from `write_next_timestep` to `_write_next_frame`.
        """
        block = self.encode_block(obj)
        self.file.writelines(block)
