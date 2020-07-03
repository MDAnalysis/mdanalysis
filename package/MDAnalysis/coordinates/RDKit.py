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

"""RDKit molecule --- :mod:`MDAnalysis.coordinates.RDKit`
================================================================

Read coordinates data from an `RDKit <https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.Mol>`_ :class:`rdkit.Chem.rdchem.Mol` with :class:`RDKitReader` 
into a MDAnalysis Universe. Convert it back to a :class:`rdkit.Chem.rdchem.Mol` with 
:class:`RDKitConverter`.


Example
-------

>>> from rdkit import Chem
>>> import MDAnalysis as mda
>>> mol = Chem.MolFromMol2File("docking_poses.mol2", removeHs=False)
>>> u = mda.Universe(mol)
>>> u
<Universe with 42 atoms>
>>> u.trajectory
<RDKitReader with 10 frames of 42 atoms>
>>> u.atoms.convert_to("RDKIT")
<rdkit.Chem.rdchem.RWMol object at 0x7fcebb958148>


Classes
-------

.. autoclass:: RDKitReader
   :members:

.. autoclass:: RDKitConverter
   :members:


"""

import warnings
import re

import numpy as np

from ..exceptions import NoDataError
from ..topology.guessers import guess_atom_element
from ..core.topologyattrs import _TOPOLOGY_ATTRS
from . import memory
from . import base

try:
    from rdkit import Chem
except ImportError:
    pass
else:
    RDBONDTYPE = {
        'AROMATIC': Chem.BondType.AROMATIC,
        'SINGLE': Chem.BondType.SINGLE,
        'DOUBLE': Chem.BondType.DOUBLE,
        'TRIPLE': Chem.BondType.TRIPLE,
    }
    RDBONDORDER = {
        1: Chem.BondType.SINGLE,
        1.5: Chem.BondType.AROMATIC,
        "ar": Chem.BondType.AROMATIC,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
    }
    # add string version of the key for each bond
    RDBONDORDER.update({str(key): value for key, value in RDBONDORDER.items()})
    RDATTRIBUTES = {
        "altLocs": "AltLoc",
        "chainIDs": "ChainId",
        "icodes": "InsertionCode",
        "names": "Name",
        "occupancies": "Occupancy",
        "resnames": "ResidueName",
        "resids": "ResidueNumber",
        "segindices": "SegmentNumber",
        "tempfactors": "TempFactor",
    }
    PERIODIC_TABLE = Chem.GetPeriodicTable()


class RDKitReader(memory.MemoryReader):
    """Coordinate reader for RDKit.

    .. versionadded:: 2.0.0
    """
    format = 'RDKIT'

    # Structure.coordinates always in Angstrom
    units = {'time': None, 'length': 'Angstrom'}

    @staticmethod
    def _format_hint(thing):
        """Can this reader read *thing*?"""
        try:
            from rdkit import Chem
        except ImportError:
            # if we can't import rdkit, it's probably not rdkit
            return False
        else:
            return isinstance(thing, Chem.Mol)

    def __init__(self, filename, **kwargs):
        """Read coordinates from an RDKit molecule.
        Each conformer in the original RDKit molecule will be read as a frame 
        in the resulting universe.

        Parameters
        ----------

        filename : rdkit.Chem.rdchem.Mol
            RDKit molecule
        """
        n_atoms = filename.GetNumAtoms()
        coordinates = np.array([
            conf.GetPositions() for conf in filename.GetConformers()],
            dtype=np.float32)
        if coordinates.size == 0:
            warnings.warn("No coordinates found in the RDKit molecule")
            coordinates = np.empty((1, n_atoms, 3), dtype=np.float32)
            coordinates[:] = np.nan
        super(RDKitReader, self).__init__(coordinates, order='fac', **kwargs)


class RDKitConverter(base.ConverterBase):
    """Convert MDAnalysis AtomGroup or Universe to `RDKit <https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.Mol>`_ :class:`rdkit.Chem.rdchem.Mol`.

    Example
    -------

    .. code-block:: python

        import MDAnalysis as mda
        from MDAnalysis.tests.datafiles import PDB_full
        u = mda.Universe(PDB_full)
        mol = u.select_atoms('resname DMS').convert_to('RDKIT')


    .. versionadded:: 2.0.0
    """

    lib = 'RDKIT'
    units = {'time': None, 'length': 'Angstrom'}

    def convert(self, obj):
        """Write selection at current trajectory frame to :class:`~rdkit.Chem.rdchem.Mol`.

        Parameters
        -----------
        obj : AtomGroup or Universe
        """
        try:
            from rdkit import Chem
        except ImportError:
            raise ImportError("RDKit is required for the RDKitConverter but "
                              "it's not installed. Try installing it with \n"
                              "conda install -c conda-forge rdkit")
        try:
            # make sure to use atoms (Issue 46)
            ag = obj.atoms
        except AttributeError:
            raise TypeError("No `atoms` attribute in object of type {}, "
                            "please use a valid AtomGroup or Universe".format(
                                type(obj))) from None

        try:
            elements = ag.elements
        except NoDataError:
            raise AttributeError(
                "The `elements` attribute is required for the RDKitConverter "
                "but is not present in this AtomGroup. Please refer to the "
                "documentation to guess elements from other attributes or "
                "type `help(mda.topology.guessers)`") from None

        # attributes accepted in PDBResidueInfo object
        pdb_attrs = {}
        for attr in RDATTRIBUTES.keys():
            if hasattr(ag, attr):
                pdb_attrs[attr] = getattr(ag, attr)
        # others
        other_attrs = {}
        for attr in ["bfactors", "charges", "segids", "types"]:
            if hasattr(ag, attr):
                other_attrs[attr] = getattr(ag, attr)

        mol = Chem.RWMol()
        atom_mapper = {}

        for i, (atom, element) in enumerate(zip(ag, elements)):
            # create atom
            rdatom = Chem.Atom(element.capitalize())
            # disable adding H to the molecule
            rdatom.SetNoImplicit(True)
            # add PDB-like properties
            mi = Chem.AtomPDBResidueInfo()
            for attr, values in pdb_attrs.items():
                _add_mda_attr_to_rdkit(attr, values[i], mi)
            rdatom.SetMonomerInfo(mi)
            # other properties
            for attr in other_attrs.keys():
                value = other_attrs[attr][i]
                attr = _TOPOLOGY_ATTRS[attr].singular
                if isinstance(value, np.float):
                    rdatom.SetDoubleProp("_MDAnalysis_%s" % attr, float(value))
                elif isinstance(value, np.int):
                    rdatom.SetIntProp("_MDAnalysis_%s" % attr, int(value))
                else:
                    rdatom.SetProp("_MDAnalysis_%s" % attr, value)
            rdatom.SetIntProp("_MDAnalysis_index", int(atom.ix))
            # add atom
            index = mol.AddAtom(rdatom)
            # map index in universe to index in mol
            atom_mapper[atom.ix] = index

        try:
            bonds = ag.bonds
            if (len(bonds) == 0) and (ag.n_atoms > 1):
                # force guessing bonds
                raise NoDataError
        except NoDataError:
            warnings.warn(
                "No `bonds` attribute in this AtomGroup. Guessing bonds based"
                "on atoms coordinates")
            ag.guess_bonds()
            bonds = ag.bonds

        # only keep bonds where both atoms belong to the AtomGroup
        bonds = bonds.atomgroup_intersection(ag, strict=True)

        for bond in bonds:
            bond_indices = [atom_mapper[i] for i in bond.indices]
            try:
                bond_type = bond.type.upper()
            except AttributeError:
                # bond type can be a tuple for PDB files
                bond_type = None
            bond_type = RDBONDTYPE.get(bond_type, RDBONDORDER.get(
                bond.order, Chem.BondType.SINGLE))
            mol.AddBond(*bond_indices, bond_type)

        # sanitization
        Chem.SanitizeMol(mol)

        return mol


def _add_mda_attr_to_rdkit(attr, value, mi):
    """Converts an MDAnalysis atom attribute into the RDKit equivalent and 
    stores it into an RDKit AtomPDBResidueInfo object.

    Parameters
    ----------

    attr : str
        Name of the atom attribute in MDAnalysis in the singular form
    value : object, np.int or np.float
        Attribute value as found in the AtomGroup
    mi : rdkit.Chem.rdchem.AtomPDBResidueInfo
        MonomerInfo object that will store the relevant atom attributes
    """
    if isinstance(value, np.generic):
        # convert numpy types to python standard types
        value = value.item()
    if attr == "names":
        # RDKit needs the name to be properly formated for a
        # PDB file (1 letter elements start at col 14)
        name = re.findall('(\D+|\d+)', value)
        if len(name) == 2:
            symbol, number = name
        else:
            symbol, number = name[0], ""
        value = "{:>2}{:<2}".format(symbol, number)
    # set attribute value in RDKit MonomerInfo
    rdattr = RDATTRIBUTES[attr]
    getattr(mi, "Set%s" % rdattr)(value)
