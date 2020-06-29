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
            rdatom = Chem.Atom(element)
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
                attr = attr[:-1] # plural to singular
                if isinstance(value, np.float):
                    rdatom.SetDoubleProp("_MDAnalysis_%s" % attr, float(value))
                elif isinstance(value, np.int):
                    rdatom.SetIntProp("_MDAnalysis_%s" % attr, int(value))
                else:
                    rdatom.SetProp("_MDAnalysis_%s" % attr, value)
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

        border_atom_indices = []
        for bond in bonds:
            try:
                bond_indices = [atom_mapper[i] for i in bond.indices]
            except KeyError:
                # one of the atoms of the bond is not part of the atomgroup
                # save the bond atom that is in the atomgroup for later
                for i in bond.indices:
                    if i in atom_mapper.keys():
                        border_atom_indices.append(atom_mapper[i])
                        break
                # skip the rest
                continue
            try:
                bond_type = bond.type.upper()
            except AttributeError:
                # bond type can be a tuple for PDB files
                bond_type = None
            bond_type = RDBONDTYPE.get(bond_type, RDBONDORDER.get(
                bond.order, Chem.BondType.SINGLE))
            mol.AddBond(*bond_indices, bond_type)

        mol.UpdatePropertyCache(strict=False)

        # infer bond orders and formal charges from the connectivity
        _infer_bo_and_charges(mol, border_atom_indices)

        # sanitize
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
    try:  # get value in MDA atom
        value = getattr(atom, attr)
    except AttributeError:
        pass
    else:
        if isinstance(value, np.generic):
            # convert numpy types to python standard types
            value = value.item()
        if attr == "name":
            # RDKit needs the name to be properly formated for a
            # PDB file (1 letter elements start at col 14)
            name = re.findall('(\D+|\d+)', value)
            if len(name) == 2:
                symbol, number = name
            else:
                symbol, number = name[0], ""
            value = "{:>2}".format(symbol) + "{:<2}".format(number)
        # set attribute value in RDKit MonomerInfo
        getattr(mi, "Set%s" % rdattr)(value)


def _infer_bo_and_charges(mol, border_atom_indices):
    """Infer bond orders and formal charges from a molecule.

    - Step 1
    Since most MD topology files don't explicitely retain informations on bond 
    orders or charges, it has to be guessed from the topology. This is done by 
    looping other each atom and comparing its expected valence to the current 
    valence, called `delta_v`. If two neighbouring atoms have a common 
    positive delta_v, the bond between them most likely has a bond order of 
    1+delta_v. If an atom doesn't share a delta_v with any of its neighbours, 
    it likely needs a formal charge of -delta_v.

    - Step 2
    Some atoms can be "mutilated" by a selection (i.e. one of their bonds is cut). The previous step is likely to assign a formal charge to such atoms even if they weren't charged in the original topology. This step converts the resulting charges to radical electrons, or in some cases to higher order bonds. This ensures the atomgroup is not artificially charged because of the previous step.

    Parameters
    ----------

    mol : rdkit.Chem.rdchem.RWMol
        The molecule is modified inplace and must have all hydrogens added

    border_atom_indices : list
        List of border atoms indices
    """
    # Step 1
    for atom in mol.GetAtoms():
        # create delta_v for each possible valence
        expected_vs = PERIODIC_TABLE.GetValenceList(atom.GetAtomicNum())
        current_v = atom.GetTotalValence()
        delta_vs = [expected_v - current_v for expected_v in expected_vs]

        # if there's only one possible valence state and the correpsonding
        # delta_v is negative, it means we can only add a positive charge to
        # the atom
        if (len(delta_vs) == 1) and (delta_vs[0] < 0):
            charge = -delta_vs[0]
            atom.SetFormalCharge(charge)
            mol.UpdatePropertyCache(strict=False)
        else:
            neighbors = atom.GetNeighbors()
            # check if one of the neighbors has a common delta_v
            for i, na in enumerate(neighbors, start=1):
                # create delta_v for the neighbor
                na_expected_vs = PERIODIC_TABLE.GetValenceList(
                    na.GetAtomicNum())
                na_current = na.GetTotalValence()
                na_delta = [
                    na_expected - na_current for na_expected in na_expected_vs]
                # smallest common delta_v, else NaN
                common_delta = min(set(delta_vs).intersection(na_delta),
                                   default=np.nan)
                # common_delta == 0 means we don't need to do anything
                if common_delta != 0:
                    # if they have no delta_v in common
                    if common_delta is np.nan:
                        # if it's the last neighbor
                        if i == len(neighbors):
                            charge = -delta_vs[0]  # negative
                            atom.SetFormalCharge(charge)
                            mol.UpdatePropertyCache(strict=False)
                    # if they both need a supplementary bond
                    else:
                        bond = mol.GetBondBetweenAtoms(
                            atom.GetIdx(), na.GetIdx())
                        bond.SetBondType(RDBONDORDER[common_delta+1])
                        mol.UpdatePropertyCache(strict=False)
                        break  # out of neighbors loop
    
    # Step 2
    for i in border_atom_indices:
        atom = mol.GetAtomWithIdx(i)
        charge = atom.GetFormalCharge()
        neighbors = atom.GetNeighbors()
        # check if a neighbor atom also bears a charge
        for i, na in enumerate(neighbors, 1):
            na_charge = na.GetFormalCharge()
            if na_charge < 0:
                # both atoms have a negative charge 
                # convert to higher order bond
                common_delta = max([charge, na_charge])
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), na.GetIdx())
                bond.SetBondType(RDBONDORDER[-common_delta+1])
                na.SetFormalCharge(na_charge - common_delta)
                atom.SetFormalCharge(0)
                atom.SetNumRadicalElectrons(common_delta - charge)
                break
            elif i == len(neighbors):
                # no neighbor shares a negative charge
                atom.SetNumRadicalElectrons(-atom.GetFormalCharge())
                atom.SetFormalCharge(0)
        mol.UpdatePropertyCache(strict=False)
