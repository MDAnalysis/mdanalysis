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

"""RDKit molecule I/O --- :mod:`MDAnalysis.coordinates.RDKit`
================================================================

Read coordinates data from an `RDKit`_ :class:`rdkit.Chem.rdchem.Mol` with
:class:`RDKitReader` into an MDAnalysis Universe. Convert it back to an
:class:`rdkit.Chem.rdchem.Mol` with :class:`RDKitConverter`.


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

.. _RDKit: https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.Mol


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
    from rdkit.Chem import AllChem
except ImportError:
    pass
else:
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
        "bfactors": "TempFactor",
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
    """Convert MDAnalysis :class:`~MDAnalysis.core.groups.AtomGroup` or
    :class:`~MDAnalysis.core.universe.Universe` to `RDKit`_
    :class:`rdkit.Chem.rdchem.Mol`

    MDanalysis attributes are stored in each RDKit atom of the resulting
    molecule in two different ways:

    * in an `AtomPDBResidueInfo` object available through the
      ``atom.GetMonomerInfo()`` method if it's an attribute that is typically
      found in a PDB file,
    * directly as an atom property available through the
      ``atom.GetPropsAsDict()`` method for the others.

    Supported attributes:

    +-----------------------+-------------------------------------------+
    | MDAnalysis attribute  | RDKit                                     |
    +=======================+===========================================+
    | altLocs               | atom.GetMonomerInfo().GetAltLoc()         |
    +-----------------------+-------------------------------------------+
    | chainIDs              | atom.GetMonomerInfo().GetChainId()        |
    +-----------------------+-------------------------------------------+
    | icodes                | atom.GetMonomerInfo().GetInsertionCode()  |
    +-----------------------+-------------------------------------------+
    | names                 | atom.GetMonomerInfo().GetName()           |
    +-----------------------+-------------------------------------------+
    | occupancies           | atom.GetMonomerInfo().GetOccupancy()      |
    +-----------------------+-------------------------------------------+
    | resnames              | atom.GetMonomerInfo().GetResidueName()    |
    +-----------------------+-------------------------------------------+
    | resids                | atom.GetMonomerInfo().GetResidueNumber()  |
    +-----------------------+-------------------------------------------+
    | segindices            | atom.GetMonomerInfo().GetSegmentNumber()  |
    +-----------------------+-------------------------------------------+
    | tempfactors           | atom.GetMonomerInfo().GetTempFactor()     |
    +-----------------------+-------------------------------------------+
    | bfactors              | atom.GetMonomerInfo().GetTempFactor()     |
    +-----------------------+-------------------------------------------+
    | charges               | atom.GetDoubleProp("_MDAnalysis_charge")  |
    +-----------------------+-------------------------------------------+
    | indices               | atom.GetIntProp("_MDAnalysis_index")      |
    +-----------------------+-------------------------------------------+
    | segids                | atom.GetProp("_MDAnalysis_segid")         |
    +-----------------------+-------------------------------------------+
    | types                 | atom.GetProp("_MDAnalysis_type")          |
    +-----------------------+-------------------------------------------+

    Example
    -------

    .. code-block:: python

        import MDAnalysis as mda
        from MDAnalysis.tests.datafiles import PDB_full
        u = mda.Universe(PDB_full)
        mol = u.select_atoms('resname DMS').convert_to('RDKIT')


    Notes
    -----

    The converter requires the :class:`~MDAnalysis.core.topologyattrs.Elements`
    attribute to be present in the topology, else it will fail.
    It also requires the `bonds` attribute, although they will be automatically
    guessed if not present.
    If both `tempfactors` and `bfactors` attributes are present, the conversion
    will fail, since only one of these should be present.


    .. versionadded:: 2.0.0

    .. _RDKit: https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.Mol
    """

    lib = 'RDKIT'
    units = {'time': None, 'length': 'Angstrom'}

    def convert(self, obj, NoImplicit=True):
        """Write selection at current trajectory frame to
        :class:`rdkit.Chem.rdchem.Mol`.

        Parameters
        -----------
        obj : AtomGroup or Universe

        NoImplicit : bool
            Prevent adding hydrogens to the molecule
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

        if "H" not in ag.elements:
            warnings.warn(
                "No hydrogen atom could be found in the topology, but the "
                "converter requires all hydrogens to be explicit. Please "
                "check carefully the output molecule as the converter is "
                "likely to add negative charges and assign incorrect bond "
                "orders to structures with implicit hydrogens."
            )

        # attributes accepted in PDBResidueInfo object
        pdb_attrs = {}
        if hasattr(ag, "bfactors") and hasattr(ag, "tempfactors"):
            raise AttributeError(
                "Both `tempfactors` and `bfactors` attributes are present but "
                "only one can be assigned to the RDKit molecule. Please "
                "delete the unnecessary one and retry."
            )
        for attr in RDATTRIBUTES.keys():
            if hasattr(ag, attr):
                pdb_attrs[attr] = getattr(ag, attr)

        other_attrs = {}
        for attr in ["charges", "segids", "types"]:
            if hasattr(ag, attr):
                other_attrs[attr] = getattr(ag, attr)

        mol = Chem.RWMol()
        # map index in universe to index in mol
        atom_mapper = {}

        for i, (atom, element) in enumerate(zip(ag, elements)):
            # create atom
            rdatom = Chem.Atom(element.capitalize())
            # disable adding H to the molecule
            rdatom.SetNoImplicit(NoImplicit)
            # add PDB-like properties
            mi = Chem.AtomPDBResidueInfo()
            for attr, values in pdb_attrs.items():
                _add_mda_attr_to_rdkit(attr, values[i], mi)
            rdatom.SetMonomerInfo(mi)
            # other properties
            for attr in other_attrs.keys():
                value = other_attrs[attr][i]
                attr = _TOPOLOGY_ATTRS[attr].singular
                _set_atom_property(rdatom, attr, value)
            _set_atom_property(rdatom, "index", int(atom.ix))
            # add atom
            index = mol.AddAtom(rdatom)
            atom_mapper[atom.ix] = index

        try:
            if (len(ag.bonds) == 0) and (ag.n_atoms > 1):
                # force guessing bonds
                raise NoDataError
        except NoDataError:
            warnings.warn(
                "No `bonds` attribute in this AtomGroup. Guessing bonds based "
                "on atoms coordinates")
            ag.guess_bonds()

        terminal_atom_indices = []
        for bond in ag.bonds:
            try:
                bond_indices = [atom_mapper[i] for i in bond.indices]
            except KeyError:
                # one of the atoms of the bond is not part of the atomgroup.
                # can happen for terminal atoms.
                # save the bond atom that is in the atomgroup for later
                terminal_atom_indices.extend([atom_mapper[i]
                                              for i in bond.indices
                                              if i in atom_mapper.keys()])
                # skip adding this bond
                continue
            bond_type = RDBONDORDER.get(bond.order, Chem.BondType.SINGLE)
            mol.AddBond(*bond_indices, bond_type)

        mol.UpdatePropertyCache(strict=False)

        # infer bond orders and formal charges from the connectivity
        _infer_bo_and_charges(mol, terminal_atom_indices)
        mol = _standardize_patterns(mol)

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
    if isinstance(value, np.generic):
        # convert numpy types to python standard types
        value = value.item()
    if attr == "names":
        # RDKit needs the name to be properly formated for a
        # PDB file (1 letter elements start at col 14)
        name = re.findall(r'(\D+|\d+)', value)
        if len(name) == 2:
            symbol, number = name
        else:
            symbol, number = name[0], ""
        value = "{:>2.2}{:<2.2}".format(symbol, number)
    # set attribute value in RDKit MonomerInfo
    rdattr = RDATTRIBUTES[attr]
    getattr(mi, "Set%s" % rdattr)(value)


def _set_atom_property(atom, attr, value):
    """Converts an MDAnalysis atom attribute into an RDKit atom property"""
    if isinstance(value, (float, np.float)):
        atom.SetDoubleProp("_MDAnalysis_%s" % attr, float(value))
    elif isinstance(value, (int, np.int)):
        atom.SetIntProp("_MDAnalysis_%s" % attr, int(value))
    else:
        atom.SetProp("_MDAnalysis_%s" % attr, value)


def _infer_bo_and_charges(mol, terminal_atom_indices=[]):
    """Infer bond orders and formal charges from a molecule.

    Since most MD topology files don't explicitely retain informations on bond
    orders or charges, it has to be guessed from the topology. This is done by
    looping other each atom and comparing its expected valence to the current
    valence to get the Number of Unpaired Electrons (NUE).
    If an atom has a negative NUE, it needs a positive formal charge (-NUE).
    If two neighbouring atoms have UEs, the bond between them most
    likely has to be increased by the value of the smallest NUE.
    If after this process, an atom still has UEs, it's either a radical
    (because one its bonds was cut when creating the AtomGroup) or it needs a
    negative formal charge of -NUE. Since these radical atoms can be detected
    when looping over the bonds of the AtomGroup, only atoms that are not part
    of this "terminal_atoms" list will be assigned a negative formal charge.

    Parameters
    ----------

    mol : rdkit.Chem.rdchem.RWMol
        The molecule is modified inplace and must have all hydrogens added

    terminal_atom_indices : list
        List of terminal atoms indices, i.e. atoms at the edges of a molecule
    """

    for atom in mol.GetAtoms():
        # get NUE for each possible valence
        expected_vs = PERIODIC_TABLE.GetValenceList(atom.GetAtomicNum())
        current_v = atom.GetTotalValence()
        nue = [v - current_v for v in expected_vs]

        # if there's only one possible valence state and the corresponding
        # NUE is negative, it means we can only add a positive charge to
        # the atom
        if (len(nue) == 1) and (nue[0] < 0):
            atom.SetFormalCharge(-nue[0])
            mol.UpdatePropertyCache(strict=False)
            continue
        else:
            neighbors = atom.GetNeighbors()
            # check if one of the neighbors has a common NUE
            for i, na in enumerate(neighbors, start=1):
                # get NUE for the neighbor
                na_expected_vs = PERIODIC_TABLE.GetValenceList(
                    na.GetAtomicNum())
                na_current_v = na.GetTotalValence()
                na_nue = [v - na_current_v for v in na_expected_vs]
                # smallest common NUE
                common_nue = min(
                    min([i for i in nue if i >= 0], default=0),
                    min([i for i in na_nue if i >= 0], default=0)
                )
                # a common NUE of 0 means we don't need to do anything
                if common_nue != 0:
                    # increase bond order
                    bond = mol.GetBondBetweenAtoms(
                        atom.GetIdx(), na.GetIdx())
                    order = common_nue + 1
                    bond.SetBondType(RDBONDORDER[order])
                    mol.UpdatePropertyCache(strict=False)
                    if i < len(neighbors):
                        # recalculate nue for atom
                        current_v = atom.GetTotalValence()
                        nue = [v - current_v for v in expected_vs]

            # if the atom still has unpaired electrons
            current_v = atom.GetTotalValence()
            nue = [v - current_v for v in expected_vs][0]
            if nue > 0:
                # keep the radical if it's a terminal atom
                # else transform it to a negative charge
                if atom.GetIdx() not in terminal_atom_indices:
                    atom.SetFormalCharge(-nue)
                    atom.SetNumRadicalElectrons(0)
                    mol.UpdatePropertyCache(strict=False)


def _standardize_patterns(mol):
    """Standardize functional groups using reactions from SMARTS patterns

    Due to the way reactions work, we first have to split the molecule by
    fragments. Then, for each fragment, we apply the standardization reactions.
    If a pattern is matched N times in the molecule, the reaction will return N
    products as an array of shape (N, 1). Only the first product will be kept
    and the same reaction will be reapplied to the product N times in total.
    Finally, the fragments are recombined.
    """

    fragments = []
    for reactant in Chem.GetMolFrags(mol, asMols=True):

        for name, reaction in [
            ("Cterm", "[C-;v3:1]=[O:2]>>[C;+0:1]=[O:2]"),
            ("Nterm", "[N-;v2;H1:1]>>[N;+0:1]"),
            ("keto-enolate", "[C-:1]-[C:2]=[O:3]>>[C;+0:1]=[C:2]-[O;-1:3]"),
            ("ARG", "[N;H1:1]-[C-;v3:2](-[N;H2:3])-[N;H2:4]"
                    ">>[N:1]-[C;+0:2](-[N:3])=[N;+1:4]"),
            ("sulfone", "[S;v4:1](-[O-;v1:2])-[O-;v1:3]"
                        ">>[S;v6:1](=[O;+0:2])=[O;+0:3]"),
            ("nitro", "[N;v3:1](-[O-;v1:2])-[O-;v1:3]"
                      ">>[N;+1:1](-[O;-1:2])=[O;+0:3]"),
            ("anion-*=*-anion", "[*-:1]-[*:2]=[*:3]-[*-:4]"
                                ">>[*;+0:1]=[*:2]-[*:3]=[*;+0:4]"),
            ("anion-*=*-*=*-anion", "[*-:1][*:2]=[*:3][*:4]=[*:5][*-:6]"
             ">>[*;+0:1]=[*:2]-[*:3]=[*:4]-[*:5]=[*;+0:6]"),
            ("anion-*=*-*=*-*=*-anion",
             "[*-:1][*:2]=[*:3][*:4]=[*:5][*:6]=[*:7][*-:8]"
             ">>[*;+0:1]=[*:2]-[*:3]=[*:4]-[*:5]=[*:6]-[*:7]=[*;+0:8]"),
        ]:
            # count how many times the reaction should be run
            pattern = Chem.MolFromSmarts(reaction.split(">>")[0])
            n_matches = len(reactant.GetSubstructMatches(pattern))

            # run the reaction for each matched pattern
            rxn = AllChem.ReactionFromSmarts(reaction)
            for n in range(n_matches):
                products = rxn.RunReactants((reactant,))
                # only keep the first product
                if products:
                    product = products[0][0]
                    product.UpdatePropertyCache(strict=False)
                    # map back atom properties from the reactant to the product
                    _reassign_props_after_reaction(reactant, product)
                    # apply the next reaction to the product
                    reactant = product
                else:
                    # exit the n_matches loop if there's no product. Example
                    # where this is needed: SO^{4}_{2-} will match the sulfone
                    # pattern 6 times but the reaction is only needed once
                    break

        fragments.append(reactant)

    # recombine fragments
    mol = fragments.pop(0)
    for fragment in fragments:
        mol = Chem.CombineMols(mol, fragment)

    return mol


def _reassign_props_after_reaction(reactant, product):
    """Maps back atomic properties from the reactant to the product.
    The product molecule is modified inplace.
    """
    for atom in product.GetAtoms():
        try:
            atom.GetIntProp("old_mapno")
        except KeyError:
            pass
        else:
            atom.ClearProp("old_mapno")
            idx = atom.GetUnsignedProp("react_atom_idx")
            old_atom = reactant.GetAtomWithIdx(idx)
            for prop, value in old_atom.GetPropsAsDict().items():
                if prop.startswith("_MDAnalysis"):
                    attr = prop.split("_")[-1]
                    _set_atom_property(atom, attr, value)
        atom.ClearProp("react_atom_idx")
