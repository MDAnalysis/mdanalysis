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

Read coordinates data from an `RDKit <https://www.rdkit.org/docs/>`__ :class:`rdkit.Chem.rdchem.Mol` with
:class:`RDKitReader` into an MDAnalysis Universe. Convert it back to a
:class:`rdkit.Chem.rdchem.Mol` with :class:`RDKitConverter`.


Example
-------

To read an RDKit molecule and then convert the AtomGroup back to an RDKit
molecule::

    >>> from rdkit import Chem
    >>> import MDAnalysis as mda
    >>> mol = Chem.MolFromMol2File("docking_poses.mol2", removeHs=False)
    >>> u = mda.Universe(mol)
    >>> u
    <Universe with 42 atoms>
    >>> u.trajectory
    <RDKitReader with 10 frames of 42 atoms>
    >>> u.atoms.convert_to("RDKIT")
    <rdkit.Chem.rdchem.Mol object at 0x7fcebb958148>


Classes
-------

.. autoclass:: RDKitReader
   :members:

.. autoclass:: RDKitConverter
   :members:

.. autofunction:: _infer_bo_and_charges

.. autofunction:: _standardize_patterns

.. autofunction:: _rebuild_conjugated_bonds

"""

import warnings
import re
import copy

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
    :class:`~MDAnalysis.core.universe.Universe` to RDKit
    :class:`~rdkit.Chem.rdchem.Mol`

    MDanalysis attributes are stored in each RDKit
    :class:`~rdkit.Chem.rdchem.Atom` of the resulting molecule in two different
    ways:

    * in an :class:`~rdkit.Chem.rdchem.AtomPDBResidueInfo` object available
      through the :meth:`~rdkit.Chem.rdchem.Atom.GetMonomerInfo` method if it's
      an attribute that is typically found in a PDB file,
    * directly as an atom property available through the
      :meth:`~rdkit.Chem.rdchem.Atom.GetProp` methods for the others.

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

    To access MDAnalysis properties::

        >>> import MDAnalysis as mda
        >>> from MDAnalysis.tests.datafiles import PDB_full
        >>> u = mda.Universe(PDB_full)
        >>> mol = u.select_atoms('resname DMS').convert_to('RDKIT')
        >>> mol.GetAtomWithIdx(0).GetMonomerInfo().GetResidueName()
        'DMS'

    To create a molecule for each frame of a trajectory::

        from MDAnalysisTests.datafiles import PSF, DCD
        from rdkit.Chem.Descriptors3D import Asphericity

        u = mda.Universe(PSF, DCD)
        elements = mda.topology.guessers.guess_types(u.atoms.names)
        u.add_TopologyAttr('elements', elements)
        ag = u.select_atoms("resid 1-10")

        for ts in u.trajectory:
            mol = ag.convert_to("RDKIT")
            x = Asphericity(mol)


    Notes
    -----

    The converter requires the :class:`~MDAnalysis.core.topologyattrs.Elements`
    attribute to be present in the topology, else it will fail.

    It also requires the `bonds` attribute, although they will be automatically
    guessed if not present.

    If both ``tempfactors`` and ``bfactors`` attributes are present, the
    conversion will fail, since only one of these should be present.
    Refer to Issue #1901 for a solution

    Hydrogens should be explicit in the topology file. If this is not the case,
    use the parameter ``NoImplicit=False`` when using the converter to allow
    implicit hydrogens and disable inferring bond orders and charges.

    Since one of the main use case of the converter is converting trajectories
    and not just a topology, creating a new molecule from scratch for every
    frame would be too slow so the converter uses a caching system. The cache
    only remembers the id of the last AtomGroup that was converted, as well
    as the arguments that were passed to the converter. This means that using
    ``u.select_atoms("protein").convert_to("RDKIT")`` will not benefit from the
    cache since the selection is deleted from memory as soon as the conversion
    is finished. Instead, users should do this in two steps by first saving the
    selection in a variable and then converting the saved AtomGroup. It also
    means that ``ag.convert_to("RDKIT")`` followed by
    ``ag.convert_to("RDKIT", NoImplicit=False)`` will not use the cache.
    Finally if you're modifying the AtomGroup in place between two conversions,
    the id of the AtomGroup won't change and thus the converter will use the
    cached molecule. For this reason, you can pass a ``cache=False`` argument
    to the converter to bypass the caching system.
    Note that the cached molecule doesn't contain the coordinates of the atoms.


    .. versionadded:: 2.0.0

    """

    lib = 'RDKIT'
    units = {'time': None, 'length': 'Angstrom'}
    _cache = dict()

    def convert(self, obj, cache=True, NoImplicit=True, max_iter=200):
        """Write selection at current trajectory frame to
        :class:`~rdkit.Chem.rdchem.Mol`.

        Parameters
        -----------
        obj : :class:`~MDAnalysis.core.groups.AtomGroup` or :class:`~MDAnalysis.core.universe.Universe`

        cache : bool
            Use a cached copy of the molecule's topology when available. To be
            used, the cached molecule and the new one have to be made from the
            same AtomGroup object (same id) and with the same arguments passed
            to the converter (with the exception of this `cache` argument)
        NoImplicit : bool
            Prevent adding hydrogens to the molecule
        max_iter : int
            Maximum number of iterations to standardize conjugated systems.
            See :func:`_rebuild_conjugated_bonds`
        """
        # parameters passed to atomgroup_to_mol and used by the cache
        kwargs = dict(NoImplicit=NoImplicit, max_iter=max_iter)

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

        if cache:
            # key used to search the cache
            key = f"<{id(ag):#x}>" + ",".join(f"{key}={value}"
                                            for key, value in kwargs.items())
            try:
                mol = self._cache[key]
            except KeyError:
                # only keep the current molecule in cache
                self._cache.clear()
                # create the topology
                self._cache[key] = mol = self.atomgroup_to_mol(ag, **kwargs)
            # continue on copy of the cached molecule
            mol = copy.deepcopy(mol)
        else:
            self._cache.clear()
            mol = self.atomgroup_to_mol(ag, **kwargs)

        # add a conformer for the current Timestep
        if hasattr(ag, "positions"):
            if np.isnan(ag.positions).any():
                warnings.warn("NaN detected in coordinates, the output "
                              "molecule will not have 3D coordinates assigned")
            else:
                # assign coordinates
                conf = Chem.Conformer(mol.GetNumAtoms())
                for atom in mol.GetAtoms():
                    idx = atom.GetIntProp("_MDAnalysis_index")
                    xyz = ag.positions[idx].astype(float)
                    conf.SetAtomPosition(atom.GetIdx(), xyz)
                mol.AddConformer(conf)
                # assign R/S to atoms and Z/E to bonds
                Chem.AssignStereochemistryFrom3D(mol)
                Chem.SetDoubleBondNeighborDirections(mol)

        return mol

    def atomgroup_to_mol(self, ag, NoImplicit=True, max_iter=200):
        """Converts an AtomGroup to an RDKit molecule.

        Parameters
        -----------
        ag : MDAnalysis.core.groups.AtomGroup
            The AtomGroup to convert
        NoImplicit : bool
            Prevent adding hydrogens to the molecule
        max_iter : int
            Maximum number of iterations to standardize conjugated systems.
            See :func:`_rebuild_conjugated_bonds`
        """
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
                "orders to structures with implicit hydrogens. Alternatively, "
                "you can use the parameter `NoImplicit=False` when using the "
                "converter to allow implicit hydrogens and disable inferring "
                "bond orders and charges."
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
            # enable/disable adding implicit H to the molecule
            rdatom.SetNoImplicit(NoImplicit)
            # add PDB-like properties
            mi = Chem.AtomPDBResidueInfo()
            for attr, values in pdb_attrs.items():
                _add_mda_attr_to_rdkit(attr, values[i], mi)
            rdatom.SetMonomerInfo(mi)
            # other properties
            for attr in other_attrs.keys():
                value = other_attrs[attr][i]
                attr = "_MDAnalysis_%s" % _TOPOLOGY_ATTRS[attr].singular
                _set_atom_property(rdatom, attr, value)
            _set_atom_property(rdatom, "_MDAnalysis_index", i)
            # add atom
            index = mol.AddAtom(rdatom)
            atom_mapper[atom.ix] = index

        try:
            ag.bonds
        except NoDataError:
            warnings.warn(
                "No `bonds` attribute in this AtomGroup. Guessing bonds based "
                "on atoms coordinates")
            ag.guess_bonds()

        for bond in ag.bonds:
            try:
                bond_indices = [atom_mapper[i] for i in bond.indices]
            except KeyError:
                continue
            bond_type = RDBONDORDER.get(bond.order, Chem.BondType.SINGLE)
            mol.AddBond(*bond_indices, bond_type)

        mol.UpdatePropertyCache(strict=False)

        if NoImplicit:
            # infer bond orders and formal charges from the connectivity
            _infer_bo_and_charges(mol)
            mol = _standardize_patterns(mol, max_iter)

        # sanitize
        Chem.SanitizeMol(mol)

        return mol


def _add_mda_attr_to_rdkit(attr, value, mi):
    """Converts an MDAnalysis atom attribute into the RDKit equivalent and
    stores it into an RDKit :class:`~rdkit.Chem.rdchem.AtomPDBResidueInfo`.

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
        # RDKit needs the name to be properly formatted for a
        # PDB file (1 letter elements start at col 14)
        name = re.findall(r'(\D+|\d+)', value)
        if len(name) == 2:
            symbol, number = name
            if len(number) > 2 and len(symbol) == 1:
                value = "{}{}".format(symbol, number)
            else:
                value = "{:>2.2}{:<2.2}".format(symbol, number)
        else:
            # no number in the name
            value = " {:<}".format(name[0])

    # set attribute value in RDKit MonomerInfo
    rdattr = RDATTRIBUTES[attr]
    getattr(mi, "Set%s" % rdattr)(value)


def _set_atom_property(atom, attr, value):
    """Saves any attribute and value into an RDKit atom property"""
    if isinstance(value, (float, np.floating)):
        atom.SetDoubleProp(attr, float(value))
    elif isinstance(value, (int, np.integer)):
        atom.SetIntProp(attr, int(value))
    else:
        atom.SetProp(attr, value)


def _infer_bo_and_charges(mol):
    """Infer bond orders and formal charges from a molecule.

    Since most MD topology files don't explicitly retain information on bond
    orders or charges, it has to be guessed from the topology. This is done by
    looping over each atom and comparing its expected valence to the current
    valence to get the Number of Unpaired Electrons (NUE).
    If an atom has a negative NUE, it needs a positive formal charge (-NUE).
    If two neighbouring atoms have UEs, the bond between them most
    likely has to be increased by the value of the smallest NUE.
    If after this process, an atom still has UEs, it needs a negative formal
    charge of -NUE.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.RWMol
        The molecule is modified inplace and must have all hydrogens added

    Notes
    -----
    This algorithm is order dependant. For example, for a carboxylate group
    R-C(-O)-O the first oxygen read will receive a double bond and the other
    one will be charged. It will also affect more complex conjugated systems.
    """

    for atom in sorted(mol.GetAtoms(), reverse=True,
                       key=lambda a: _get_nb_unpaired_electrons(a)[0]):
        # get NUE for each possible valence
        nue = _get_nb_unpaired_electrons(atom)
        # if there's only one possible valence state and the corresponding
        # NUE is negative, it means we can only add a positive charge to
        # the atom
        if (len(nue) == 1) and (nue[0] < 0):
            atom.SetFormalCharge(-nue[0])
            mol.UpdatePropertyCache(strict=False)
        # go to next atom if above case or atom has no unpaired electron
        if (len(nue) == 1) and (nue[0] <= 0):
            continue
        else:
            neighbors = sorted(atom.GetNeighbors(), reverse=True,
                               key=lambda a: _get_nb_unpaired_electrons(a)[0])
            # check if one of the neighbors has a common NUE
            for i, na in enumerate(neighbors, start=1):
                # get NUE for the neighbor
                na_nue = _get_nb_unpaired_electrons(na)
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
                        nue = _get_nb_unpaired_electrons(atom)

            # if the atom still has unpaired electrons
            nue = _get_nb_unpaired_electrons(atom)[0]
            if nue > 0:
                # transform it to a negative charge
                atom.SetFormalCharge(-nue)
                atom.SetNumRadicalElectrons(0)
                mol.UpdatePropertyCache(strict=False)


def _get_nb_unpaired_electrons(atom):
    """Calculate the number of unpaired electrons (NUE) of an atom

    Parameters
    ----------
    atom: rdkit.Chem.rdchem.Atom
        The atom for which the NUE will be computed

    Returns
    -------
    nue : list
        The NUE for each possible valence of the atom
    """
    expected_vs = PERIODIC_TABLE.GetValenceList(atom.GetAtomicNum())
    current_v = atom.GetTotalValence() - atom.GetFormalCharge()
    return [v - current_v for v in expected_vs]


def _standardize_patterns(mol, max_iter=200):
    """Standardizes functional groups

    Uses :func:`_rebuild_conjugated_bonds` to standardize conjugated systems,
    and SMARTS reactions for other functional groups.
    Due to the way reactions work, we first have to split the molecule by
    fragments. Then, for each fragment, we apply the standardization reactions.
    Finally, the fragments are recombined.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.RWMol
        The molecule to standardize
    max_iter : int
        Maximum number of iterations to standardize conjugated systems

    Returns
    -------
    mol : rdkit.Chem.rdchem.Mol
        The standardized molecule

    Notes
    -----
    The following functional groups are transformed in this order:

    +---------------+------------------------------------------------------------------------------+
    | Name          | Reaction                                                                     |
    +===============+==============================================================================+
    | conjugated    | ``[*-;!O:1]-[*:2]=[*:3]-[*-:4]>>[*+0:1]=[*:2]-[*:3]=[*+0:4]``                |
    +---------------+------------------------------------------------------------------------------+
    | conjugated-N+ | ``[N;X3;v3:1]-[*:2]=[*:3]-[*-:4]>>[N+:1]=[*:2]-[*:3]=[*+0:4]``               |
    +---------------+------------------------------------------------------------------------------+
    | conjugated-O- | ``[O:1]=[#6:2]-[*:3]=[*:4]-[*-:5]>>[O-:1]-[*:2]=[*:3]-[*:4]=[*+0:5]``        |
    +---------------+------------------------------------------------------------------------------+
    | Cterm         | ``[C-;X2:1]=[O:2]>>[C+0:1]=[O:2]``                                           |
    +---------------+------------------------------------------------------------------------------+
    | Nterm         | ``[N-;X2;H1:1]>>[N+0:1]``                                                    |
    +---------------+------------------------------------------------------------------------------+
    | keto-enolate  | ``[#6-:1]-[#6:2]=[O:3]>>[#6+0:1]=[#6:2]-[O-:3]``                             |
    +---------------+------------------------------------------------------------------------------+
    | arginine      | ``[N;H1:1]-[C-;X3;H0:2](-[N;H2:3])-[N;H2:4]>>[N:1]-[C+0:2](-[N:3])=[N+:4]``  |
    +---------------+------------------------------------------------------------------------------+
    | histidine     | ``[#6+0;H0:1]=[#6+0:2]-[#7;X3:3]-[#6-;X3:4]>>[#6:1]=[#6:2]-[#7+:3]=[#6+0:4]``|
    +---------------+------------------------------------------------------------------------------+
    | sulfone       | ``[S;X4;v4:1](-[O-;X1:2])-[O-;X1:3]>>[S:1](=[O+0:2])=[O+0:3]``               |
    +---------------+------------------------------------------------------------------------------+
    | nitro         | ``[N;X3;v3:1](-[O-;X1:2])-[O-;X1:3]>>[N+:1](-[O-:2])=[O+0:3]``               |
    +---------------+------------------------------------------------------------------------------+
 
    """

    # standardize conjugated systems
    _rebuild_conjugated_bonds(mol, max_iter)

    fragments = []
    for reactant in Chem.GetMolFrags(mol, asMols=True):

        for name, reaction in [
            ("Cterm", "[C-;X2:1]=[O:2]>>[C+0:1]=[O:2]"),
            ("Nterm", "[N-;X2;H1:1]>>[N+0:1]"),
            ("keto-enolate", "[#6-:1]-[#6:2]=[O:3]>>[#6+0:1]=[#6:2]-[O-:3]"),
            ("ARG", "[N;H1:1]-[C-;X3;H0:2](-[N;H2:3])-[N;H2:4]"
                    ">>[N:1]-[C+0:2](-[N:3])=[N+:4]"),
            ("HIP", "[#6+0;H0:1]=[#6+0:2]-[#7;X3:3]-[#6-;X3:4]"
                    ">>[#6:1]=[#6:2]-[#7+:3]=[#6+0:4]"),
            ("sulfone", "[S;X4;v4:1](-[O-;X1:2])-[O-;X1:3]"
                        ">>[S:1](=[O+0:2])=[O+0:3]"),
            ("nitro", "[N;X3;v3:1](-[O-;X1:2])-[O-;X1:3]"
                      ">>[N+:1](-[O-:2])=[O+0:3]"),
        ]:
            reactant.UpdatePropertyCache(strict=False)
            Chem.Kekulize(reactant)
            reactant = _run_reaction(reaction, reactant)

        fragments.append(reactant)

    # recombine fragments
    mol = fragments.pop(0)
    for fragment in fragments:
        mol = Chem.CombineMols(mol, fragment)

    return mol


def _run_reaction(reaction, reactant):
    """Runs a reaction until all reactants are transformed

    If a pattern is matched N times in the molecule, the reaction will return N
    products as an array of shape (N, 1). Only the first product will be kept
    and the same reaction will be reapplied to the product N times in total.

    Parameters
    ----------
    reaction : str
        SMARTS reaction
    reactant : rdkit.Chem.rdchem.Mol
        The molecule to transform

    Returns
    -------
    product : rdkit.Chem.rdchem.Mol
        The final product of the reaction
    """
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
            # map back atom properties from the reactant to the product
            _reassign_props_after_reaction(reactant, product)
            # apply the next reaction to the product
            product.UpdatePropertyCache(strict=False)
            reactant = product
        else:
            # exit the n_matches loop if there's no product. Example
            # where this is needed: SO^{4}_{2-} will match the sulfone
            # pattern 6 times but the reaction is only needed once
            break
    return reactant


def _rebuild_conjugated_bonds(mol, max_iter=200):
    """Rebuild conjugated bonds without negatively charged atoms at the
    beginning and end of the conjugated system

    Depending on the order in which atoms are read during the conversion, the
    :func:`_infer_bo_and_charges` function might write conjugated systems with
    a double bond less and both edges of the system as anions instead of the
    usual alternating single and double bonds. This function corrects this
    behaviour by using an iterative procedure.
    The problematic molecules always follow the same pattern:
    ``anion[-*=*]n-anion`` instead of ``*=[*-*=]n*``, where ``n`` is the number
    of successive single and double bonds. The goal of the iterative procedure
    is to make ``n`` as small as possible by consecutively transforming
    ``anion-*=*`` into ``*=*-anion`` until it reaches the smallest pattern with
    ``n=1``. This last pattern is then transformed from ``anion-*=*-anion`` to
    ``*=*-*=*``.
    Since ``anion-*=*`` is the same as ``*=*-anion`` in terms of SMARTS, we can
    control that we don't transform the same triplet of atoms back and forth by
    adding their indices to a list.
    The molecule needs to be kekulized first to also cover systems
    with aromatic rings.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.RWMol
        The molecule to transform, modified inplace
    max_iter : int
        Maximum number of iterations performed by the function
    """
    mol.UpdatePropertyCache(strict=False)
    Chem.Kekulize(mol)
    pattern = Chem.MolFromSmarts("[*-;!O]-[*+0]=[*+0]")
    # number of unique matches with the pattern
    n_matches = len(set([match[0]
                         for match in mol.GetSubstructMatches(pattern)]))
    if n_matches == 0:
        # nothing to standardize
        return
    # check if there's an even number of anion-*=* patterns
    elif n_matches % 2 == 0:
        end_pattern = Chem.MolFromSmarts("[*-;!O]-[*+0]=[*+0]-[*-]")
    else:
        # as a last resort, the only way to standardize is to find a nitrogen
        # that can accept a double bond and a positive charge
        # or a carbonyl that will become an enolate
        end_pattern = Chem.MolFromSmarts(
            "[*-;!O]-[*+0]=[*+0]-[$([#7;X3;v3]),$([#6+0]=O)]")
    backtrack = []
    for _ in range(max_iter):
        # simplest case where n=1
        end_match = mol.GetSubstructMatch(end_pattern)
        if end_match:
            # index of each atom
            anion1, a1, a2, anion2 = end_match
            term_atom = mol.GetAtomWithIdx(anion2)
            # [*-]-*=*-C=O
            if term_atom.GetAtomicNum() == 6 and term_atom.GetFormalCharge() == 0:
                for neighbor in term_atom.GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(anion2, neighbor.GetIdx())
                    if neighbor.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2:
                        bond.SetBondType(Chem.BondType.SINGLE)
                        neighbor.SetFormalCharge(-1)
            else:
                # [*-]-*=*-N
                if term_atom.GetAtomicNum() == 7 and term_atom.GetFormalCharge() == 0:
                    end_charge = 1
                # [*-]-*=*-[*-]
                else:
                    end_charge = 0
                mol.GetAtomWithIdx(anion2).SetFormalCharge(end_charge)
            # common part of the conjugated systems: [*-]-*=*
            mol.GetAtomWithIdx(anion1).SetFormalCharge(0)
            mol.GetBondBetweenAtoms(anion1, a1).SetBondType(
                Chem.BondType.DOUBLE)
            mol.GetBondBetweenAtoms(a1, a2).SetBondType(Chem.BondType.SINGLE)
            mol.GetBondBetweenAtoms(a2, anion2).SetBondType(
                Chem.BondType.DOUBLE)
            mol.UpdatePropertyCache(strict=False)

        # shorten the anion-anion pattern from n to n-1
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            # check if we haven't already transformed this triplet
            for match in matches:
                # sort the indices for the comparison
                g = tuple(sorted(match))
                if g in backtrack:
                    # already transformed
                    continue
                else:
                    # take the first one that hasn't been tried
                    anion, a1, a2 = match
                    backtrack.append(g)
                    break
            else:
                anion, a1, a2 = matches[0]
            # charges
            mol.GetAtomWithIdx(anion).SetFormalCharge(0)
            mol.GetAtomWithIdx(a2).SetFormalCharge(-1)
            # bonds
            mol.GetBondBetweenAtoms(anion, a1).SetBondType(
                Chem.BondType.DOUBLE)
            mol.GetBondBetweenAtoms(a1, a2).SetBondType(Chem.BondType.SINGLE)
            mol.UpdatePropertyCache(strict=False)
            # start new iteration
            continue

        # no more changes to apply
        return

    # reached max_iter
    warnings.warn("The standardization could not be completed within a "
                  "reasonable number of iterations")


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
                _set_atom_property(atom, prop, value)
            # fix bonds with "crossed" stereo
            for bond in atom.GetBonds():
                if bond.GetStereo() == Chem.BondStereo.STEREOANY:
                    bond.SetStereo(Chem.BondStereo.STEREONONE)
        atom.ClearProp("react_atom_idx")
