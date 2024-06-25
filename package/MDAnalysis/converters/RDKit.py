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

"""RDKit molecule I/O --- :mod:`MDAnalysis.converters.RDKit`
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


.. warning::
    The RDKit converter is currently *experimental* and may not work as
    expected for all molecules. Currently the converter accurately
    infers the structures of approximately 99% of the `ChEMBL27`_ dataset.
    Work is currently ongoing on further improving this and updates to the
    converter are expected in future releases of MDAnalysis.
    Please see `Issue #3339`_ and the `RDKitConverter benchmark`_ for more
    details.



Classes
-------

.. autoclass:: RDKitReader
   :members:

.. autoclass:: RDKitConverter
   :members:

.. autofunction:: _infer_bo_and_charges

.. autofunction:: _standardize_patterns

.. autofunction:: _rebuild_conjugated_bonds


.. Links

.. _`ChEMBL27`: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_27/
.. _`Issue #3339`: https://github.com/MDAnalysis/mdanalysis/issues/3339
.. _`RDKitConverter benchmark`: https://github.com/MDAnalysis/RDKitConverter-benchmark
"""

import copy
import warnings
from functools import lru_cache
from io import StringIO

import numpy as np
from numpy.lib import NumpyVersion

from . import base
from ..coordinates import memory
from ..coordinates.PDB import PDBWriter
from ..core.topologyattrs import _TOPOLOGY_ATTRS
from ..exceptions import NoDataError

try:
    # TODO: remove this guard when RDKit has a release
    # that supports NumPy 2
    if NumpyVersion(np.__version__) < "2.0.0":
        from rdkit import Chem
        from rdkit.Chem import AllChem
    else:
        raise ImportError
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
    }
    PERIODIC_TABLE = Chem.GetPeriodicTable()


# charges that should be assigned to monatomic cations
# structure --> atomic number : formal charge
# anion charges are directly handled by the code using the typical valence
# of the atom
MONATOMIC_CATION_CHARGES = {
    3: 1, 11: 1, 19: 1, 37: 1, 47: 1, 55: 1,
    12: 2, 20: 2, 29: 2, 30: 2, 38: 2, 56: 2,
    26: 2,  # Fe could also be +3
    13: 3,
}
# reactions uses by _standardize_patterns to fix challenging cases
# must have single reactant and product, and cannot add any atom
STANDARDIZATION_REACTIONS = [
    "[C-;X2;H0:1]=[O:2]>>[C+0:1]=[O:2]",  # Cterm
    "[N-;X2;H1;$(N-[*^3]):1]>>[N+0:1]",  # Nterm
    "[#6-:1]-[#6:2]=[O:3]>>[#6+0:1]=[#6:2]-[O-:3]",  # keto-enolate
    "[C-;v3:1]-[#7+0;v3;H2:2]>>[#6+0:1]=[#7+:2]",  # ARG
    "[#6+0;H0:1]=[#6+0:2]-[#7;X3:3]-[#6-;X3:4]"
    ">>[#6:1]=[#6:2]-[#7+:3]=[#6+0:4]",  # HIP
    "[S;D4;!v6:1]-[*-:2]>>[S;v6:1]=[*+0:2]",  # sulfone
    "[#7+0;X3:1]-[*-:2]>>[#7+:1]=[*+0:2]",  # charged-N
]
_deduce_PDB_atom_name = PDBWriter(StringIO())._deduce_PDB_atom_name


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
    |                       | atom.GetProp("_MDAnalysis_name")          |
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

    Hydrogens should be explicit in the topology file. If this is not the case,
    use the parameter ``NoImplicit=False`` when using the converter to allow
    implicit hydrogens and disable inferring bond orders and charges.

    Since one of the main use case of the converter is converting trajectories
    and not just a topology, creating a new molecule from scratch for every
    frame would be too slow so the converter uses a caching system. The cache
    only stores the 2 most recent AtomGroups that were converted, and is
    sensitive to the arguments that were passed to the converter. The number of
    objects cached can be changed with the function
    :func:`set_converter_cache_size`. However, ``ag.convert_to("RDKIT")``
    followed by ``ag.convert_to("RDKIT", NoImplicit=False)`` will not use the
    cache since the arguments given are different. You can pass a
    ``cache=False`` argument to the converter to bypass the caching system.

    The ``_MDAnalysis_index`` property of the resulting molecule corresponds
    to the index of the specific :class:`~MDAnalysis.core.groups.AtomGroup`
    that was converted, which may not always match the ``index`` property.

    To get a better understanding of how the converter works under the hood,
    please refer to the following RDKit UGM presentation:

    * `Video (4:55 to 8:05) <https://youtu.be/5b5wYmK4URU>`__
    * `Slides <https://github.com/rdkit/UGM_2020/blob/master/Presentations/C%C3%A9dricBouysset_From_RDKit_to_the_Universe.pdf>`__

    There are some molecules containing specific patterns that the converter
    cannot currently tackle correctly. See
    `Issue #3339 <https://github.com/MDAnalysis/mdanalysis/issues/3339>`__ for
    more info.

    .. versionadded:: 2.0.0

    .. versionchanged:: 2.2.0
        Improved the accuracy of the converter. Atoms in the resulting molecule
        now follow the same order as in the AtomGroup. The output of
        ``atom.GetMonomerInfo().GetName()`` now follows the guidelines for PDB
        files while the original name is still available through
        ``atom.GetProp("_MDAnalysis_name")``

    """

    lib = 'RDKIT'
    units = {'time': None, 'length': 'Angstrom'}

    def convert(self, obj, cache=True, NoImplicit=True, max_iter=200,
                force=False):
        """Write selection at current trajectory frame to
        :class:`~rdkit.Chem.rdchem.Mol`.

        Parameters
        -----------
        obj : :class:`~MDAnalysis.core.groups.AtomGroup` or :class:`~MDAnalysis.core.universe.Universe`
        cache : bool
            Use a cached copy of the molecule's topology when available. To be
            used, the cached molecule and the new one have to be made from the
            same AtomGroup selection and with the same arguments passed
            to the converter
        NoImplicit : bool
            Prevent adding hydrogens to the molecule
        max_iter : int
            Maximum number of iterations to standardize conjugated systems.
            See :func:`_rebuild_conjugated_bonds`
        force : bool
            Force the conversion when no hydrogens were detected but
            ``NoImplicit=True``. Useful for inorganic molecules mostly.
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

        # parameters passed to atomgroup_to_mol
        kwargs = dict(NoImplicit=NoImplicit, max_iter=max_iter, force=force)
        if cache:
            mol = atomgroup_to_mol(ag, **kwargs)
            mol = copy.deepcopy(mol)
        else:
            mol = atomgroup_to_mol.__wrapped__(ag, **kwargs)

        # add a conformer for the current Timestep
        if hasattr(ag, "positions"):
            if np.isnan(ag.positions).any():
                warnings.warn("NaN detected in coordinates, the output "
                              "molecule will not have 3D coordinates assigned")
            else:
                # assign coordinates
                conf = Chem.Conformer(mol.GetNumAtoms())
                for atom in mol.GetAtoms():
                    idx = atom.GetIdx()
                    xyz = ag.positions[idx].astype(float)
                    conf.SetAtomPosition(idx, xyz)
                mol.AddConformer(conf)
                # assign R/S to atoms and Z/E to bonds
                Chem.AssignStereochemistryFrom3D(mol)
                Chem.SetDoubleBondNeighborDirections(mol)

        return mol


@lru_cache(maxsize=2)
def atomgroup_to_mol(ag, NoImplicit=True, max_iter=200, force=False):
    """Converts an AtomGroup to an RDKit molecule without coordinates.

    Parameters
    -----------
    ag : MDAnalysis.core.groups.AtomGroup
        The AtomGroup to convert
    NoImplicit : bool
        Prevent adding hydrogens to the molecule and allow bond orders and
        formal charges to be guessed from the valence of each atom.
    max_iter : int
        Maximum number of iterations to standardize conjugated systems.
        See :func:`_rebuild_conjugated_bonds`
    force : bool
        Force the conversion when no hydrogens were detected but
        ``NoImplicit=True``. Mostly useful for inorganic molecules.
    """
    try:
        elements = ag.elements
    except NoDataError:
        raise AttributeError(
            "The `elements` attribute is required for the RDKitConverter "
            "but is not present in this AtomGroup. Please refer to the "
            "documentation to guess elements from other attributes or "
            "type `help(MDAnalysis.topology.guessers)`") from None

    if "H" not in ag.elements:
        if force:
            warnings.warn(
                "No hydrogen atom found in the topology. "
                "Forcing to continue the conversion."
            )
        elif NoImplicit:
            raise AttributeError(
                "No hydrogen atom could be found in the topology, but the "
                "converter requires all hydrogens to be explicit. You can use "
                "the parameter ``NoImplicit=False`` when using the converter "
                "to allow implicit hydrogens and disable inferring bond "
                "orders and charges. You can also use ``force=True`` to "
                "ignore this error.")

    # attributes accepted in PDBResidueInfo object
    pdb_attrs = {}
    for attr in RDATTRIBUTES.keys():
        if hasattr(ag, attr):
            pdb_attrs[attr] = getattr(ag, attr)
    resnames = pdb_attrs.get("resnames", None)
    if resnames is None:
        def get_resname(idx):
            return ""
    else:
        def get_resname(idx):
            return resnames[idx]

    other_attrs = {}
    for attr in ["charges", "segids", "types", "names"]:
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
            _add_mda_attr_to_rdkit(attr, values[i], mi, get_resname(i))
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
        # reorder atoms to match MDAnalysis, since the reactions from
        # _standardize_patterns will mess up the original order
        order = np.argsort([atom.GetIntProp("_MDAnalysis_index")
                            for atom in mol.GetAtoms()])
        mol = Chem.RenumberAtoms(mol, order.astype(int).tolist())

    # sanitize if possible
    err = Chem.SanitizeMol(mol, catchErrors=True)
    if err:
        warnings.warn("Could not sanitize molecule: "
                      f"failed during step {err!r}")

    return mol


def set_converter_cache_size(maxsize):
    """Set the maximum cache size of the RDKit converter

    Parameters
    ----------
    maxsize : int or None
        If int, the cache will only keep the ``maxsize`` most recent
        conversions in memory. Using ``maxsize=None`` will remove all limits
        to the cache size, i.e. everything is cached.
    """
    global atomgroup_to_mol   # pylint: disable=global-statement
    atomgroup_to_mol = lru_cache(maxsize=maxsize)(atomgroup_to_mol.__wrapped__)


def _add_mda_attr_to_rdkit(attr, value, mi, resname=""):
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
    resname : str
        Residue name of the atom, if available
    """
    if isinstance(value, np.generic):
        # convert numpy types to python standard types
        value = value.item()
    if attr == "names":
        # RDKit needs the name to be properly formatted for a PDB file
        value = _deduce_PDB_atom_name(value, resname)

    # set attribute value in RDKit MonomerInfo
    rdattr = RDATTRIBUTES[attr]
    getattr(mi, "Set%s" % rdattr)(value)


def _set_str_prop(atom, attr, value):
    atom.SetProp(attr, value)


def _set_float_prop(atom, attr, value):
    atom.SetDoubleProp(attr, value)


def _set_np_float_prop(atom, attr, value):
    atom.SetDoubleProp(attr, float(value))


def _set_int_prop(atom, attr, value):
    atom.SetIntProp(attr, value)


def _set_np_int_prop(atom, attr, value):
    atom.SetIntProp(attr, int(value))


def _ignore_prop(atom, attr, value):
    pass


_atom_property_dispatcher = {
    str: _set_str_prop,
    float: _set_float_prop,
    np.float32: _set_np_float_prop,
    np.float64: _set_np_float_prop,
    int: _set_int_prop,
    np.int8: _set_np_int_prop,
    np.int16: _set_np_int_prop,
    np.int32: _set_np_int_prop,
    np.int64: _set_np_int_prop,
    np.uint8: _set_np_int_prop,
    np.uint16: _set_np_int_prop,
    np.uint32: _set_np_int_prop,
    np.uint64: _set_np_int_prop,
}


def _set_atom_property(atom, attr, value):
    """Saves any attribute and value into an RDKit atom property"""
    _atom_property_dispatcher.get(type(value), _ignore_prop)(atom, attr, value)


def _atom_sorter(atom):
    """Sorts atoms in the molecule in a way that makes it easy for the bond
    order and charge infering code to get the correct state on the first
    try. Currently sorts by number of unpaired electrons, then by number of
    heavy atom neighbors (i.e. atoms at the edge first)."""
    num_heavy_neighbors = len([
        neighbor for neighbor in atom.GetNeighbors()
        if neighbor.GetAtomicNum() > 1]
    )
    return (-_get_nb_unpaired_electrons(atom)[0], num_heavy_neighbors)


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
    # heavy atoms sorted by number of heavy atom neighbors (lower first) then
    # NUE (higher first)
    atoms = sorted([atom for atom in mol.GetAtoms()
                    if atom.GetAtomicNum() > 1],
                   key=_atom_sorter)

    for atom in atoms:
        # monatomic ions
        if atom.GetDegree() == 0:
            atom.SetFormalCharge(MONATOMIC_CATION_CHARGES.get(
                                 atom.GetAtomicNum(),
                                 -_get_nb_unpaired_electrons(atom)[0]))
            mol.UpdatePropertyCache(strict=False)
            continue
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
            for na in neighbors:
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
                    # go to next atom if one of the valences is complete
                    nue = _get_nb_unpaired_electrons(atom)
                    if any([n == 0 for n in nue]):
                        break

            # if atom valence is still not filled
            nue = _get_nb_unpaired_electrons(atom)
            if not any([n == 0 for n in nue]):
                # transform nue to charge
                atom.SetFormalCharge(-nue[0])
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
    | conjugated    | ``[*-:1]-[*:2]=[*:3]-[*-:4]>>[*+0:1]=[*:2]-[*:3]=[*+0:4]``                   |
    +---------------+------------------------------------------------------------------------------+
    | conjugated N+ | ``[N;X3;v3:1]-[*:2]=[*:3]-[*-:4]>>[N+:1]=[*:2]-[*:3]=[*+0:4]``               |
    +---------------+------------------------------------------------------------------------------+
    | conjugated O- | ``[O:1]=[#6+0,#7+:2]-[*:3]=[*:4]-[*-:5]>>[O-:1]-[*:2]=[*:3]-[*:4]=[*+0:5]``  |
    +---------------+------------------------------------------------------------------------------+
    | conjug. S=O   | ``[O-:1]-[S;D4;v4:2]-[*:3]=[*:4]-[*-:5]>>[O+0:1]=[*:2]=[*:3]-[*:4]=[*+0:5]`` |
    +---------------+------------------------------------------------------------------------------+
    | Cterm         | ``[C-;X2;H0:1]=[O:2]>>[C+0:1]=[O:2]``                                        |
    +---------------+------------------------------------------------------------------------------+
    | Nterm         | ``[N-;X2;H1;$(N-[*^3]):1]>>[N+0:1]``                                         |
    +---------------+------------------------------------------------------------------------------+
    | keto-enolate  | ``[#6-:1]-[#6:2]=[O:3]>>[#6+0:1]=[#6:2]-[O-:3]``                             |
    +---------------+------------------------------------------------------------------------------+
    | arginine      | ``[C-;v3:1]-[#7+0;v3;H2:2]>>[#6+0:1]=[#7+:2]``                               |
    +---------------+------------------------------------------------------------------------------+
    | histidine     | ``[#6+0;H0:1]=[#6+0:2]-[#7;X3:3]-[#6-;X3:4]>>[#6:1]=[#6:2]-[#7+:3]=[#6+0:4]``|
    +---------------+------------------------------------------------------------------------------+
    | sulfone       | ``[S;D4;!v6:1]-[*-:2]>>[S;v6:1]=[*+0:2]``                                    |
    +---------------+------------------------------------------------------------------------------+
    | charged N     | ``[#7+0;X3:1]-[*-:2]>>[#7+:1]=[*+0:2]``                                      |
    +---------------+------------------------------------------------------------------------------+

    """

    # standardize conjugated systems
    _rebuild_conjugated_bonds(mol, max_iter)

    # list of sanitized reactions
    reactions = []
    for rxn in STANDARDIZATION_REACTIONS:
        reaction = AllChem.ReactionFromSmarts(rxn)
        reactions.append(reaction)

    # fragment mol (reactions must have single reactant and product)
    fragments = list(Chem.GetMolFrags(mol, asMols=True))
    for reactant in fragments:
        _apply_reactions(reactions, reactant)

    # recombine fragments
    newmol = fragments.pop(0)
    for fragment in fragments:
        newmol = Chem.CombineMols(newmol, fragment)

    return newmol


def _apply_reactions(reactions, reactant):
    """Applies a series of unimolecular reactions to a molecule. The reactions
    must not add any atom to the product. The molecule is modified inplace.

    Parameters
    ----------
    reactions : List[rdkit.Chem.rdChemReactions.ChemicalReaction]
        Reactions from SMARTS. Each reaction is applied until no more
        transformations can be made.
    reactant : rdkit.Chem.rdchem.Mol
        The molecule to transform inplace

    """
    reactant.UpdatePropertyCache(strict=False)
    Chem.Kekulize(reactant)
    for reaction in reactions:
        while reaction.RunReactantInPlace(reactant):
            reactant.UpdatePropertyCache(strict=False)
        reactant.UpdatePropertyCache(strict=False)
        Chem.Kekulize(reactant)


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
    This function also handles conjugated systems with triple bonds.
    The molecule needs to be kekulized first to also cover systems
    with aromatic rings.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.RWMol
        The molecule to transform, modified inplace
    max_iter : int
        Maximum number of iterations performed by the function

    Notes
    -----
    The molecule is modified inplace
    """
    mol.UpdatePropertyCache(strict=False)
    Chem.Kekulize(mol)
    # pattern used to find problematic conjugated bonds
    # there's usually an even number of matches for this
    pattern = Chem.MolFromSmarts("[*-{1-2}]-,=[*+0]=,#[*+0]")
    # pattern used to finish fixing a series of conjugated bonds
    base_end_pattern = Chem.MolFromSmarts(
        "[*-{1-2}]-,=[*+0]=,#[*+0]-,=[*-{1-2}]")
    # used when there's an odd number of matches for `pattern`
    odd_end_pattern = Chem.MolFromSmarts(
        "[*-]-[*+0]=[*+0]-[*-,$([#7;X3;v3]),$([#6+0,#7+1]=O),"
        "$([S;D4;v4]-[O-])]")
    # number of unique matches with the pattern
    n_matches = len(set([match[0]
                         for match in mol.GetSubstructMatches(pattern)]))
    # nothing to standardize
    if n_matches == 0:
        return
    # single match (unusual)
    elif n_matches == 1:
        # as a last resort, the only way to standardize is to find a nitrogen
        # that can accept a double bond and a positive charge
        # or a carbonyl that will become an enolate
        end_pattern = odd_end_pattern
    # at least 2 matches
    else:
        # give priority to base end pattern and then deal with the odd one if
        # necessary
        end_pattern = base_end_pattern
    backtrack = []
    backtrack_cycles = 0
    for _ in range(max_iter):
        # check for ending pattern
        end_match = mol.GetSubstructMatch(end_pattern)
        if end_match:
            # index of each atom
            anion1, a1, a2, anion2 = end_match
            term_atom = mol.GetAtomWithIdx(anion2)
            # edge-case 1: C-[O-] or [N+]-[O-]
            # [*-]-*=*-[C,N+]=O --> *=*-*=[C,N+]-[O-]
            # transform the =O to -[O-]
            if (
                term_atom.GetAtomicNum() == 6
                and term_atom.GetFormalCharge() == 0
            ) or (
                term_atom.GetAtomicNum() == 7
                and term_atom.GetFormalCharge() == 1
            ):
                for neighbor in term_atom.GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(anion2, neighbor.GetIdx())
                    if (neighbor.GetAtomicNum() == 8 and
                            bond.GetBondTypeAsDouble() == 2):
                        bond.SetBondType(Chem.BondType.SINGLE)
                        neighbor.SetFormalCharge(-1)
                        break
            # edge-case 2: S=O
            # [*-]-*=*-[Sv4]-[O-] --> *=*-*=[Sv6]=O
            # transform -[O-] to =O
            elif (term_atom.GetAtomicNum() == 16 and
                  term_atom.GetFormalCharge() == 0):
                for neighbor in term_atom.GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(anion2, neighbor.GetIdx())
                    if (neighbor.GetAtomicNum() == 8 and
                        neighbor.GetFormalCharge() == -1 and
                            bond.GetBondTypeAsDouble() == 1):
                        bond.SetBondType(Chem.BondType.DOUBLE)
                        neighbor.SetFormalCharge(0)
                        break
            # not an edge case:
            # increment the charge
            else:
                term_atom.SetFormalCharge(term_atom.GetFormalCharge() + 1)
            # common to all cases:
            # [*-]-*=*-[R] --> *=*-*=[R]
            # increment the charge and switch single and double bonds
            a = mol.GetAtomWithIdx(anion1)
            a.SetFormalCharge(a.GetFormalCharge() + 1)
            b = mol.GetBondBetweenAtoms(anion1, a1)
            b.SetBondType(RDBONDORDER[b.GetBondTypeAsDouble() + 1])
            b = mol.GetBondBetweenAtoms(a1, a2)
            b.SetBondType(RDBONDORDER[b.GetBondTypeAsDouble() - 1])
            b = mol.GetBondBetweenAtoms(a2, anion2)
            b.SetBondType(RDBONDORDER[b.GetBondTypeAsDouble() + 1])
            mol.UpdatePropertyCache(strict=False)
            continue

        # switch the position of the charge and the double bond
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            # check if we haven't already transformed this triplet
            for match in matches:
                # store order-independent atom indices
                g = set(match)
                # already transformed --> try the next one
                if g in backtrack:
                    continue
                # add to backtracking and start the switch
                else:
                    anion, a1, a2 = match
                    backtrack.append(g)
                    break
            # already performed all changes
            else:
                if backtrack_cycles == 1:
                    # might be stuck because there were 2 "odd" end patterns
                    # misqualified as a single base one
                    end_pattern = odd_end_pattern
                elif backtrack_cycles > 1:
                    # likely stuck changing the same bonds over and over so
                    # let's stop here
                    mol.UpdatePropertyCache(strict=False)
                    return
                match = matches[0]
                anion, a1, a2 = match
                backtrack = [set(match)]
                backtrack_cycles += 1

            # switch charges
            a = mol.GetAtomWithIdx(anion)
            a.SetFormalCharge(a.GetFormalCharge() + 1)
            a = mol.GetAtomWithIdx(a2)
            a.SetFormalCharge(a.GetFormalCharge() - 1)
            # switch bond orders
            b = mol.GetBondBetweenAtoms(anion, a1)
            b.SetBondType(RDBONDORDER[b.GetBondTypeAsDouble() + 1])
            b = mol.GetBondBetweenAtoms(a1, a2)
            b.SetBondType(RDBONDORDER[b.GetBondTypeAsDouble() - 1])
            mol.UpdatePropertyCache(strict=False)
            # update number of matches for the end pattern
            n_matches = len(set([match[0] for match in matches]))
            if n_matches == 1:
                end_pattern = odd_end_pattern
            # start new iteration
            continue

        # no more changes to apply
        mol.UpdatePropertyCache(strict=False)
        return

    # reached max_iter
    warnings.warn("The standardization could not be completed within a "
                  "reasonable number of iterations")
