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

Read coordinates data from an `RDKit <https://www.rdkit.org/docs/>`__
:class:`rdkit.Chem.rdchem.Mol` with
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

Instead of using the default bond-order and charge inferring algorithm that
relies on the topology alone and the presence of explicit hydrogen atoms, you
can also use alternative builtin algorithms (see the
:mod:`~MDAnalysis.converters.RDKitInferring` module) or even your own function
to modify the RDKit molecule::

    >>> template = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    >>> inferer = mda.converters.RDKitInferring.TemplateInferer(template=template)
    >>> u.atoms.convert_to.rdkit(inferer=inferer)
    <rdkit.Chem.rdchem.Mol at 0x7f70ee6f3ca0>
    >>> def dummy_inferer(mol):
    ...     # assign bond orders and do any other modification here
    ...     return mol
    >>> u.atoms.convert_to.rdkit(inferer=dummy_inferer)
    <rdkit.Chem.rdchem.Mol at 0x7f70ee6f2b90>


Classes
-------

.. autoclass:: RDKitReader
   :members:

.. autoclass:: RDKitConverter
   :members:


.. Links

.. _`ChEMBL27`: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_27/
.. _`Issue #3339`: https://github.com/MDAnalysis/mdanalysis/issues/3339
.. _`RDKitConverter benchmark`: https://github.com/MDAnalysis/RDKitConverter-benchmark
"""

import copy
import warnings
from contextlib import suppress
from functools import lru_cache
from io import StringIO

import numpy as np

from ..coordinates import memory
from ..coordinates.PDB import PDBWriter
from ..core.topologyattrs import _TOPOLOGY_ATTRS
from ..exceptions import NoDataError
from . import base
from .RDKitInferring import RDBONDORDER, MDAnalysisInferer

with suppress(ImportError):
    from rdkit import Chem

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
DEFAULT_INFERER = MDAnalysisInferer()
# only here for backwards compatibility
_infer_bo_and_charges = DEFAULT_INFERER._infer_bo_and_charges
_standardize_patterns = DEFAULT_INFERER._standardize_patterns
MONATOMIC_CATION_CHARGES = MDAnalysisInferer.MONATOMIC_CATION_CHARGES
STANDARDIZATION_REACTIONS = MDAnalysisInferer.STANDARDIZATION_REACTIONS

_deduce_PDB_atom_name = PDBWriter(StringIO())._deduce_PDB_atom_name


class RDKitReader(memory.MemoryReader):
    """Coordinate reader for RDKit.

    .. versionadded:: 2.0.0
    """

    format = "RDKIT"

    # Structure.coordinates always in Angstrom
    units = {"time": None, "length": "Angstrom"}

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
        coordinates = np.array(
            [conf.GetPositions() for conf in filename.GetConformers()],
            dtype=np.float32,
        )
        if coordinates.size == 0:
            warnings.warn("No coordinates found in the RDKit molecule")
            coordinates = np.empty((1, n_atoms, 3), dtype=np.float32)
            coordinates[:] = np.nan
        super(RDKitReader, self).__init__(coordinates, order="fac", **kwargs)


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
    use the parameter ``implicit_hydrogens=True`` when using the converter to
    allow implicit hydrogens, and ``inferer=None`` to disable inferring bond
    orders and charges. You can also pass your own callable function to
    ``inferer`` to assign bond orders and charges as you see fit::

        >>> from rdkit import Chem
        >>> from rdkit.Chem.AllChem import AssignBondOrdersFromTemplate
        >>> template = Chem.MolFromSmiles("NC(Cc1ccccc1)C(=O)O")
        >>> def infer_from_template(mol: Chem.Mol) -> Chem.Mol:
        ...     return AssignBondOrdersFromTemplate(template, mol)
        >>> mol = u.atoms.convert_to.rdkit(inferer=infer_from_template)


    Note that all builtin inferer functions can be found in the
    :mod:`RDKitInferring` module.

    Since one of the main use case of the converter is converting trajectories
    and not just a topology, creating a new molecule from scratch for every
    frame would be too slow so the converter uses a caching system. The cache
    only stores the 2 most recent AtomGroups that were converted, and is
    sensitive to the arguments that were passed to the converter. The number of
    objects cached can be changed with the function
    :func:`set_converter_cache_size`. However, ``ag.convert_to("RDKIT")``
    followed by ``ag.convert_to("RDKIT", implicit_hydrogens=False)`` will not
    use the cache since the arguments given are different. You can pass a
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

    .. versionchanged:: 2.8.0
        Deprecated ``max_iter`` (moved to the inferer class
        :class:`~MDAnalysis.converters.RDKitInferring.MDAnalysisInferer`) and
        deprecated ``NoImplicit`` in favor of ``implicit_hydrogens``. Added
        ``inferer`` to specify a callable that can transform the molecule (this
        operation is cached).

    """

    lib = "RDKIT"
    units = {"time": None, "length": "Angstrom"}

    def convert(
        self,
        obj,
        cache=True,
        implicit_hydrogens=False,
        force=False,
        inferer=DEFAULT_INFERER,
        **kwargs,
    ):
        """Write selection at current trajectory frame to
        :class:`~rdkit.Chem.rdchem.Mol`.

        Parameters
        -----------
        obj : :class:`~MDAnalysis.core.groups.AtomGroup` or
            :class:`~MDAnalysis.core.universe.Universe`
        cache : bool
            Use a cached copy of the molecule's topology when available. To be
            used, the cached molecule and the new one have to be made from the
            same AtomGroup selection and with the same arguments passed
            to the converter
        inferer : Optional[Callable[[Chem.Mol], Chem.Mol]]
            A callable to infer bond orders and charges for the RDKit molecule
            created by the converter. If ``None``, inferring is skipped.
        implicit_hydrogens : bool
           Whether to allow implicit hydrogens on the molecule or not.
        force : bool
            Force the conversion when no hydrogens were detected but
            ``inferer`` is not ``None``. Useful for inorganic molecules mostly.
        """

        try:
            from rdkit import Chem
        except ImportError:
            raise ImportError(
                "RDKit is required for the RDKitConverter but "
                "it's not installed. Try installing it with \n"
                "`conda install -c conda-forge rdkit` or `pip install rdkit`"
            ) from None
        try:
            # make sure to use atoms (Issue 46)
            ag = obj.atoms
        except AttributeError:
            raise TypeError(
                f"No `atoms` attribute in object of type {type(obj)!r}, "
                "please use a valid AtomGroup or Universe"
            ) from None

        if (max_iter := kwargs.get("max_iter")) is not None:
            warnings.warn(
                "Using `max_iter` is deprecated, use `MDAnalysisInferer"
                "(max_iter=...)` instead",
                DeprecationWarning,
            )
            if isinstance(inferer, MDAnalysisInferer):
                inferer = MDAnalysisInferer(max_iter=max_iter)

        if (NoImplicit := kwargs.get("NoImplicit")) is not None:
            warnings.warn(
                "Using `NoImplicit` is deprecated, use `implicit_hydrogens` "
                " instead. To disable bond order and formal charge inferring, "
                "use `inferer=None`",
                DeprecationWarning,
            )
            implicit_hydrogens = not NoImplicit
            # backwards compatibility
            if implicit_hydrogens:
                inferer = None

        # parameters passed to atomgroup_to_mol
        params = dict(
            implicit_hydrogens=implicit_hydrogens, force=force, inferer=inferer
        )
        if cache:
            mol = atomgroup_to_mol(ag, **params)
            mol = copy.deepcopy(mol)
        else:
            mol = atomgroup_to_mol.__wrapped__(ag, **params)

        # add a conformer for the current Timestep
        if hasattr(ag, "positions"):
            if np.isnan(ag.positions).any():
                warnings.warn(
                    "NaN detected in coordinates, the output "
                    "molecule will not have 3D coordinates assigned"
                )
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
def atomgroup_to_mol(
    ag,
    implicit_hydrogens=False,
    force=False,
    inferer=DEFAULT_INFERER,
    **kwargs,
):
    """Converts an AtomGroup to an RDKit molecule without coordinates.

    Parameters
    -----------
    ag : MDAnalysis.core.groups.AtomGroup
        The AtomGroup to convert
    implicit_hydrogens : bool
        Whether to allow implicit hydrogens on the molecule or not.
    force : bool
        Force the conversion when no hydrogens were detected but ``inferer``
        is not ``None``. Useful for inorganic molecules mostly.
    inferer : Optional[Callable[[rdkit.Chem.rdchem.Mol], rdkit.Chem.rdchem.Mol]]
        A callable to infer bond orders and charges for the RDKit molecule
        created by the converter. If ``None``, inferring is skipped.


    .. versionchanged:: 2.8.0
        Deprecated ``NoImplicit`` in favor of ``implicit_hydrogens``. Added
        ``inferer`` to specify a callable that can transform the molecule (this
        operation is cached).
    """
    try:
        elements = ag.elements
    except NoDataError:
        raise AttributeError(
            "The `elements` attribute is required for the RDKitConverter "
            "but is not present in this AtomGroup. Please refer to the "
            "documentation to guess elements from other attributes or "
            "type `help(MDAnalysis.topology.guessers)`"
        ) from None

    if "H" not in ag.elements:
        if force:
            warnings.warn(
                "No hydrogen atom found in the topology. "
                "Forcing to continue the conversion."
            )
        elif inferer is not None:
            raise AttributeError(
                "No hydrogen atom could be found in the topology, but the "
                "converter requires all hydrogens to be explicit. You can use "
                "the parameter ``inferer=None`` when using the converter "
                "to disable inferring bond orders and charges. You can also "
                "use ``force=True`` to ignore this error."
            )

    if (NoImplicit := kwargs.pop("NoImplicit", None)) is not None:
        warnings.warn(
            "Using `NoImplicit` is deprecated, use `implicit_hydrogens` "
            "instead. To disable bond order and formal charge inferring, use "
            "`inferer=None`",
            DeprecationWarning,
        )
        implicit_hydrogens = not NoImplicit
        # backwards compatibility
        if implicit_hydrogens:
            inferer = None
    if kwargs:
        raise ValueError(f"Found unexpected arguments: {kwargs}.")

    # attributes accepted in PDBResidueInfo object
    pdb_attrs = {}
    for attr in RDATTRIBUTES:
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
    NoImplicit = not implicit_hydrogens
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
        for attr, values in other_attrs.items():
            attr = f"_MDAnalysis_{_TOPOLOGY_ATTRS[attr].singular}"
            _set_atom_property(rdatom, attr, values[i])
        _set_atom_property(rdatom, "_MDAnalysis_index", i)
        # add atom
        index = mol.AddAtom(rdatom)
        atom_mapper[atom.ix] = index

    try:
        ag.bonds
    except NoDataError:
        warnings.warn(
            "No `bonds` attribute in this AtomGroup. Guessing bonds based "
            "on atoms coordinates"
        )
        ag.guess_bonds()

    for bond in ag.bonds:
        try:
            bond_indices = [atom_mapper[i] for i in bond.indices]
        except KeyError:
            continue
        bond_type = RDBONDORDER.get(bond.order, Chem.BondType.SINGLE)
        mol.AddBond(*bond_indices, bond_type)

    mol.UpdatePropertyCache(strict=False)

    if inferer is not None:
        # infer bond orders and formal charges
        mol = inferer(mol)

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
    global atomgroup_to_mol  # pylint: disable=global-statement
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
    getattr(mi, f"Set{rdattr}")(value)


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
