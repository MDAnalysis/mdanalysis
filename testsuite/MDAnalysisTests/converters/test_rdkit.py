# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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

import copy
import warnings
from contextlib import suppress
from io import StringIO

import MDAnalysis as mda
import numpy as np
import pytest
from MDAnalysis.topology.guessers import guess_atom_element
from MDAnalysisTests.datafiles import GRO, PDB_full, PDB_helix, mol2_molecule
from MDAnalysisTests.util import import_not_available
from numpy.testing import assert_allclose, assert_equal

with suppress(ImportError):
    from MDAnalysis.converters.RDKit import (
        RDATTRIBUTES,
        _add_mda_attr_to_rdkit,
        _set_atom_property,
        atomgroup_to_mol,
    )
    from MDAnalysis.converters.RDKitInferring import MDAnalysisInferer, TemplateInferer
    from rdkit import Chem
    from rdkit.Chem import AllChem


requires_rdkit = pytest.mark.skipif(
    import_not_available("rdkit"), reason="requires RDKit"
)


@pytest.mark.skipif(
    not import_not_available("rdkit"), reason="only for min dependencies build"
)
class TestRequiresRDKit(object):
    def test_converter_requires_rdkit(self):
        u = mda.Universe(PDB_full)
        with pytest.raises(
            ImportError, match="RDKit is required for the RDKitConverter"
        ):
            u.atoms.convert_to("RDKIT")


@requires_rdkit
class MolFactory:
    @staticmethod
    def mol2_mol():
        return Chem.MolFromMol2File(mol2_molecule, removeHs=False)

    @staticmethod
    def smiles_mol():
        mol = Chem.MolFromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMultipleConfs(mol, numConfs=3)
        return mol

    @staticmethod
    def dummy_product():
        mol = Chem.RWMol()
        atom = Chem.Atom(1)
        atom.SetIntProp("old_mapno", 0)
        atom.SetUnsignedProp("react_atom_idx", 0)
        mol.AddAtom(atom)
        return mol

    @staticmethod
    def dummy_reactant():
        mol = Chem.RWMol()
        atom = Chem.Atom(1)
        atom.SetIntProp("_MDAnalysis_index", 1)
        mol.AddAtom(atom)
        return mol


@pytest.fixture(scope="function")
def rdmol(request):
    return getattr(MolFactory, request.param)()


@pytest.fixture(scope="function")
def product(request):
    return getattr(MolFactory, request.param)()


def is_isomorphic(mol, ref, useChirality=False):
    return mol.HasSubstructMatch(
        ref, useChirality=useChirality
    ) and ref.HasSubstructMatch(mol, useChirality=useChirality)


@requires_rdkit
class TestRDKitReader(object):
    @pytest.mark.parametrize(
        "rdmol, n_frames",
        [
            ("mol2_mol", 1),
            ("smiles_mol", 3),
        ],
        indirect=["rdmol"],
    )
    def test_coordinates(self, rdmol, n_frames):
        universe = mda.Universe(rdmol)
        assert universe.trajectory.n_frames == n_frames
        expected = np.array(
            [conf.GetPositions() for conf in rdmol.GetConformers()], dtype=np.float32
        )
        assert_equal(expected, universe.trajectory.coordinate_array)

    def test_no_coordinates(self):
        with pytest.warns(UserWarning, match="No coordinates found"):
            u = mda.Universe.from_smiles("CCO", generate_coordinates=False)
        expected = np.empty((1, u.atoms.n_atoms, 3), dtype=np.float32)
        expected[:] = np.nan
        assert_equal(u.trajectory.coordinate_array, expected)

    def test_compare_mol2reader(self):
        universe = mda.Universe(MolFactory.mol2_mol())
        mol2 = mda.Universe(mol2_molecule)
        assert universe.trajectory.n_frames == mol2.trajectory.n_frames
        assert_equal(universe.trajectory.ts.positions, mol2.trajectory.ts.positions)


@requires_rdkit
class TestRDKitConverter(object):
    @pytest.fixture
    def pdb(self):
        return mda.Universe(PDB_full)

    @pytest.fixture
    def mol2(self):
        u = mda.Universe(mol2_molecule)
        # add elements
        elements = np.array(
            [guess_atom_element(x) for x in u.atoms.types], dtype=object
        )
        u.add_TopologyAttr("elements", elements)
        return u

    @pytest.fixture
    def peptide(self):
        u = mda.Universe(GRO)
        elements = mda.topology.guessers.guess_types(u.atoms.names)
        u.add_TopologyAttr("elements", elements)
        return u.select_atoms("resid 2-12")

    @pytest.fixture
    def uo2(self):
        return mda.Universe.from_smiles("O=O")

    def test_no_atoms_attr(self):
        rdkit_converter = mda._CONVERTERS["RDKIT"]().convert
        with pytest.raises(TypeError, match=""):
            rdkit_converter("foo")

    @pytest.mark.parametrize("smi", ["[H]", "C", "O", "[He]"])
    def test_single_atom_mol(self, smi):
        u = mda.Universe.from_smiles(smi, addHs=False, generate_coordinates=False)
        mol = u.atoms.convert_to.rdkit(inferer=None)
        assert mol.GetNumAtoms() == 1
        assert mol.GetAtomWithIdx(0).GetSymbol() == smi.strip("[]")

    @pytest.mark.parametrize(
        "resname, n_atoms, n_fragments",
        [
            ("PRO", 14, 1),
            ("ILE", 38, 1),
            ("ALA", 20, 2),
            ("GLY", 21, 3),
        ],
    )
    def test_mol_from_selection(self, peptide, resname, n_atoms, n_fragments):
        mol = peptide.select_atoms("resname %s" % resname).convert_to("RDKIT")
        assert n_atoms == mol.GetNumAtoms()
        assert n_fragments == len(Chem.GetMolFrags(mol))

    @pytest.mark.parametrize(
        "sel_str, atom_index",
        [
            ("resid 1", 0),
            ("resname LYS and name NZ", 1),
            ("resid 34 and altloc B", 2),
        ],
    )
    def test_monomer_info(self, pdb, sel_str, atom_index):
        sel = pdb.select_atoms(sel_str)
        mda_atom = sel.atoms[atom_index]
        umol = sel.convert_to.rdkit(inferer=None)
        atom = umol.GetAtomWithIdx(atom_index)
        mi = atom.GetMonomerInfo()
        assert mda_atom.altLoc == mi.GetAltLoc()
        assert mda_atom.chainID == mi.GetChainId()
        assert mda_atom.icode == mi.GetInsertionCode()
        assert mda_atom.name == mi.GetName().strip()
        assert mda_atom.occupancy == mi.GetOccupancy()
        assert mda_atom.resname == mi.GetResidueName()
        assert mda_atom.resid == mi.GetResidueNumber()
        assert mda_atom.segindex == mi.GetSegmentNumber()
        assert mda_atom.tempfactor == mi.GetTempFactor()

    @pytest.mark.parametrize("rdmol", ["mol2_mol", "smiles_mol"], indirect=True)
    def test_identical_topology(self, rdmol):
        u = mda.Universe(rdmol)
        umol = u.atoms.convert_to("RDKIT")
        assert is_isomorphic(rdmol, umol)
        u2 = mda.Universe(umol)
        assert_equal(u.atoms.bonds, u2.atoms.bonds)
        assert_equal(u.atoms.elements, u2.atoms.elements)
        assert_equal(u.atoms.names, u2.atoms.names)
        assert_allclose(u.atoms.positions, u2.atoms.positions, rtol=0, atol=1e-7)

    def test_raise_requires_elements(self):
        u = mda.Universe(mol2_molecule)

        # Delete topology attribute (PR #3069)
        u.del_TopologyAttr("elements")

        with pytest.raises(
            AttributeError,
            match="`elements` attribute is required for the RDKitConverter",
        ):
            u.atoms.convert_to("RDKIT")

    def test_warn_guess_bonds(self):
        u = mda.Universe(PDB_helix)
        with pytest.warns(UserWarning, match="No `bonds` attribute in this AtomGroup"):
            u.atoms.convert_to("RDKIT")

    def test_bonds_outside_sel(self):
        u = mda.Universe(Chem.MolFromSmiles("CC(=O)C"))
        ag = u.select_atoms("index 1")
        ag.convert_to.rdkit(inferer=None)

    def test_error_no_hydrogen(self, uo2):
        with pytest.raises(
            AttributeError,
            match="the converter requires all hydrogens to be " "explicit",
        ):
            uo2.atoms.convert_to("RDKIT")

    def test_error_no_hydrogen_skip_inferring(self, uo2):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            uo2.atoms.convert_to.rdkit(inferer=None)

    def test_warning_no_hydrogen_force(self, uo2):
        with pytest.warns(UserWarning, match="Forcing to continue the conversion"):
            uo2.atoms.convert_to.rdkit(force=True)

    @pytest.mark.parametrize(
        "attr, value, expected",
        [
            ("names", "N", " N  "),
            ("names", "CA", " CA "),
            ("names", "CAT", " CAT"),
            ("names", "N1", " N1 "),
            ("names", "CE2", " CE2"),
            ("names", "C12", " C12"),
            ("names", "HD12", "HD12"),
            ("names", "C123", "C123"),
            ("altLocs", "A", "A"),
            ("chainIDs", "B", "B"),
            ("icodes", "C", "C"),
            ("occupancies", 0.5, 0.5),
            ("resnames", "LIG", "LIG"),
            ("resids", 123, 123),
            ("segindices", 1, 1),
            ("tempfactors", 0.8, 0.8),
        ],
    )
    def test_add_mda_attr_to_rdkit(self, attr, value, expected):
        mi = Chem.AtomPDBResidueInfo()
        _add_mda_attr_to_rdkit(attr, value, mi)
        rdvalue = getattr(mi, "Get%s" % RDATTRIBUTES[attr])()
        assert rdvalue == expected

    @pytest.mark.parametrize("idx", [0, 10, 42])
    def test_other_attributes(self, mol2, idx):
        mol = mol2.atoms.convert_to("RDKIT")
        rdatom = mol.GetAtomWithIdx(idx)
        mda_atom = mol2.atoms[idx]
        assert mda_atom.charge == rdatom.GetDoubleProp("_MDAnalysis_charge")
        assert mda_atom.segid == rdatom.GetProp("_MDAnalysis_segid")
        assert mda_atom.type == rdatom.GetProp("_MDAnalysis_type")

    @pytest.mark.parametrize(
        "sel_str",
        [
            "resname ALA",
            "resname PRO and segid A",
        ],
    )
    def test_index_property(self, pdb, sel_str):
        ag = pdb.select_atoms(sel_str)
        mol = ag.convert_to.rdkit(inferer=None)
        expected = [i for i in range(len(ag))]
        indices = [a.GetIntProp("_MDAnalysis_index") for a in mol.GetAtoms()]
        assert_equal(indices, expected)

    def test_assign_coordinates(self, pdb):
        mol = pdb.atoms.convert_to.rdkit(inferer=None)
        positions = mol.GetConformer().GetPositions()
        assert_allclose(positions, pdb.atoms.positions)

    def test_assign_stereochemistry(self, mol2):
        umol = mol2.atoms.convert_to("RDKIT")
        rdmol = Chem.MolFromMol2File(mol2_molecule, removeHs=False)
        assert is_isomorphic(rdmol, umol, useChirality=True)

    def test_trajectory_coords(self):
        u = mda.Universe.from_smiles(
            "CCO", numConfs=3, rdkit_kwargs=dict(randomSeed=42)
        )
        for ts in u.trajectory:
            mol = u.atoms.convert_to("RDKIT")
            positions = mol.GetConformer().GetPositions()
            assert_allclose(positions, ts.positions)

    def test_nan_coords(self):
        u = mda.Universe.from_smiles("CCO")
        xyz = u.atoms.positions
        xyz[0][2] = np.nan
        u.atoms.positions = xyz
        with pytest.warns(UserWarning, match="NaN detected"):
            mol = u.atoms.convert_to("RDKIT")
        with pytest.raises(ValueError, match="Bad Conformer Id"):
            mol.GetConformer()

    def test_cache(self):
        """5 tests are performed here:

        * while iterating on timesteps in a trajectory, the number of cached
          objects should not change
        * the cache should remember about 2 atomgroups by default
        * the cache size can be increased
        * the cache is sensitive to arguments passed to the converter as they
          might change the output molecule
        * the cache can be ignored
        """
        cached_func = mda.converters.RDKit.atomgroup_to_mol
        # create universes
        utraj = mda.Universe.from_smiles("CCO", numConfs=5)
        uc = mda.Universe.from_smiles("C")
        ucc = mda.Universe.from_smiles("CC")
        uccc = mda.Universe.from_smiles("CCC")
        # test (1): iterating over frames in a trajectory
        previous_cache = None
        for ts in utraj.trajectory:
            utraj.atoms.convert_to("RDKIT")
            cache = cached_func.cache_info()
            if previous_cache:
                # the cache shouldn't change when iterating on timesteps
                assert cache.currsize == previous_cache.currsize
                previous_cache = copy.deepcopy(cache)
        # test (2): only 2 molecules should be cached by default
        cached_func.cache_clear()
        uc.atoms.convert_to("RDKIT")
        ucc.atoms.convert_to("RDKIT")
        uccc.atoms.convert_to("RDKIT")
        cache = cached_func.cache_info()
        assert cache.currsize == 2
        assert cache.misses == 3
        ucc.atoms.convert_to("RDKIT")  # should be inside of the cache
        assert cached_func.cache_info().hits == 1
        uc.atoms.convert_to("RDKIT")  # outside of the cache
        assert cached_func.cache_info().hits == 1
        assert cached_func.cache_info().misses == 4
        # test (3): increase cache size
        mda.converters.RDKit.set_converter_cache_size(3)
        cached_func = mda.converters.RDKit.atomgroup_to_mol
        assert cached_func.cache_info().maxsize == 3
        # test (4): caching is sensitive to converter arguments
        previous_cache = cached_func.cache_info()
        uc.atoms.convert_to.rdkit(implicit_hydrogens=True)
        cache = cached_func.cache_info()
        assert cache.misses == previous_cache.misses + 1
        # test (5): skip cache
        uc.atoms.convert_to.rdkit(cache=False)
        new_cache = cached_func.cache_info()
        assert cache == new_cache

    def test_sanitize_fail_warning(self):
        params = Chem.SmilesParserParams()
        params.removeHs = False
        params.sanitize = False
        mol = Chem.MolFromSmiles("[H]-N(-[H])(-[H])-[H]", params)
        u = mda.Universe(mol)
        with pytest.warns() as record:
            u.atoms.convert_to.rdkit(inferer=None)
        assert "Could not sanitize molecule" not in "\n".join(
            [str(r.message) for r in record]
        )

    def test_deprecation_maw_iter(self, mol2):
        with pytest.warns(DeprecationWarning, match="Using `max_iter` is deprecated"):
            mol2.atoms.convert_to.rdkit(max_iter=2)

    def test_deprecation_NoImplicit(self, mol2):
        with pytest.warns(DeprecationWarning, match="Using `NoImplicit` is deprecated"):
            mol2.atoms.convert_to.rdkit(NoImplicit=True)

    def test_deprecation_atomgroup_to_mol_NoImplicit(self, mol2):
        with pytest.warns(DeprecationWarning, match="Using `NoImplicit` is deprecated"):
            atomgroup_to_mol(mol2.atoms, NoImplicit=False)

    def test_atomgroup_to_mol_unexpected_kwargs(self, mol2):
        with pytest.raises(
            ValueError, match="Found unexpected arguments: {'foo': 'bar'}"
        ):
            atomgroup_to_mol(mol2.atoms, foo="bar")

    def test_custom_callable_inferer(self, mol2):
        def potatoe(mol):
            new = Chem.Mol(mol)
            new.SetProp("_Name", "🥔")
            return new

        mol = mol2.atoms.convert_to.rdkit(inferer=potatoe)
        assert mol.GetProp("_Name") == "🥔"


@requires_rdkit
class TestRDKitMDAnalysisInferer(object):
    def add_Hs_remove_bo_and_charges(self, mol):
        """Add hydrogens and remove bond orders and charges from a molecule"""
        mH = Chem.AddHs(mol)
        for atom in mH.GetAtoms():
            atom.SetIsAromatic(False)
            atom.SetFormalCharge(0)
            atom.SetNoImplicit(True)
        for bond in mH.GetBonds():
            bond.SetIsAromatic(False)
            bond.SetBondType(Chem.BondType.SINGLE)
        mH.UpdatePropertyCache(strict=False)
        return mH

    def enumerate_reordered_mol(self, mol):
        """Enumerates all possible starting atoms for a given molecule"""
        # go through each possible starting atom
        for root_atom in mol.GetAtoms():
            smi = Chem.MolToSmiles(mol, rootedAtAtom=root_atom.GetIdx())
            reordered_mol = Chem.MolFromSmiles(smi, sanitize=False)
            for atom in reordered_mol.GetAtoms():
                atom.SetNoImplicit(True)
            reordered_mol.UpdatePropertyCache(strict=False)
            yield reordered_mol

    def assign_bond_orders_and_charges(self, mol):
        """Returns a sanitized molecule with infered bond orders and charges"""
        inferer = MDAnalysisInferer()
        inferer._infer_bo_and_charges(mol)
        mol = inferer._standardize_patterns(mol)
        Chem.SanitizeMol(mol)
        return mol

    def assert_isomorphic_resonance_structure(self, mol, ref):
        """Checks if 2 molecules are isomorphic using their resonance
        structures
        """
        isomorphic = mol.HasSubstructMatch(ref)
        if not isomorphic:
            isomorphic = bool(Chem.ResonanceMolSupplier(mol).GetSubstructMatch(ref))
        assert isomorphic, f"{Chem.MolToSmiles(ref)} != {Chem.MolToSmiles(mol)}"

    @pytest.mark.parametrize(
        "smi, out",
        [
            ("C(-[H])(-[H])(-[H])-[H]", "C"),
            ("[C](-[H])(-[H])-[C](-[H])-[H]", "C=C"),
            (
                "[C]1(-[H])-[C](-[H])-[C](-[H])-[C](-[H])-[C](-[H])-[C]1(-[H])",
                "c1ccccc1",
            ),
            ("C-[C](-[H])-[O]", "C(=O)C"),
            ("[H]-[C](-[O])-[N](-[H])-[H]", "C(=O)N"),
            ("[N]-[C]-[H]", "N#C"),
            ("C-[C](-[O]-[H])-[O]", "CC(=O)O"),
            ("[P](-[O]-[H])(-[O]-[H])(-[O]-[H])-[O]", "P(O)(O)(O)=O"),
            ("[P](-[O]-[H])(-[O]-[H])(-[O])-[O]", "P([O-])(O)(O)=O"),
            ("[P](-[O]-[H])(-[O])(-[O])-[O]", "P([O-])([O-])(O)=O"),
            ("[P](-[O])(-[O])(-[O])-[O]", "P([O-])([O-])([O-])=O"),
            ("[H]-[O]-[N]-[O]", "ON=O"),
            ("[N]-[C]-[O]", "N#C[O-]"),
        ],
    )
    def test_infer_bond_orders(self, smi, out):
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        inferer = MDAnalysisInferer()
        inferer._infer_bo_and_charges(mol)
        Chem.SanitizeMol(mol)
        mol = Chem.RemoveHs(mol)
        molref = Chem.MolFromSmiles(out)
        assert is_isomorphic(mol, molref), "{} != {}".format(Chem.MolToSmiles(mol), out)

    @pytest.mark.parametrize(
        "smi, atom_idx, charge",
        [
            ("[C](-[H])(-[H])(-[H])-[O]", 4, -1),
            ("[N]-[C]-[O]", 2, -1),
            ("[N](-[H])(-[H])(-[H])-[H]", 0, 1),
            ("C-[C](-[O])-[O]", 3, -1),
        ],
    )
    def test_infer_charges(self, smi, atom_idx, charge):
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        inferer = MDAnalysisInferer()
        inferer._infer_bo_and_charges(mol)
        Chem.SanitizeMol(mol)
        assert mol.GetAtomWithIdx(atom_idx).GetFormalCharge() == charge

    @pytest.mark.parametrize(
        "smi, out",
        [
            ("[S](-[O]-[H])(-[O]-[H])(-[O])-[O]", "S(=O)(=O)(O)O"),
            ("[S](-[O]-[H])(-[O])(-[O])-[O]", "S(=O)(=O)([O-])O"),
            ("[S](-[O])(-[O])(-[O])-[O]", "S(=O)(=O)([O-])[O-]"),
            ("C-[N](-[H])-[C](-[N](-[H])-[H])-[N](-[H])-[H]", "CNC(N)=[N+](-[H])-[H]"),
            ("[O]-[C](-[H])-[C](-[H])-[H]", "C([O-])=C"),
            ("C-[N](-[O])-[O]", "C[N+](=O)[O-]"),
            ("C(-[N](-[O])-[O])-[N](-[O])-[O]", "C([N+](=O)[O-])[N+](=O)[O-]"),
            ("C-[N](-[O])-[O].C-[N](-[O])-[O]", "C[N+](=O)[O-].C[N+](=O)[O-]"),
            ("[C-](=O)-C", "[C](=O)-C"),
            ("[H]-[N-]-C", "[H]-[N]-C"),
            (
                "[O]-[C]1-[C](-[H])-[C](-[H])-[C](-[H])-[C](-[H])-[C](-[H])1",
                "[O-]c1ccccc1",
            ),
            ("[O]-[C]1-[C](-[H])-[C](-[H])-[C](-[H])-[C]1-[O]", "[O-]C1=CC=CC1=O"),
            ("[H]-[C]-[C]-[C](-[H])-[C](-[H])-[H]", "C#CC=C"),
            ("[H]-[C]-[C]-[C]-[C]-[H]", "C#CC#C"),
        ],
    )
    def test_standardize_patterns(self, smi, out):
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        mol = self.assign_bond_orders_and_charges(mol)
        mol = Chem.RemoveHs(mol)
        molref = Chem.MolFromSmiles(out)
        assert is_isomorphic(mol, molref), "{} != {}".format(Chem.MolToSmiles(mol), out)

    @pytest.mark.parametrize(
        "attr, value, getter",
        [
            ("index", 42, "GetIntProp"),
            ("index", np.int8(42), "GetIntProp"),
            ("index", np.int16(42), "GetIntProp"),
            ("index", np.int32(42), "GetIntProp"),
            ("index", np.int64(42), "GetIntProp"),
            ("index", np.uint8(42), "GetIntProp"),
            ("index", np.uint16(42), "GetIntProp"),
            ("index", np.uint32(42), "GetIntProp"),
            ("index", np.uint64(42), "GetIntProp"),
            ("charge", 4.2, "GetDoubleProp"),
            ("charge", np.float32(4.2), "GetDoubleProp"),
            ("charge", np.float64(4.2), "GetDoubleProp"),
            ("type", "C.3", "GetProp"),
        ],
    )
    def test_set_atom_property(self, attr, value, getter):
        atom = Chem.Atom(1)
        prop = "_MDAnalysis_%s" % attr
        _set_atom_property(atom, prop, value)
        assert getattr(atom, getter)(prop) == value

    def test_ignore_prop(self):
        atom = Chem.Atom(1)
        _set_atom_property(atom, "foo", {"bar": "baz"})
        assert "foo" not in atom.GetPropsAsDict().items()

    @pytest.mark.parametrize(
        "smi",
        [
            "c1ccc(cc1)-c1ccccc1-c1ccccc1",
            "c1cc[nH]c1",
            "O=C([C@H](CC1=C[NH1+]=CN1)[NH3+])[O-]",
        ],
    )
    def test_transfer_properties(self, smi):
        mol = Chem.MolFromSmiles(smi)
        mol = self.add_Hs_remove_bo_and_charges(mol)
        old = {}
        for atom in mol.GetAtoms():
            ix = atom.GetIdx()
            atom.SetIntProp("_MDAnalysis_index", ix)
            atom.SetProp("dummy", f"foo_{ix}")
            old[ix] = {"_MDAnalysis_index": ix, "dummy": f"foo_{ix}"}
        newmol = self.assign_bond_orders_and_charges(mol)
        new = {}
        for a in newmol.GetAtoms():
            ix = a.GetIntProp("_MDAnalysis_index")
            new[ix] = {"_MDAnalysis_index": ix, "dummy": a.GetProp("dummy")}
            props = a.GetPropsAsDict().keys()
            assert "old_mapno" not in props
            assert "react_atom_idx" not in props
        assert new == old

    @pytest.mark.parametrize(
        "smi",
        [
            "c1ccc(cc1)-c1ccccc1-c1ccccc1",
            "c1cc[nH]c1",
            "c1ccc(cc1)-c1ccc(-c2ccccc2)c(-c2ccccc2)c1-c1ccccc1",
            "c1ccc2c(c1)c1ccccc1c1ccccc1c1ccccc1c1ccccc21",
            "c1csc(c1)-c1ccoc1-c1cc[nH]c1",
            "C1=C2C(=NC=N1)N=CN2",
            "CN1C=NC(=C1SC2=NC=NC3=C2NC=N3)[N+](=O)[O-]",
            "c1c[nH]c(c1)-c1ccc(s1)-c1ccoc1-c1c[nH]cc1-c1ccccc1",
            "C=CC=CC=CC=CC=CC=C",
            "NCCCCC([NH3+])C(=O)[O-]",
            "CC(C=CC1=C(C)CCCC1(C)C)=CC=CC(C)=CC=[NH+]C",
            "C#CC=C",
            # HID HIE HIP residues, see PR #2941
            "O=C([C@H](CC1=CNC=N1)N)O",
            "O=C([C@H](CC1=CN=CN1)N)O",
            "O=C([C@H](CC1=C[NH1+]=CN1)[NH3+])[O-]",
            # fixes from PR #3044
            "CCOC(=O)c1cc2cc(C(=O)O)ccc2[nH]1",
            "[O-][n+]1cccnc1",
            "C[n+]1ccccc1",
            "[PH4+]",
            "c1nc[nH]n1",
            "CC(=O)C=C(C)N",
            "CC(C)=CC=C[O-]",
            "O=S(C)(C)=NC",
            "Cc1ccc2c3ncn(Cc4ccco4)c(O)c-3nc2c1",
            "CCCC/C=C/C#CC#CCCCCCCCC(=O)O",
            "c1c2c(=O)n3cccc(C)c3nc2n(C)c(=N)c1C(=O)NCc1cnccc1",
            "N#Cc1c[nH]c(C(=O)NC(=O)c2cc[n+]([O-])cc2)n1",
            "C[C@@H](Oc1cc(F)ccc1Nc1ncnc2cc(N=S3(=O)CCC3)cc(F)c12)C(=O)NCC#N",
            "[O-][n+]1onc2ccccc21",
            "Cc1cc[n+](CC[n+]2ccc(C)cc2)cc1",
            "[O-]c1ccccc1",
            "[O-]C=CC=CCC=CC=[N+](C)C",
            "C=[N+](-[O-])-C",
            "C-[N-]-C(=O)C",
            # amino acids
            "C[C@H](N)C(=O)O",  # A
            "NCC(=O)O",  # G
            "CC[C@H](C)[C@H](N)C(=O)O",  # I
            "CC(C)C[C@H](N)C(=O)O",  # L
            "O=C(O)[C@@H]1CCCN1",  # P
            "CC(C)[C@H](N)C(=O)O",  # V
            "N[C@@H](Cc1ccccc1)C(=O)O",  # F
            "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O",  # W
            "N[C@@H](Cc1ccc(O)cc1)C(=O)O",  # Y
            "N[C@@H](CC(=O)O)C(=O)O",  # D
            "N[C@@H](CCC(=O)O)C(=O)O",  # E
            "N=C(N)NCCC[C@H](N)C(=O)O",  # R
            "N[C@@H](Cc1c[nH]cn1)C(=O)O",  # H
            "NCCCC[C@H](N)C(=O)O",  # K
            "N[C@@H](CO)C(=O)O",  # S
            "C[C@@H](O)[C@H](N)C(=O)O",  # T
            "N[C@@H](CS)C(=O)O",  # C
            "CSCC[C@H](N)C(=O)O",  # M
            "NC(=O)C[C@H](N)C(=O)O",  # N
            "NC(=O)CC[C@H](N)C(=O)O",  # Q
            # nucleic acids
            "Nc1ncnc2c1ncn2[C@H]1C[C@H](OP(=O)(O)O)[C@@H](COP(=O)(O)O)O1",  # A
            "Cc1cn([C@H]2C[C@H](OP(=O)(O)O)[C@@H](COP(=O)(O)O)O2)c(=O)[nH]c1=O",  # T
            "O=c1ccn([C@H]2C[C@H](OP(=O)(O)O)[C@@H](COP(=O)(O)O)O2)c(=O)[nH]1",  # U
            "Nc1ccn([C@H]2C[C@H](OP(=O)(O)O)[C@@H](COP(=O)(O)O)O2)c(=O)n1",  # C
            "Nc1nc2c(ncn2[C@H]2C[C@H](OP(=O)(O)O)[C@@H](COP(=O)(O)O)O2)c(=O)[nH]1",  # G
            # RDKitConverter benchmark 2022-05-23, fixed by better sorting
            "CC(C)NS(C)(=O)=Nc1ccc(-c2ccc(-c3nc4[nH]c(O[C@@H]5CO[C@@H]6[C@H](O)CO[C@H]56)nc4cc3Cl)cc2)cc1 CHEMBL4111464",
            "C[n+]1c2ccccc2c(C(=O)O)c2[nH]c3ccc(Br)cc3c21 CHEMBL511591",
            "C[n+]1ccccc1-c1ccc(NC(=O)c2ccc(C(=O)Nc3ccc(-c4cccc[n+]4C)cc3)cc2)cc1 CHEMBL3230729",
            "Cc1ccc(/C(O)=C(/C([S-])=NCc2ccccc2)[n+]2cccc(CO)c2)cc1[N+](=O)[O-] CHEMBL1596426",
            "CC(C)(O)c1cc2c(cc1C(C)(C)O)-c1ccccc1-2 CHEMBL1990966",
            "[O-]/N=C1/c2ccccc2-c2nc3ccc(C(F)(F)F)cc3nc21 CHEMBL4557878",
            "N#Cc1c[nH]c(C(=O)Nc2ccc(-c3cccc[n+]3[O-])cc2C2=CCCCC2)n1 CHEMBL1172116",
        ],
    )
    def test_order_independant(self, smi):
        # generate mol with hydrogens but without bond orders
        ref = Chem.MolFromSmiles(smi)
        stripped_mol = self.add_Hs_remove_bo_and_charges(ref)
        # enumerate mols with different starting atoms
        for m in self.enumerate_reordered_mol(stripped_mol):
            m = self.assign_bond_orders_and_charges(m)
            m = Chem.RemoveHs(m)
            self.assert_isomorphic_resonance_structure(m, ref)

    @pytest.mark.xfail(reason="Not currently tackled by the RDKitConverter")
    @pytest.mark.parametrize(
        "smi",
        [
            "C-[N+]#N",
            "C-N=[N+]=[N-]",
            "C-[O+]=C",
            "C-[N+]#[C-]",
        ],
    )
    def test_order_independant_issue_3339(self, smi):
        self.test_order_independant(smi)

    def test_warn_conjugated_max_iter(self):
        smi = "[C-]C=CC=CC=CC=CC=CC=C[C-]"
        mol = Chem.MolFromSmiles(smi)
        inferer = MDAnalysisInferer(max_iter=2)
        with pytest.warns(UserWarning, match="reasonable number of iterations"):
            inferer._rebuild_conjugated_bonds(mol)

    def test_deprecation_warning_max_iter(self):
        smi = "c1ccccc1"
        mol = Chem.MolFromSmiles(smi)
        inferer = MDAnalysisInferer()
        with pytest.warns(
            DeprecationWarning,
            match="Specifying `max_iter` is deprecated and will be removed",
        ):
            inferer._standardize_patterns(mol, max_iter=1)

    @pytest.mark.parametrize(
        "smi",
        [
            "[Li+]",
            "[Na+]",
            "[K+]",
            "[Rb+]",
            "[Ag+]",
            "[Cs+]",
            "[Mg+2]",
            "[Ca+2]",
            "[Cu+2]",
            "[Zn+2]",
            "[Sr+2]",
            "[Ba+2]",
            "[Al+3]",
            "[Fe+2]",
            "[Cl-]",
            "[O-2]",
            "[Na+].[Cl-]",
        ],
    )
    def test_ions(self, smi):
        ref = Chem.MolFromSmiles(smi)
        stripped_mol = self.add_Hs_remove_bo_and_charges(ref)
        mol = self.assign_bond_orders_and_charges(stripped_mol)
        assert is_isomorphic(mol, ref)

    @pytest.mark.parametrize(
        "smi",
        [
            "O=C([C@H](CC1=C[NH1+]=CN1)[NH3+])[O-]",
            "O=S(C)(C)=NC",
        ],
    )
    def test_reorder_atoms(self, smi):
        mol = Chem.MolFromSmiles(smi)
        mol = Chem.AddHs(mol)
        # remove bond order and charges info
        pdb = Chem.MolToPDBBlock(mol)
        u = mda.Universe(StringIO(pdb), format="PDB")
        # atoms are reordered during infering, and will be reordered back to
        # their original order with Chem.RenumberAtoms
        mu = u.atoms.convert_to.rdkit()
        values = [a.GetSymbol() for a in mu.GetAtoms()]
        expected = [a.GetSymbol() for a in mol.GetAtoms()]
        assert values == expected

    def test_pdb_names(self):
        u = mda.Universe(PDB_helix)
        mol = u.atoms.convert_to.rdkit()
        names = u.atoms.names
        rd_names = np.array([a.GetProp("_MDAnalysis_name") for a in mol.GetAtoms()])
        assert (names == rd_names).all()

    @pytest.mark.parametrize(
        "smi",
        [
            r"F/C(Br)=C(Cl)/I",
            r"F\C(Br)=C(Cl)\I",
            "F-C(Br)=C(Cl)-I",
        ],
    )
    def test_bond_stereo_not_stereoany(self, smi):
        u = mda.Universe.from_smiles(smi)
        mol = u.atoms.convert_to.rdkit(force=True)
        for bond in mol.GetBonds():
            if bond.GetBondTypeAsDouble() == 2:
                assert bond.GetStereo() != Chem.BondStereo.STEREOANY

    def test_atom_sorter(self):
        mol = Chem.MolFromSmiles("[H]-[C](-[H])-[C](-[H])-[C]-[C]-[H]", sanitize=False)
        # corresponding mol: C=C-C#C
        # atom indices:      1 3 5 6
        mol.UpdatePropertyCache()
        sorted_atoms = sorted(
            [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1],
            key=MDAnalysisInferer._atom_sorter,
        )
        sorted_indices = [atom.GetIdx() for atom in sorted_atoms]
        assert sorted_indices == [6, 5, 1, 3]


@requires_rdkit
class TestRDKitTemplateInferer:
    @pytest.fixture(scope="function")
    def template(self):
        params = Chem.SmilesParserParams()
        params.removeHs = False
        return Chem.MolFromSmiles("[H]-[N+](-[H])(-[H])-[H]", params)

    def test_template_inferring(self, template):
        params = Chem.SmilesParserParams()
        params.removeHs = False
        params.sanitize = False
        mol = Chem.MolFromSmiles("[H]-N(-[H])(-[H])-[H]", params)
        assert mol.GetAtomWithIdx(1).GetFormalCharge() == 0

        inferer = TemplateInferer(template=template)
        mol = inferer(mol)
        assert mol.GetAtomWithIdx(1).GetFormalCharge() == 1

    def test_from_convert_to_method(self, template):
        u = mda.Universe(template)
        mol = u.atoms.convert_to.rdkit(inferer=TemplateInferer(template=template))
        assert mol.GetAtomWithIdx(1).GetFormalCharge() == 1
