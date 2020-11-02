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
import pytest
import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element
import numpy as np
from numpy.testing import (assert_equal,
                           assert_almost_equal)

from MDAnalysisTests.datafiles import mol2_molecule, PDB_full, GRO, PDB_helix
from MDAnalysisTests.util import import_not_available


try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from MDAnalysis.coordinates.RDKit import (
        RDATTRIBUTES,
        _add_mda_attr_to_rdkit,
        _infer_bo_and_charges,
        _standardize_patterns,
        _rebuild_conjugated_bonds,
        _set_atom_property,
        _reassign_props_after_reaction,
    )
except ImportError:
    pass


requires_rdkit = pytest.mark.skipif(import_not_available("rdkit"),
                                    reason="requires RDKit")


@pytest.mark.skipif(not import_not_available("rdkit"),
                    reason="only for min dependencies build")
class TestRequiresRDKit(object):
    def test_converter_requires_rdkit(self):
        u = mda.Universe(PDB_full)
        with pytest.raises(ImportError,
                           match="RDKit is required for the RDKitConverter"):
            u.atoms.convert_to("RDKIT")


@requires_rdkit
class MolFactory:
    def mol2_mol():
        return Chem.MolFromMol2File(mol2_molecule, removeHs=False)

    def smiles_mol():
        mol = Chem.MolFromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
        mol = Chem.AddHs(mol)
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=3)
        return mol

    def dummy_product():
        mol = Chem.RWMol()
        atom = Chem.Atom(1)
        atom.SetIntProp("old_mapno", 0)
        atom.SetUnsignedProp("react_atom_idx", 0)
        mol.AddAtom(atom)
        return mol

    def dummy_product_nomap():
        mol = Chem.RWMol()
        atom = Chem.Atom(1)
        atom.SetUnsignedProp("react_atom_idx", 0)
        mol.AddAtom(atom)
        return mol

    def dummy_reactant_noprops():
        mol = Chem.RWMol()
        atom = Chem.Atom(1)
        mol.AddAtom(atom)
        return mol

    def dummy_reactant():
        mol = Chem.RWMol()
        atom = Chem.Atom(1)
        atom.SetProp("foo", "bar")
        atom.SetIntProp("_MDAnalysis_index", 1)
        atom.SetDoubleProp("_MDAnalysis_charge", 4.2)
        atom.SetProp("_MDAnalysis_type", "C.3")
        mol.AddAtom(atom)
        return mol


@pytest.fixture(scope="function")
def rdmol(request):
    return getattr(MolFactory, request.param)()


@pytest.fixture(scope="function")
def product(request):
    return getattr(MolFactory, request.param)()


@requires_rdkit
class TestRDKitReader(object):
    @pytest.mark.parametrize("rdmol, n_frames", [
        ("mol2_mol", 1),
        ("smiles_mol", 3),
    ], indirect=["rdmol"])
    def test_coordinates(self, rdmol, n_frames):
        universe = mda.Universe(rdmol)
        assert universe.trajectory.n_frames == n_frames
        expected = np.array([
            conf.GetPositions() for conf in rdmol.GetConformers()],
            dtype=np.float32)
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
        assert_equal(universe.trajectory.ts.positions,
                     mol2.trajectory.ts.positions)


@requires_rdkit
class TestRDKitConverter(object):
    @pytest.fixture
    def pdb(self):
        return mda.Universe(PDB_full)

    @pytest.fixture
    def mol2(self):
        u = mda.Universe(mol2_molecule)
        # add elements
        elements = np.array([guess_atom_element(x) for x in u.atoms.types],
                            dtype=object)
        u.add_TopologyAttr('elements', elements)
        return u

    @pytest.fixture
    def peptide(self):
        u = mda.Universe(GRO)
        elements = mda.topology.guessers.guess_types(u.atoms.names)
        u.add_TopologyAttr('elements', elements)
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
        u = mda.Universe.from_smiles(smi, addHs=False,
                                     generate_coordinates=False)
        rdkit_converter = mda._CONVERTERS["RDKIT"]().convert
        mol = rdkit_converter(u.atoms, NoImplicit=False)
        assert mol.GetNumAtoms() == 1
        assert mol.GetAtomWithIdx(0).GetSymbol() == smi.strip("[]")

    @pytest.mark.parametrize("resname, n_atoms, n_fragments", [
        ("PRO", 14, 1),
        ("ILE", 38, 1),
        ("ALA", 20, 2),
        ("GLY", 21, 3),
    ])
    def test_mol_from_selection(self, peptide, resname, n_atoms, n_fragments):
        mol = peptide.select_atoms("resname %s" % resname).convert_to("RDKIT")
        assert n_atoms == mol.GetNumAtoms()
        assert n_fragments == len(Chem.GetMolFrags(mol))

    @pytest.mark.parametrize("sel_str, atom_index", [
        ("resid 1", 0),
        ("resname LYS and name NZ", 1),
        ("resid 34 and altloc B", 2),
    ])
    def test_monomer_info(self, pdb, sel_str, atom_index):
        rdkit_converter = mda._CONVERTERS["RDKIT"]().convert
        sel = pdb.select_atoms(sel_str)
        umol = rdkit_converter(sel, NoImplicit=False)
        atom = umol.GetAtomWithIdx(atom_index)
        mda_index = atom.GetIntProp("_MDAnalysis_index")
        mi = atom.GetMonomerInfo()

        for mda_attr, rd_attr in RDATTRIBUTES.items():
            rd_value = getattr(mi, "Get%s" % rd_attr)()
            if hasattr(sel, mda_attr):
                mda_value = getattr(sel, mda_attr)[mda_index]
                if mda_attr == "names":
                    rd_value = rd_value.strip()
                assert rd_value == mda_value

    @pytest.mark.parametrize("rdmol", ["mol2_mol", "smiles_mol"],
                             indirect=True)
    def test_identical_topology(self, rdmol):
        u = mda.Universe(rdmol)
        umol = u.atoms.convert_to("RDKIT")
        assert rdmol.HasSubstructMatch(umol) and umol.HasSubstructMatch(rdmol)
        u2 = mda.Universe(umol)
        assert_equal(u.atoms.bonds, u2.atoms.bonds)
        assert_equal(u.atoms.elements, u2.atoms.elements)
        assert_equal(u.atoms.names, u2.atoms.names)
        assert_almost_equal(u.atoms.positions, u2.atoms.positions, decimal=7)

    def test_raise_requires_elements(self):
        u = mda.Universe(mol2_molecule)
        with pytest.raises(
            AttributeError,
            match="`elements` attribute is required for the RDKitConverter"
        ):
            u.atoms.convert_to("RDKIT")

    def test_warn_guess_bonds(self):
        u = mda.Universe(PDB_helix)
        with pytest.warns(UserWarning,
                          match="No `bonds` attribute in this AtomGroup"):
            u.atoms.convert_to("RDKIT")

    def test_bonds_outside_sel(self):
        u = mda.Universe(Chem.MolFromSmiles("CC(=O)C"))
        ag = u.select_atoms("index 1")
        rdkit_converter = mda._CONVERTERS["RDKIT"]().convert
        rdkit_converter(ag, NoImplicit=False)

    def test_error_no_hydrogen(self, uo2):
        with pytest.raises(AttributeError,
                           match="the converter requires all hydrogens to be "
                                 "explicit"):
            uo2.atoms.convert_to("RDKIT")

    def test_error_no_hydrogen_implicit(self, uo2):
        rdkit_converter = mda._CONVERTERS["RDKIT"]().convert
        with pytest.warns(None) as record:
            rdkit_converter(uo2.atoms, NoImplicit=False)
        assert len(record) == 0

    def test_warning_no_hydrogen_force(self, uo2):
        rdkit_converter = mda._CONVERTERS["RDKIT"]().convert
        with pytest.warns(UserWarning,
                          match="Forcing to continue the conversion"):
            rdkit_converter(uo2.atoms, NoImplicit=False, force=True)

    @pytest.mark.parametrize("attr, value, expected", [
        ("names", "C1", " C1 "),
        ("names", "C12", " C12"),
        ("names", "Cl1", "Cl1 "),
        ("altLocs", "A", "A"),
        ("chainIDs", "B", "B"),
        ("icodes", "C", "C"),
        ("occupancies", 0.5, 0.5),
        ("resnames", "LIG", "LIG"),
        ("resids", 123, 123),
        ("segindices", 1, 1),
        ("tempfactors", 0.8, 0.8),
        ("bfactors", 0.8, 0.8),
    ])
    def test_add_mda_attr_to_rdkit(self, attr, value, expected):
        mi = Chem.AtomPDBResidueInfo()
        _add_mda_attr_to_rdkit(attr, value, mi)
        rdvalue = getattr(mi, "Get%s" % RDATTRIBUTES[attr])()
        assert rdvalue == expected

    def test_bfactors_tempfactors_raises_error(self):
        u = mda.Universe.from_smiles("C")
        bfactors = np.array(u.atoms.n_atoms*[1.0], dtype=np.float32)
        u.add_TopologyAttr('bfactors', bfactors)
        u.add_TopologyAttr('tempfactors', bfactors)
        with pytest.raises(
            AttributeError,
            match="Both `tempfactors` and `bfactors` attributes are present"
        ):
            u.atoms.convert_to("RDKIT")

    @pytest.mark.parametrize("idx", [0, 10, 42])
    def test_other_attributes(self, mol2, idx):
        mol = mol2.atoms.convert_to("RDKIT")
        rdatom = mol.GetAtomWithIdx(idx)
        rdprops = rdatom.GetPropsAsDict()
        mda_idx = int(rdprops["_MDAnalysis_index"])
        for prop in ["charge", "segid", "type"]:
            rdprop = rdprops["_MDAnalysis_%s" % prop]
            mdaprop = getattr(mol2.atoms[mda_idx], prop)
            assert rdprop == mdaprop

    @pytest.mark.parametrize("sel_str", [
        "resname ALA",
        "resname PRO and segid A",
    ])
    def test_index_property(self, pdb, sel_str):
        rdkit_converter = mda._CONVERTERS["RDKIT"]().convert
        ag = pdb.select_atoms(sel_str)
        mol = rdkit_converter(ag, NoImplicit=False)
        expected = [i for i in range(len(ag))]
        indices = sorted([a.GetIntProp("_MDAnalysis_index")
                          for a in mol.GetAtoms()])
        assert_equal(indices, expected)

    def test_assign_coordinates(self, pdb):
        rdkit_converter = mda._CONVERTERS["RDKIT"]().convert
        mol = rdkit_converter(pdb.atoms, NoImplicit=False)
        positions = mol.GetConformer().GetPositions()
        indices = sorted(mol.GetAtoms(),
                         key=lambda a: a.GetIntProp("_MDAnalysis_index"))
        indices = [a.GetIdx() for a in indices]
        assert_almost_equal(positions[indices], pdb.atoms.positions)

    def test_assign_stereochemistry(self, mol2):
        umol = mol2.atoms.convert_to("RDKIT")
        rdmol = Chem.MolFromMol2File(mol2_molecule, removeHs=False)
        assert rdmol.HasSubstructMatch(
            umol, useChirality=True) and umol.HasSubstructMatch(
            rdmol, useChirality=True)

    def test_trajectory_coords(self):
        u = mda.Universe.from_smiles(
            "CCO", numConfs=3, rdkit_kwargs=dict(randomSeed=42))
        for ts in u.trajectory:
            mol = u.atoms.convert_to("RDKIT")
            positions = mol.GetConformer().GetPositions()
            indices = sorted(mol.GetAtoms(),
                             key=lambda a: a.GetIntProp("_MDAnalysis_index"))
            indices = [a.GetIdx() for a in indices]
            assert_almost_equal(positions[indices], ts.positions)

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
        cached_func = mda.coordinates.RDKit.atomgroup_to_mol
        # create universes
        utraj = mda.Universe.from_smiles("CCO", numConfs=5)
        uc = mda.Universe.from_smiles("C")
        ucc = mda.Universe.from_smiles("CC")
        uccc = mda.Universe.from_smiles("CCC")
        # test (1): iterating over frames in a trajectory
        previous_cache = None
        for ts in utraj.trajectory:
            mol = utraj.atoms.convert_to("RDKIT")
            cache = cached_func.cache_info()
            if previous_cache:
                # the cache shouldn't change when iterating on timesteps
                assert cache.currsize == previous_cache.currsize
                previous_cache = copy.deepcopy(cache)
        # test (2): only 2 molecules should be cached by default
        cached_func.cache_clear()
        mol = uc.atoms.convert_to("RDKIT")
        mol = ucc.atoms.convert_to("RDKIT")
        mol = uccc.atoms.convert_to("RDKIT")
        cache = cached_func.cache_info()
        assert cache.currsize == 2
        assert cache.misses == 3
        mol = ucc.atoms.convert_to("RDKIT")  # should be inside of the cache
        assert cached_func.cache_info().hits == 1
        mol = uc.atoms.convert_to("RDKIT")  # outside of the cache
        assert cached_func.cache_info().hits == 1
        assert cached_func.cache_info().misses == 4
        # test (3): increase cache size
        mda.coordinates.RDKit.set_converter_cache_size(3)
        cached_func = mda.coordinates.RDKit.atomgroup_to_mol
        assert cached_func.cache_info().maxsize == 3
        # test (4): caching is sensitive to converter arguments
        previous_cache = cached_func.cache_info()
        rdkit_converter = mda._CONVERTERS["RDKIT"]().convert
        mol = rdkit_converter(uc.atoms, NoImplicit=False)
        cache = cached_func.cache_info()
        assert cache.misses == previous_cache.misses + 1
        # test (5): skip cache
        mol = rdkit_converter(uc.atoms, cache=False)
        new_cache = cached_func.cache_info()
        assert cache == new_cache


@requires_rdkit
class TestRDKitFunctions(object):
    @pytest.mark.parametrize("smi, out", [
        ("C(-[H])(-[H])(-[H])-[H]", "C"),
        ("[C](-[H])(-[H])-[C](-[H])-[H]", "C=C"),
        ("[C]1(-[H])-[C](-[H])-[C](-[H])-[C](-[H])-[C](-[H])-[C]1(-[H])",
         "c1ccccc1"),
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
    ])
    def test_infer_bond_orders(self, smi, out):
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        _infer_bo_and_charges(mol)
        Chem.SanitizeMol(mol)
        mol = Chem.RemoveHs(mol)
        molref = Chem.MolFromSmiles(out)
        assert mol.HasSubstructMatch(molref) and molref.HasSubstructMatch(
            mol), "{} != {}".format(Chem.MolToSmiles(mol), out)

    @pytest.mark.parametrize("smi, atom_idx, charge", [
        ("[C](-[H])(-[H])(-[H])-[O]", 4, -1),
        ("[N]-[C]-[O]", 2, -1),
        ("[N](-[H])(-[H])(-[H])-[H]", 0, 1),
        ("C-[C](-[O])-[O]", 3, -1),
    ])
    def test_infer_charges(self, smi, atom_idx, charge):
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        _infer_bo_and_charges(mol)
        Chem.SanitizeMol(mol)
        assert mol.GetAtomWithIdx(atom_idx).GetFormalCharge() == charge

    @pytest.mark.parametrize("smi, out", [
        ("[S](-[O]-[H])(-[O]-[H])(-[O])-[O]", "S(=O)(=O)(O)O"),
        ("[S](-[O]-[H])(-[O])(-[O])-[O]", "S(=O)(=O)([O-])O"),
        ("[S](-[O])(-[O])(-[O])-[O]", "S(=O)(=O)([O-])[O-]"),
        ("C-[N](-[H])-[C](-[N](-[H])-[H])-[N](-[H])-[H]",
         "CNC(N)=[N+](-[H])-[H]"),
        ("[O]-[C](-[H])-[C](-[H])-[H]", "C([O-])=C"),
        ("C-[N](-[O])-[O]", "C[N+](=O)[O-]"),
        ("C(-[N](-[O])-[O])-[N](-[O])-[O]", "C([N+](=O)[O-])[N+](=O)[O-]"),
        ("C-[N](-[O])-[O].C-[N](-[O])-[O]", "C[N+](=O)[O-].C[N+](=O)[O-]"),
        ("[C-](=O)-C", "[C](=O)-C"),
        ("[H]-[N-]-C", "[H]-[N]-C"),
        ("[O]-[C]1-[C](-[H])-[C](-[H])-[C](-[H])-[C](-[H])-[C](-[H])1", 
         "[O-]c1ccccc1"),
        ("[O]-[C]1-[C](-[H])-[C](-[H])-[C](-[H])-[C]1-[O]", 
         "[O-]C1=CC=CC1=O"),
        ("[H]-[C]-[C]-[C](-[H])-[C](-[H])-[H]", "C#CC=C"),
        ("[H]-[C]-[C]-[C]-[C]-[H]", "C#CC#C"),
    ])
    def test_standardize_patterns(self, smi, out):
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        _infer_bo_and_charges(mol)
        mol = _standardize_patterns(mol)
        Chem.SanitizeMol(mol)
        mol = Chem.RemoveHs(mol)
        molref = Chem.MolFromSmiles(out)
        assert mol.HasSubstructMatch(molref) and molref.HasSubstructMatch(
            mol), "{} != {}".format(Chem.MolToSmiles(mol), out)

    @pytest.mark.parametrize("attr, value, getter", [
        ("index", 42, "GetIntProp"),
        ("index", np.int(42), "GetIntProp"),
        ("charge", 4.2, "GetDoubleProp"),
        ("charge", np.float(4.2), "GetDoubleProp"),
        ("type", "C.3", "GetProp"),
    ])
    def test_set_atom_property(self, attr, value, getter):
        atom = Chem.Atom(1)
        prop = "_MDAnalysis_%s" % attr
        _set_atom_property(atom, prop, value)
        assert getattr(atom, getter)(prop) == value

    @pytest.mark.parametrize("rdmol, product, name", [
        ("dummy_reactant", "dummy_product", "props"),
        ("dummy_reactant_noprops", "dummy_product", "noprops"),
        ("dummy_reactant", "dummy_product_nomap", "nomap"),
    ], indirect=["rdmol", "product"])
    def test_reassign_props_after_reaction(self, rdmol, product, name):
        _reassign_props_after_reaction(rdmol, product)
        atom = product.GetAtomWithIdx(0)
        if name == "props":
            assert atom.GetProp("foo") == "bar"
            assert atom.GetIntProp("_MDAnalysis_index") == 1
            assert atom.GetDoubleProp("_MDAnalysis_charge") == 4.2
            assert atom.GetProp("_MDAnalysis_type") == "C.3"
            with pytest.raises(KeyError, match="old_mapno"):
                atom.GetIntProp("old_mapno")
            with pytest.raises(KeyError, match="react_atom_idx"):
                atom.GetUnsignedProp("react_atom_idx")
        elif name == "noprops":
            with pytest.raises(KeyError, match="old_mapno"):
                atom.GetIntProp("old_mapno")
            with pytest.raises(KeyError, match="react_atom_idx"):
                atom.GetUnsignedProp("react_atom_idx")
        elif name == "nomap":
            with pytest.raises(KeyError, match="react_atom_idx"):
                atom.GetUnsignedProp("react_atom_idx")
            with pytest.raises(KeyError, match="_MDAnalysis_index"):
                atom.GetIntProp("_MDAnalysis_index")

    @pytest.mark.parametrize("smi_in", [
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
    ])
    def test_order_independant(self, smi_in):
        # generate mol with hydrogens but without bond orders
        ref = Chem.MolFromSmiles(smi_in)
        template = Chem.AddHs(ref)
        for atom in template.GetAtoms():
            atom.SetIsAromatic(False)
            atom.SetFormalCharge(0)
        for bond in template.GetBonds():
            bond.SetIsAromatic(False)
            bond.SetBondType(Chem.BondType.SINGLE)

        # go through each possible starting atom
        for a in template.GetAtoms():
            smi = Chem.MolToSmiles(template, rootedAtAtom=a.GetIdx())
            m = Chem.MolFromSmiles(smi, sanitize=False)
            for atom in m.GetAtoms():
                atom.SetFormalCharge(0)
                atom.SetNoImplicit(True)
            m.UpdatePropertyCache(strict=False)
            _infer_bo_and_charges(m)
            m = _standardize_patterns(m)
            Chem.SanitizeMol(m)
            m = Chem.RemoveHs(m)
            assert m.HasSubstructMatch(ref) and ref.HasSubstructMatch(
                m), (f"(input) {Chem.MolToSmiles(ref)} != "
                     f"{Chem.MolToSmiles(m)} (output) root atom {a.GetIdx()}")

    def test_warn_conjugated_max_iter(self):
        smi = "[C-]C=CC=CC=CC=CC=CC=C[C-]"
        mol = Chem.MolFromSmiles(smi)
        with pytest.warns(UserWarning,
                          match="reasonable number of iterations"):
            _rebuild_conjugated_bonds(mol, 2)
