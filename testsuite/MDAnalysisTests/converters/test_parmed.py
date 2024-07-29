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
import pytest
import MDAnalysis as mda

import numpy as np
from numpy.testing import (assert_allclose, assert_equal)
from numpy.lib import NumpyVersion

from MDAnalysisTests.coordinates.base import _SingleFrameReader
from MDAnalysisTests.coordinates.reference import RefAdKSmall
from MDAnalysis.converters.ParmEd import ParmEdConverter

from MDAnalysisTests.datafiles import (
    GRO,
    PDB,
    PSF,
    PSF_NAMD,
    PSF_cmap,
    PDB_small,
    PRM,
    PRM12,
    PRM_UreyBradley,
)

# TODO: remove this guard when parmed has a release
# that support NumPy 2
if NumpyVersion(np.__version__) < "2.0.0":
    pmd = pytest.importorskip('parmed')
else:
    pmd = pytest.importorskip('parmed_skip_with_numpy2')



class TestParmEdReaderGRO:
    ref_filename = GRO

    universe = mda.Universe(pmd.load_file(GRO))
    ref = mda.Universe(GRO)

    def test_dimensions(self):
        assert_allclose(
            self.universe.trajectory.ts.dimensions, 
            self.ref.trajectory.ts.dimensions,
            rtol=0,
            atol=1e-3,
            err_msg=("ParmEdReader failed to get unitcell dimensions "
                     "from ParmEd"))
    
    def test_coordinates(self):
        up = self.universe.atoms.positions
        rp = self.ref.atoms.positions
        assert_allclose(up, rp, rtol=0, atol=1e-3)
    

class BaseTestParmEdReader(_SingleFrameReader):
    def setUp(self):
        self.universe = mda.Universe(pmd.load_file(self.ref_filename))
        self.ref = mda.Universe(self.ref_filename)

    def test_dimensions(self):
        assert_allclose(
            self.universe.trajectory.ts.dimensions, 
            self.ref.trajectory.ts.dimensions,
            rtol=0,
            atol=1e-3,
            err_msg=("ParmEdReader failed to get unitcell dimensions "
                     "from ParmEd"))
    
    def test_coordinates(self):
        up = self.universe.atoms.positions
        rp = self.ref.atoms.positions
        assert_allclose(up, rp, rtol=0, atol=1e-3)
    

class TestParmEdReaderPDB(BaseTestParmEdReader):

    ref_filename = RefAdKSmall.filename
    
    def test_uses_ParmEdReader(self):
        from MDAnalysis.coordinates.ParmEd import ParmEdReader

        assert isinstance(self.universe.trajectory, ParmEdReader), "failed to choose ParmEdReader"


def _parmed_param_eq(a, b):
    a_idx = [a.atom1.idx, a.atom2.idx]
    b_idx = [b.atom1.idx, b.atom2.idx]

    for i in (3, 4, 5):
        atom = 'atom{}'.format(i)
        if hasattr(a, atom):
            if not hasattr(b, atom):
                return False
            a_idx.append(getattr(a, atom).idx)
            b_idx.append(getattr(b, atom).idx)
    
    atoms = a_idx == b_idx or a_idx == b_idx[::-1]
    return atoms and a.type == b.type


class BaseTestParmEdConverter:

    equal_atom_attrs = ('name', 'altloc')
    almost_equal_atom_attrs = ('mass', 'charge', 'occupancy')
    expected_attrs = ('atoms', 'bonds', 'angles', 'dihedrals', 'impropers',
                     'cmaps', 'urey_bradleys')


    @pytest.fixture(scope='class')
    def ref(self):
        # skip_bonds controls whether to search for bonds if it's not in the file
        return pmd.load_file(self.ref_filename, skip_bonds=True)  

    @pytest.fixture(scope='class')
    def universe(self, ref):
        return mda.Universe(self.ref_filename)

    @pytest.fixture(scope='class')
    def output(self, universe):
        return universe.atoms.convert_to('PARMED')

    @pytest.fixture(scope='class')
    def roundtrip(self, ref):
        u = mda.Universe(ref)
        return u.atoms.convert_to('PARMED')

    def test_equivalent_connectivity_counts(self, universe, output):
        for attr in ('atoms', 'bonds', 'angles', 'dihedrals', 'impropers',
                     'cmaps', 'urey_bradleys'):
            u = getattr(universe, attr, [])
            o = getattr(output, attr)
            assert len(u) == len(o)
    
    @pytest.mark.parametrize('attr', ('bonds', 'angles', 'dihedrals', 'impropers',
                     'cmaps'))
    def test_equivalent_connectivity_values(self, universe, output, attr):
        u = getattr(universe._topology, attr, [])
        vals = u.values if u else []
        o = getattr(output, attr)
        for param in o:
            ix = [param.atom1.idx, param.atom2.idx]
            try:
                ix.append(param.atom3.idx)
            except:
                pass
            try:
                ix.append(param.atom4.idx)
            except:
                pass
            try:
                ix.append(param.atom5.idx)
            except:
                pass
            ix = tuple(ix)
            assert ix in vals or ix[::-1] in vals

    def test_equivalent_ureybradley_values(self, universe, output):
        try:
            vals = universe._topology.ureybradleys.values
        except AttributeError:
            vals = []
        o = output.urey_bradleys
        for param in o:
            ix = (param.atom1.idx, param.atom2.idx)
            assert ix in vals or ix[::-1] in vals

    def test_equivalent_atoms(self, ref, output):
        for r, o in zip(ref.atoms, output.atoms):
            for attr in self.equal_atom_attrs:
                ra = getattr(r, attr)
                oa = getattr(o, attr)
                assert ra == oa, 'atom {} not equal for atoms {} and {}'.format(attr, r, o)
            
            for attr in self.almost_equal_atom_attrs:
                ra = getattr(r, attr)
                oa = getattr(o, attr)
                assert_allclose(ra, oa, rtol=0, atol=1e-2,
                    err_msg=(f'atom {attr} not almost equal for atoms '
                             f'{r} and {o}'))

    @pytest.mark.parametrize('attr', ('bonds', 'angles', 'impropers',
                     'cmaps'))
    def test_equivalent_connectivity_types(self, ref, roundtrip, attr):
        original = getattr(ref, attr)
        for p in getattr(roundtrip, attr):
            for q in original:
                _parmed_param_eq(p, q)
            assert any(_parmed_param_eq(p, q) for q in original)

    def test_equivalent_dihedrals(self, ref, roundtrip):
        original = ref.dihedrals
        for p in roundtrip.dihedrals:
            assert any((_parmed_param_eq(p, q) and
                        p.improper == q.improper and
                        p.ignore_end == q.ignore_end) for q in original)

    def test_missing_attr(self):
        n_atoms = 10
        u = mda.Universe.empty(n_atoms)
        u.add_TopologyAttr("resid", [1])
        u.add_TopologyAttr("segid", ["DUM"])
        u.add_TopologyAttr("mass", [1] * n_atoms)
        with pytest.warns(UserWarning,
                          match="Supplied AtomGroup was missing the following "
                                "attributes"):
            # should miss names and resnames
            u.atoms.convert_to("PARMED")


class BaseTestParmEdConverterSubset(BaseTestParmEdConverter):

    start_i = 0
    end_i = 0
    skip_i = 1

    @pytest.fixture(scope='class')
    def ref(self):
        # skip_bonds controls whether to search for bonds if it's not in the file
        struct = pmd.load_file(self.ref_filename, skip_bonds=True)
        return struct[self.start_i:self.end_i:self.skip_i]

    @pytest.fixture(scope='class')
    def universe(self):
        u = mda.Universe(self.ref_filename)
        return mda.Merge(u.atoms[self.start_i:self.end_i:self.skip_i])


class BaseTestParmEdConverterFromParmed(BaseTestParmEdConverter):

    equal_atom_attrs = ('name', 'number', 'altloc')

    @pytest.fixture(scope='class')
    def universe(self, ref):
        return mda.Universe(ref)

    def test_equivalent_connectivity_counts(self, ref, output):
        for attr in ('atoms', 'bonds', 'angles', 'dihedrals', 'impropers',
                     'cmaps', 'urey_bradleys'):
            r = getattr(ref, attr)
            o = getattr(output, attr)
            assert len(r) == len(o)


class TestParmEdConverterPRM(BaseTestParmEdConverter):
    ref_filename = PRM


class TestParmEdConverterParmedPRM(BaseTestParmEdConverterFromParmed):
    ref_filename = PRM_UreyBradley


# class TestParmEdConverterPDBSubset(BaseTestParmEdConverterSubset):
#     ref_filename = PDB
#     start_i = 101
#     end_i = 8023
#     skip_i = 3

# TODO: Add Subset test for PRMs when mda.Merge accepts Universes without positions


class TestParmEdConverterParmedPSF(BaseTestParmEdConverterFromParmed):
    ref_filename = PSF_cmap


class TestParmEdConverterPSF(BaseTestParmEdConverter):
    ref_filename = PSF_NAMD


class TestParmEdConverterGROSubset(BaseTestParmEdConverterSubset):
    ref_filename = GRO
    start_i = 5
    end_i = 100

# TODO: Add Subset test for PRMs when mda.Merge accepts Universes without positions

class TestParmEdConverterPDB(BaseTestParmEdConverter):
    ref_filename = PDB_small

    # Neither MDAnalysis nor ParmEd read the mass column 
    # of PDBs and are liable to guess wrong
    almost_equal_atom_attrs = ('charge', 'occupancy')

    def test_equivalent_coordinates(self, ref, output):
        assert_allclose(ref.coordinates, output.coordinates, rtol=0, atol=1e-3)


def test_incorrect_object_passed_typeerror():
    err = "No atoms found in obj argument"
    with pytest.raises(TypeError, match=err):
        c = ParmEdConverter()
        c.convert("we still don't support emojis :(")


def test_old_import_warning():
    wmsg = "Please import the ParmEd classes from MDAnalysis.converters"
    with pytest.warns(DeprecationWarning, match=wmsg):
        from MDAnalysis.coordinates.ParmEd import ParmEdConverter
