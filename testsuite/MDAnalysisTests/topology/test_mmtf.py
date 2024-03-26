import pytest
from numpy.testing import assert_equal, assert_allclose
import mmtf
from unittest import mock

import MDAnalysis as mda
from MDAnalysis.core.groups import AtomGroup

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import MMTF, MMTF_gz, MMTF_skinny, MMTF_skinny2


class MMTFBase(ParserBase):
    expected_attrs = [
        'ids', 'names', 'types', 'altLocs', 'tempfactors', 'occupancies',
        'charges', 'names', 'resnames', 'resids', 'resnums', 'icodes',
        'segids', 'bonds', 'models'
    ]


class TestMMTFParser(MMTFBase):
    parser = mda.topology.MMTFParser.MMTFParser
    ref_filename = MMTF
    guessed_attrs = ['masses']
    expected_n_atoms = 512
    expected_n_residues = 124
    expected_n_segments = 8


class TestMMTFParser_gz(TestMMTFParser):
    ref_filename = MMTF_gz
    expected_n_atoms = 1140
    expected_n_residues = 36
    expected_n_segments = 4


class TestMMTFSkinny(MMTFBase):
    parser = mda.topology.MMTFParser.MMTFParser
    ref_filename = MMTF_skinny
    # for all attributes often in MMTF,
    # check that we get expected error on access
    # (sort so pytest gets reliable order)
    guessed_attrs = ['ids', 'masses', 'segids']
    expected_n_atoms = 660
    expected_n_residues = 134
    expected_n_segments = 2


class TestMMTFSkinny2(MMTFBase):
    parser = mda.topology.MMTFParser.MMTFParser
    ref_filename = MMTF_skinny2
    guessed_attrs = ['ids', 'masses', 'segids']
    expected_n_atoms = 169
    expected_n_residues = 44
    expected_n_segments = 2


class TestMMTFUniverse(object):
    @pytest.fixture()
    def u(self):
        return mda.Universe(MMTF)

    def test_bonds(self, u):
        assert len(u.bonds) == 458

    def test_ids(self, u):
        assert_equal(u.atoms.ids[:3], [1, 2, 3])

    def test_names(self, u):
        assert_equal(u.atoms.names[:3], ["O5'", "C5'", "C4'"])

    def test_resnames(self, u):
        assert_equal(u.residues.resnames[:3], ['DG', 'DA', 'DA'])

    def test_segids(self, u):
        assert_equal(u.segments[:3].segids, ['A', 'B', 'C'])

    def test_resids(self, u):
        assert_equal(u.residues.resids[-3:], [2008, 2009, 2010])

    def test_occupancies(self, u):
        # pylint: disable=unsubscriptable-object
        assert_equal(u.atoms.occupancies[:3], [1.0, 1.0, 1.0])

    def test_bfactors(self, u):
        assert_equal(u.atoms.bfactors[:3], [9.48, 10.88, 10.88])

    def test_types(self, u):
        assert_equal(u.atoms.types[:3], ['O', 'C', 'C'])

    def test_models(self, u):
        assert all(u.atoms.models == 0)

    def test_icodes(self, u):
        assert all(u.atoms.icodes == '')

    def test_altlocs(self, u):
        assert all(u.atoms.altLocs[:3] == '')

    def test_guessed_masses(self, u):
        expected = [15.999, 12.011, 12.011, 15.999, 12.011, 15.999, 12.011]
        assert_allclose(u.atoms.masses[:7], expected)


class TestMMTFUniverseFromDecoder(TestMMTFUniverse):
    @pytest.fixture()
    def u(self):
        return mda.Universe(mmtf.parse(MMTF))


class TestMMTFgzUniverse(object):
    @pytest.fixture()
    def u(self):
        return mda.Universe(MMTF_gz)

    def test_models(self, u):
        # has 2 models
        assert_equal(u.segments[:2].models, [0, 0])
        assert_equal(u.segments[2:].models, [1, 1])

    def test_universe_models(self, u):
        assert len(u.models) == 2
        for m in u.models:
            assert isinstance(m, AtomGroup)
            assert len(m) == 570

    def test_guessed_masses(self, u):
        expected = [15.999, 12.011, 12.011, 15.999, 12.011, 15.999, 12.011]
        assert_allclose(u.atoms.masses[:7], expected)


class TestMMTFgzUniverseFromDecoder(TestMMTFgzUniverse):
    @pytest.fixture()
    def u(self):
        top = mmtf.parse_gzip(MMTF_gz)
        return mda.Universe(top)


class TestMMTFFetch(TestMMTFUniverse):
    @pytest.fixture()
    def u(self):
        top = mmtf.parse(MMTF)
        with mock.patch('mmtf.fetch') as mock_fetch:
            mock_fetch.return_value = top

            return mda.fetch_mmtf('173D')  # string is irrelevant


class TestSelectModels(object):
    # tests for 'model' keyword in select_atoms
    @pytest.fixture()
    def u(self):
        return mda.Universe(MMTF_gz)

    def test_model_selection(self, u):
        m1 = u.select_atoms('model 0')
        m2 = u.select_atoms('model 1')

        assert len(m1) == 570
        assert len(m2) == 570

    def test_model_multiple(self, u):
        m2plus = u.select_atoms('model 1-10')

        assert len(m2plus) == 570

    def test_model_multiple_2(self, u):
        m2plus = u.select_atoms('model 1:10')

        assert len(m2plus) == 570

    def test_model_multiple_3(self, u):
        m1and2 = u.select_atoms('model 0-1')

        assert len(m1and2) == 1140

    def test_model_multiple_4(self, u):
        m1and2 = u.select_atoms('model 0:1')

        assert len(m1and2) == 1140

    def test_model_multiple_5(self, u):
        m1and2 = u.select_atoms('model 0 1')

        assert len(m1and2) == 1140
