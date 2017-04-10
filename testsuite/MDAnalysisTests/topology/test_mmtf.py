from __future__ import absolute_import
from numpy.testing import (
    assert_,
    assert_array_equal,
)
import mmtf
import mock

import MDAnalysis as mda
from MDAnalysis.core.groups import AtomGroup

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import MMTF, MMTF_gz


class TestMMTFParser(ParserBase):
    parser = mda.topology.MMTFParser.MMTFParser
    filename = MMTF
    expected_attrs = ['ids', 'names', 'types', 'altLocs',
                      'bfactors', 'occupancies', 'charges', 'names',
                      'resnames', 'resids', 'resnums', 'icodes',
                      'segids', 'bonds', 'models']
    guessed_attrs = ['masses']
    expected_n_atoms = 512
    expected_n_residues = 124
    expected_n_segments = 8


class TestMMTFParser_gz(TestMMTFParser):
    filename = MMTF_gz
    expected_n_atoms = 1140
    expected_n_residues = 36
    expected_n_segments = 4


class TestMMTFUniverse(object):
    def setUp(self):
        self.u = mda.Universe(MMTF)

    def tearDown(self):
        del self.u

    def test_bonds(self):
        assert_(len(self.u.bonds) == 458)

    def test_ids(self):
        assert_array_equal(self.u.atoms.ids[:3], [1, 2, 3])

    def test_names(self):
        assert_array_equal(self.u.atoms.names[:3], ["O5'", "C5'", "C4'"])

    def test_altlocs(self):
        assert_array_equal(self.u.atoms.altLocs[:3], [' ', ' ', ' '])

    def test_resnames(self):
        assert_array_equal(self.u.residues.resnames[:3], ['DG', 'DA', 'DA'])

    def test_segids(self):
        assert_array_equal(self.u.segments[:3].segids, ['A', 'B', 'C'])

    def test_resids(self):
        assert_array_equal(self.u.residues.resids[-3:], [2008, 2009, 2010])

    def test_occupancies(self):
        # pylint: disable=unsubscriptable-object
        assert_array_equal(self.u.atoms.occupancies[:3], [1.0, 1.0, 1.0])

    def test_bfactors(self):
        assert_array_equal(self.u.atoms.bfactors[:3], [9.48, 10.88, 10.88])

    def test_types(self):
        assert_array_equal(self.u.atoms.types[:3], ['O', 'C', 'C'])

    def test_models(self):
        assert_(all(self.u.atoms.models == 0))

    def test_icodes(self):
        assert_(all(self.u.atoms.icodes == ''))

    def test_altlocs(self):
        assert_(all(self.u.atoms.altLocs[:3] == ''))


class TestMMTFUniverseFromDecoder(TestMMTFUniverse):
    def setUp(self):
        top = mmtf.parse(MMTF)
        self.u = mda.Universe(top)


class TestMMTFgzUniverse(object):
    def setUp(self):
        self.u = mda.Universe(MMTF_gz)

    def tearDown(self):
        del self.u

    def test_models(self):
        # has 2 models
        assert_array_equal(self.u.segments[:2].models, [0, 0])
        assert_array_equal(self.u.segments[2:].models, [1, 1])

    def test_universe_models(self):
        u = self.u
        assert_(len(u.models) == 2)
        for m in u.models:
            assert_(isinstance(m, AtomGroup))
            assert_(len(m) == 570)


class TestMMTFgzUniverseFromDecoder(TestMMTFgzUniverse):
    def setUp(self):
        top = mmtf.parse_gzip(MMTF_gz)
        self.u = mda.Universe(top)


class TestMMTFFetch(TestMMTFUniverse):
    @mock.patch('mmtf.fetch')
    def setUp(self, mock_fetch):
        top = mmtf.parse(MMTF)
        mock_fetch.return_value = top
        self.u = mda.fetch_mmtf('173D')  # string is irrelevant
        

class TestSelectModels(object):
    # tests for 'model' keyword in select_atoms
    def setUp(self):
        self.u = mda.Universe(MMTF_gz)

    def tearDown(self):
        del self.u

    def test_model_selection(self):
        m1 = self.u.select_atoms('model 0')
        m2 = self.u.select_atoms('model 1')

        assert_(len(m1) == 570)
        assert_(len(m2) == 570)

    def test_model_multiple(self):
        m2plus = self.u.select_atoms('model 1-10')

        assert_(len(m2plus) == 570)

    def test_model_multiple_2(self):
        m2plus = self.u.select_atoms('model 1:10')

        assert_(len(m2plus) == 570)

    def test_model_multiple_3(self):
        m1and2 = self.u.select_atoms('model 0-1')

        assert_(len(m1and2) == 1140)

    def test_model_multiple_4(self):
        m1and2 = self.u.select_atoms('model 0:1')

        assert_(len(m1and2) == 1140)

    def test_model_multiple_5(self):
        m1and2 = self.u.select_atoms('model 0 1')

        assert_(len(m1and2) == 1140)
