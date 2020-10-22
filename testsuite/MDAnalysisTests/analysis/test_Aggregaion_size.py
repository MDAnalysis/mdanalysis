from __future__ import absolute_import
from six.moves import zip

import pytest
from numpy.testing import assert_equal, assert_almost_equal
import numpy as np

import MDAnalysis as mda

from MDAnalysis.analysis.Aggregation_size import AggregationSize

from MDAnalysisTests.datafiles import Asph_gro, Asph_xtc


@pytest.fixture()
def atomgroup():
    return mda.Universe(Asph_gro, Asph_xtc).select_atoms('resname AsphC')


def test_monomer_on_off(atomgroup):
    ag_size_no_monomer_off = AggregationSize(atomgroup, number_of_mol=50,
                                             aggregation_criteria="Closest",
                                             no_monomer=False).run()
    assert_almost_equal(ag_size_no_monomer_off.results[:, 1],
                        [1.92308, 3.125, 7.1429, 8.3333, 6.25, 7.1429],
                        decimal=4)


def test_agregation_closest(atomgroup):
    ag_size = AggregationSize(atomgroup, number_of_mol=50,
                              aggregation_criteria="Closest",
                              no_monomer=True).run()
    assert_almost_equal(ag_size.results[:, 1],
                        [3, 3.2667, 7.1429, 9.8, 7, 7.1429], decimal=4)


def test_agregation_com(atomgroup):
    ag_size = AggregationSize(atomgroup, number_of_mol=50,
                              aggregation_criteria="COM",
                              no_monomer=True).run()
    assert_almost_equal(ag_size.results[:, 2], [2, 2.95, 3.5909, 3.1304, 3.7917, 4.6522],
                        decimal=4)


def test_agregation_specific(atomgroup):
    ag_size = AggregationSize(atomgroup, number_of_mol=50,
                              aggregation_criteria="Specific",
                              no_monomer=True, atom_1="C1",
                              atom_2="C1").run()
    assert_almost_equal(ag_size.results[:, 3],
                        [2, 3.6972, 5.1977, 3.0840, 5.7573, 8.5312], decimal=4)


def test_gyr_on_off(atomgroup):
    ag_size_gyr_on = AggregationSize(atomgroup, number_of_mol=50,
                                     aggregation_criteria="Closest",
                                     no_monomer=True,
                                     ave_gyr_calc=True).run()
    assert_almost_equal(ag_size_gyr_on.results[:, 4],
                        [9.7073, 8.7964, 11.8574, 12.1666, 11.5719, 10.9102],
                        decimal=4)
