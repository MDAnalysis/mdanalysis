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
    ag_size_no_monomer_off = AggregationSize(atomgroup, number_of_mol=30,
                                             aggregation_criteria="Closest",
                                             no_monomer=False).run()
    assert_almost_equal(ag_size_no_monomer_off.results[:, 1],
                        [1.4286, 2.7273, 2.7273, 3.75, 3.3333, 4.2857],
                        decimal=4)


def test_agregation_closest(atomgroup):
    ag_size = AggregationSize(atomgroup, number_of_mol=30,
                              aggregation_criteria="Closest",
                              no_monomer=True).run()
    assert_almost_equal(ag_size.results[:, 1],
                        [2.2857, 3.7143, 4.1667, 5.4, 4, 4.2857], decimal=4)


def test_agregation_com(atomgroup):
    ag_size = AggregationSize(atomgroup, number_of_mol=30,
                              aggregation_criteria="COM",
                              no_monomer=True).run()
    assert_almost_equal(ag_size.results[:, 2], [2, 3.4211, 3.1, 3, 3, 3],
                        decimal=4)


def test_agregation_specific(atomgroup):
    ag_size = AggregationSize(atomgroup, number_of_mol=30,
                              aggregation_criteria="Specific",
                              no_monomer=True, atom_1="C1",
                              atom_2="C1").run()
    assert_almost_equal(ag_size.results[:, 3],
                        [2, 3.871, 3.2807, 3.3226, 4, 3.8052], decimal=4)


def test_gyr_on_off(atomgroup):
    ag_size_gyr_on = AggregationSize(atomgroup, number_of_mol=30,
                                     aggregation_criteria="Closest",
                                     no_monomer=True,
                                     ave_gyr_calc=True).run()
    assert_almost_equal(ag_size_gyr_on.results[:, 4],
                        [8.1545, 10.2749, 10.7043, 11.229, 10.4413, 10.1946],
                        decimal=4)
