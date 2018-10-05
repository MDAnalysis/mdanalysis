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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import print_function, absolute_import


import MDAnalysis
import MDAnalysis.analysis.hbonds
import itertools
import pytest
from MDAnalysis import SelectionError

from numpy.testing import (
    assert_equal, assert_array_equal, assert_almost_equal,
    assert_array_almost_equal, assert_allclose,
)
import numpy as np


from six import StringIO

from MDAnalysisTests.datafiles import PDB_helix, GRO, XTC, waterPSF, waterDCD

# For type guessing:
from MDAnalysis.topology.core import guess_atom_type
from MDAnalysis.core.topologyattrs import Atomtypes


def guess_types(names):
    """GRO doesn't supply types, this returns an Attr"""
    return Atomtypes(np.array([guess_atom_type(name) for name in names], dtype=object))


class TestHydrogenBondAnalysis(object):
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(PDB_helix)

    @staticmethod
    @pytest.fixture(scope='class')
    def values(universe):
        return {
            'num_bb_hbonds': universe.atoms.n_residues - universe.select_atoms('resname PRO').n_residues - 4,
            'donor_resid': np.array([5, 6, 8, 9, 10, 11, 12, 13]),
            'acceptor_resnm': np.array(['ALA', 'ALA', 'ALA', 'ALA', 'ALA', 'PRO', 'ALA', 'ALA'], dtype='U4'),
        }

    kwargs = {
        'selection1': 'protein',
        'selection2': 'protein',
        'detect_hydrogens': "distance",
        'distance': 3.0,
        'angle': 150.0,
    }

    @pytest.fixture(scope='class')
    def h(self, universe):
        kw = self.kwargs.copy()
        # kw.update(kwargs)
        h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(universe, **kw)
        # remove in 1.0
        if kw['detect_hydrogens'] == 'heuristic':
            with pytest.warns(DeprecationWarning):
                h.run(verbose=False)
        else:
            h.run(verbose=False)
        return h

    def test_helix_backbone(self, values, h):

        assert len(h.timeseries[0]) == values['num_bb_hbonds'], "wrong number of backbone hydrogen bonds"
        assert h.timesteps, [0.0]

    def test_generate_table(self, values, h):

        h.generate_table()
        assert len(h.table) == values['num_bb_hbonds'], "wrong number of backbone hydrogen bonds in table"
        assert isinstance(h.table, np.core.records.recarray)
        assert_array_equal(h.table.donor_resid, values['donor_resid'])
        assert_array_equal(h.table.acceptor_resnm, values['acceptor_resnm'])

    def test_atoms_too_far(self):
        pdb = '''
ATOM      1  N   LEU     1      32.310  13.778  14.372  1.00  0.00      SYST N 0
ATOM      2  OW  SOL     2       3.024   4.456   4.147  1.00  0.00      SYST H 0'''

        u = MDAnalysis.Universe(StringIO(pdb), format="pdb")
        h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u, 'resname SOL', 'protein')
        h.run(verbose=False)
        assert h.timeseries == [[]]

    def test_acceptor_OC1_OC2(self):
        gro = '''test
3
    1ALA    OC1    1   0.000   0.000   0.000
    2ALA      N    2   0.300   0.000   0.000
    2ALA     H1    3   0.200   0.000   0.000
7.29748 7.66094 9.82962'''

        u = MDAnalysis.Universe(StringIO(gro), format="gro")
        h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u)
        h.run(verbose=False)
        assert h.timeseries[0][0][2] == 'ALA2:H1'

    def test_true_traj(self):
        u = MDAnalysis.Universe(GRO, XTC)
        u.add_TopologyAttr(guess_types(u.atoms.names))
        h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u,'protein','resname ASP', distance=3.0, angle=120.0)
        h.run()
        assert len(h.timeseries) == 10

    def test_count_by_time(self, values, h):

        c = h.count_by_time()
        assert c.tolist(), [(0.0, values['num_bb_hbonds'])]

    def test_count_by_type(self, values, h):

        c = h.count_by_type()
        assert_equal(c.frequency, values['num_bb_hbonds'] * [1.0])

    # TODO: Rename
    def test_count_by_type2(self, values, h):

        t = h.timesteps_by_type()
        assert_equal(t.time, values['num_bb_hbonds'] * [0.0])

class TestHydrogenBondAnalysisPBC(TestHydrogenBondAnalysis):
    # This system is identical to above class
    # But has a box introduced, and atoms moved into neighbouring images
    # The results however should remain identical if PBC is used
    # If pbc:True in kwargs is changed, these tests should fail
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        u = MDAnalysis.Universe(PDB_helix)
        # transfer to memory to changes to coordinates are reset
        u.transfer_to_memory()
        # place in huge oversized box
        # real system is about 30A wide at most
        boxsize = 150.
        box = np.array([boxsize, boxsize, boxsize, 90., 90., 90.])
        u.dimensions = box

        # then scatter the atoms throughout various images of this box
        u.atoms[::4].translate([boxsize * 2, 0, 0])
        u.atoms[1::4].translate([0, boxsize *4, 0])
        u.atoms[2::4].translate([-boxsize * 5, 0, -boxsize * 2])

        return u

    kwargs = {
        'selection1': 'protein',
        'selection2': 'protein',
        'detect_hydrogens': "distance",
        'distance': 3.0,
        'angle': 150.0,
        'pbc': True,
    }



class TestHydrogenBondAnalysisHeuristic(TestHydrogenBondAnalysis):
        kwargs = {
            'selection1': 'protein',
            'selection2': 'protein',
            'detect_hydrogens': "heuristic",
            'distance': 3.0,
            'angle': 150.0,
        }


class TestHydrogenBondAnalysisHeavy(TestHydrogenBondAnalysis):

    kwargs = {
        'selection1': 'protein',
        'selection2': 'protein',
        'detect_hydrogens': "heuristic",
        'distance': 3.5,
        'angle': 150.0,
        'distance_type': 'heavy',
    }


class TestHydrogenBondAnalysisHeavyFail(TestHydrogenBondAnalysisHeavy):
    kwargs = {
        'selection1': 'protein',
        'selection2': 'protein',
        'detect_hydrogens': "heuristic",
        'distance': 3.0,
        'angle': 150.0,
        'distance_type': 'heavy',
    }

    @staticmethod
    @pytest.fixture()
    def values(universe):
        return {
            'num_bb_hbonds': 0,  # no H-bonds with a D-A distance < 3.0 A (they start at 3.05 A)
            'donor_resid': np.array([]),
            'acceptor_resnm': np.array([], dtype="<U3"),
        }


class TestHydrogenBondAnalysisChecking(object):
    def _setUp(self):
        self.universe = u = MDAnalysis.Universe(PDB_helix)
        self.kwargs = {
            'selection1': 'protein',
            'selection2': 'protein',
            'detect_hydrogens': "distance",
            'distance': 3.0,
            'angle': 150.0,
        }

    # ignore SelectionWarning
    @pytest.mark.filterwarnings('ignore')
    def _run(self, **kwargs):
        kw = kwargs.copy()
        kw.update(kwargs)
        h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(self.universe, **kw)
        h.run(verbose=False)
        return h

    def run_HBA_no_update(self, s1, s2, s1type):
        """test that HydrogenBondAnalysis() raises SelectionError for missing donors/acceptors"""
        # no donors/acceptors; only raises error if no updates
        return self._run(
            selection1=s1,
            selection2=s2,
            update_selection1=False,
            update_selection2=False,
            selection1_type=s1type,
        )

    def runOK(self, s1, s2, s1type):
        """test that HydrogenBondAnalysis() works for protein/protein"""
        try:
            h = self.run_HBA_no_update(s1, s2, s1type)
        except:
            raise pytest.fail("HydrogenBondAnalysis protein/protein failed")
        else:
            return True

    @pytest.mark.parametrize('s1, s2, s1type', itertools.product(
        ("protein", "resname ALA and not backbone"),
        ("protein", "resname ALA and not backbone"),
        ("donor", "acceptor", "both")
    ))
    def test_check_static_selections_pt(self, s1, s2, s1type):
        self._setUp()
        protein = "protein"

        if s1 == s2 == protein:
            self.runOK(s1, s2, s1type)
        else:
            with pytest.raises(SelectionError):
                self.run_HBA_no_update(s1, s2, s1type)

    def run_HBA_update(self, s1, s2, s1type):
        # no donors/acceptors; should not raise error because updates=True
        return self._run(
            selection1=s1,
            selection2=s2,
            update_selection1=True,
            update_selection2=True,
            selection1_type=s1type,
        )

    def run_HBA_dynamic_selections(self, *args):
        try:
            h = self.run_HBA_update(*args)
        except:
            raise pytest.fail("HydrogenBondAnalysis with update=True failed")
        else:
            return True

    @pytest.mark.parametrize('s1, s2, s1type', itertools.product(
        ("protein", "resname ALA and not backbone"),
        ("protein", "resname ALA and not backbone"),
        ("donor", "acceptor", "both")
    ))
    def test_run_empty_selections_pt(self, s1, s2, s1type):
        self._setUp()

        self.run_HBA_dynamic_selections(s1, s2, s1type)


class TestHydrogenBondAnalysisTIP3P(object):
    @staticmethod
    @pytest.fixture()
    def universe():
        return MDAnalysis.Universe(waterPSF, waterDCD)

    kwargs = {
        'selection1': 'all',
        'selection2': 'all',
        'detect_hydrogens': "distance",
        'distance': 3.0,
        'angle': 120.0,
    }

    @pytest.fixture()
    def h(self, universe):
        h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(universe, **self.kwargs)
        h.run(verbose=False)
        h.generate_table()
        return h

    @pytest.fixture()
    def normalized_timeseries(self, h):
        # timeseries in normalized form: (t, d_indx1, a_indx1, d_index0, a_index0, donor, acceptor, dist, angle)
        #                   array index:  0     1        2        3         4        5      6        7      8
        timeseries = [[t] + item
                      for t, hframe in zip(h.timesteps, h.timeseries)
                      for item in hframe]
        return timeseries

    # keys are the names in the h.table
    reference = {
        'distance': {'mean': 2.0208776, 'std': 0.31740859},
        'angle': {'mean': 155.13521, 'std': 12.98955},
    }
    
    @pytest.fixture()
    def reference_table(self, normalized_timeseries):
        # reference values for the table only
        return {
            'donor_resnm': ["TIP3"] * len(normalized_timeseries),
            'acceptor_resnm': ["TIP3"] * len(normalized_timeseries),
        }

    # index into timeseries (ADJUST ONCE donor_idx and acceptor_ndx are removed)
    # with keys being field names in h.table
    columns = {
        'time': 0,
        'donor_index': 1,
        'acceptor_index': 2,
        'distance': 5,
        'angle': 6,
    }

    # hackish way to allow looping over self.reference and generating tests
    _functions = {
        'mean': np.mean,
        'std': np.std,
    }

    def test_timeseries(self, h, normalized_timeseries):
        h = h
        assert len(h.timeseries) == 10
        assert len(normalized_timeseries) == 29

        for observable in self.reference:
            idx = self.columns[observable]
            for quantity, reference in self.reference[observable].items():
                func = self._functions[quantity]
                assert_allclose(
                    func([item[idx] for item in normalized_timeseries]), reference,
                    rtol=1e-5, atol=0,
                    err_msg="{quantity}({observable}) does not match reference".format(**vars())
                )

    def test_table_atoms(self, h, normalized_timeseries, reference_table):
        h = h
        table = h.table

        assert len(h.table) == len(normalized_timeseries)

        # test that timeseries and table agree on index data and
        # hydrogen bond information at atom level
        for name, idx in self.columns.items():
            assert_array_almost_equal(table.field(name), [data[idx] for data in normalized_timeseries],
                                      err_msg="table[{name}] and timeseries[{idx} do not agree".format(**vars()))

        # test at residue level (issue #801
        # https://github.com/MDAnalysis/mdanalysis/issues/801)
        for name, ref in reference_table.items():
            assert_array_equal(h.table.field(name), ref, err_msg="resname for {0} do not match (Issue #801)")
