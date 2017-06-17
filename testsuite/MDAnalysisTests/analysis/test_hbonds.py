# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
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
from MDAnalysis import SelectionError, SelectionWarning

from numpy.testing import (assert_, assert_equal, assert_array_equal,
                           assert_almost_equal, assert_array_almost_equal,
                           assert_raises, dec)
import numpy as np

import itertools
import warnings
from six import StringIO

from MDAnalysisTests import parser_not_found
from MDAnalysisTests.datafiles import PDB_helix, GRO, XTC, waterPSF, waterDCD

# For type guessing:
from MDAnalysis.topology.core import guess_atom_type
from MDAnalysis.core.topologyattrs import Atomtypes

def guess_types(names):
    """GRO doesn't supply types, this returns an Attr"""
    return Atomtypes(np.array([guess_atom_type(name) for name in names], dtype=object))


class TestHydrogenBondAnalysis(object):
    def setUp(self):
        self.universe = u = MDAnalysis.Universe(PDB_helix)
        self.kwargs = {
            'selection1': 'protein',
            'selection2': 'protein',
            'detect_hydrogens': "distance",
            'distance': 3.0,
            'angle': 150.0,
        }
        # ideal helix with 1 proline:
        self.values = {
            'num_bb_hbonds':  u.atoms.n_residues - u.select_atoms('resname PRO').n_residues - 4,
            'donor_resid': np.array([5,  6,  8,  9, 10, 11, 12, 13]),
            'acceptor_resnm': np.array(['ALA', 'ALA', 'ALA', 'ALA', 'ALA', 'PRO', 'ALA', 'ALA'], dtype='U4'),
            }

    def _run(self, **kwargs):
        kw = self.kwargs.copy()
        kw.update(kwargs)
        h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(self.universe, **kw)
        h.run(verbose=False)
        return h

    def test_helix_backbone(self):
        h = self._run()
        assert_equal(len(h.timeseries[0]),
                     self.values['num_bb_hbonds'], "wrong number of backbone hydrogen bonds")
        assert_equal(h.timesteps, [0.0])

    def test_zero_vs_1based(self):
        h = self._run()
        if h.timeseries[0]:
            assert_equal((int(h.timeseries[0][0][0])-int(h.timeseries[0][0][2])),1)
            assert_equal((int(h.timeseries[0][0][1])-int(h.timeseries[0][0][3])),1)

    def test_generate_table(self):
        h = self._run()
        h.generate_table()
        assert_equal(len(h.table),
                     self.values['num_bb_hbonds'], "wrong number of backbone hydrogen bonds in table")
        assert_(isinstance(h.table, np.core.records.recarray))
        assert_array_equal(h.table.donor_resid, self.values['donor_resid'])
        assert_array_equal(h.table.acceptor_resnm, self.values['acceptor_resnm'])

    @staticmethod
    def test_atoms_too_far():
        pdb = '''
ATOM      1  N   LEU     1      32.310  13.778  14.372  1.00  0.00      SYST N 0
ATOM      2  OW  SOL     2       3.024   4.456   4.147  1.00  0.00      SYST H 0'''

        u = MDAnalysis.Universe(StringIO(pdb), format="pdb")
        h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u, 'resname SOL', 'protein')
        h.run(verbose=False)
        assert_equal(h.timeseries, [[]])

    @staticmethod
    def test_acceptor_OC1_OC2():
        gro = '''test
3
    1ALA    OC1    1   0.000   0.000   0.000
    2ALA      N    2   0.300   0.000   0.000
    2ALA     H1    3   0.200   0.000   0.000
7.29748 7.66094 9.82962'''

        u = MDAnalysis.Universe(StringIO(gro), format="gro")
        h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u)
        h.run(verbose=False)
        assert_equal(h.timeseries[0][0][4], 'ALA2:H1')

    @staticmethod
    def test_true_traj():
        u = MDAnalysis.Universe(GRO, XTC)
        u.add_TopologyAttr(guess_types(u.atoms.names))
        h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u,'protein','resname ASP', distance=3.0, angle=120.0)
        h.run()
        assert_equal(len(h.timeseries), 10)

    def test_count_by_time(self):
        h = self._run()
        c = h.count_by_time()
        assert_equal(c.tolist(), [(0.0, self.values['num_bb_hbonds'])])

    def test_count_by_type(self):
        h = self._run()
        c = h.count_by_type()
        assert_equal(c.frequency, self.values['num_bb_hbonds'] * [1.0])

    def test_count_by_type(self):
        h = self._run()
        t = h.timesteps_by_type()
        assert_equal(t.time, self.values['num_bb_hbonds'] * [0.0])

    def tearDown(self):
        del self.universe


class TestHydrogenBondAnalysisHeuristic(TestHydrogenBondAnalysis):
    def setUp(self):
        super(TestHydrogenBondAnalysisHeuristic, self).setUp()
        self.kwargs['detect_hydrogens'] = "heuristic"


class TestHydrogenBondAnalysisHeavy(TestHydrogenBondAnalysis):
    def setUp(self):
        super(TestHydrogenBondAnalysisHeavy, self).setUp()
        self.kwargs['distance_type'] = "heavy"
        self.kwargs["distance"] = 3.5


class TestHydrogenBondAnalysisHeavyFail(TestHydrogenBondAnalysisHeavy):
    def setUp(self):
        super(TestHydrogenBondAnalysisHeavyFail, self).setUp()
        self.kwargs["distance"] = 3.0
        self.values['num_bb_hbonds'] = 0  # no H-bonds with a D-A distance < 3.0 A (they start at 3.05 A)
        self.values['donor_resid'] = np.array([])
        self.values['acceptor_resnm'] = np.array([], dtype="<U3")


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

    def _tearDown(self):
        del self.universe

    def _run(self, **kwargs):
        kw = self.kwargs.copy()
        kw.update(kwargs)
        with warnings.catch_warnings():
            # ignore SelectionWarning
            warnings.simplefilter("ignore")
            h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(self.universe, **kw)
            h.run(verbose=False)
        return h

    def test_check_static_selections(self):
        self._setUp()
        try:
            def run_HBA(s1, s2, s1type):
                """test that HydrogenBondAnalysis() raises SelectionError for missing donors/acceptors"""
                # no donors/acceptors; only raises error if no updates
                return self._run(selection1=s1, selection2=s2,
                                 update_selection1=False, update_selection2=False,
                                 selection1_type=s1type,
                )
            protein = "protein"
            nothing = "resname ALA and not backbone"
            for s1, s2, s1type in itertools.product((protein, nothing),
                                                    (protein, nothing),
                                                    ("donor", "acceptor", "both")):
                if s1 == s2 == protein:
                    def runOK():
                        """test that HydrogenBondAnalysis() works for protein/protein"""
                        try:
                            h = run_HBA(s1, s2, s1type)
                        except:
                            raise AssertionError("HydrogenBondAnalysis protein/protein failed")
                        else:
                            return True
                    yield runOK
                else:
                    yield assert_raises, SelectionError, run_HBA, s1, s2, s1type
        finally:
            self._tearDown()

    def test_run_empty_selections(self):
        self._setUp()
        try:
            def run_HBA(s1, s2, s1type):
                # no donors/acceptors; should not raise error because updates=True
                return self._run(selection1=s1, selection2=s2,
                                 update_selection1=True, update_selection2=True,
                                 selection1_type=s1type,
                )
            protein = "protein"
            nothing = "resname ALA and not backbone"
            for s1, s2, s1type in itertools.product((protein, nothing),
                                                    (protein, nothing),
                                                    ("donor", "acceptor", "both")):
                def run_HBA_dynamic_selections(*args):
                    try:
                        h = run_HBA(*args)
                    except:
                        raise AssertionError("HydrogenBondAnalysis with update=True failed")
                    else:
                        return True
                yield run_HBA_dynamic_selections, s1, s2, s1type
        finally:
            self._tearDown()


class TestHydrogenBondAnalysisTIP3P(object):
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        self.universe = u = MDAnalysis.Universe(waterPSF, waterDCD)
        self.kwargs = {
            'selection1': 'all',
            'selection2': 'all',
            'detect_hydrogens': "distance",
            'distance': 3.0,
            'angle': 120.0,
        }
        self.h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(self.universe, **self.kwargs)
        self.h.run(verbose=False)
        self.h.generate_table()
        self.normalized_timeseries = self._normalize_timeseries()

        # keys are the names in the h.table
        self.reference = {
            'distance': {'mean': 2.0208776, 'std': 0.31740859},
            'angle': {'mean': 155.13521, 'std': 12.98955},
        }

        # reference values for the table only
        self.reference_table = {
            'donor_resnm': ["TIP3"] * len(self.normalized_timeseries),
            'acceptor_resnm': ["TIP3"] * len(self.normalized_timeseries),
        }

        # index into timeseries (ADJUST ONCE donor_idx and acceptor_ndx are removed)
        # with keys being field names in h.table
        self.columns = {
            'time': 0,
            'donor_idx': 1,
            'acceptor_idx': 2,
            'donor_index': 3,
            'acceptor_index': 4,
            'distance': 7,
            'angle': 8,
        }

        # hackish way to allow looping over self.reference and generating tests
        self._functions = {
            'mean': np.mean,
            'std': np.std,
        }

    def _normalize_timeseries(self):
        # timeseries in normalized form: (t, d_indx1, a_indx1, d_index0, a_index0, donor, acceptor, dist, angle)
        #                   array index:  0     1        2        3         4        5      6        7      8
        timeseries = [[t] + item
                      for t, hframe in zip(self.h.timesteps, self.h.timeseries)
                      for item in hframe]
        return timeseries

    def test_timeseries(self):
        h = self.h
        assert_equal(len(h.timeseries), 10)
        assert_equal(len(self.normalized_timeseries), 29)

        for observable in self.reference:
            idx = self.columns[observable]
            for quantity, reference in self.reference[observable].items():
                func = self._functions[quantity]
                assert_almost_equal(
                    func([item[idx] for item in self.normalized_timeseries]), reference,
                    decimal=5,
                    err_msg="{quantity}({observable}) does not match reference".format(**vars()))

    def test_table_atoms(self):
        h = self.h
        table = h.table

        assert_equal(len(h.table), len(self.normalized_timeseries))

        # test that timeseries and table agree on index data and
        # hydrogen bond information at atom level
        for name, idx in self.columns.items():
            assert_array_almost_equal(table.field(name), [data[idx] for data in self.normalized_timeseries],
                                      err_msg="table[{name}] and timeseries[{idx} do not agree".format(**vars()))

        # test at residue level (issue #801
        # https://github.com/MDAnalysis/mdanalysis/issues/801)
        for name, ref in self.reference_table.items():
            assert_array_equal(h.table.field(name), ref,
                               err_msg="resname for {0} do not match (Issue #801)")

