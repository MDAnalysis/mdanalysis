# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import print_function

import MDAnalysis
import MDAnalysis.analysis.hbonds
from MDAnalysis import SelectionError, SelectionWarning

from numpy.testing import (TestCase, assert_, assert_equal, assert_array_equal,
                           assert_raises)
import numpy as np

import itertools
import warnings

from MDAnalysisTests.datafiles import PDB_helix, GRO, XTC


class TestHydrogenBondAnalysis(TestCase):
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
            'num_bb_hbonds':  u.atoms.n_residues - u.SYSTEM.PRO.n_residues - 4,
            'donor_resid': np.array([5,  6,  8,  9, 10, 11, 12, 13]),
            'acceptor_resnm': np.array(['ALA', 'ALA', 'ALA', 'ALA', 'ALA', 'PRO', 'ALA', 'ALA'], dtype='U4'),
            }

    def _run(self, **kwargs):
        kw = self.kwargs.copy()
        kw.update(kwargs)
        h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(self.universe, **kw)
        h.run(quiet=True)
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
    def test_true_traj():
        u = MDAnalysis.Universe(GRO, XTC)
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
            h.run(quiet=True)
        return h

    def test_check_static_selections(self):
        self._setUp()  # manually set up (because with yield cannot use TestCase)
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
            self._tearDown() # manually tear down (because with yield cannot use TestCase)


    def test_run_empty_selections(self):
        self._setUp()  # manually set up (because with yield cannot use TestCase)
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
            self._tearDown() # manually tear down (because with yield cannot use TestCase)
