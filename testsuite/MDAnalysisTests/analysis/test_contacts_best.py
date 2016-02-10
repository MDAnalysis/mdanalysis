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

from numpy.testing import (
    assert_almost_equal,
)

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.contacts import BestHummerContacts, calculate_contacts

from MDAnalysisTests.datafiles import (
    contacts_villin_folded, 
    contacts_villin_unfolded,
    contacts_file,
)


class TestInterRDF(object):
    def setUp(self):
        self.folded = mda.Universe(contacts_villin_folded)
        self.unfolded = mda.Universe(contacts_villin_unfolded)

    def tearDown(self):
        del self.folded, self.unfolded

    def test_math_folded(self):

        u = self.folded

        # read the text files
        data = [l.split() for l in open(contacts_file).readlines()]
        # convert to 0-based indexing
        data = [ (int(i)-1, int(j)-1, float(d))  for i, j, d in data]
        # get r and r0
        data = [ (np.linalg.norm(u.atoms[i].pos - u.atoms[j].pos), d)  for i, j, d in data]
        data = np.array(data)

        r = data[:,0]
        r0 = data[:,1]

        beta = 5.0
        lambda_constant = 1.8

        Q = 1/(1 + np.exp(beta*(r - lambda_constant * r0)))

        assert_almost_equal(Q.mean(),  1.0, decimal=3)

    def test_math_unfolded(self):

        u = self.unfolded

        # read the text files
        data = [l.split() for l in open(contacts_file).readlines()]
        # convert to 0-based indexing
        data = [ (int(i)-1, int(j)-1, float(d))  for i, j, d in data]
        # get r and r0
        data = [ (np.linalg.norm(u.atoms[i].pos - u.atoms[j].pos), d)  for i, j, d in data]
        data = np.array(data)

        r = data[:,0]
        r0 = data[:,1]

        beta = 5.0
        lambda_constant = 1.8

        Q = 1/(1 + np.exp(beta*(r - lambda_constant * r0)))

        assert_almost_equal(Q.mean(),  0.0, decimal=1)

    def test_math_folded(self):

        # one folded, one unfolded
        f = mda.Universe(contacts_villin_folded)
        u = mda.Universe(contacts_villin_unfolded)
        sel = "protein and not name H*"

        grF = f.select_atoms(sel)
        grU = u.select_atoms(sel)

        q = BestHummerContacts(grU, grU, grF, grF)
        q.run()
        
        results = zip(*calculate_contacts(f, u, sel, sel))[1]
        assert_almost_equal(q.results, results)

    def test_math_folded(self):

        # both folded
        f = mda.Universe(contacts_villin_folded)
        u = mda.Universe(contacts_villin_folded)
        sel = "protein and not name H*"

        grF = f.select_atoms(sel)
        grU = u.select_atoms(sel)

        q = BestHummerContacts(grU, grU, grF, grF)
        q.run()
        
        results = zip(*calculate_contacts(f, u, sel, sel))[1]
        assert_almost_equal(q.results, results)        