# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
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
    assert_,
    assert_raises,
)
from MDAnalysisTests.core.groupbase import make_Universe

import MDAnalysis as mda
from MDAnalysis.topology.base import TopologyReader


class IOErrorParser(TopologyReader):
    def parse(self):
        raise IOError("Useful information")



class TestUniverseCreation(object):
    # tests concerning Universe creation and errors encountered
    @staticmethod
    def test_make_universe_no_args():
        # universe creation without args should work
        u = mda.Universe()

        assert_(isinstance(u, mda.Universe))
        assert_(u.atoms == None)

    @staticmethod
    def test_Universe_no_trajectory_AE():
        # querying trajectory without a trajectory loaded (only topology)
        u = make_Universe()

        assert_raises(AttributeError, getattr, u, 'trajectory')

    @staticmethod
    def test_Universe_topology_IE():
        assert_raises(IOError,
                      mda.Universe, 'thisfile', topology_format=IOErrorParser)

    @staticmethod
    def test_Universe_topology_IE_msg():
        # should get the original error, as well as Universe error
        try:
            mda.Universe('thisfile', topology_format=IOErrorParser)
        except IOError as e:
            assert_('Failed to load from the topology file' in e.args[0])
            assert_('Useful information' in e.args[0])
        else:
            raise AssertionError

    @staticmethod
    def test_load_new_VE():
        u = mda.Universe()
        
        assert_raises(TypeError,
                      u.load_new, 'thisfile', format='soup')
