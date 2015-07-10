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

import MDAnalysisTests
from MDAnalysisTests.plugins.knownfailure import knownfailure
from MDAnalysisTests.plugins.memleak import MemleakError
from numpy.testing import *

def test_import():
    try:
        import MDAnalysis
    except ImportError:
        raise AssertionError('Failed to import module MDAnalysis. Install MDAnalysis'
                'first to run the tests, e.g. "pip install mdanalysis"')

def test_matching_versions():
    import MDAnalysis
    assert_(MDAnalysis.__version__ == MDAnalysisTests.__version__,
            "MDAnalysis release {0} must be installed to have meaningful tests, not {1}".format(
            MDAnalysisTests.__version__, MDAnalysis.__version__))

## The following allow testing of the memleak tester plugin.
## Keep commented out unless you suspect the plugin
## might be misbehaving."""
#class A():
#    """This is a small leaky class that won't break anything."""
#    def __init__(self):
#        self.self_ref = self
#    def __del__(self):
#        pass
#
#def test_that_memleaks():
#    """Test that memleaks (Issue 323)"""
#    a = A()
#
#class TestML1(TestCase):
#    def test_that_memleaks(self):
#        """Test that memleaks (Issue 323)"""
#        self.a = A()
#
#class TestML2(TestCase):
#    def setUp(self):
#        a = A()
#    def test_that_memleaks(self):
#        """Test that memleaks (Issue 323)"""
#        pass
