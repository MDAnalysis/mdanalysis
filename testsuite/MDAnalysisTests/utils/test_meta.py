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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
import re

import MDAnalysisTests

def test_import():
    try:
        import MDAnalysis
    except ImportError:
        raise AssertionError('Failed to import module MDAnalysis. Install MDAnalysis'
                'first to run the tests, e.g. "pip install mdanalysis"')


def test_matching_versions():
    import MDAnalysis.version
    assert MDAnalysis.version.__version__ == MDAnalysisTests.__version__, \
        "MDAnalysis release {0} must be installed to have meaningful tests, not {1}".format(
            MDAnalysisTests.__version__, MDAnalysis.__version__)


def test_version_format(version=None):
    if version is None:
        import MDAnalysis.version
        version = MDAnalysis.version.__version__
    # see https://github.com/MDAnalysis/mdanalysis/wiki/SemanticVersioning for format definition
    m = re.match(r'(?P<MAJOR>\d+)\.(?P<MINOR>\d+)\.(?P<PATCH>\d+)(-(?P<suffix>\w+))?$',
                 version)
    assert m, "version {0} does not match the MAJOR.MINOR.PATCH(-suffix) format".format(version)


def test_version_at_packagelevel():
    import MDAnalysis
    try:
        version = MDAnalysis.__version__
    except:
        raise AssertionError("MDAnalysis.__version__ missing")
    return test_version_format(version)


# The following allow testing of the memleak tester plugin.
# Keep commented out unless you suspect the plugin
# might be misbehaving. Apparently python3 is immune to these leaks!"""
#from numpy.testing import TestCase
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
