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
from __future__ import absolute_import

import re

import MDAnalysisTests

def test_import():
    try:
        import MDAnalysis
    except ImportError:
        raise AssertionError('Failed to import module MDAnalysis. Install MDAnalysis'
                'first to run the tests, e.g. "pip install mdanalysis"')


def test_matching_versions():
    import MDAnalysis
    assert MDAnalysis.__version__ == MDAnalysisTests.__version__, \
        "MDAnalysis release {0} must be installed to have meaningful tests, not {1}".format(
            MDAnalysisTests.__version__, MDAnalysis.__version__)


def is_pep440_compliant(version):
    """Check that version is PEP440 compliant

    Including local component (after '+', which is used by versioneer
    for builds after a release tag).

    Source: https://www.python.org/dev/peps/pep-0440/#appendix-b-parsing-version-strings-with-regular-expressions

    """
    VERSION_PATTERN = r"""
        v?
        (?:
            (?:(?P<epoch>[0-9]+)!)?                           # epoch
            (?P<release>[0-9]+(?:\.[0-9]+)*)                  # release segment
            (?P<pre>                                          # pre-release
                [-_\.]?
                (?P<pre_l>(a|b|c|rc|alpha|beta|pre|preview))
                [-_\.]?
                (?P<pre_n>[0-9]+)?
            )?
            (?P<post>                                         # post release
                (?:-(?P<post_n1>[0-9]+))
                |
                (?:
                    [-_\.]?
                    (?P<post_l>post|rev|r)
                    [-_\.]?
                    (?P<post_n2>[0-9]+)?
                )
            )?
            (?P<dev>                                          # dev release
                [-_\.]?
                (?P<dev_l>dev)
                [-_\.]?
                (?P<dev_n>[0-9]+)?
            )?
        )
        (?:\+(?P<local>[a-z0-9]+(?:[-_\.][a-z0-9]+)*))?       # local version
    """
    _regex = re.compile(r"^\s*" + VERSION_PATTERN + r"\s*$",
                        re.VERBOSE | re.IGNORECASE,)
    return _regex.match(version) is not None


def test_version_format(version=None):
    if version is None:
        import MDAnalysis._version
        version = MDAnalysis._version.get_versions()['version']
    # see http://wiki.mdanalysis.org/SemanticVersioning for format definition
    # and PEP440
    assert is_pep440_compliant(version), \
        "version {0} does not match the MAJOR.MINOR.PATCH(+suffix) PEP440 format".format(version)

def test_version_at_packagelevel():
    import MDAnalysis
    try:
        version = MDAnalysis.__version__
    except:
        raise AssertionError("MDAnalysis.__version__ missing")
    return test_version_format(version)

def test_versioneer_version():
    import MDAnalysis._version
    versions = MDAnalysis._version.get_versions()
    assert isinstance(versions, dict)
    assert versions['version']
    assert versions['full-revisionid']
    assert versions['error'] is None

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
