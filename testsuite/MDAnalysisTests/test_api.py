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
"""
Test the user facing API is as we expect...
"""

import MDAnalysis as mda
import pytest


def test_Universe():
    assert mda.Universe is mda.core.universe.Universe


def test_fetch_mmtf():
    assert mda.fetch_mmtf is mda.coordinates.MMTF.fetch_mmtf


def test_Writer():
    assert mda.Writer is mda.coordinates.core.writer


def test_AtomGroup():
    assert mda.AtomGroup is mda.core.groups.AtomGroup


def test_ResidueGroup():
    assert mda.ResidueGroup is mda.core.groups.ResidueGroup


def test_SegmentGroup():
    assert mda.SegmentGroup is mda.core.groups.SegmentGroup


@pytest.mark.parametrize('module',
                         ('',
                          '.analysis',
                          '.analysis.data.filenames',
                          '.analysis.distances',
                          '.analysis.encore',
                          '.analysis.encore.dimensionality_reduction',
                          '.analysis.hbonds',
                          '.analysis.hydrogenbonds',
                          '.coordinates',
                          '.core',
                          '.lib._augment',
                          '.lib._cutil',
                          '.lib.formats',
                          '.lib.pkdtree',
                          '.topology',
                          '.topology.tpr',
                          '.visualization'
                          )
                         )
def test_import_star(module):
    try:
        # This is really bad and ugly but I haven't found a way to run this
        # using `importlib.import_module...`
        exec("from MDAnalysis{} import *".format(module))
    except (ImportError, AttributeError) as error:
        # We want the test to fail, not to error
        raise AssertionError(error)
