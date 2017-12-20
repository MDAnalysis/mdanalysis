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

# test deprecated code
# - in particular stubs introduced in 0.11.0 (and which
#   will be removed in 1.0)
from __future__ import absolute_import

import pytest


@pytest.mark.parametrize('package', [
    "MDAnalysis.core.units",
    "MDAnalysis.core.util",
    "MDAnalysis.core.log",
    "MDAnalysis.core.distances",
    "MDAnalysis.core.transformations",
    "MDAnalysis.core.qcprot",
    "MDAnalysis.KDTree",
    "MDAnalysis.analysis.x3dna"

])
def test_import(package):
    try:
        with pytest.warns(DeprecationWarning):
            __import__(package)
    except ImportError:
        pytest.fail('Failed to import {0}'.format(package))


def test_analysis_x3dna():
    try:
        from MDAnalysis.analysis.x3dna import X3DNA
    except ImportError:
        raise AssertionError("MDAnalysis.analysis.x3dna not available")
