# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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

# Water Bridge Analysis
r"""Water Bridge analysis --- :mod:`MDAnalysis.analysis.hbonds.WaterBridgeAnalysis`
===============================================================================

:Author: Zhiyi Wu
:Year: 2017-2018
:Copyright: GNU Public License v3
:Maintainer: Zhiyi Wu <zhiyi.wu@gtc.ox.ac.uk>,  `@xiki-tempula`_ on GitHub


.. _`@xiki-tempula`: https://github.com/xiki-tempula

This module is moved to
:class:`MDAnalysis.analysis.hydrogenbonds.WaterBridgeAnalysis`."""

import warnings

from ..hydrogenbonds.wbridge_analysis import (
    WaterBridgeAnalysis,)

warnings.warn(
            "This module has been moved to MDAnalysis.analysis.hydrogenbonds"
            "It will be removed in MDAnalysis version 2.0"
            "Please use MDAnalysis.analysis.hydrogenbonds.wbridge_analysis instead.",
            category=DeprecationWarning
        )
