# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2020 The MDAnalysis Development Team and contributors
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

from __future__ import absolute_import
from ...due import due, Doi
from .hole import hole, HoleAnalysis
from .utils import create_vmd_surface

due.cite(Doi("10.1016/S0006-3495(93)81293-1"),
         description="HOLE program",
         path="MDAnalysis.analysis.hole",
         cite_module=True)
due.cite(Doi("10.1016/S0263-7855(97)00009-X"),
         description="HOLE program",
         path="MDAnalysis.analysis.hole",
         cite_module=True)
due.cite(Doi("10.1016/j.jmb.2013.10.024"),
         description="HOLE trajectory analysis with orderparameters",
         path="MDAnalysis.analysis.hole",
         cite_module=True)
del Doi
