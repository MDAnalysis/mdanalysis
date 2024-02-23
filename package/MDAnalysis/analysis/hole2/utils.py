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


"""
HOLE Analysis --- :mod:`MDAnalysis.analysis.hole2.helper`
=========================================================

:Author: Lily Wang
:Year: 2020
:Copyright: GNU Public License v3

.. versionadded:: 1.0

Helper functions used in :mod:`MDAnalysis.analysis.hole2.hole`
"""
import warnings

from mdahole2.analysis.utils import (
    write_simplerad2,
    check_and_fix_long_filename,
    set_up_hole_input,
    run_hole,
    collect_hole,
    create_vmd_surface
)

wmsg = ("Deprecation in version 2.8.0\n"
        "MDAnalysis.analysis.hole2 is deprecated in favour of the "
        "MDAKit madahole2 (https://www.mdanalysis.org/mdahole2/) "
        "and will be removed in MDAnalysis version 3.0.0")
warnings.warn(wmsg, category=DeprecationWarning)
