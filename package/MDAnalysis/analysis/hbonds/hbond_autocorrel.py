# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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
"""Hydrogen bond autocorrelation --- :mod:`MDAnalysis.analysis.hbonds.hbond_autocorrel` (deprecated)
====================================================================================================

:Author: Richard J. Gowers
:Year: 2014
:Copyright: Lesser GNU Public License v2.1+

.. versionadded:: 0.9.0

.. deprecated:: 2.0.0

   This module was moved to
   :mod:`MDAnalysis.analysis.hydrogenbonds.hbond_autocorrel` and access to this
   module through :mod:`MDAnalysis.analysis.hbonds.hbond_autocorrel` will be
   removed in 3.0.0.


See Also
--------
:mod:`MDAnalysis.analysis.hydrogenbonds.hbond_autocorrel`

"""
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("always", DeprecationWarning)
    wmsg = (
        "This module was moved to "
        "MDAnalysis.analysis.hydrogenbonds.hbond_autocorrel; "
        "hbonds.hbond_autocorrel will be removed in 3.0.0."
    )
    warnings.warn(wmsg, category=DeprecationWarning)

from MDAnalysis.lib.util import deprecate

from ..hydrogenbonds import hbond_autocorrel


find_hydrogen_donors = deprecate(
    hbond_autocorrel.find_hydrogen_donors,
    release="2.0.0",
    remove="3.0.0",
    message="The function was moved to "
    "MDAnalysis.analysis.hbonds.hbond_autocorrel.",
)

HydrogenBondAutoCorrel = deprecate(
    hbond_autocorrel.HydrogenBondAutoCorrel,
    release="2.0.0",
    remove="3.0.0",
    message="The class was moved to "
    "MDAnalysis.analysis.hbonds.hbond_autocorrel.",
)
