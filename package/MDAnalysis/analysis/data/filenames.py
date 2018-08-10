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
"""
Analysis data files
===================

:mod:`MDAnalysis.analysis.data` contains data files that are used as part of
analysis. These can be experimental or theoretical data. Files are stored
inside the package and made accessible via variables in
:mod:`MDAnalysis.analysis.data.filenames`. These variables are documented
below, including references to the literature and where they are used
inside :mod:`MDAnalysis.analysis`.

Data files
----------

.. data:: Rama_ref

   Reference Ramachandran histogram for :class:`MDAnalysis.analysis.dihedrals.Ramachandran`.
   The data were calculated on a data set of 500 PDB structures taken from [Lovell2003]_.

.. data:: Janin_ref

   Reference Ramachandran histogram for :class:`MDAnalysis.analysis.dihedrals.Ramachandran`.
   The data were calculated on a data set of 500 PDB structures taken from [Lovell2003]_.

"""
from __future__ import absolute_import

__all__ = [
    "Rama_ref", "Janin_ref" # reference plots for Ramachandran and Janin classes
]

from pkg_resources import resource_filename

Rama_ref = resource_filename(__name__, 'rama_ref_data.npy')
Janin_ref = resource_filename(__name__, 'janin_ref_data.npy')

# This should be the last line: clean up namespace
del resource_filename
