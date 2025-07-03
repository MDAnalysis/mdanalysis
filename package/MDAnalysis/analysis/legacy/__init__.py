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

""":mod:`MDAnalysis.analysis.legacy` --- Legacy analysis code
==========================================================

.. versionadded:: 0.16.0

The :mod:`MDAnalysis.analysis.legacy` package contains analysis
modules that are not or only incompletely tested and not regularly
maintained. They nevertheless still provide useful and sometimes
unique analysis capabilities and are therefore provided **as is**.

.. warning::

   Code in this module is not regularly maintained. Please use it very
   carefully.

If you want to use code from this module then you will have to import
it explicitly. For example, ::

   from MDAnalysis.analysis.legacy import x3dna

(For further discussion, see `Issue 743`_.)


.. _Issue 743: https://github.com/MDAnalysis/mdanalysis/issues/743

"""
