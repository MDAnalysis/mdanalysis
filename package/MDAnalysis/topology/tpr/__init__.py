# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""
TPR support
===========

The :mod:`MDAnalysis.topology.tpr` module is required for the
:mod:`MDAnalysis.topology.TPRParser` module.

.. autodata:: SUPPORTED_VERSIONS

.. rubric:: Sub-modules

* :mod:`MDAnalysis.topology.tpr.setting`
* :mod:`MDAnalysis.topology.tpr.obj`
* :mod:`MDAnalysis.topology.tpr.utils`

"""
from __future__ import absolute_import

from .setting import SUPPORTED_VERSIONS

__all__ = ["obj", "setting", "utils"]
