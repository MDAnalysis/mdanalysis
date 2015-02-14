# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#

# Based on Biopython's KDTree module and used under the Biopython license:
# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# Changes to the original __init__.py:
# 2008-08-24 Oliver Beckstein <orbeckst@gmail.com>
# Added the modified NeighborSearch module to the KDTree package

"""
.. KDTree distance searches --- :mod:`MDAnalysis.KDTree`
.. =====================================================

:Author: Thomas Hamelryck
:Year:   2002
:Licence: Biopython

The KD tree data structure can be used for all kinds of searches that
involve N-dimensional vectors. For example, neighbor searches (find all points
within a radius of a given point) or finding all point pairs in a set
that are within a certain radius of each other. See "Computational Geometry:
Algorithms and Applications" (Mark de Berg, Marc van Kreveld, Mark Overmars,
Otfried Schwarzkopf) [deBerg2000].
"""
__all__ = ['KDTree', 'NeighborSearch']

from KDTree import KDTree
