# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
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
from __future__ import absolute_import
import warnings

from .groups import (Atom, AtomGroup, Residue, ResidueGroup, Segment,
                     SegmentGroup)
from . import universe


def deprecate_class(class_new, message):
    """utility to deprecate a class"""

    class new_class(class_new):
        def __init__(self, *args, **kwargs):
            super(new_class, self).__init__(*args, **kwargs)
            warnings.warn(message, DeprecationWarning)

    return new_class


Universe = deprecate_class(
    universe.Universe,
    "MDAnalysis.core.AtomGroup.Universe has been removed."
    "Please use MDAnalysis.Universe."
    "This stub will be removed in 1.0")

_group_message = ("MDAnalysis.core.AtomGroup.{0} has been removed."
                  "Please use MDAnalysis.groups.{0}"
                  "This stub will be removed in 1.0")

Atom = deprecate_class(Atom, message=_group_message.format('Atom'))
AtomGroup = deprecate_class(
    AtomGroup, message=_group_message.format('AtomGroup'))

Residue = deprecate_class(Residue, message=_group_message.format('Residue'))
ResidueGroup = deprecate_class(
    ResidueGroup, message=_group_message.format('ResidueGroup'))

Segment = deprecate_class(Segment, message=_group_message.format('Segment'))
SegmentGroup = deprecate_class(
    SegmentGroup, message=_group_message.format('SegmentGroup'))

__all__ = [
    'Universe', 'Atom', 'AtomGroup', 'Residue', 'ResidueGroup', 'Segment',
    'SegmentGroup'
]
