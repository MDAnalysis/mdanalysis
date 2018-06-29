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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""\
Wrap/Unwrap trajectory --- :mod:`MDAnalysis.transformations.translate`
====================================================================

###
    
"""
from __future__ import absolute_import

import numpy as np

from MDAnalysis.lib.mdamath import make_whole

def unwrap(atomgroup):
    """Makes all the molecules in the given atomgroup, that have been broken over
    periodic boundary conditions, whole again. This is useful for visualisation and
    some calculations such as the center of mass.
    
    This function is most useful when used in conjunction with other transformations such 
    as translate and center. 
    
    Parameters
    ----------
    atomgroup:
        atomgroup containing all the molecules that are to be made whole again.
    
    Returns
    -------
    :class:`~MDAnalysis.coordinates.base.Timestep` object
    
    Note
    ----
    If the AtomGroup provided has molecules whose bonds are not described in the
    AtomGroup.bonds property will cause the transformation to fail, but won't be
    made whole again.
    """
    

    try:
        atomgroup.positions
    except AttributeError:
        raise AttributeError("{} is not an AtomGroup".format(atomgroup))
    
    def wrapped(ts):
        make_whole(atomgroup)
            
        return ts
    
    return wrapped
