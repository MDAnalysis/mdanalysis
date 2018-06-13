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
Translate trajectory --- :mod:`MDAnalysis.transformations.translate`
====================================================================

Translate the coordinates of a given trajectory by a given vector.
This is a wrapped function so usage is as the following example:

    ts = MDAnalysis.transformations.translate(vector)(timestep)
    
"""
from __future__ import absolute_import

import numpy as np

def translate(vector):
    """
    Translates the coordinates of a given :class:`~MDAnalysis.coordinates.base.Timestep`
    instance by a given vector.
    
    Example
    -------
        ts = MDAnalysis.transformations.translate([1,2,3])(ts)    
    
    Parameters
    ----------
    vector: list
        coordinates of the vector to which the coordinates will be translated
    ts: Timestep
        frame that will be transformed
     
    Returns
    -------
    :class:`~MDAnalysis.coordinates.base.Timestep` object
    
    """       
    def wrapped(ts):
        if len(vector)>2:
            v = np.float32(vector)
            ts.positions += v
        else:
            raise ValueError("{} vector in too short".format(vector))
        
        return ts
    
    return wrapped
