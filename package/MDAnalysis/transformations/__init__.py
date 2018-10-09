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
Trajectory transformations --- :mod:`MDAnalysis.transformations`
================================================================

The transformations submodule contains a collection of functions to modify the 
trajectory. Coordinate transformations, such as PBC corrections and molecule fitting
are often required for some analyses and visualization, and the functions in this 
module allow transformations to be applied on-the-fly. 
These transformation functions can be called by the user for any given
timestep of the trajectory, added as a workflow using :meth:`add_transformations`
of the :mod:`~MDAnalysis.coordinates.base` module, or upon Universe creation using
the keyword argument `transformations`. Note that in the two latter cases, the 
workflow cannot be changed after being defined.

In addition to the specific arguments that each transformation can take, they also 
contain a wrapped function that takes a `Timestep` object as argument. 
So, a transformation can be roughly defined as follows:

.. code-block:: python

    def transformations(*args,**kwargs):
        # do some things
            def wrapped(ts):
                # apply changes to the Timestep object
                return ts
            
            return wrapped


See `MDAnalysis.transformations.translate` for a simple example.

Currently implemented transformations are:
    
    - translate: translate the coordinates of a given trajectory frame by a given vector.
    - center_in_box: translate the coordinates of a given trajectory frame so that a given
      AtomGroup is centered in the unit cell
    - rotateby: rotates the coordinates by a given angle arround an axis formed by a direction 
      and a point    

Examples
--------

e.g. translate the coordinates of a frame:

.. code-block:: python

    u = MDAnalysis.Universe(topology, trajectory)
    new_ts = MDAnalysis.transformations.translate([1,1,1])(u.trajectory.ts)

e.g. create a workflow and adding it to the trajectory:

.. code-block:: python

    u = MDAnalysis.Universe(topology, trajectory)
    workflow = [MDAnalysis.transformations.translate([1,1,1]), 
                MDAnalysis.transformations.translate([1,2,3])]
    u.trajetory.add_transformations(*workflow)

e.g. giving a workflow as a keyword argument when defining the universe:

.. code-block:: python
    
    workflow = [MDAnalysis.transformations.translate([1,1,1]), 
                MDAnalysis.transformations.translate([1,2,3])]
    u = MDAnalysis.Universe(topology, trajectory, transformations = *workflow)
    
    
"""

from __future__ import absolute_import

from .translate import translate, center_in_box
from .rotate import rotateby
