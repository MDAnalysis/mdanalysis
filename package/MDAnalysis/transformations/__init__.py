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


"""
Trajectory transformations --- :mod:`MDAnalysis.transformations`
================================================================

The transformations submodule contains a collection of function-like classes to
modify the trajectory.
Coordinate transformations, such as PBC corrections and molecule fitting
are often required for some analyses and visualization, and the functions in
this module allow transformations to be applied on-the-fly.

A typical transformation class looks like this (note that we keep its name
lowercase because we will treat it as a function, thanks to the ``__call__``
method):

.. code-blocks:: python

    class transformation(object):
        def __init__(self, *args, **kwargs):
            #  do some things
            #  save needed args as attributes.
            self.needed_var = args[0]

        def __call__(self, ts):
            #  apply changes to the Timestep,
            #  or modify an AtomGroup and return Timestep

            return ts

As a concrete example we will write a transformation that rotates a group of
atoms around the z-axis through the center of geometry by a fixed increment
for every time step. We will use
:meth:`MDAnalysis.core.groups.AtomGroup.rotateby`
and simply increment the rotation angle every time the
transformation is called ::

    class spin_atoms(object):
        def __init__(self, atoms, dphi):
            #  Rotate atoms by dphi degrees for every ts around the z axis
            self.atoms = atoms
            self.dphi = dphi
            self.axis = np.array([0, 0, 1])

        def __call__(self, ts):
            phi = self.dphi * ts.frame
            self.atoms.rotateby(phi, self.axis)
            return ts

This transformation can be used as ::

   u = mda.Universe(PSF, DCD)
   u.trajectory.add_transformations(spin_atoms(u.select_atoms("protein"), 1.0))

Also see :mod:`MDAnalysis.transformations.translate` for a simple example.

These transformation functions can be called by the user for any given timestep
of the trajectory, added as a workflow using :meth:`add_transformations`
of the :mod:`~MDAnalysis.coordinates.base`, or upon Universe creation using
the keyword argument `transformations`. Note that in the two latter cases, the
workflow cannot be changed after being defined. for example:

.. code-block:: python

    u = mda.Universe(GRO, XTC)
    trans = transformation(args)
    u.trajectory.add_transformations(trans)

    #  it is equivalent to applying this transforamtion to each Timestep by
    ts = u.trajectory[0]
    ts_trans = trans(ts)

Transformations can also be created as a closure/nested function.
In addition to the specific arguments that each transformation can take, they
also contain a wrapped function that takes a `Timestep` object as argument.
So, a closure-style transformation can be roughly defined as follows:

.. code-block:: python

    def transformation(*args,**kwargs):
        # do some things

            def wrapped(ts):
                #  apply changes to the Timestep,
                #  or modify an AtomGroup and return Timestep

                return ts

            return wrapped

.. Note::
   Although functions (closures) work as transformations, they are not used in
   in MDAnalysis from release 2.0.0 onwards because they cannot be reliably
   serialized and thus a :class:`Universe` with such transformations cannot be
   used with common parallelization schemes (e.g., ones based on
   :mod:`multiprocessing`).
   For detailed descriptions about how to write a closure-style transformation,
   please refer to MDAnalysis 1.x documentation.


.. versionchanged:: 2.0.0
    Transformations should now be created as classes with a :meth:`__call__`
    method instead of being written as a function/closure.
"""

from .base import TransformationBase
from .translate import translate, center_in_box
from .rotate import rotateby
from .positionaveraging import PositionAverager
from .nojump import NoJump
from .fit import fit_translation, fit_rot_trans
from .wrap import wrap, unwrap
from .boxdimensions import set_dimensions
