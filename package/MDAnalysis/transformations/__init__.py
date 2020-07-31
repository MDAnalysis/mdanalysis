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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""\
Trajectory transformations --- :mod:`MDAnalysis.transformations`
================================================================

The transformations submodule contains a collection of function-like classes to
modify the trajectory.
Coordinate transformations, such as PBC corrections and molecule fitting
are often required for some analyses and visualization, and the functions in
this module allow transformations to be applied on-the-fly.

A typical transformation class looks like this:

.. code-blocks:: python

    class transfomration(object):
        def __init__(self, *args, **kwargs):
            #  do some things
            #  save needed args as attributes.
            self.needed_var = args[0]

        def __call__(self, ts):
            #  apply changes to the Timestep object
            return ts

See `MDAnalysis.transformations.translate` for a simple example.

These transformation functions can be called by the user for any given timestep
of the trajectory, added as a workflow using :meth:`add_transformations`
of the :mod:`~MDAnalysis.coordinates.base`, or upon Universe creation using
the keyword argument `transformations`. Note that in the two latter cases, the
workflow cannot be changed after being defined. for example:

.. code-block:: python

    u = mda.Universe(GRO, XTC)
    ts = u.trajectory[0]
    trans = transformation(args)
    ts = trans(ts)

    #  or add as a workflow
    u.trajectory.add_transformations(trans)

Transformations can also be created as a closure/nested function.
In addition to the specific arguments that each transformation can take, they
also contain a wrapped function that takes a `Timestep` object as argument.
So, a closure-style transformation can be roughly defined as follows:

.. code-block:: python

    def transformation(*args,**kwargs):
        # do some things
            def wrapped(ts):
                # apply changes to the Timestep object
                return ts

            return wrapped

Note, to meet the need of serialization of universe, only transformation class
are used after MDAnlaysis 2.0.0. One can still write functions (closures) as in
MDA 1.x, but that these cannot be serialized and thus will not work with all
forms of parallel analysis. For detailed descriptions about how to write a
closure-style transformation, read the code in MDA 1.x as a reference
or read MDAnalysis UserGuide.


.. versionchanged:: 2.0.0
    Transformations should now be created as classes with a :meth:`__call__`
    method instead of being written as a function/closure.
"""

from .translate import translate, center_in_box
from .rotate import rotateby
from .positionaveraging import PositionAverager
from .fit import fit_translation, fit_rot_trans
from .wrap import wrap, unwrap
