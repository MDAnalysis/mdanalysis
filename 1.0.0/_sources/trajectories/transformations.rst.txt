.. -*- coding: utf-8 -*-
.. _transformations:

On-the-fly transformations
==========================

An on-the-fly transformation is a function that silently modifies the dynamic data contained in a trajectory :class:`~MDAnalysis.coordinates.base.Timestep` (typically coordinates) as it is loaded into memory. It is called for each current time step to transform data into your desired representation. A transformation function must also return the current :class:`~MDAnalysis.coordinates.base.Timestep`, as transformations are often chained together.

The :mod:`MDAnalysis.transformations` module contains a collection of transformations. For example, :func:`~MDAnalysis.transformations.fit.fit_rot_trans` can perform a mass-weighted alignment on an :class:`~MDAnalysis.core.groups.AtomGroup` to a reference.

.. ipython:: python

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import TPR, XTC
    from MDAnalysis import transformations as trans

    u = mda.Universe(TPR, XTC)
    protein = u.select_atoms('protein')
    align_transform = trans.fit_rot_trans(protein, protein, weights='mass')
    u.trajectory.add_transformations(align_transform)

Other implemented transformations include functions to :mod:`~MDAnalysis.transformations.translate`, :mod:`~MDAnalysis.transformations.rotate`, :mod:`~MDAnalysis.transformations.fit` an :class:`~MDAnalysis.core.groups.AtomGroup` to a reference, and :mod:`~MDAnalysis.transformations.wrap` or unwrap groups in the unit cell. 

Although you can only call :meth:`~MDAnalysis.coordinates.base.ProtoReader.add_transformations` *once*, you can pass in multiple transformations in a list, which will be executed in order. For example, the below workflow:

* makes all molecules whole (unwraps them over periodic boundary conditions)
* centers the protein in the center of the box
* wraps water back into the box

.. ipython:: python

    # create new Universe for new transformations
    u = mda.Universe(TPR, XTC)
    protein = u.select_atoms('protein')
    water = u.select_atoms('resname SOL')
    workflow = [trans.unwrap(u.atoms),
                trans.center_in_box(protein, center='geometry'),
                trans.wrap(water, compound='residues')]
    u.trajectory.add_transformations(*workflow)

If your transformation does not depend on something within the :class:`~MDAnalysis.core.universe.Universe` (e.g. a chosen :class:`~MDAnalysis.core.groups.AtomGroup`), you can also create a :class:`~MDAnalysis.core.universe.Universe` directly with transformations. The code below translates coordinates 1 angstrom up on the z-axis:

.. ipython:: python

    u = mda.Universe(TPR, XTC, transformations=[trans.translate([0, 0, 1])])

If you need a different transformation, it is easy to implement your own.

----------------------
Custom transformations
----------------------

At its core, a transformation function must only take a :class:`~MDAnalysis.coordinates.base.Timestep` as its input and return the :class:`~MDAnalysis.coordinates.base.Timestep` as the output.

.. ipython:: python

    def up_by_2(ts):
        """Translates atoms up by 2 angstrom"""
        ts.positions += np.array([0.0, 0.0, 0.2])
        return ts
    
    u = mda.Universe(TPR, XTC, transformations=[up_by_2])


If your transformation needs other arguments, you will need to wrap your core transformation with a wrapper function that can accept the other arguments.

.. ipython:: python

    def up_by_x(x):
        """Translates atoms up by x angstrom"""
        def wrapped(ts):
            """Handles the actual Timestep"""
            ts.positions += np.array([0.0, 0.0, float(x)])
            return ts
        return wrapped
    
    # load Universe with transformations that move it up by 7 angstrom
    u = mda.Universe(TPR, XTC, transformations=[up_by_x(5), up_by_x(2)])

    
Alternatively, you can use :func:`functools.partial` to substitute the other arguments.

.. ipython:: python

    import functools

    def up_by_x(ts, x):
        ts.positions += np.array([0.0, 0.0, float(x)])
        return ts
    
    up_by_5 = functools.partial(up_by_x, x=5)
    u = mda.Universe(TPR, XTC, transformations=[up_by_5])

On-the-fly transformation functions can be applied to any property of a Timestep, not just the atom positions. For example, to give each frame of a trajectory a box:

.. ipython:: python
    
    def set_box(ts):
        # creates box of length 10 on x-axis, 20 on y-axis, 30 on z-axis
        # angles are all 90 degrees
        ts.dimensions = [10, 20, 30, 90, 90, 90]
        return ts
    
    u = mda.Universe(TPR, XTC, transformations=[set_box])

