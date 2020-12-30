.. Contains the formatted docstrings for the transformations located in 'mdanalysis/MDAnalysis/transformations'
.. _transformations:

*********************************************************
Trajectory transformations ("on-the-fly" transformations)
*********************************************************

.. module:: MDAnalysis.transformations

In MDAnalysis, a *transformation* is a function that modifies the data
for the current :class:`Timestep` and returns the
:class:`Timestep`. For instance, coordinate transformations, such as
PBC corrections and molecule fitting are often required for some
analyses and visualization. Transformation functions
(``transformation_1`` and ``transformation_2`` in the following
example) can be called by the user for any given :class:`Timestep` of
the trajectory,

.. code-block:: python

    u = MDAnalysis.Universe(topology, trajectory)

    for ts in u.trajectory:
       ts = transformation_2(transformation_1(ts))

where they change the coordinates of the timestep ``ts`` in
place. There is nothing special about these transformations except
that they have to be written in such a way that they change the
:class:`Timestep` in place.

As described under :ref:`workflows`, multiple transformations can be
grouped together and associated with a trajectory so that the
trajectory is **transformed on-the-fly**, i.e., the data read from the
trajectory file will be changed before it is made available in, say,
the :attr:`AtomGroup.positions` attribute.

The  submodule :mod:`MDAnalysis.transformations` contains a
collection of transformations (see :ref:`transformations-module`) that
can be immediately used but one can always write custom
transformations (see :ref:`custom-transformations`).


.. _workflows:

Workflows
---------

Instead of manually applying transformations, it is much more
convenient to associate a whole *workflow* of transformations with a
trajectory and have the transformations be called automatically.

A workflow is a sequence (tuple or list) of transformation functions
that will be applied in this order. For example,

.. code-block:: python

    workflow = [transformation_1, transformation_2]
    
would effectively result in

.. code-block:: python
		
     ts = transformation_2(transformation_1(ts))

for every time step in the trajectory.

One can add a workflow using the
:meth:`Universe.trajectory.add_transformations
<MDAnalysis.coordinates.base.ReaderBase.add_transformations>` method
of a trajectory (where the list ``workflow`` is taken from the example
above),

.. code-block:: python

    u.trajectory.add_transformations(*workflow)

or upon :class:`Universe <MDAnalysis.core.universe.Universe>`
creation using the keyword argument `transformations`:

.. code-block:: python
    
    u = MDAnalysis.Universe(topology, trajectory, transformations=workflow)

Note that in these two cases, the workflow cannot be changed after having
being added.


.. _custom-transformations:

Creating transformations
------------------------

A *transformation* is a function that takes a
:class:`~MDAnalysis.coordinates.base.Timestep` as input, modifies it, and
returns it.

A simple transformation that takes no other arguments but a :class:`Timestep`
can be defined as the following example:

.. code-block:: python

    def up_by_2(ts):
    	"""
    	Translate all coordinates by 2 angstroms up along the Z dimension.
    	"""
    	ts.positions = ts.positions + np.array([0, 0, 2], dtype=np.float32)
    	return ts


If the transformation requires other arguments besides the :class:`Timestep`,
the transformation takes these arguments, while a wrapped function takes the
:class:`Timestep` object as argument.  So, a transformation can be roughly
defined as follows:

.. code-block:: python

    def up_by_x(distance):
    	"""
    	Creates a transformation that will translate all coordinates by a given amount along the Z dimension.
    	"""
    	def wrapped(ts):
        	ts.positions = ts.positions + np.array([0, 0, distance], dtype=np.float32)
        return ts
    return wrapped
    
    
An alternative to using a wrapped function is using partials from :mod:`functools`. The
above function can be written as:

.. code-block:: python

    import functools

    def up_by_x(ts, distance):
    	ts.positions = ts.positions + np.array([0, 0, distance], dtype=np.float32)
    	return ts

    up_by_2 = functools.partial(up_by_x, distance=2)


See :func:`MDAnalysis.transformations.translate` for a simple
example of such a type of function.


.. _transformations-module:

Transformations in MDAnalysis
-----------------------------

The module :mod:`MDAnalysis.transformations` contains transformations that can
be immediately used in your own :ref:`workflows<workflows>`. In order to use
any of these transformations, the module must first be imported:

.. code-block:: python

   import MDAnalysis.transformations

A workflow can then be added to a trajectory as described above.

See :ref:`implemented-transformations` for more on the existing
transformations in :mod:`MDAnalysis.transformations`.
    

How to transformations
----------------------

Translating the coordinates of a single frame (although one would normally add
the transformation to a :ref:`workflow<workflows>`, as shown in the subsequent
examples):

.. code-block:: python

    u = MDAnalysis.Universe(topology, trajectory)
    new_ts = MDAnalysis.transformations.translate([1,1,1])(u.trajectory.ts)
    
Create a workflow and add it to the trajectory:

.. code-block:: python

    u = MDAnalysis.Universe(topology, trajectory)
    workflow = [MDAnalysis.transformations.translate([1,1,1]), 
                MDAnalysis.transformations.translate([1,2,3])]
    u.trajectory.add_transformations(*workflow)

Giving a workflow as a keyword argument when defining the universe:

.. code-block:: python
    
    workflow = [MDAnalysis.transformations.translate([1,1,1]), 
                MDAnalysis.transformations.translate([1,2,3])]
    u = MDAnalysis.Universe(topology, trajectory, transformations=workflow)
    

.. _implemented-transformations:

Currently implemented transformations
-------------------------------------

.. toctree::
   
   ./transformations/translate
   ./transformations/rotate
   ./transformations/positionaveraging
   ./transformations/fit
