.. Contains the formatted docstrings for the transformations located in 'mdanalysis/MDAnalysis/transformations'
.. _transformations:

*********************************************************
Trajectory transformations ("on-the-fly" transformations)
*********************************************************

.. module:: MDAnalysis.transformations

In MDAnalysis, a *transformation* is a function/function-like class
that modifies the data for the current :class:`Timestep` and returns the
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

The submodule :mod:`MDAnalysis.transformations` contains a
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

A simple *transformation* can also be a function that takes a
:class:`~MDAnalysis.coordinates.base.Timestep` as input, modifies it, and
returns it. If it takes no other arguments but a :class:`Timestep`
can be defined as the following example:

.. code-block:: python

    def up_by_2(ts):
        """
        Translate all coordinates by 2 angstroms up along the Z dimension.
        """
        ts.positions = ts.positions + np.array([0, 0, 2], dtype=np.float32)
        return ts

If the transformation requires other arguments besides the :class:`Timestep`,
the following two methods can be used to create such transformation:


.. _custom-transformations-class:

Creating complex transformation classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

It is implemented by inheriting from
:class:`MDAnalysis.transformations.base.TransformationBase`,
which defines :func:`__call__` for the transformation class
and can be applied directly to a :class:`Timestep`. :func:`_transform` has to
be defined and include the operations on the :class:`MDAnalysis.coordinates.base.Timestep`.

So, a transformation class can be roughly defined as follows:

.. code-block:: python

    from MDAnalysis.transformations import TransformationBase

    class up_by_x_class(TransformationBase):
        def __init__(self, distance):
            self.distance = distance

        def _transform(self, ts):
            ts.positions = ts.positions + np.array([0, 0, self.distance], dtype=np.float32)
            return ts

It is the default construction method in :mod:`MDAnalysis.transformations`
from release 2.0.0 onwards because it can be reliably serialized.
See :class:`MDAnalysis.transformations.translate` for a simple example.


.. _custom-transformations-closure:

Creating complex transformation closure functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Transformation can also be a wrapped function takes the :class:`Timestep` object as argument. 
So in this case, a transformation function (closure) can be roughly defined as follows:

.. code-block:: python

    def up_by_x_func(distance):
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

Although functions (closures) work as transformations, they are not used in
in MDAnalysis from release 2.0.0 onwards because they cannot be reliably
serialized and thus a :class:`Universe` with such transformations cannot be
used with common parallelization schemes (e.g., ones based on
:mod:`multiprocessing`).
For detailed descriptions about how to write a closure-style transformation,
please refer to MDAnalysis 1.x documentation.

.. _transformations-module:

Transformations in MDAnalysis
-----------------------------

The module :mod:`MDAnalysis.transformations` contains transformations that can
be immediately used in your own :ref:`workflows<workflows>`. In order to use
any of these transformations, the module must first be imported:

.. code-block:: python

   import MDAnalysis.transformations

A workflow can then be added to a trajectory as described above. Notably,
the parameter `max_threads` can be defined when creating a transformation
instance to limit the maximum threads.
(See :class:`MDAnalysis.transformations.base.TransformationBase` for more details) 
Whether a specific transformation can be used along with parallel analysis
can be assessed by checking its 
:attr:`~MDAnalysis.transformations.base.TransformationBase.parallelizable`
attribute.

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
    

.. _building-block-transformation:

Building blocks for Transformation Classes
------------------------------------------
Transformations normally ultilize the power of NumPy to get better performance
on array operations. However, when it comes to parallelism, NumPy will sometimes
oversubscribe the threads, either by hyper threading (when it uses OpenBlas backend),
or by working with other parallel engines (e.g. Dask).

In MDAnalysis, we use `threadpoolctl <https://github.com/joblib/threadpoolctl>`_
inside :class:`~MDAnalysis.transformations.base.TransformationBase` to control the maximum threads for transformations.

It is also possible to apply a global thread limit by setting the external environmental
varibale, e.g. :code:`OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
BLIS_NUM_THREADS=1 python script.py`. Read more about parallelism and resource management
in `scikit-learn documentations <https://scikit-learn.org/dev/computing/parallelism.html>`_.

Users are advised to benchmark code because interaction between different
libraries can lead to sub-optimal performance with defaults.

.. toctree::

   ./transformations/base

.. _implemented-transformations:

Currently implemented transformations
-------------------------------------

.. toctree::
   
   ./transformations/translate
   ./transformations/rotate
   ./transformations/positionaveraging
   ./transformations/fit
   ./transformations/wrap
   ./transformations/nojump
   ./transformations/boxdimensions

