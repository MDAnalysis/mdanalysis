.. Contains the formatted docstrings for the transformations located in 'mdanalysis/MDAnalysis/transformations'
.. _transformations:

**************************
Trajectory transformations
**************************

.. module:: MDAnalysis.transformations

The transformations submodule contains a collection of functions to modify the 
trajectory. Coordinate transformations, such as PBC corrections and molecule fitting
are often required for some analyses and visualization, and the functions in this 
module allow transformations to be applied on-the-fly. 
These transformation functions can be called by the user for any given
timestep of the trajectory, added as a workflow using :meth:`add_transformations`
of the :mod:`~MDAnalysis.coordinates.base` module, or upon Universe creation using
the keyword argument `transformations`. Note that in the two latter cases, the 
workflow cannot be changed after being defined.

A simple transformation that takes no other arguments but the `Timestep` can be defined
as the following example:

.. code-block:: python

    def up_by_2(ts):
    	"""
    	Translate all coordinates by 2 angstroms up along the Z dimension.
    	"""
    	ts.positions = ts.positions + np.array([0, 0, 2], dtype=np.float32)
    	return ts


If the transformation requires other arguments besides the `Timestep`, the transformation 
takes these arguments, while a wrapped function takes the `Timestep` object as 
argument. 
So, a transformation can be roughly defined as follows:

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


See :func:`MDAnalysis.transformations.translate` for a simple example.    
    

.. rubric:: Examples

Translating the coordinates of a frame:

.. code-block:: python

    u = MDAnalysis.Universe(topology, trajectory)
    new_ts = MDAnalysis.transformations.translate([1,1,1])(u.trajectory.ts)

e.g. create a workflow and adding it to the trajectory:

.. code-block:: python

    u = MDAnalysis.Universe(topology, trajectory)
    workflow = [MDAnalysis.transformations.translate([1,1,1]), 
                MDAnalysis.transformations.translate([1,2,3])]
    u.trajcetory.add_transformations(*workflow)

e.g. giving a workflow as a keyword argument when defining the universe:

.. code-block:: python
    
    workflow = [MDAnalysis.transformations.translate([1,1,1]), 
                MDAnalysis.transformations.translate([1,2,3])]
    u = MDAnalysis.Universe(topology, trajectory, transformations = *workflow)
    
    
.. rubric:: Currently implemented transformations

.. toctree::
   
   ./transformations/translate
   ./transformations/rotate