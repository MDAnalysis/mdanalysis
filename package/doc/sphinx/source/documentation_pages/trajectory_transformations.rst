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


See :func:`MDAnalysis.transformations.translate` for a simple example.    
    

.. rubric:: Examples

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
    
    
.. rubric:: Currently implemented transformations

.. toctree::
   
   ./transformations/translate
