.. Contains the formatted docstrings from the analysis modules located in 'mdanalysis/MDAnalysis/analysis', although in some cases the documentation imports functions and docstrings from other files which are also curated to reStructuredText markup.

**************************
Specific analysis modules
**************************

These modules are located in the mdanalysis/MDAnalysis/analysis/ directory and abstract functionality for specific analysis tasks (for example, bilayer :mod:`leaflet` selection). In many cases it will be preferable to use these provided tools over manually-coded alternatives because some components may be written in C for time-critical routines. 

.. automodule:: align
   :members:

.. automodule:: contacts
   :members:

.. automodule:: distances
   :members:

.. automodule:: leaflet
   
	.. autoclass:: LeafletFinder

	.. autofunction:: optimize_cutoff	  







