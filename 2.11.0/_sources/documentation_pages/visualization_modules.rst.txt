.. Contains the formatted docstrings from the visualization modules located in 'mdanalysis/MDAnalysis/visualization'

*********************
Visualization modules
*********************

The :mod:`MDAnalysis.visualization` namespace contains code to carry out
analyses which return data that is specifically-tailored for visualization.

Please see the individual module documentation for additional references and
citation information.

These modules are not imported by default; in order to use them one has to
import them from :mod:`MDAnalysis.visualization`, for instance: ::

    import MDAnalysis.visualization.streamlines

.. Note::

  Some of the modules require additional Python packages such as matplotlib_ or
  scipy_.

.. _matplotlib: http://matplotlib.org
.. _scipy: https://scipy.org/


Visualization of Lipid Flow
===========================

.. toctree::
   :maxdepth: 1

   visualization/streamlines
   visualization/streamlines_3D



