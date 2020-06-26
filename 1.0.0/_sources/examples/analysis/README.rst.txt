.. -*- coding: utf-8 -*-
.. _analysis:

========
Analysis
========

The :mod:`~MDAnalysis.analysis` module of MDAnalysis provides the tools needed to analyse your data. 
Several analyses are included with the package. These range from standard algorithms 
(e.g. :ref:`calculating root mean squared quantities <alignment-and-rms>`) to unique algorithms such as 
the `path similarity analysis <examples/analysis/trajectory_similarity/psa.ipynb>`__.

Generally these bundled analyses are contributed by various researchers who use the code for their own work. 
Please refer to the individual module documentation or relevant user guide tutorials for additional 
references and citation information.

If you need functionality that is not already provided in MDAnalysis, there are 
`several ways to write your own analysis <examples/analysis/custom_trajectory_analysis.ipynb>`__.



Imports and dependencies
========================
 
Analysis modules are not imported by default. In order to use them, you will need to import 
them separately, e.g.::

    from MDAnalysis.analysis import align

.. note::

    Several modules in :mod:`MDAnalysis.analysis` require additional Python packages. 
    For example, :mod:`~MDAnalysis.analysis.encore` makes use of `scikit-learn <http://scikit-learn.org/>`__. 
    The Python packages are not automatically installed with `pip`, although they are with `conda`.

    Other modules require external programs. For example, :mod:`~MDAnalysis.analysis.hole` requires 
    the `HOLE <http://www.holeprogram.org/>`_ programs. You will need to install these yourself.


.. include:: /examples/analysis/alignment_and_rms/README.rst

.. include:: /examples/analysis/distances_and_contacts/README.rst

.. include:: /examples/analysis/trajectory_similarity/README.rst

.. include:: /examples/analysis/structure/README.rst

.. include:: /examples/analysis/volumetric/README.rst

.. include:: /examples/analysis/reduced_dimensions/README.rst

.. include:: /examples/analysis/polymers_and_membranes/README.rst