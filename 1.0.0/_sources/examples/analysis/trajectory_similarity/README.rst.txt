.. -*- coding: utf-8 -*-
.. _trajectory-similarity:

======================
Trajectory similarity
======================

A molecular dynamics trajectory with :math:`N` atoms can be considered through a path through :math:`3N`-dimensional molecular configuration space. MDAnalysis contains a number of algorithms to compare the conformational ensembles of different trajectories. Most of these are in the :mod:`MDAnalysis.analysis.encore` module (:cite:`tiberti_encore_2015`) and compare estimated probability distributions to measure similarity. The `path similarity analysis </examples/analysis/trajectory_similarity/psa.html>`_ compares the RMSD between pairs of structures in conformation transition paths. :mod:`MDAnalysis.analysis.encore` also contains functions for evaluating the conformational convergence of a trajectory using the `similarity over conformation clusters </examples/analysis/trajectory_similarity/clustering_ensemble_similarity.html>`_ or `similarity in a reduced dimensional space </examples/analysis/trajectory_similarity/dimension_reduction_ensemble_similarity.html>`_.


.. toctree::
   :maxdepth: 1

   /examples/analysis/trajectory_similarity/psa
   /examples/analysis/trajectory_similarity/harmonic_ensemble_similarity
   /examples/analysis/trajectory_similarity/clustering_ensemble_similarity
   /examples/analysis/trajectory_similarity/dimension_reduction_ensemble_similarity
   /examples/analysis/trajectory_similarity/convergence