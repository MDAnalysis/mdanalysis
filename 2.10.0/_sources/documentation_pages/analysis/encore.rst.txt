===============================================================================
 ENCORE Ensemble Similarity Calculations --- :mod:`MDAnalysis.analysis.encore`
===============================================================================

:Author: Matteo Tiberti, Wouter Boomsma, Tone Bengtsen
:Year: 2015-2017
:Copyright: Lesser GNU Public License v2.1+
:Maintainer: Matteo Tiberti <matteo.tiberti@gmail.com>, mtiberti on github

.. versionadded:: 0.16.0


.. deprecated:: 2.8.0
   This module is deprecated in favour of the 
   MDAKit `mdaencore <https://mdanalysis.org/mdaencore/>`_ and will be removed
   in MDAnalysis 3.0.0.

The module contains implementations of similarity measures between protein
ensembles described in :footcite:p:`LindorffLarsen2009`. The implementation and examples
are described in :footcite:p:`Tiberti2015`.

The module includes facilities for handling ensembles and trajectories through
the :class:`Universe` class, performing clustering or dimensionality reduction
of the ensemble space, estimating multivariate probability distributions from
the input data, and more. ENCORE can be used to compare experimental and
simulation-derived ensembles, as well as estimate the convergence of
trajectories from time-dependent simulations.

ENCORE includes three different methods for calculations of similarity measures
between ensembles implemented in individual functions:


+ **Harmonic Ensemble Similarity** : :func:`~MDAnalysis.analysis.encore.similarity.hes`
+ **Clustering Ensemble Similarity** : :func:`~MDAnalysis.analysis.encore.similarity.ces`
+ **Dimensional Reduction Ensemble Similarity** : :func:`~MDAnalysis.analysis.encore.similarity.dres`

as well as two methods to evaluate the convergence of trajectories:

+ **Clustering based convergence evaluation** : :func:`~MDAnalysis.analysis.encore.similarity.ces_convergence`
+ **Dimensionality-reduction based convergence evaluation** : :func:`~MDAnalysis.analysis.encore.similarity.dres_convergence`


When using this module in published work please cite :footcite:p:`Tiberti2015`.


Modules
-------

.. toctree::
   :maxdepth: 1

   ./encore/similarity
   ./encore/clustering
   ./encore/dimensionality_reduction
   ./encore/confdistmatrix
   ./encore/covariance
   ./encore/bootstrap
   ./encore/utils
   

References
----------

.. footbibliography::
