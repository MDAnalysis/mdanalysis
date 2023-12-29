# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

r"""
Calculating path similarity --- :mod:`MDAnalysis.analysis.psa`
==========================================================================

:Author: Sean Seyler
:Year: 2015
:Copyright: GNU Public License v3

.. versionadded:: 0.10.0

.. deprecated:: 2.8.0

  This module is deprecated in favour of the mdakit
  `pathsimanalysis <https://github.com/MDAnalysis/PathSimAnalysis>`_ and
  will be removed in MDAnalysis 3.0.0.


The module contains code to calculate the geometric similarity of trajectories
using path metrics such as the Hausdorff or Fréchet distances
:footcite:p:`Seyler2015`. The path metrics are functions of two paths and return a
nonnegative number, i.e., a distance. Two paths are identical if their distance
is zero, and large distances indicate dissimilarity. Each path metric is a
function of the individual points (e.g., coordinate snapshots) that comprise
each path and, loosely speaking, identify the two points, one per path of a
pair of paths, where the paths deviate the most.  The distance between these
points of maximal deviation is measured by the root mean square deviation
(RMSD), i.e., to compute structural similarity.

One typically computes the pairwise similarity for an ensemble of paths to
produce a symmetric distance matrix, which can be clustered to, at a glance,
identify patterns in the trajectory data. To properly analyze a path ensemble,
one must select a suitable reference structure to which all paths (each
conformer in each path) will be universally aligned using the rotations
determined by the best-fit rmsds. Distances between paths and their structures
are then computed directly with no further alignment. This pre-processing step
is necessary to preserve the metric properties of the Hausdorff and Fréchet
metrics; using the best-fit rmsd on a pairwise basis does not generally
preserve the triangle inequality.

Note
----
The `PSAnalysisTutorial`_ outlines a typical application of PSA to
a set of trajectories, including doing proper alignment,
performing distance comparisons, and generating heat
map-dendrogram plots from hierarchical clustering.

.. _`PSAnalysisTutorial`: https://github.com/Becksteinlab/PSAnalysisTutorial


Helper functions and variables
------------------------------
The following convenience functions are used by other functions in this module.

.. autofunction:: sqnorm
.. autofunction:: get_msd_matrix
.. autofunction:: get_coord_axes


Classes, methods, and functions
-------------------------------

.. autofunction:: get_path_metric_func
.. autofunction:: hausdorff
.. autofunction:: hausdorff_wavg
.. autofunction:: hausdorff_avg
.. autofunction:: hausdorff_neighbors
.. autofunction:: discrete_frechet
.. autofunction:: dist_mat_to_vec

.. autoclass:: Path
   :members:

   .. attribute:: u_original

      :class:`MDAnalysis.Universe` object with a trajectory

   .. attribute:: u_reference

      :class:`MDAnalysis.Universe` object containing a reference structure

   .. attribute:: select

      string, selection for
      :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` to select frame
      from :attr:`Path.u_reference`

   .. attribute:: path_select

      string, selection for
      :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` to select atoms
      to compose :attr:`Path.path`

   .. attribute:: ref_frame

      int, frame index to select frame from :attr:`Path.u_reference`

   .. attribute:: u_fitted

      :class:`MDAnalysis.Universe` object with the fitted trajectory

   .. attribute:: path

      :class:`numpy.ndarray` object representation of the fitted trajectory

.. autoclass:: PSAPair

   .. attribute:: npaths

      int, total number of paths in the comparison in which *this*
      :class:`PSAPair` was generated

   .. attribute:: matrix_id

      (int, int), (row, column) indices of the location of *this*
      :class:`PSAPair` in the corresponding pairwise distance matrix

   .. attribute:: pair_id

      int, ID of *this* :class:`PSAPair` (the pair_id:math:`^\text{th}`
      comparison) in the distance vector corresponding to the pairwise distance
      matrix

   .. attribute:: nearest_neighbors

      dict, contains the nearest neighbors by frame index and the
      nearest neighbor distances for each path in *this* :class:`PSAPair`

   .. attribute:: hausdorff_pair

      dict, contains the frame indices of the Hausdorff pair for each path in
      *this* :class:`PSAPair` and the corresponding (Hausdorff) distance

.. autoclass:: PSAnalysis
   :members:

   .. attribute:: universes

      list of :class:`MDAnalysis.Universe` objects containing trajectories

   .. attribute:: u_reference

      :class:`MDAnalysis.Universe` object containing a reference structure

   .. attribute:: select

      string, selection for
      :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` to select frame
      from :attr:`PSAnalysis.u_reference`

   .. attribute:: path_select

      string, selection for
      :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` to select atoms
      to compose :attr:`Path.path`

   .. attribute:: ref_frame

      int, frame index to select frame from :attr:`Path.u_reference`

   .. attribute:: paths

      list of :class:`numpy.ndarray` objects representing the set/ensemble of
      fitted trajectories

   .. attribute:: D

      :class:`numpy.ndarray` which stores the calculated distance matrix


See Also
--------
:mod:`pathsimanalysis.psa`


.. Markup definitions
.. ------------------
..
.. |3Dp| replace:: :math:`N_p \times N \times 3`
.. |2Dp| replace:: :math:`N_p \times (3N)`
.. |3Dq| replace:: :math:`N_q \times N \times 3`
.. |2Dq| replace:: :math:`N_q \times (3N)`
.. |3D| replace:: :math:`N_p\times N\times 3`
.. |2D| replace:: :math:`N_p\times 3N`
.. |Np| replace:: :math:`N_p`


.. Rubric:: References

.. footbibliography::
"""

import warnings

from pathsimanalysis import (
    get_path_metric_func,
    sqnorm,
    get_msd_matrix,
    reshaper,
    get_coord_axes,
    hausdorff,
    hausdorff_wavg,
    hausdorff_avg,
    hausdorff_neighbors,
    discrete_frechet,
    dist_mat_to_vec,
    Path,
    PSAPair,
    PSAnalysis,
)
