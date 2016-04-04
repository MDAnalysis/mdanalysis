# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
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

The module contains code to calculate the geometric similarity of trajectories
using path metrics such as the Hausdorff or Fréchet distances
[Seyler2015]_. The path metrics are functions of two paths and return a
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

.. SeeAlso:: The `PSAnalysisTutorial`_ outlines a typical application of PSA to
             a set of trajectories, including doing proper alignment,
             performing distance comparisons, and generating heat
             map-dendrogram plots from hierarchical clustering.


.. Rubric:: References


.. [Seyler2015] Seyler SL, Kumar A, Thorpe MF, Beckstein O (2015)
                Path Similarity Analysis: A Method for Quantifying
                Macromolecular Pathways. PLoS Comput Biol 11(10): e1004568.
                doi: `10.1371/journal.pcbi.1004568`_

.. _`10.1371/journal.pcbi.1004568`: http://dx.doi.org/10.1371/journal.pcbi.1004568
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

.. autoclass:: PDBToBinaryTraj
   :members:

   .. attribute:: universe

      :class:`MDAnalysis.Universe` object with a trajectory

   .. attribute:: frames

      :attr:`MDAnalysis.Universe.trajectory`

   .. attribute:: newname

      string, filename for converted trajectory, including file extension

.. autoclass:: Path
   :members:

   .. attribute:: u_original

      :class:`MDAnalysis.Universe` object with a trajectory

   .. attribute:: u_reference

      :class:`MDAnalysis.Universe` object containing a reference structure

   .. attribute:: ref_select

      string, selection for
      :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.select_atoms` to select frame
      from :attr:`Path.u_reference`

   .. attribute:: path_select

      string, selection for
      :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.select_atoms` to select atoms
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

   .. attribute:: ref_select

      string, selection for
      :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.select_atoms` to select frame
      from :attr:`PSAnalysis.u_reference`

   .. attribute:: path_select

      string, selection for
      :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.select_atoms` to select atoms
      to compose :attr:`Path.path`

   .. attribute:: ref_frame

      int, frame index to select frame from :attr:`Path.u_reference`

   .. attribute:: filename

      string, name of file to store calculated distance matrix
      (:attr:`PSAnalysis.D`)

   .. attribute:: paths

      list of :class:`numpy.ndarray` objects representing the set/ensemble of
      fitted trajectories

   .. attribute:: D

      string, name of file to store calculated distance matrix
      (:attr:`PSAnalysis.D`)


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

"""
import six
from six.moves import range, cPickle

import numpy as np

import MDAnalysis
import MDAnalysis.analysis.align
from MDAnalysis import NoDataError

import os

import logging
logger = logging.getLogger('MDAnalysis.analysis.psa')


def get_path_metric_func(name):
    """Selects a path metric function by name.

    :Arguments:
      *name*
         string, name of path metric

    :Returns:
      The path metric function specified by *name* (if found).
    """
    path_metrics = {
            'hausdorff' : hausdorff,
            'weighted_average_hausdorff' : hausdorff_wavg,
            'average_hausdorff' : hausdorff_avg,
            'hausdorff_neighbors' : hausdorff_neighbors,
            'discrete_frechet' : discrete_frechet
    }
    try:
        return path_metrics[name]
    except KeyError as key:
        print("Path metric {0} not found. Valid selections: ".format(key))
        for name in path_metrics.keys(): print("  \"{0}\"".format(name))


def sqnorm(v, axis=None):
    """Compute the sum of squares of elements along specified axes.

    :Arguments:
      *v*
         :class:`numpy.ndarray` of coordinates
      *axes*
         None or int or tuple of ints, optional
         Axes or axes along which a sum is performed. The default (*axes* =
         ``None``) performs a sum over all the dimensions of the input array.
         The value of *axes* may be negative, in which case it counts from the
         last axis to the zeroth axis.

    :Returns:
      float, the sum of the squares of the elements of *v* along *axes*
    """
    return np.sum(v*v, axis=axis)


def get_msd_matrix(P, Q, axis=None):
    """Generate the matrix of pairwise mean-squared deviations (MSDs) between
    all pairs of points in *P* and *Q*, each pair having a point from *P* and a
    point from *Q*.

    *P* (*Q*) is a :class:`numpy.ndarray` of :math:`N_p` (:math:`N_q`) time
    steps, :math:`N` atoms, and :math:`3N` coordinates (e.g.,
    :meth:`MDAnalysis.core.AtomGroup.AtomGroup.coordinates`). The pairwise MSD
    matrix has dimensions :math:`N_p` by :math:`N_q`.

    :Arguments:
      *P*
         :class:`numpy.ndarray` representing a path
      *Q*
         :class:`numpy.ndarray` representing a path

    :Returns:
      :class:`numpy.ndarray` of pairwise MSDs between points in *P* and points
      in *Q*
    """
    return np.asarray([sqnorm(p - Q, axis=axis) for p in P])


def get_coord_axes(path):
    """Return the number of atoms and the axes (*axis*) corresponding to atoms
    and coordinates for a given path.

    The *path* is assumed to be a :class:`numpy.ndarray` where the 0th axis
    corresponds to a frame (a snapshot of coordinates). The :math:`3N`
    (Cartesian) coordinates are assumed to be either:
        (1) all in the 1st axis, starting with the x,y,z coordinates of the
            first atom, followed by the *x*,*y*,*z* coordinates of the 2nd, etc.
        (2) in the 1st *and* 2nd axis, where the 1st axis indexes the atom
            number and the 2nd axis contains the *x*,*y*,*z* coordinates of each
            atom.

    :Arguments:
      *path*
         :class:`numpy.ndarray` representing a path

    :Returns:
      (int, (int, ...)), the number of atoms and the axes containing coordinates
    """
    path_dimensions = len(path.shape)
    if path_dimensions == 3:
        N = path.shape[1]
        axis = (1,2) # 1st axis: atoms, 2nd axis: x,y,z coords
    elif path_dimensions == 2:
        # can use mod to check if total # coords divisible by 3
        N = path.shape[1] / 3
        axis = (1,) # 1st axis: 3N structural coords (x1,y1,z1,...,xN,xN,zN)
    else:
        err_str = "Path must have 2 or 3 dimensions; the first dimensions (axis"\
                + " 0) must correspond to frames, axis 1 (and axis 2, if"       \
                + " present) must contain atomic coordinates.".format(N)
        raise ValueError(err_str)
    return N, axis


def hausdorff(P, Q):
    r"""Calculate the Hausdorff distance between two paths.

    *P* (*Q*) is a :class:`numpy.ndarray` of :math:`N_p` (:math:`N_q`) time
    steps, :math:`N` atoms, and :math:`3N` coordinates (e.g.,
    :meth:`MDAnalysis.core.AtomGroup.AtomGroup.coordinates`). *P* (*Q*) has
    either shape |3Dp| (|3Dq|), or |2Dp| (|2Dq|) in flattened form.

    :Arguments:
      *P*
         :class:`numpy.ndarray` representing a path
      *Q*
         :class:`numpy.ndarray` representing a path

    :Returns:
      float, the Hausdorff distance between paths *P* and *Q*

    Example::
     >>> from MDAnalysis.tests.datafiles import PSF, DCD
     >>> u = Universe(PSF,DCD)
     >>> mid = len(u.trajectory)/2
     >>> ca = u.select_atoms('name CA')
     >>> P = numpy.array([
     ...                ca.positions for _ in u.trajectory[:mid:]
     ...              ]) # first half of trajectory
     >>> Q = numpy.array([
     ...                ca.positions for _ in u.trajectory[mid::]
     ...              ]) # second half of trajectory
     >>> hausdorff(P,Q)
     4.7786639840135905
     >>> hausdorff(P,Q[::-1]) # hausdorff distance w/ reversed 2nd trajectory
     4.7786639840135905
    """
    N, axis = get_coord_axes(P)
    d = get_msd_matrix(P, Q, axis=axis)
    return ( max( np.amax(np.amin(d, axis=0)),                                  \
                  np.amax(np.amin(d, axis=1)) ) / N  )**0.5


def hausdorff_wavg(P, Q):
    r"""Calculate the weighted average Hausdorff distance between two paths.

    *P* (*Q*) is a :class:`numpy.ndarray` of :math:`N_p` (:math:`N_q`) time
    steps, :math:`N` atoms, and :math:`3N` coordinates (e.g.,
    :meth:`MDAnalysis.core.AtomGroup.AtomGroup.coordinates`). *P* (*Q*) has
    either shape |3Dp| (|3Dq|), or |2Dp| (|2Dq|) in flattened form. The nearest
    neighbor distances for *P* (to *Q*) and those of *Q* (to *P*) are averaged
    individually to get the average nearest neighbor distance for *P* and
    likewise for *Q*. These averages are then summed and divided by 2 to get a
    measure that gives equal weight to *P* and *Q*.

    :Arguments:
      *P*
         :class:`numpy.ndarray` representing a path
      *Q*
         :class:`numpy.ndarray` representing a path

    :Returns:
      float, the weighted average Hausdorff distance between paths *P* and *Q*

    Example::
     >>> from MDAnalysis import Universe
     >>> from MDAnalysis.tests.datafiles import PSF, DCD
     >>> u = Universe(PSF,DCD)
     >>> mid = len(u.trajectory)/2
     >>> ca = u.select_atoms('name CA')
     >>> P = numpy.array([
     ...                ca.positions for _ in u.trajectory[:mid:]
     ...              ]) # first half of trajectory
     >>> Q = numpy.array([
     ...                ca.positions for _ in u.trajectory[mid::]
     ...              ]) # second half of trajectory
     >>> hausdorff_wavg(P,Q)
     2.5669644353703447
     >>> hausdorff_wavg(P,Q[::-1]) # weighted avg hausdorff dist w/ Q reversed
     2.5669644353703447
    """
    N, axis = get_coord_axes(P)
    d = get_msd_matrix(P, Q, axis=axis)
    out = 0.5*( np.mean(np.amin(d,axis=0)) + np.mean(np.amin(d,axis=1)) )
    return ( out / N )**0.5


def hausdorff_avg(P, Q):
    r"""Calculate the average Hausdorff distance between two paths.

    *P* (*Q*) is a :class:`numpy.ndarray` of :math:`N_p` (:math:`N_q`) time
    steps, :math:`N` atoms, and :math:`3N` coordinates (e.g.,
    :meth:`MDAnalysis.core.AtomGroup.AtomGroup.coordinates`). *P* (*Q*) has
    either shape |3Dp| (|3Dq|), or |2Dp| (|2Dq|) in flattened form. The nearest
    neighbor distances for *P* (to *Q*) and those of *Q* (to *P*) are all
    averaged together to get a mean nearest neighbor distance. This measure
    biases the average toward the path that has more snapshots, whereas weighted
    average Hausdorff gives equal weight to both paths.

    :Arguments:
      *P*
         :class:`numpy.ndarray` representing a path
      *Q*
         :class:`numpy.ndarray` representing a path

    :Returns:
      float, the average Hausdorff distance between paths *P* and *Q*

    Example::
     >>> from MDAnalysis.tests.datafiles import PSF, DCD
     >>> u = Universe(PSF,DCD)
     >>> mid = len(u.trajectory)/2
     >>> ca = u.select_atoms('name CA')
     >>> P = numpy.array([
     ...                ca.positions for _ in u.trajectory[:mid:]
     ...              ]) # first half of trajectory
     >>> Q = numpy.array([
     ...                ca.positions for _ in u.trajectory[mid::]
     ...              ]) # second half of trajectory
     >>> hausdorff_avg(P,Q)
     2.5669646575869005
     >>> hausdorff_avg(P,Q[::-1]) # hausdorff distance w/ reversed 2nd trajectory
     2.5669646575869005
    """
    N, axis = get_coord_axes(P)
    d = get_msd_matrix(P, Q, axis=axis)
    out = np.mean( np.append( np.amin(d,axis=0), np.amin(d,axis=1) ) )
    return ( out / N )**0.5


def hausdorff_neighbors(P, Q):
    r"""Calculate the Hausdorff distance between two paths.

    .. |3Dp| replace:: :math:`N_p \times N \times 3`
    .. |2Dp| replace:: :math:`N_p \times (3N)`
    .. |3Dq| replace:: :math:`N_q \times N \times 3`
    .. |2Dq| replace:: :math:`N_q \times (3N)`

    *P* (*Q*) is a :class:`numpy.ndarray` of :math:`N_p` (:math:`N_q`) time
    steps, :math:`N` atoms, and :math:`3N` coordinates (e.g.,
    :meth:`MDAnalysis.core.AtomGroup.AtomGroup.coordinates`). *P* (*Q*) has
    either shape |3Dp| (|3Dq|), or |2Dp| (|2Dq|) in flattened form.

    :Arguments:
      *P*
         :class:`numpy.ndarray` representing a path
      *Q*
         :class:`numpy.ndarray` representing a path

    :Returns:
      dictionary of two pairs of numpy arrays, the first pair containing the
      indices of (Hausdorff) nearest neighbors for *P* and *Q*, respectively, the
      second containing (corresponding) nearest neighbor distances for *P* and
      *Q*, respectively
    """
    N, axis = get_coord_axes(P)
    d = get_msd_matrix(P, Q, axis=axis)
    nearest_neighbors =                                                         \
            {'frames' :                                                         \
                (np.argmin(d, axis=1), np.argmin(d, axis=0)),                   \
             'distances' :                                                      \
                ((np.amin(d,axis=1)/N)**0.5, (np.amin(d, axis=0)/N)**0.5),      \
            }
    return nearest_neighbors


def discrete_frechet(P, Q):
    r"""Calculate the discrete Frechet distance between two paths.

    *P* (*Q*) is a :class:`numpy.ndarray` of :math:`N_p` (:math:`N_q`) time
    steps, :math:`N` atoms, and :math:`3N` coordinates (e.g.,
    :meth:`MDAnalysis.core.AtomGroup.AtomGroup.coordinates`). *P* (*Q*) has
    either shape |3Dp| (|3Dq|), or :|2Dp| (|2Dq|) in flattened form.

    :Arguments:
      *P*
         :class:`numpy.ndarray` representing a path
      *Q*
         :class:`numpy.ndarray` representing a path

    :Returns:
      float, the discrete Frechet distance between paths *P* and *Q*

    Example::
     >>> u = Universe(PSF,DCD)
     >>> mid = len(u.trajectory)/2
     >>> ca = u.select_atoms('name CA')
     >>> P = np.array([
     ...                ca.positions for _ in u.trajectory[:mid:]
     ...              ]) # first half of trajectory
     >>> Q = np.array([
     ...                ca.positions for _ in u.trajectory[mid::]
     ...              ]) # second half of trajectory
     >>> discrete_frechet(P,Q)
     4.7786639840135905
     >>> discrete_frechet(P,Q[::-1]) # frechet distance w/ 2nd trj reversed 2nd
     6.8429011177113832
    """
    N, axis = get_coord_axes(P)
    Np, Nq = len(P), len(Q)
    d = get_msd_matrix(P, Q, axis=axis)
    ca = -np.ones((Np, Nq))

    def c(i, j):
        """Compute the coupling distance for two partial paths formed by *P* and
        *Q*, where both begin at frame 0 and end (inclusive) at the respective
        frame indices :math:`i-1` and :math:`j-1`. The partial path of *P* (*Q*)
        up to frame *i* (*j*) is formed by the slicing ``P[0:i]`` (``Q[0:j]``).

        :func:`c` is called recursively to compute the coupling distance
        between the two full paths *P* and *Q*  (i.e., the discrete Frechet
        distance) in terms of coupling distances between their partial paths.

        :Arguments:
          *i*
             int, partial path of *P* through final frame *i-1*
          *j*
             int, partial path of *Q* through final frame *j-1*

        :Returns:
          float, the coupling distance between partial paths ``P[0:i]`` and
          ``Q[0:j]``
        """
        if ca[i,j] != -1 : return ca[i,j]
        if i > 0:
            if j > 0: ca[i,j] = max( min(c(i-1,j),c(i,j-1),c(i-1,j-1)), d[i,j] )
            else:     ca[i,j] = max( c(i-1,0), d[i,0] )
        elif j > 0:   ca[i,j] = max( c(0,j-1), d[0,j] )
        else:         ca[i,j] = d[0,0]
        return        ca[i,j]

    return ( c(Np-1, Nq-1) / N )**0.5


def dist_mat_to_vec(N, i, j):
    """Convert distance matrix indices (in the upper triangle) to the index of
    the corresponding distance vector.

    This is a convenience function to locate distance matrix elements (and the
    pair generating it) in the corresponding distance vector. The row index *j*
    should be greater than *i+1*, corresponding to the upper triangle of the
    distance matrix.

    :Arguments:
      *N*
         int, size of the distance matrix (of shape *N*-by-*N*)
      *i*
         int, row index (starting at 0) of the distance matrix
      *j*
         int, column index (starting at 0) of the distance matrix

    :Returns:
      int, index (of the matrix element) in the corresponding distance vector
    """
    if i > N or j > N:
        err_str = "Matrix indices are out of range; i and j must be less than"  \
                + " N = {0:d}".format(N)
        raise ValueError(err_str)
    if j > i:
        return (N*i) + j - (i+2)*(i+1)/2
    elif j < i:
        warn_str = "Column index entered (j = {:d} is smaller than row index"   \
                 + " (i = {:d}). Using symmetric element in upper triangle of"  \
                 + " distance matrix instead: i --> j, j --> i"
        print(warn_str)
        return (N*j) + i - (j+2)*(j+1)/2
    else:
        err_str = "Error in processing matrix indices; i and j must be integers"\
                + " less than integer N = {0:d} such that j >= i+1.".format(N)
        raise ValueError(err_str)


class PDBToBinaryTraj(object):

    def __init__(self, universe, output_type='.dcd', infix=''):
        self.universe = universe
        self.universe.atoms.write('new_top.pdb') # write first frame as topology

        self.frames = self.universe.trajectory
        base, ext = os.path.splitext(self.frames.filename)
        path, name = os.path.split(base)
        self.newname = name + infix + output_type

    def convert(self):
        w = MDAnalysis.Writer(self.newname, self.frames.numatoms)
        for ts in self.frames:
            w.write(ts)
        w.close_trajectory()


class Path(object):
    """Pre-process a :class:`MDAnalysis.Universe` object: (1) fit the
    trajectory to a reference structure, (2) convert fitted time series to a
    :class:`numpy.ndarray` representation of :attr:`Path.path`.

    The analysis is performed with :meth:`PSAnalysis.run` and stores the result
    in the :class:`numpy.ndarray` distance matrix :attr:`PSAnalysis.D`.
    :meth:`PSAnalysis.run` also generates a fitted trajectory and path from
    alignment of the original trajectories to a reference structure.

    .. versionadded:: 0.9.1
    """

    def __init__(self, universe, reference, ref_select='name CA',
                 path_select='all', ref_frame=0):
        """Setting up trajectory alignment and fitted path generation.

        :Arguments:
          *universe*
             :class:`MDAnalysis.Universe` object containing a trajectory
          *reference*
             reference structure; :class:`MDAnalysis.Universe` object; if
             ``None`` then *traj* is used (uses the current time step of the
             object) [``None``]
          *ref_select*
             The selection to operate on for rms fitting; can be one of:

             1. any valid selection string for
                :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.select_atoms` that
                produces identical selections in *mobile* and *reference*; or
             2. a dictionary ``{'mobile':sel1, 'reference':sel2}`` (the
                :func:`MDAnalysis.analysis.align.fasta2select` function returns
                such a dictionary based on a ClustalW_ or STAMP_ sequence
                alignment); or
             3. a tuple ``(sel1, sel2)``

             When using 2. or 3. with *sel1* and *sel2* then these selections
             can also each be a list of selection strings (to generate an
             AtomGroup with defined atom order as described under
             :ref:`ordered-selections-label`).
          *path_select*
             atom selection composing coordinates of (fitted) path; if ``None``
             then *path_select* is set to *ref_select* [``None``]

        .. _ClustalW: http://www.clustal.org/
        .. _STAMP: http://www.compbio.dundee.ac.uk/manuals/stamp.4.2/

        """
        self.u_original = universe
        self.u_reference = reference
        self.ref_select = ref_select
        self.ref_frame = ref_frame
        self.path_select = path_select

        self.top_name = self.u_original.filename
        self.trj_name = self.u_original.trajectory.filename
        self.newtrj_name = None
        self.u_fitted = None
        self.path = None
        self.natoms = None


    def fit_to_reference(self, filename=None, prefix='', postfix='_fit',
                         rmsdfile=None, targetdir=os.path.curdir,
                         mass_weighted=False, tol_mass=0.1):
        """Align each trajectory frame to the reference structure with
        :func:`MDAnalysis.analysis.align.rms_fit_trj`.

        :Arguments:
          *filename*
             file name for the RMS-fitted trajectory or pdb; defaults to the
             original trajectory filename (from :attr:`Path.u_original`) with
             *prefix* prepended
          *prefix*
             prefix for auto-generating the new output filename
          *rmsdfile*
             file name for writing the RMSD time series [``None``]
          *mass_weighted*
             do a mass-weighted RMSD fit
          *tol_mass*
             Reject match if the atomic masses for matched atoms differ by more
             than *tol_mass* [0.1]

        :Returns:
          :class:`MDAnalysis.Universe` object containing a fitted trajectory
        """
        head, tail = os.path.split(self.trj_name)
        oldname, ext = os.path.splitext(tail)
        filename = filename or oldname
        self.newtrj_name = os.path.join(targetdir, filename + postfix + ext)
        self.u_reference.trajectory[self.ref_frame] # select frame from ref traj
        MDAnalysis.analysis.align.rms_fit_trj(self.u_original, self.u_reference,\
                select=self.ref_select, filename=self.newtrj_name,              \
                rmsdfile=rmsdfile, prefix=prefix, mass_weighted=mass_weighted,  \
                tol_mass=tol_mass)
        return MDAnalysis.Universe(self.top_name, self.newtrj_name)


    def to_path(self, fitted=False, select=None, flat=False):
        r"""Generates a coordinate time series from the fitted universe
        trajectory.

        Given a selection of *N* atoms from *select*, the atomic positions for
        each frame in the fitted universe (:attr:`Path.u_fitted`) trajectory
        (with |Np| total frames) are appended sequentially to form a 3D or 2D
        (if *flat* is ``True``) :class:`numpy.ndarray` representation of the
        fitted trajectory (with dimensions |3D| or |2D|, respectively).

        :Arguments:
          *fitted*
             construct a :attr:`Path.path` from the :attr:`Path.u_fitted`
             trajectory; if ``False`` then :attr:`Path.path` is generated with
             the trajectory from :attr:`Path.u_original` [``False``]
          *select*
             the selection for constructing the coordinates of each frame in
             :attr:`Path.path`; if ``None`` then :attr:`Path.path_select`
             is used, else it is overridden by *select* [``None``]
          *flat*
             represent :attr:`Path.path` as a 2D (|2D|) :class:`numpy.ndarray`;
             if ``False`` then :attr:`Path.path` is a 3D (|3D|)
             :class:`numpy.ndarray` [``False``]

        :Returns:
          :class:`numpy.ndarray` representing a time series of atomic positions
          of an :class:`MDAnalysis.core.AtomGroup.AtomGroup` selection from
          :attr:`Path.u_fitted.trajectory`
        """
        select = select if select is not None else self.path_select
        if fitted:
            if not isinstance(self.u_fitted, MDAnalysis.Universe):
                raise TypeError("Fitted universe not found. Generate a fitted " +
                        "universe with fit_to_reference() first, or explicitly "+
                        "set argument \"fitted\" to \"False\" to generate a "   +
                        "path from the original universe.")
            u = self.u_fitted
        else:
            u = self.u_original
        frames = u.trajectory
        atoms = u.select_atoms(select)
        self.natoms = len(atoms)
        frames.rewind()
        if flat:
            return np.array([atoms.positions.flatten() for _ in frames])
        else:
            return np.array([atoms.positions for _ in frames])


    def run(self, align=False, filename=None, postfix='_fit', rmsdfile=None,
            targetdir=os.path.curdir, mass_weighted=False, tol_mass=0.1,
            flat=False):
        r"""Generate a path from a trajectory and reference structure, aligning
        to a reference structure if specified.

        This is a convenience method to generate a fitted trajectory from an
        inputted universe (:attr:`Path.u_original`) and reference structure
        (:attr:`Path.u_reference`). :meth:`Path.fit_to_reference` and
        :meth:`Path.to_path` are used consecutively to generate a new universe
        (:attr:`Path.u_fitted`) containing the fitted trajectory along with the
        corresponding :attr:`Path.path` represented as an
        :class:`numpy.ndarray`. The method returns a tuple of the topology name
        and new trajectory name, which can be fed directly into an
        :class:`MDAnalysis.Universe` object after unpacking the tuple using the
        ``*`` operator, as in
        ``MDAnalysis.Universe(*(top_name, newtraj_name))``.

        :Arguments:
          *align*
             Align trajectory to atom selection :attr:`Path.ref_select` of
             :attr:`Path.u_reference`. If ``True``, a universe containing an
             aligned trajectory is produced with :meth:`Path.fit_to_reference`
             [``False``]
          *filename*
             filename for the RMS-fitted trajectory or pdb; defaults to the
             original trajectory filename (from :attr:`Path.u_original`) with
             *prefix* prepended
          *prefix*
             prefix for auto-generating the new output filename
          *rmsdfile*
             file name for writing the RMSD time series [``None``]
          *mass_weighted*
             do a mass-weighted RMSD fit
          *tol_mass*
             Reject match if the atomic masses for matched atoms differ by more
             than *tol_mass* [0.1]
          *flat*
             represent :attr:`Path.path` with 2D (|2D|) :class:`numpy.ndarray`;
             if ``False`` then :attr:`Path.path` is a 3D (|3D|)
             :class:`numpy.ndarray` [``False``]

        :Returns:
          A tuple of the topology name and new trajectory name.
        """
        if align:
            self.u_fitted = self.fit_to_reference(                              \
                                filename=filename, postfix=postfix,             \
                                rmsdfile=rmsdfile, targetdir=targetdir,         \
                                mass_weighted=False, tol_mass=0.1)
        self.path = self.to_path(fitted=align, flat=flat)
        return self.top_name, self.newtrj_name


    def get_num_atoms(self):
        """Return the number of atoms used to construct the :class:`Path`.

        Must run :method:`Path.to_path()` prior to calling this method.

        :Returns:
          int, the number of atoms in the :class:`Path`
        """
        if self.natoms is None:
            err_str = "No path data; do 'Path.to_path()' first."
            raise ValueError(err_str)
        return self.natoms


class PSAPair(object):
    """Generate nearest neighbor and Hausdorff pair information between a pair
    of paths from an all-pairs comparison generated by :class:`PSA`.

    The nearest neighbors for each path of a pair of paths is generated by
    :meth:`PSAPair.compute_nearest_neighbors` and stores the result
    in a dictionary (:attr:`nearest_neighbors`): each path has a
    :class:`numpy.ndarray` of the frames of its nearest neighbors, and a
    :class:`numpy.ndarray` of its nearest neighbor distances
    :attr:`PSAnalysis.D`. For example, *nearest_neighbors['frames']* is a pair
    of :class:`numpy.ndarray`, the first being the frames of the nearest
    neighbors of the first path, *i*, the second being those of the second path,
    *j*.

    The Hausdorff pair for the pair of paths is found by calling
    :meth:`find_hausdorff_pair` (locates the nearest neighbor pair having the
    largest overall distance separating them), which stores the result in a
    dictionary (:attr:`hausdorff_pair`) containing the frames (indices) of the
    pair along with the corresponding (Hausdorff) distance.
    *hausdorff_pair['frame']* contains a pair of frames in the first path, *i*,
    and the second path, *j*, respectively, that correspond to the Hausdorff
    distance between them.

    .. versionadded:: 0.11
    """

    def __init__(self, npaths, i, j):
        """Set up a :class:`PSAPair` for a pair of paths that are part of a
        :class:`PSA` comparison of *npaths* total paths.

        Each unique pair of paths compared using :class:`PSA` is related by
        their nearest neighbors (and corresponding distances) and the Hausdorff
        pair and distance. :class:`PSAPair` is a convenience class for
        calculating and encapsulating nearest neighbor and Hausdorff pair
        information for one pair of paths.

        Given *npaths*, :class:`PSA` performs and all-pairs comparison among all
        paths for a total of :math:`\text{npaths}*(\text{npaths}-1)/2` unique
        comparisons. If distances between paths are computed, the all-pairs
        comparison can be summarized in a symmetric distance matrix whose upper
        triangle can be mapped to a corresponding distance vector form in a
        one-to-one manner. A particular comparison of a pair of paths in a
        given instance of :class:`PSAPair` is thus unique identified by the row
        and column indices in the distance matrix representation (whether or not
        distances are actually computed), or a single ID (index) in the
        corresponding distance vector.

        :Arguments:
          *npaths*
             int, total number of paths in :class:`PSA` used to generate *this*
             :class:`PSAPair`
          *i*
             int, row index (starting at 0) of the distance matrix
          *j*
             int, column index (starting at 0) of the distance matrix
        """
        self.npaths = npaths
        self.matrix_idx = (i,j)
        self.pair_idx = self._dvec_idx(i,j)

        # Set by calling hausdorff_nn
        self.nearest_neighbors = {'frames' : None, 'distances' : None}

        # Set by self.getHausdorffPair
        self.hausdorff_pair = {'frames' : (None, None), 'distance' : None}


    def _dvec_idx(self, i, j):
        """Convert distance matrix indices (in the upper triangle) to the index
        of the corresponding distance vector.

        This is a convenience function to locate distance matrix elements (and
        the pair generating it) in the corresponding distance vector. The row
        index *j* should be greater than *i+1*, corresponding to the upper
        triangle of the distance matrix.

        :Arguments:
          *i*
             int, row index (starting at 0) of the distance matrix
          *j*
             int, column index (starting at 0) of the distance matrix

        :Returns:
          int, (matrix element) index in the corresponding distance vector
        """
        return (self.npaths*i) + j - (i+2)*(i+1)/2


    def compute_nearest_neighbors(self, P,Q, N=None):
        """Generates Hausdorff nearest neighbor lists of *frames* (by index) and
        *distances* for *this* pair of paths corresponding to distance matrix
        indices (*i*,*j*).

        :method:`PSAPair.compute_nearest_neighbors` calls
        :func:`hausdorff_neighbors` to populate the dictionary of the nearest
        neighbor lists of frames (by index) and distances
        (:attr:`PSAPair.nearest_neighbors`). This method must explicitly take as
        arguments a pair of paths, *P* and *Q*, where *P* is the
        :math:`i^\text{th}` path and *Q* is the :math:`j^\text{th}` path among
        the set of *N* total paths in the comparison.

        :Arguments:
          *P*
             :class:`numpy.ndarray` representing a path
          *Q*
             :class:`numpy.ndarray` representing a path
          *N*
             int, size of the distance matrix (of shape *N*-by-*N*) [``None``]
        """
        hn = hausdorff_neighbors(P, Q)
        self.nearest_neighbors['frames'] = hn['frames']
        self.nearest_neighbors['distances'] = hn['distances']


    def find_hausdorff_pair(self):
        r"""Find the Hausdorff pair (of frames) for *this* pair of paths.

        :method:`PSAPair.find_hausdorff_pair` requires that
        `:method:`PSAPair.compute_nearest_neighbors` be called first to
        generate the nearest neighbors (and corresponding distances) for each
        path in *this* :class:`PSAPair`. The Hausdorff pair is the nearest
        neighbor pair (of snapshots/frames), one in the first path and one in
        the second, with the largest separation distance.
        """
        if self.nearest_neighbors['distances'] is None:
            err_str = "Nearest neighbors have not been calculated yet;"         \
                    + " run compute_nearest_neighbors() first."
            raise NoDataError(err_str)

        nn_idx_P, nn_idx_Q = self.nearest_neighbors['frames']
        nn_dist_P, nn_dist_Q = self.nearest_neighbors['distances']
        max_nn_dist_P = max(nn_dist_P)
        max_nn_dist_Q = max(nn_dist_Q)
        if max_nn_dist_P > max_nn_dist_Q:
            max_nn_idx_P = np.argmax(nn_dist_P)
            self.hausdorff_pair['frames'] = max_nn_idx_P, nn_idx_P[max_nn_idx_P]
            self.hausdorff_pair['distance']  = max_nn_dist_P
        else:
            max_nn_idx_Q = np.argmax(nn_dist_Q)
            self.hausdorff_pair['frames'] = nn_idx_Q[max_nn_idx_Q], max_nn_idx_Q
            self.hausdorff_pair['distance'] = max_nn_dist_Q


    def get_nearest_neighbors(self, frames=True, distances=True):
        """Returns the nearest neighbor frame indices, distances, or both, for
        each path in *this* :class:`PSAPair`.

        :method:`PSAPair.get_nearest_neighbors` requires that the nearest
        neighbors (:attr:`nearest_neighbors`) be initially computed by first
        calling :method:`compute_nearest_neighbors`. At least one of *frames*
        or *distances* must be ``True``, or else a ``NoDataError`` is raised.

        :Arguments:
          *frames*
             boolean; if ``True``, return nearest neighbor frame indices
             [``True``]
          *distances*
             boolean; if ``True``, return nearest neighbor distances [``True``]

        :Returns:
          If both *frames* and *distances* are ``True``, return the entire
          dictionary (:attr:`nearest_neighbors`); if only *frames* is
          ``True``, return a pair of :class:`numpy.ndarray` containing the
          indices of the frames (for the pair of paths) of the nearest
          neighbors; if only *distances* is ``True``, return a pair of
          :class:`numpy.ndarray` of the nearest neighbor distances (for the pair
          of paths).
        """
        if self.nearest_neighbors['distances'] is None:
            err_str = "Nearest neighbors have not been calculated yet;"         \
                    + " run compute_nearest_neighbors() first."
            raise NoDataError(err_str)

        if frames:
            if distances:
                return self.nearest_neighbors
            else:
                return self.nearest_neighbors['frames']
        elif distances:
            return self.nearest_neighbors['distances']
        else:
            err_str = "Need to select Hausdorff pair \"frames\" or"             \
                + " \"distances\" or both. \"frames\" and \"distances\" cannot" \
                + " both be set to False."
            raise NoDataError(err_str)


    def get_hausdorff_pair(self, frames=True, distance=True):
        """Returns the Hausdorff pair of frames indices, the Hausdorff distance,
        or both, for the paths in *this* :class:`PSAPair`.

        :method:`PSAPair.get_hausdorff_pair` requires that the Hausdorff pair
        (and distance) be initially found by first calling
        :method:`find_hausdorff_pair`. At least one of *frames* or *distance*
        must be ``True``, or else a ``NoDataError`` is raised.

        :Arguments:
          *frames*
             boolean; if ``True``, return the indices of the frames
             of the Hausdorff pair [``True``]
          *distances*
             boolean; if ``True``, return Hausdorff distance [``True``]

        :Returns:
          If both *frames* and *distance* are ``True``, return the entire
          dictionary (:attr:`hausdorff_pair`); if only *frames* is
          ``True``, return a pair of ``int`` containing the indices of the
          frames (one index per path) of the Hausdorff pair; if only *distance*
          is ``True``, return the Hausdorff distance for this path pair.
        """
        if self.hausdorff_pair['distance'] is None:
            err_str = "Hausdorff pair has not been calculated yet;"             \
                    + " run find_hausdorff_pair() first."
            raise NoDataError(err_str)

        if frames:
            if distance:
                return self.hausdorff_pair
            else:
                return self.hausdorff_pair['frames']
        elif distance:
            return self.hausdorff_pair['distance']
        else:
            err_str = "Need to select Hausdorff pair \"frames\" or"             \
                + " \"distance\" or both. \"frames\" and \"distance\" cannot" \
                + " both be set to False."
            raise NoDataError(err_str)


class PSAnalysis(object):
    """Perform Path Similarity Analysis (PSA) on a set of trajectories.

    The analysis is performed with :meth:`PSAnalysis.run` and stores the result
    in the :class:`numpy.ndarray` distance matrix :attr:`PSAnalysis.D`.
    :meth:`PSAnalysis.run` also generates a fitted trajectory and path from
    alignment of the original trajectories to a reference structure.

    .. versionadded:: 0.8
    """
    def __init__(self, universes, reference=None, ref_select='name CA',
                 ref_frame=0, path_select=None, labels=None,
                 targetdir=os.path.curdir):
        """Setting up Path Similarity Analysis.

        The mutual similarity between all unique pairs of trajectories
        are computed using a selected path metric.

        :Arguments:
          *universes*
             a list of universes (:class:`MDAnalysis.Universe` object), each
             containing a trajectory
          *reference*
             reference coordinates; :class:`MDAnalysis.Universe` object; if
             ``None`` the first time step of the first item in *trajs* is used
             [``None``]
          *ref_select*
             The selection to operate on; can be one of:

             1. any valid selection string for
                :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.select_atoms` that
                produces identical selections in *mobile* and *reference*; or
             2. a dictionary ``{'mobile':sel1, 'reference':sel2}`` (the
                :func:`MDAnalysis.analysis.align.fasta2select` function returns
                such a dictionary based on a ClustalW_ or STAMP_ sequence
                alignment); or
             3. a tuple ``(sel1, sel2)``

             When using 2. or 3. with *sel1* and *sel2* then these selections
             can also each be a list of selection strings (to generate an
             AtomGroup with defined atom order as described under
             :ref:`ordered-selections-label`).
          *mass_weighted*
             do a mass-weighted RMSD fit [``False``]
          *tol_mass*
             Reject match if the atomic masses for matched atoms differ by more
             than *tol_mass* [0.1]
          *ref_frame*
             frame index to select frame from *reference* [0]
          *path_select*
             atom selection composing coordinates of (fitted) path; if ``None``
             then *path_select* is set to *ref_select* [``None``]
          *targetdir*
            output files are saved there [.]
          *labels*
             list of strings, names of trajectories to be analyzed
             (:class:`MDAnalysis.Universe`); if ``None``, defaults to trajectory
             names [``None``]
        """
        self.universes = universes
        self.u_reference = self.universes[0] if reference is None else reference
        self.ref_select = ref_select
        self.ref_frame = ref_frame
        self.path_select = self.ref_select if path_select is None else path_select
        if targetdir is None:
            try:
                targetdir = os.path.join(os.path.curdir, 'psadata')
                os.makedirs(targetdir)
            except OSError:
                if not os.path.isdir(targetdir):
                    raise
        self.targetdir = os.path.realpath(targetdir)

        # Set default directory names for storing topology/reference structures,
        # fitted trajectories, paths, distance matrices, and plots
        self.datadirs = {'fitted_trajs' : '/fitted_trajs',
                         'paths' : '/paths',
                         'distance_matrices' : '/distance_matrices',
                         'plots' : '/plots'}
        for dir_name, directory in six.iteritems(self.datadirs):
            try:
                full_dir_name = os.path.join(self.targetdir, dir_name)
                os.makedirs(full_dir_name)
            except OSError:
                if not os.path.isdir(full_dir_name):
                    raise

        # Keep track of topology, trajectory, and related files
        trj_names = []
        for i, u in enumerate(self.universes):
            head, tail = os.path.split(u.trajectory.filename)
            filename, ext = os.path.splitext(tail)
            trj_names.append(filename)
        self.trj_names = trj_names
        self.fit_trj_names = None
        self.path_names = None
        self.top_name = self.universes[0].filename if len(universes) != 0 else None
        self.labels = labels or self.trj_names

        # Names of persistence (pickle) files where topology and trajectory
        # filenames are stored--should not be modified by user
        self._top_pkl = os.path.join(self.targetdir, "psa_top-name.pkl")
        self._trjs_pkl = os.path.join(self.targetdir, "psa_orig-traj-names.pkl")
        self._fit_trjs_pkl = os.path.join(self.targetdir, "psa_fitted-traj-names.pkl")
        self._paths_pkl = os.path.join(self.targetdir, "psa_path-names.pkl")
        self._labels_pkl = os.path.join(self.targetdir, "psa_labels.pkl")
        # Pickle topology and trajectory filenames for this analysis to curdir
        with open(self._top_pkl, 'wb') as output:
            cPickle.dump(self.top_name, output)
        with open(self._trjs_pkl, 'wb') as output:
            cPickle.dump(self.trj_names, output)
        with open(self._labels_pkl, 'wb') as output:
            cPickle.dump(self.labels, output)

        self.natoms = None
        self.npaths = None
        self.paths = None
        self.D = None   # pairwise distances
        self._HP = None # (distance vector order) list of all Hausdorff pairs
        self._NN = None # (distance vector order) list of all nearest neighbors
        self._psa_pairs = None # (distance vector order) list of all PSAPairs


    def generate_paths(self, **kwargs):
        """Generate paths, aligning each to reference structure if necessary.

        :Keywords:
          *align*
             Align trajectories to atom selection :attr:`PSAnalysis.ref_select`
             of :attr:`PSAnalysis.u_reference` [``False``]
          *filename*
             strings representing base filename for fitted trajectories and
             paths [``None``]
          *infix*
             additional tag string that is inserted into the output filename of
             the fitted trajectory files ['']
          *mass_weighted*
             do a mass-weighted RMSD fit
          *tol_mass*
             Reject match if the atomic masses for matched atoms differ by more
             than *tol_mass*
          *ref_frame*
             frame index to select frame from *reference*
          *flat*
             represent :attr:`Path.path` as a 2D (|2D|) :class:`numpy.ndarray`;
             if ``False`` then :attr:`Path.path` is a 3D (|3D|)
             :class:`numpy.ndarray` [``False``]
          *save*
             boolean; if ``True``, pickle list of names for fitted trajectories
             [``True``]
          *store*
             boolean; if ``True`` then writes each path (:class:`numpy.ndarray`)
             in :attr:`PSAnalysis.paths` to compressed npz (numpy) files
             [``False``]

        The fitted trajectories are written to new files in the
        "/trj_fit" subdirectory in :attr:`PSAnalysis.targetdir` named
        "filename(*trajectory*)XXX*infix*_psa", where "XXX" is a number between
        000 and 999; the extension of each file is the same as its original.
        Optionally, the trajectories can also be saved in numpy compressed npz
        format in the "/paths" subdirectory in :attr:`PSAnalysis.targetdir` for
        persistence and can be accessed as the attribute
        :attr:`PSAnalysis.paths`.
        """
        align = kwargs.pop('align', False)
        filename = kwargs.pop('filename', 'fitted')
        infix = kwargs.pop('infix', '')
        mass_weighted = kwargs.pop('mass_weighted', False)
        tol_mass = kwargs.pop('tol_mass', False)
        ref_frame = kwargs.pop('ref_frame', self.ref_frame)
        flat = kwargs.pop('flat', False)
        save = kwargs.pop('save', True)
        store = kwargs.pop('store', False)

        paths = []
        fit_trj_names = []
        for i, u in enumerate(self.universes):
            p = Path(u, self.u_reference, ref_select=self.ref_select,           \
                     path_select=self.path_select, ref_frame=ref_frame)
            trj_dir = self.targetdir + self.datadirs['fitted_trajs']
            postfix = '{0}{1}{2:03n}'.format(infix, '_psa', i+1)
            top_name, fit_trj_name = p.run(align=align, filename=filename,      \
                                           postfix=postfix,                     \
                                           targetdir=trj_dir,                   \
                                           mass_weighted=mass_weighted,         \
                                           tol_mass=tol_mass, flat=flat)
            paths.append(p.path)
            fit_trj_names.append(fit_trj_name)
        self.natoms, axis = get_coord_axes(paths[0])
        self.paths = paths
        self.npaths = len(paths)
        self.fit_trj_names = fit_trj_names
        if save:
            with open(self._fit_trjs_pkl, 'wb') as output:
                cPickle.dump(self.fit_trj_names, output)
        if store:
            filename = kwargs.pop('filename', None)
            self.save_paths(filename=filename)


    def run(self, **kwargs):
        """Perform path similarity analysis on the trajectories to compute
        the distance matrix.

        A number of parameters can be changed from the defaults. The
        result is stored as the array :attr:`PSAnalysis.D`.

        :Keywords:
          *metric*
             selection string specifying the path metric to measure pairwise
             distances among :attr:`PSAnalysis.paths` [``'hausdorff'``]
          *start*, *stop*, *step*
             start and stop frame index with step size: analyze
             ``trajectory[start:stop:step]`` [``None``]
          *store*
             boolean; if ``True`` then writes :attr:`PSAnalysis.D` to text and
             compressed npz (numpy) files [``True``]
          *filename*
             string, filename to save :attr:`PSAnalysis.D`
        """
        metric = kwargs.pop('metric', 'hausdorff')
        start = kwargs.pop('start', None)
        stop = kwargs.pop('stop', None)
        step = kwargs.pop('step', None)
        store = kwargs.pop('store', True)

        if type(metric) is str:
            metric_func = get_path_metric_func(metric)
        else:
            metric_func = metric
        numpaths = self.npaths
        D = np.zeros((numpaths,numpaths))

        for i in range(0, numpaths-1):
            for j in range(i+1, numpaths):
                P = self.paths[i][start:stop:step]
                Q = self.paths[j][start:stop:step]
                D[i,j] = metric_func(P, Q)
                D[j,i] = D[i,j]
        self.D = D
        if store:
            filename = kwargs.pop('filename', str(metric))
            self.save_result(filename=filename)


    def run_pairs_analysis(self, **kwargs):
        """Perform PSA Hausdorff (nearest neighbor) pairs analysis on all unique
        pairs of paths in :attr:`PSAnalysis.paths`.

        Partial results can be stored in separate lists, where each list is
        indexed according to distance vector convention (i.e., element *(i,j)*
        in distance matrix representation corresponds to element
        :math:`s=N*i+j-(i+1)*(i+2)` in distance vector representation, which is
        the :math:`s^\text{th}` comparison). For each unique pair of paths, the
        nearest neighbors for that pair can be stored in :attr:`NN` and the
        Hausdorff pair in :attr:`HP`. :attr:`PP` stores the full information
        of Hausdorff pairs analysis that is available for each pair of path,
        including nearest neighbors lists and the Hausdorff pairs.

        :Keywords:
          *start*, *stop*, *step*
             start and stop frame index with step size: analyze
             ``trajectory[start:stop:step]`` [``None``]
          *neighbors*
             boolean; if ``True``, then stores dictionary of nearest neighbor
             frames/distances in :attr:`PSAnalysis.NN` [``False``]
          *hausdorff_pairs*
             boolean; if ``True``, then stores dictionary of Hausdorff pair
             frames/distances in :attr:`PSAnalysis.HP` [``False``]
        """
        start = kwargs.pop('start', None)
        stop = kwargs.pop('stop', None)
        step = kwargs.pop('step', None)
        neighbors = kwargs.pop('neighbors', False)
        hausdorff_pairs = kwargs.pop('hausdorff_pairs', False)

        numpaths = self.npaths
        self._NN = [] # list of nearest neighbors pairs
        self._HP = [] # list of Hausdorff pairs
        self._psa_pairs = [] # list of PSAPairs

        for i in range(0, numpaths-1):
            for j in range(i+1, numpaths):
                pp = PSAPair(i, j, numpaths)
                P = self.paths[i][start:stop:step]
                Q = self.paths[j][start:stop:step]
                pp.compute_nearest_neighbors(P, Q, self.natoms)
                pp.find_hausdorff_pair()
                self._psa_pairs.append(pp)
                if neighbors:
                    self._NN.append(pp.get_nearest_neighbors())
                if hausdorff_pairs:
                    self._HP.append(pp.get_hausdorff_pair())


    def save_result(self, filename=None):
        """Save distance matrix :attr:`PSAnalysis.D` to a numpy compressed npz
        file and text file.

        :Arguments:
          *filename*
             string, specifies filename [``None``]

        The data are saved with :func:`numpy.savez_compressed` and
        :func:`numpy.savetxt` in the directory specified by
        :attr:`PSAnalysis.targetdir`.
        """
        filename = filename or 'psa_distances'
        head = self.targetdir + self.datadirs['distance_matrices']
        outfile = os.path.join(head, filename)
        if self.D is None:
            raise NoDataError("Distance matrix has not been calculated yet")
        np.save(outfile + '.npy', self.D)
        np.savetxt(outfile + '.dat', self.D)
        logger.info("Wrote distance matrix to file %r.npz", outfile)
        logger.info("Wrote distance matrix to file %r.dat", outfile)
        return filename


    def save_paths(self, filename=None):
        """Save fitted :attr:`PSAnalysis.paths` to numpy compressed npz files.

        :Arguments:
          *filename*
             string, specifies filename [``None``]

        The data are saved with :func:`numpy.savez_compressed` in the directory
        specified by :attr:`PSAnalysis.targetdir`.
        """
        filename = filename or 'path_psa'
        head = self.targetdir + self.datadirs['paths']
        outfile = os.path.join(head, filename)
        if self.paths is None:
            raise NoDataError("Paths have not been calculated yet")
        path_names = []
        for i, path in enumerate(self.paths):
            current_outfile = "{0}{1:03n}.npy".format(outfile, i+1)
            np.save(current_outfile, self.paths[i])
            path_names.append(current_outfile)
            logger.info("Wrote path to file %r", current_outfile)
        self.path_names = path_names
        with open(self._paths_pkl, 'wb') as output:
            cPickle.dump(self.path_names, output)
        return filename


    def load(self):
        """Load fitted paths specified by 'psa_path-names.pkl' in
        :attr:`PSAnalysis.targetdir`.
        """
        if not os.path.exists(self._paths_pkl):
            raise NoDataError("Fitted trajectories cannot be loaded; save file" +
                              "{0} does not exist.".format(self._paths_pkl))
        self.path_names = np.load(self._paths_pkl)
        self.paths = [np.load(pname) for pname in self.path_names]
        if os.path.exists(self._labels_pkl):
            self.labels = np.load(self._labels_pkl)
        print("Loaded paths from " + self._paths_pkl)


    def plot(self, filename=None, linkage='ward', count_sort=False,
             distance_sort=False, figsize=4.5, labelsize=12):
        """Plot a clustered distance matrix using method *linkage* along with
        the corresponding dendrogram. Rows (and columns) are identified using
        the list of strings specified by :attr:`PSAnalysis.labels`.

        :Arguments:
          *filename*
             string, save figure to *filename* [``None``]
          *linkage*
             string, name of linkage criterion for clustering [``'ward'``]
          *count_sort*
            boolean, see scipy.cluster.hierarchy.dendrogram [``False``]
          *distance_sort*
            boolean, see scipy.cluster.hierarchy.dendrogram [``False``]
          *figsize*
             set the vertical size of plot in inches [``4.5``]
          *labelsize*
             set the font size for colorbar labels; font size for path labels on
             dendrogram default to 3 points smaller [``12``]

        If *filename* is supplied then the figure is also written to file (the
        suffix determines the file type, e.g. pdf, png, eps, ...). All other
        keyword arguments are passed on to :func:`pylab.imshow`.
        """
        from matplotlib.pyplot import figure, colorbar, cm, savefig, clf

        if self.D is None:
            err_str = "No distance data; do 'PSAnalysis.run(store=True)' first."
            raise ValueError(err_str)
        npaths = len(self.D)
        dist_matrix = self.D

        dgram_loc, hmap_loc, cbar_loc = self._get_plot_obj_locs()
        aspect_ratio = 1.25
        clf()
        fig = figure(figsize=(figsize*aspect_ratio, figsize))
        ax_hmap = fig.add_axes(hmap_loc)
        ax_dgram = fig.add_axes(dgram_loc)

        Z, dgram = self.cluster(dist_matrix,                                    \
                                method=linkage,                                 \
                                count_sort=count_sort,                          \
                                distance_sort=distance_sort)
        rowidx = colidx = dgram['leaves'] # get row-wise ordering from clustering
        ax_dgram.invert_yaxis() # Place origin at up left (from low left)

        minDist, maxDist = 0, np.max(dist_matrix)
        dist_matrix_clus = dist_matrix[rowidx,:]
        dist_matrix_clus = dist_matrix_clus[:,colidx]
        im = ax_hmap.matshow(dist_matrix_clus, aspect='auto', origin='lower',   \
                    cmap=cm.YlGn, vmin=minDist, vmax=maxDist)
        ax_hmap.invert_yaxis() # Place origin at upper left (from lower left)
        ax_hmap.locator_params(nbins=npaths)
        ax_hmap.set_xticks(np.arange(npaths), minor=True)
        ax_hmap.set_yticks(np.arange(npaths), minor=True)
        ax_hmap.tick_params(axis='x', which='both', labelleft='off',            \
                        labelright='off', labeltop='on', labelsize=0)
        ax_hmap.tick_params(axis='y', which='both', labelleft='on',             \
                labelright='off', labeltop='off', labelsize=0)
        rowlabels = [self.labels[i] for i in rowidx]
        collabels = [self.labels[i] for i in colidx]
        ax_hmap.set_xticklabels(collabels, rotation='vertical',                 \
                size=(labelsize-4), multialignment='center', minor=True)
        ax_hmap.set_yticklabels(rowlabels, rotation='horizontal',               \
                size=(labelsize-4), multialignment='left', ha='right',          \
                minor=True)

        ax_color = fig.add_axes(cbar_loc)
        colorbar(im, cax=ax_color, ticks=np.linspace(minDist, maxDist, 10),  \
                format="%0.1f")
        ax_color.tick_params(labelsize=labelsize)

        # Remove major ticks from both heat map axes
        for tic in ax_hmap.xaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False
            tic.label1On = tic.label2On = False
        for tic in ax_hmap.yaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False
            tic.label1On = tic.label2On = False
        # Remove minor ticks from both heat map axes
        for tic in ax_hmap.xaxis.get_minor_ticks():
            tic.tick1On = tic.tick2On = False
        for tic in ax_hmap.yaxis.get_minor_ticks():
            tic.tick1On = tic.tick2On = False
        # Remove tickmarks from colorbar
        for tic in ax_color.yaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False

        if filename is not None:
            head = self.targetdir + self.datadirs['plots']
            outfile = os.path.join(head, filename)
            savefig(outfile, dpi=300, bbox_inches='tight')

        return Z, dgram, dist_matrix_clus


    def plot_annotated_heatmap(self, filename=None, linkage='ward',             \
                               count_sort=False, distance_sort=False,           \
                               figsize=8, annot_size=6.5):
        """Plot a clustered distance matrix using method *linkage* with
        annotated distances in the matrix. Rows (and columns) are identified
        using the list of strings specified by :attr:`PSAnalysis.labels`.

        :Arguments:
          *filename*
             string, save figure to *filename* [``None``]
          *count_sort*
            boolean, see scipy.cluster.hierarchy.dendrogram [``False``]
          *distance_sort*
            boolean, see scipy.cluster.hierarchy.dendrogram [``False``]
          *linkage*
             string, name of linkage criterion for clustering [``'ward'``]
          *figsize*
             set the vertical size of plot in inches [``8``]
          *annot_size*
             float, font size of annotation labels on heat map [``6.5``]

        If *filename* is supplied then the figure is also written to file (the
        suffix determines the file type, e.g. pdf, png, eps, ...). All other
        keyword arguments are passed on to :func:`pylab.imshow`.
        """
        from matplotlib.pylab import figure, colorbar, cm, savefig, clf

        try:
            import seaborn.apionly as sns
        except ImportError:
            raise ImportError(
                """ERROR --- The seaborn package cannot be found!

                The seaborn API could not be imported. Please install it first.
                You can try installing with pip directly from the
                internet:

                  pip install seaborn

                Alternatively, download the package from

                  http://pypi.python.org/pypi/seaborn/

                and install in the usual manner.
                """
            )

        if self.D is None:
            err_str = "No distance data; do 'PSAnalysis.run(store=True)' first."
            raise ValueError(err_str)
        dist_matrix = self.D

        Z, dgram = self.cluster(dist_matrix,                                    \
                                method=linkage,                                 \
                                count_sort=count_sort,                          \
                                distance_sort=distance_sort,                    \
                                no_plot=True)
        rowidx = colidx = dgram['leaves'] # get row-wise ordering from clustering
        dist_matrix_clus = dist_matrix[rowidx,:]
        dist_matrix_clus = dist_matrix_clus[:,colidx]

        clf()
        aspect_ratio = 1.25
        fig = figure(figsize=(figsize*aspect_ratio, figsize))
        ax_hmap = fig.add_subplot(111)
        ax_hmap = sns.heatmap(dist_matrix_clus,                                 \
                         linewidths=0.25, cmap=cm.YlGn, annot=True, fmt='3.1f', \
                         square=True, xticklabels=rowidx, yticklabels=colidx,   \
                         annot_kws={"size": 7}, ax=ax_hmap)

        # Remove major ticks from both heat map axes
        for tic in ax_hmap.xaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False
            tic.label1On = tic.label2On = False
        for tic in ax_hmap.yaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False
            tic.label1On = tic.label2On = False
        # Remove minor ticks from both heat map axes
        for tic in ax_hmap.xaxis.get_minor_ticks():
            tic.tick1On = tic.tick2On = False
        for tic in ax_hmap.yaxis.get_minor_ticks():
            tic.tick1On = tic.tick2On = False

        if filename is not None:
            head = self.targetdir + self.datadirs['plots']
            outfile = os.path.join(head, filename)
            savefig(outfile, dpi=600, bbox_inches='tight')

        return Z, dgram, dist_matrix_clus


    def plot_nearest_neighbors(self, filename=None, idx=0,                      \
                               labels=('Path 1', 'Path 2'), figsize=4.5,        \
                               multiplot=False, aspect_ratio=1.75,              \
                               labelsize=12):
        """Plot nearest neighbor distances as a function of normalized frame
        number (mapped to the interval *[0,1]*).

        :Arguments:
          *filename*
             string, save figure to *filename* [``None``]
          *idx*
             integer, index of path (pair) comparison to plot [``0``]
          *labels*
             (string, string), pair of names to label nearest neighbor distance
             curves [``('Path 1', 'Path 2')``]
          *figsize*
             float, set the vertical size of plot in inches [``4.5``]
          *multiplot*
             boolean, set to ``True`` to enable plotting multiple nearest
             neighbor distances on the same figure [``False``]
          *aspect_ratio*
             float, set the ratio of width to height of the plot [``1.75``]
          *labelsize*
             set the font size for colorbar labels; font size for path labels on
             dendrogram default to 3 points smaller [``12``]

        If *filename* is supplied then the figure is also written to file (the
        suffix determines the file type, e.g. pdf, png, eps, ...). All other
        keyword arguments are passed on to :func:`pylab.imshow`.
        """
        from matplotlib.pyplot import figure, savefig, tight_layout, clf, show
        try:
            import seaborn.apionly as sns
        except ImportError:
            raise ImportError(
                """ERROR --- The seaborn package cannot be found!

                The seaborn API could not be imported. Please install it first.
                You can try installing with pip directly from the
                internet:

                  pip install seaborn

                Alternatively, download the package from

                  http://pypi.python.org/pypi/seaborn/

                and install in the usual manner.
                """
            )

        colors = sns.xkcd_palette(["cherry", "windows blue"])

        if self._NN is None:
            err_str = ("No nearest neighbor data; run "
                       "'PSAnalysis.run_nearest_neighbors()' first.")
            raise ValueError(err_str)

        sns.set_style('whitegrid')

        if not multiplot:
            clf()
        fig = figure(figsize=(figsize*aspect_ratio, figsize))
        ax = fig.add_subplot(111)

        nn_dist_P, nn_dist_Q = self._NN[idx]['distances']
        frames_P = len(nn_dist_P)
        frames_Q = len(nn_dist_Q)
        progress_P = np.asarray(range(frames_P))/(1.0*frames_P)
        progress_Q = np.asarray(range(frames_Q))/(1.0*frames_Q)

        ax.plot(progress_P, nn_dist_P, color=colors[0], lw=1.5, label=labels[0])
        ax.plot(progress_Q, nn_dist_Q, color=colors[1], lw=1.5, label=labels[1])

        ax.legend()
        ax.set_xlabel(r'(normalized) progress by frame number', fontsize=12)
        ax.set_ylabel(r'nearest neighbor rmsd ($\AA$)', fontsize=12)
        ax.tick_params(axis='both', which='major', labelsize=12, pad=4)

        sns.despine(bottom=True, left=True, ax=ax)
        tight_layout()

        if filename is not None:
            head = self.targetdir + self.datadirs['plots']
            outfile = os.path.join(head, filename)
            savefig(outfile, dpi=300, bbox_inches='tight')
        show()


    def cluster(self, distArray, method='ward', count_sort=False,               \
                distance_sort=False, no_plot=False, no_labels=True,             \
                color_threshold=4):
        """Cluster trajectories and optionally plot the dendrogram.

        :Arguments:
          *method*
            string, name of linkage criterion for clustering [``'ward'``]
          *no_plot*
            boolean, if ``True``, do not render the dendrogram [``False``]
          *no_labels*
            boolean, if ``True`` then do not label dendrogram [``True``]
          *color_threshold*
            For brevity, let t be the color_threshold. Colors all the
            descendent links below a cluster node k the same color if k is
            the first node below the cut threshold t. All links connecting
            nodes with distances greater than or equal to the threshold are
            colored blue. If t is less than or equal to zero, all nodes are
            colored blue. If color_threshold is None or ‘default’,
            corresponding with MATLAB(TM) behavior, the threshold is set to
            0.7*max(Z[:,2]). [``4``]]
        :Returns:
          list of indices representing the row-wise order of the objects
          after clustering
        """
        import matplotlib
        from scipy.cluster.hierarchy import linkage, dendrogram
        from brewer2mpl import get_map

        color_list = get_map('Set1', 'qualitative', 9).mpl_colors
        matplotlib.rcParams['lines.linewidth'] = 0.5

        Z = linkage(distArray, method=method)
        dgram = dendrogram(Z, no_labels=no_labels, orientation='right',         \
                           count_sort=count_sort, distance_sort=distance_sort,  \
                           no_plot=no_plot, color_threshold=color_threshold,    \
                           color_list=color_list)
        return Z, dgram


    def _get_plot_obj_locs(self):
        """Find and return coordinates for dendrogram, heat map, and colorbar.

        :Returns:
          tuple of coordinates for placing the dendrogram, heat map, and
          colorbar in the plot.
        """
        plot_xstart = 0.04
        plot_ystart = 0.04
        label_margin = 0.155

        dgram_height = 0.2 # dendrogram heights(s)
        hmap_xstart = plot_xstart + dgram_height + label_margin

        # Set locations for dendrogram(s), matrix, and colorbar
        hmap_height = 0.8
        hmap_width = 0.6
        dgram_loc = [plot_xstart, plot_ystart, dgram_height, hmap_height]
        cbar_width = 0.02
        cbar_xstart = hmap_xstart + hmap_width + 0.01
        cbar_loc = [cbar_xstart, plot_ystart, cbar_width, hmap_height]
        hmap_loc =  [hmap_xstart, plot_ystart, hmap_width, hmap_height]

        return dgram_loc, hmap_loc, cbar_loc


    def get_num_atoms(self):
        """Return the number of atoms used to construct the :class:`Path`s in
        :class:`PSA`.

        .. note::

        Must run :method:`PSAnalysis.generate_paths()` prior to calling this
        method.

        :Returns:
          int, the number of atoms in :class:`PSA`'s :class:`Path`s'
        """
        if self.natoms is None:
            err_str = "No path data; do 'PSAnalysis.generate_paths()' first."
            raise ValueError(err_str)
        return self.natoms


    def get_num_paths(self):
        """Return the number of paths in :class:`PSA`.

        .. note::

        Must run :method:`PSAnalysis.generate_paths()` prior to calling this
        method.

        :Returns:
          int, the number of paths in :class:`PSA`
        """
        if self.npaths is None:
            err_str = "No path data; do 'PSAnalysis.generate_paths()' first."
            raise ValueError(err_str)
        return self.npaths


    def get_paths(self):
        """Return the paths in :class:`PSA`.

        .. note::

        Must run :method:`PSAnalysis.generate_paths()` prior to calling this
        method.

        :Returns:
          list of :class:`numpy.ndarray` representations of paths in
          :class:`PSA`
        """
        if self.paths is None:
            err_str = "No path data; do 'PSAnalysis.generate_paths()' first."
            raise ValueError(err_str)
        return self.paths


    def get_pairwise_distances(self, vectorform=False):
        """Return the distance matrix (or vector) of pairwise path distances.

        .. note::

        Must run :method:`PSAnalysis.run(store=True)` prior to calling this
        method.

        :Arguments:
          *vectorform*
            boolean, if ``True``, return the distance vector instead [``False``]

        :Returns:
          :class:`numpy.ndarray` representation of the distance matrix (or
          vector)
        """
        if self.D is None:
            err_str = "No distance data; do 'PSAnalysis.run(store=True)' first."
            raise ValueError(err_str)
        if vectorform:
            from scipy.spatial.distance import squareform
            return squareform(self.D)
        else:
            return self.D


    @property
    def psa_pairs(self):
        """Get the list of :class:`PSAPair`s for each pair of paths.

        :method:`psa_pairs` is a list of :class:`PSAPair` whose elements are
        pairs of paths that have been compared using
        :method:`PSAnalysis.run_pairs_analysis()`. Each :class:`PSAPair`
        contains nearest neighbor and Hausdorff pair information specific to a
        pair of paths. The nearest neighbor frames and distances for a
        :class:`PSAPair` can be accessed in the nearest neighbor dictionary
        using the keys 'frames' and 'distances', respectively. E.g.,
        :attr:`PSAPair.nearest_neighbors['distances']` returns a *pair* of
        :class:`numpy.ndarray` corresponding to the nearest neighbor distances
        for each path. Similarly, Hausdorff pair information can be accessed
        using :attr:`PSAPair.hausdorff_pair` with the keys 'frames' and
        'distance'.

        .. note::

        Must run :method:`PSAnalysis.run_pairs_analysis()` prior to calling this
        method.

        :Returns:
          list of all :class:`PSAPair` objects (in distance vector order)
        """
        if self._psa_pairs is None:
            err_str = "No nearest neighbors data; do"                           \
                    + " 'PSAnalysis.run_pairs_analysis()' first."
            raise ValueError(err_str)
        return self._psa_pairs


    @property
    def hausdorff_pairs(self):
        """Get the Hausdorff pair for each (unique) pairs of paths.

        This method returns a list of Hausdorff pair information, where each
        element is a dictionary containing the pair of frames and the
        (Hausdorff) distance between a pair of paths. See
        :method:`PSAnalysis.psa_pairs` and :attr:`PSAPair.hausdorff_pair` for
        more information about accessing Hausdorff pair data.

        .. note::

        Must run :method:`PSAnalysis.run_pairs_analysis(hausdorff_pairs=True)`
        prior to calling this method.

        :Returns:
          list of all Hausdorff pairs (in distance vector order)
        """
        if self._HP is None:
            err_str = "No Hausdorff pairs data; do "                            \
                    + "'PSAnalysis.run_pairs_analysis(hausdorff_pairs=True)' "  \
                    + "first."
            raise ValueError(err_str)
        return self._HP

    @property
    def nearest_neighbors(self):
        """Get the nearest neighbors for each (unique) pair of paths.

        This method returns a list of nearest neighbor information, where each
        element is a dictionary containing the nearest neighbor frames and
        distances between a pair of paths. See :method:`PSAnalysis.psa_pairs`
        and :attr:`PSAPair.nearest_neighbors` for more information about
        accessing nearest neighbor data.

        .. note::

        Must run :method:`PSAnalysis.run_pairs_analysis(neighbors=True)` prior
        to calling this method.

        :Returns:
          list of all nearest neighbors (in distance vector order)
        """
        if self._NN is None:
            err_str = "No nearest neighbors data; do"                           \
                    + " 'PSAnalysis.run_pairs_analysis(neighbors=True)' first."
            raise ValueError(err_str)
        return self._NN
