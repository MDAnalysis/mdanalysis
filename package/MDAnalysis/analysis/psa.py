# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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

.. [Seyler2015] Sean L. Seyler, Avishek Kumar, Michael F. Thorpe, Oliver Beckstein.
   *Path Similarity Analysis: a Method for Quantifying Macromolecular
   Pathways.* `arXiv:1505.04807`_ (2015).

.. _`arXiv:1505.04807`: http://arxiv.org/abs/1505.04807
.. _`PSAnalysisTutorial`: https://github.com/Becksteinlab/PSAnalysisTutorial


Helper functions and variables
------------------------------

The following global variables are used by the functions and classes in this
module.

.. data::   hausdorff_names
.. data::   frechet_names

The following functions are used by the other functions in this module.

.. autofunction:: sqnorm


Classes, methods, and functions
-------------------------------

.. autofunction:: get_path_metric_func
.. autofunction:: hausdorff
.. autofunction:: discrete_frechet

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


.. autoclass:: PSA
   :members:

   .. attribute:: universes

      list of :class:`MDAnalysis.Universe` objects containing trajectories

   .. attribute:: u_reference

      :class:`MDAnalysis.Universe` object containing a reference structure

   .. attribute:: ref_select

      string, selection for
      :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.select_atoms` to select frame
      from :attr:`PSA.u_reference`

   .. attribute:: path_select

      string, selection for
      :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.select_atoms` to select atoms
      to compose :attr:`Path.path`

   .. attribute:: ref_frame

      int, frame index to select frame from :attr:`Path.u_reference`

   .. attribute:: filename

      string, name of file to store calculated distance matrix (:attr:`PSA.D`)

   .. attribute:: paths

      list of :class:`numpy.ndarray` objects representing the set/ensemble of
      fitted trajectories

   .. attribute:: D

      string, name of file to store calculated distance matrix (:attr:`PSA.D`)


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

import numpy as np

import MDAnalysis
import MDAnalysis.analysis.align
from MDAnalysis import NoDataError

import cPickle as pickle
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
            'discrete_frechet' : discrete_frechet
    }
    try:
        return path_metrics[name]
    except KeyError as key:
        print("Path metric {} not found. Valid selections: ".format(key))
        for name in path_metrics.keys(): print("  \"{}\"".format(name))


def hausdorff(P,Q, N=None):
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
      *N*
         int, number of atoms; if ``None`` then *P* and *Q* are each assumed to
         be a 3D :class:`numpy.ndarray`; else, they are 2D (flattened)

    :Returns:
      float, the Hausdorff distance between paths *P* and *Q*

    Example::
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
    if N is None:
        N = P.shape[1] # number of atoms from 2nd dim of P
        axis = (1,2)
    else:
        axis = 1
    d = np.array([sqnorm(p - Q, axis) for p in P])
    return ( max( np.amax(np.amin(d, axis=0)),             \
                  np.amax(np.amin(d, axis=1)) ) / N  )**0.5


def discrete_frechet(P,Q, N=None):
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
      *N*
         int, number of atoms; if ``None`` then *P* and *Q* are each assumed to
         be a 3D :class:`numpy.ndarray`; else, they are 2D (flattened)

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
    if N is None:
        N = P.shape[1] # number of atoms from 2nd dim of P
        axis = (1,2)
    else:
        axis = 1
    Np, Nq = len(P), len(Q)
    d = np.array([sqnorm(p - Q, axis) for p in P])
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


def sqnorm(v, axis=None):
    """Compute the sum of squares of elements along specified axes.

    :Arguments:
      *v*
         :class:`numpy.ndarray` of coordinates
      *axis*
         None or int or tuple of ints, optional
         Axis or axes along which a sum is performed. The default (axis = None)
         is perform a sum over all the dimensions of the input array. axis may
         be negative, in which case it counts from the last to the first axis.

    :Returns:
      float, the sum of the squares of the elements of *v* along *axis*
    """
    return np.sum(v*v, axis=axis)


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

    The analysis is performed with :meth:`PSA.run` and stores the result
    in the :class:`numpy.ndarray` distance matrix :attr:`PSA.D`. :meth:`PSA.run`
    also generates a fitted trajectory and path from alignment of the original
    trajectories to a reference structure.

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
        MDAnalysis.analysis.align.rms_fit_trj(self.u_original, self.u_reference,
                select=self.ref_select, filename=self.newtrj_name,
                rmsdfile=rmsdfile, prefix=prefix, mass_weighted=mass_weighted,
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
             represent :attr:`Path.path` as a 2D (|2D|) :class:`numpy.ndarray`; if
             ``False`` then :attr:`Path.path` is a 3D (|3D|)
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
             represent :attr:`Path.path` with a 2D (|2D|) :class:`numpy.ndarray`;
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


class PSA(object):
    """Perform Path Similarity Analysis (PSA) on a set of trajectories.

    The analysis is performed with :meth:`PSA.run` and stores the result
    in the :class:`numpy.ndarray` distance matrix :attr:`PSA.D`. :meth:`PSA.run`
    also generates a fitted trajectory and path from alignment of the original
    trajectories to a reference structure.

    .. versionadded:: 0.8
    """
    def __init__(self, universes, reference=None, ref_select='name CA',
                 ref_frame=0, path_select=None, labels=None,
                 targetdir=None):
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
        for dir_name, directory in self.datadirs.iteritems():
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
            pickle.dump(self.top_name, output)
        with open(self._trjs_pkl, 'wb') as output:
            pickle.dump(self.trj_names, output)
        with open(self._labels_pkl, 'wb') as output:
            pickle.dump(self.labels, output)

        self.paths = None
        self.D = None


    def generate_paths(self, **kwargs):
        r"""Generate paths, aligning each to reference structure if necessary.

        :Keywords:
          *align*
             Align trajectories to atom selection :attr:`PSA.ref_select` of
             :attr:`PSA.u_reference` [``False``]
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
             represent :attr:`Path.path` as a 2D (|2D|) :class:`numpy.ndarray`; if
             ``False`` then :attr:`Path.path` is a 3D (|3D|)
             :class:`numpy.ndarray` [``False``]
          *save*
             boolean; if ``True``, pickle list of names for fitted trajectories
             [``True``]
          *store*
             boolean; if ``True`` then writes each path (:class:`numpy.ndarray`)
             in :attr:`PSA.paths` to compressed npz (numpy) files [``False``]

        The fitted trajectories are written to new files in the
        "/trj_fit" subdirectory in :attr:`PSA.targetdir` named
        "filename(*trajectory*)XXX*infix*_psa", where "XXX" is a number between
        000 and 999; the extension of each file is the same as its original.
        Optionally, the trajectories can also be saved in numpy compressed npz
        format in the "/paths" subdirectory in :attr:`PSA.targetdir` for
        persistence and can be accessed as the attribute :attr:`PSA.paths`.
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
            postfix = '{}{}{:03n}'.format(infix, '_psa', i+1)
            top_name, fit_trj_name = p.run(align=align, filename=filename,      \
                                           postfix=postfix,                     \
                                           targetdir=trj_dir,                   \
                                           mass_weighted=mass_weighted,         \
                                           tol_mass=tol_mass, flat=flat)
            paths.append(p.path)
            fit_trj_names.append(fit_trj_name)
        self.paths = paths
        self.fit_trj_names = fit_trj_names
        if save:
            with open(self._fit_trjs_pkl, 'wb') as output:
                pickle.dump(self.fit_trj_names, output)
        if store:
            filename = kwargs.pop('filename', None)
            self.save_paths(filename=filename)


    def run(self, **kwargs):
        """Perform path similarity analysis on the trajectories to compute
        the distance matrix.

        A number of parameters can be changed from the defaults. The
        result is stored as the array :attr:`PSA.D`.

        :Keywords:
          *metric*
             selection string specifying the path metric to measure pairwise
             distances among :attr:`PSA.paths` [``'hausdorff'``]
          *start*, *stop*, *step*
             start and stop frame index with step size: analyze
             ``trajectory[start:stop:step]`` [``None``]
          *store*
             boolean; if ``True`` then writes :attr:`PSA.D` to text and
             compressed npz (numpy) files [``True``]
          *filename*
             string, filename to save :attr:`PSA.D`
        """
        metric = kwargs.pop('metric', 'hausdorff')
        start = kwargs.pop('start', None)
        stop = kwargs.pop('stop', None)
        step = kwargs.pop('step', None)
        store = kwargs.pop('store', True)

        metric_func = get_path_metric_func(metric)
        npaths = len(self.paths)
        D = np.zeros((npaths,npaths))
        for i in xrange(0, npaths-1):
            for j in xrange(i+1, npaths):
                P = self.paths[i][start:stop:step]
                Q = self.paths[j][start:stop:step]
                D[i,j] = D[j,i] = metric_func(P, Q)
        self.D = D
        if store:
            filename = kwargs.pop('filename', metric)
            self.save_result(filename=filename)


    def save_result(self, filename=None):
        """Save distance matrix :attr:`PSA.D` to a numpy compressed npz file and
        text file.

        :Arguments:
          *filename*
             string, specifies filename [``None``]

        The data are saved with :func:`numpy.savez_compressed` and
        :func:`numpy.savetxt` in the directory specified by
        :attr:`PSA.targetdir`.
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
        """Save fitted :attr:`PSA.paths` to numpy compressed npz files.

        :Arguments:
          *filename*
             string, specifies filename [``None``]

        The data are saved with :func:`numpy.savez_compressed` in the directory
        specified by :attr:`PSA.targetdir`.
        """
        filename = filename or 'path_psa'
        head = self.targetdir + self.datadirs['paths']
        outfile = os.path.join(head, filename)
        if self.paths is None:
            raise NoDataError("Paths have not been calculated yet")
        path_names = []
        for i, path in enumerate(self.paths):
            current_outfile = "{}{:03n}.npy".format(outfile, i+1)
            np.save(current_outfile, self.paths[i])
            path_names.append(current_outfile)
            logger.info("Wrote path to file %r", current_outfile)
        self.path_names = path_names
        with open(self._paths_pkl, 'wb') as output:
            pickle.dump(self.path_names, output)
        return filename


    def load(self):
        """Load fitted paths specified by 'psa_path-names.pkl' in
        :attr:`PSA.targetdir`.
        """
        if not os.path.exists(self._paths_pkl):
            raise NoDataError("Fitted trajectories cannot be loaded; save file" +
                              "{} does not exist.".format(self._paths_pkl))
        self.path_names = np.load(self._paths_pkl)
        self.paths = [np.load(pname) for pname in self.path_names]
        if os.path.exists(self._labels_pkl):
            self.labels = np.load(self._labels_pkl)
        print("Loaded paths from " + self._paths_pkl)


    def plot(self, filename=None, linkage='ward', count_sort=False,
             distance_sort=False, figsize=4.5, labelsize=12):
        """Plot a clustered distance matrix using method *linkage* along with
        the corresponding dendrogram. Rows (and columns) are identified using
        the list of strings specified by :attr:`PSA.labels`.

        :Arguments:
          *filename*
             string, save figure to *filename* [``None``]
          *linkage*
             string, name of linkage criterion for clustering [``'ward'``]
          *figsize*
             set the vertical size of plot in inches [``6``]
          *labelsize*
             set the font size for colorbar labels; font size for path labels on
             dendrogram default to 3 points smaller [``12``]

        If *filename* is supplied then the figure is also written to file (the
        suffix determines the file type, e.g. pdf, png, eps, ...). All other
        keyword arguments are passed on to :func:`pylab.imshow`.
        """
        from pylab import cm, savefig
        from matplotlib.pylab import figure, colorbar
        from scipy.spatial.distance import squareform

        if self.D is None:
            raise ValueError("No distance data; do 'PSA.run(store=True)' first.")
        npaths = len(self.D)
        distMatrix = self.D

        dgram_loc, hmap_loc, cbar_loc = self._get_plot_obj_locs()
        aspect_ratio = 1.2
        fig = figure(figsize=(figsize*aspect_ratio,figsize))
        ax_hmap = fig.add_axes(hmap_loc)
        ax_dgram = fig.add_axes(dgram_loc)

        Z, dgram = self.cluster(distMatrix, method=linkage,
                            count_sort=count_sort, distance_sort=distance_sort)
        rowidx = colidx = dgram['leaves'] # get row-wise ordering from clustering
        ax_dgram.invert_yaxis() # Place origin at up left (from low left)

        minDist, maxDist = 0, np.max(distMatrix)
        clustMatrix = distMatrix[rowidx,:]
        clustMatrix = clustMatrix[:,colidx]
        im = ax_hmap.matshow(clustMatrix, aspect='auto', origin='lower',        \
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
                size=(labelsize-4), multialignment='left', ha='right',         \
                minor=True)

        ax_color = fig.add_axes(cbar_loc)
        colorbar(im, cax=ax_color, ticks=np.linspace(minDist, maxDist, 10),  \
                format="%0.2f")
        ax_color.tick_params(labelsize=labelsize)

        # Remove tickmarks from border of heatmap and colorbar
        for tic in ax_hmap.xaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False
            tic.label1On = tic.label2On = False
        for tic in ax_hmap.xaxis.get_minor_ticks():
            tic.tick1On = tic.tick2On = False
        for tic in ax_hmap.yaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False
            tic.label1On = tic.label2On = False
        for tic in ax_hmap.yaxis.get_minor_ticks():
            tic.tick1On = tic.tick2On = False
        for tic in ax_color.yaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False

        if filename is not None:
            head = self.targetdir + self.datadirs['plots']
            outfile = os.path.join(head, filename)
            savefig(outfile, dpi=300, bbox_inches='tight')

        return Z, dgram


    def cluster(self, distArray, method='ward', count_sort=False,
                distance_sort=False, no_labels=True, color_threshold=4):
        """
            Cluster trajectories.

            :Arguments:
              *method*
                string, name of linkage criterion for clustering [``'ward'``]
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
        dgram = dendrogram(Z, no_labels=no_labels, orientation='right',
                           count_sort=count_sort, distance_sort=distance_sort,
                           color_threshold=color_threshold,
                           color_list=color_list)
        return Z, dgram


    def _get_plot_obj_locs(self):
        """
            Return the coordinates for dendrogram, heat map, and colorbar.

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
