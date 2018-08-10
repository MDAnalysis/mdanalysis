# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
#

"""Fast distance array computation --- :mod:`MDAnalysis.lib.distances`
===================================================================

Fast C-routines to calculate distance arrays from coordinate
arrays. Many of the functions also exist in parallel versions, that
typically provide higher performance than the serial code.
The boolean attribute MDAnalysis.lib.distances.USED_OPENMP can be
checked to see if OpenMP was used in the compilation of MDAnalysis.

Selection of acceleration ("backend")
-------------------------------------

All functions take the optional keyword *backend*, which determines
the type of acceleration. Currently, the following choices are
implemented (*backend* is case-insensitive):

.. Table:: Available *backends* for accelerated distance functions.

   ========== ========================= ======================================
   *backend*  module                    description
   ========== ========================= ======================================
   "serial"   :mod:`c_distances`        serial implementation in C/Cython

   "OpenMP"   :mod:`c_distances_openmp` parallel implementation in C/Cython
                                        with OpenMP
   ========== ========================= ======================================

.. versionadded:: 0.13.0

Functions
---------

.. autofunction:: distance_array(reference, configuration [, box [, result [, backend]]])
.. autofunction:: self_distance_array(reference [, box [,result [, backend]]])
.. autofunction:: calc_bonds(atom1, atom2 [, box, [, result [, backend]]])
.. autofunction:: calc_angles(atom1, atom2, atom3 [,box [, result [, backend]]])
.. autofunction:: calc_dihedrals(atom1, atom2, atom3, atom4 [,box [, result [, backend]]])
.. autofunction:: apply_PBC(coordinates, box [, backend])
.. autofunction:: capped_distance(reference, configuration, max_cutoff [, min_cutoff [, box [, method]]])
.. autofunction:: self_capped_distance(reference, max_cutoff, [, min_cutoff [, box [, method]]])
.. autofunction:: transform_RtoS(coordinates, box [, backend])
.. autofunction:: transform_StoR(coordinates, box [,backend])
.. autofunction:: MDAnalysis.lib._augment.augment_coordinates(coordinates, box, radius)
.. autofunction:: MDAnalysis.lib._augment.undo_augment(indices, translation, nreal)

"""
from __future__ import division, absolute_import
from six.moves import range

import numpy as np
from numpy.lib.utils import deprecate

from .mdamath import triclinic_vectors, triclinic_box
from ._augment import augment_coordinates, undo_augment



# hack to select backend with backend=<backend> kwarg. Note that
# the cython parallel code (prange) in parallel.distances is
# independent from the OpenMP code
import importlib
_distances = {}
_distances['serial'] = importlib.import_module(".c_distances",
                                         package="MDAnalysis.lib")
try:
    _distances['openmp'] = importlib.import_module(".c_distances_openmp",
                                          package="MDAnalysis.lib")
except ImportError:
    pass
del importlib

def _run(funcname, args=None, kwargs=None, backend="serial"):
    """Helper function to select a backend function *funcname*."""
    args = args if args is not None else tuple()
    kwargs = kwargs if kwargs is not None else dict()
    backend = backend.lower()
    try:
        func = getattr(_distances[backend], funcname)
    except KeyError:
        raise ValueError("Function {0} not available with backend {1}; try one of: {2}".format(
            funcname, backend, ", ".join(_distances.keys())))
    return func(*args, **kwargs)

# serial versions are always available (and are typically used within
# the core and topology modules)
from .c_distances import (calc_distance_array,
                          calc_distance_array_ortho,
                          calc_distance_array_triclinic,
                          calc_self_distance_array,
                          calc_self_distance_array_ortho,
                          calc_self_distance_array_triclinic,
                          coord_transform,
                          calc_bond_distance,
                          calc_bond_distance_ortho,
                          calc_bond_distance_triclinic,
                          calc_angle,
                          calc_angle_ortho,
                          calc_angle_triclinic,
                          calc_dihedral,
                          calc_dihedral_ortho,
                          calc_dihedral_triclinic,
                          ortho_pbc,
                          triclinic_pbc)

from .c_distances_openmp import OPENMP_ENABLED as USED_OPENMP


def _box_check(box):
    """Take a box input and deduce what type of system it represents based
    on the shape of the array and whether all angles are 90.

    Parameters
    ----------
    box : array
        Box information of unknown format.

    Returns
    -------
    boxtype : str
        * ``ortho`` orthogonal box
        * ``tri_vecs`` triclinic box vectors
        * ``tri_box`` triclinic box lengths and angles

    Raises
    ------
    TypeError
        If box is not float32.
    ValueError
        If box type not detected.
    """
    if box.dtype != np.float32:
        raise TypeError("Box must be of type float32")

    boxtype = 'unknown'

    if box.shape == (3,):
        boxtype = 'ortho'
    elif box.shape == (3, 3):
        if np.all([box[0][1] == 0.0,  # Checks that tri box is properly formatted
                      box[0][2] == 0.0,
                      box[1][2] == 0.0]):
            boxtype = 'tri_vecs'
        else:
            boxtype = 'tri_vecs_bad'
    elif box.shape == (6,):
        if np.all(box[3:] == 90.):
            boxtype = 'ortho'
        else:
            boxtype = 'tri_box'

    if boxtype == 'unknown':
        raise ValueError("box input not recognised"
                         ", must be an array of box dimensions")

    return boxtype


def _check_array(coords, desc):
    """Check an array is a valid array of coordinates

    Must be:
       (n,3) in shape
       float32 data
    """
    if (coords.ndim != 2 or coords.shape[1] != 3):
        raise ValueError("{0} must be a sequence of 3 dimensional coordinates"
                         "".format(desc))
    _check_array_dtype(coords, desc)


def _check_array_dtype(coords, desc):
    """Check whether an array contains values of dtype: np.float32 or not"""
    if coords.dtype != np.float32:
        raise TypeError("{0} must be of type float32".format(desc))


def _check_results_array(results, size):
    """Check the results array is ok to use

    Must be:
      same shape as size
      float64
    """
    if results.shape != size:
        raise ValueError("Result array has incorrect size,"
                         "should be {0}, got {1}".format(size, results.shape))
    if results.dtype != np.float64:
        raise TypeError("Results array must be of type float64")


def _check_lengths_match(*arrays):
    """Check all arrays are same shape"""
    ref = arrays[0].shape

    if not all( a.shape == ref for a in arrays):
        raise ValueError("Input arrays must all be same shape"
                         "Got {0}".format([a.shape for a in arrays]))

def distance_array(reference, configuration, box=None, result=None, backend="serial"):
    """Calculate all distances between a reference set and another configuration.

    If there are *i* positions in reference, and *j* positions in configuration,
    will calculate a *i* x *j* array of distances
    If an *box* is supplied then a minimum image convention is used when
    calculating distances.

    If a 2D numpy array of dtype ``numpy.float64`` with the shape ``(len(reference),
    len(configuration))`` is provided in *result* then this preallocated array is
    filled. This can speed up calculations.

    Parameters
    ----------
    reference : numpy.ndarray
        Reference coordinate array of shape ``(n, 3)`` (``dtype`` is arbitrary,
        will be converted to ``dtype=numpy.float32`` internally)
    configuration : numpy.ndarray
        Configuration coordinate array of shape ``(m, 3)`` (``dtype`` is
        arbitrary, will be converted to ``dtype=numpy.float32`` internally)
    box : numpy.ndarray or None
        Dimensions of the cell; if provided, the minimum image convention is
        applied. The dimensions must be provided in the same format as returned
        by by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`: ``[lx,
        ly, lz, alpha, beta, gamma]``.
    result : numpy.ndarray(dtype=numpy.float64), optional
        Preallocated result array which must have the
        shape ``(len(ref), len(conf))`` and ``dtype=numpy.float64``.
        Avoids creating the array which saves time when the function
        is called repeatedly. [``None``]
    backend : str
        Select the type of acceleration; "serial" is always available. Other
        possibilities are "OpenMP" (OpenMP).

    Returns
    -------
    d : numpy.ndarray
        ``(len(reference),len(configuration))`` numpy array with the distances
        ``d[i,j]`` between reference coordinates `i` and configuration
        coordinates `j`.

    Note
    ----
    This method is slower than it could be because internally we need to make
    copies of the ref and conf arrays.


    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    .. versionchanged:: 0.19.0
       Internal dtype conversion of input coordinates to ``numpy.float32``.
    """
    ref = reference.astype(np.float32, order='C', copy=True)
    conf = configuration.astype(np.float32, order='C', copy=True)

    _check_array(conf, 'conf')
    _check_array(ref, 'ref')

    if box is not None:
        boxtype = _box_check(box)
        # Convert [A,B,C,alpha,beta,gamma] to [[A],[B],[C]]
        if (boxtype == 'tri_box'):
            box = triclinic_vectors(box)
        if (boxtype == 'tri_vecs_bad'):
            box = triclinic_vectors(triclinic_box(box[0], box[1], box[2]))

    confnum = conf.shape[0]
    refnum = ref.shape[0]

    if result is not None:
        _check_results_array(result, (refnum, confnum))
        distances = np.asarray(result)
    else:
        distances = np.zeros((refnum, confnum), np.float64)

    if box is not None:
        if boxtype == 'ortho':
            _run("calc_distance_array_ortho",
                   args=(ref, conf, box, distances),
                   backend=backend)
        else:
            _run("calc_distance_array_triclinic",
                   args=(ref, conf, box, distances),
                   backend=backend)
    else:
        _run("calc_distance_array",
               args=(ref, conf, distances),
               backend=backend)

    return distances


def self_distance_array(reference, box=None, result=None, backend="serial"):
    """Calculate all distances within a configuration *reference*.

    If a *box* is supplied then a minimum image convention is used before
    calculating distances.

    If a 1D numpy array of dtype ``numpy.float64`` with the shape
    ``(N*(N-1)/2)`` is provided in *result* then this preallocated array
    is filled. This can speed up calculations.

    Parameters
    ----------
    reference : numpy.ndarray
        Reference coordinate array with ``N=len(ref)`` coordinates (``dtype`` is
        arbitrary, will be converted to ``dtype=numpy.float32`` internally)
    box : numpy.ndarray or None
        Dimensions of the cell; if provided, the minimum image convention is
        applied. The dimensions must be provided in the same format as returned
        by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`: ``[lx,
        ly, lz, alpha, beta, gamma]``.
    result : numpy.ndarray(dtype=numpy.float64), optional
        Preallocated result array which must have the shape ``(N*(N-1)/2,)`` and
        dtype ``numpy.float64``. Avoids creating the array which saves time when
        the function is called repeatedly. [``None``]
    backend : str
        Select the type of acceleration; "serial" is always available. Other
        possibilities are "OpenMP" (OpenMP).

    Returns
    -------
    d : numpy.ndarray
        ``N*(N-1)/2`` numpy 1D array with the distances dist[i,j] between ref
        coordinates i and j at position d[k]. Loop through d:

        .. code-block:: python

            for i in range(N):
                for j in range(i+1, N):
                    k += 1
                    dist[i,j] = d[k]

    Note
    ----
    This method is slower than it could be because internally we need to make
    copies of the coordinate array.


    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    .. versionchanged:: 0.19.0
       Internal dtype conversion of input coordinates to ``numpy.float32``.
    """
    ref = reference.astype(np.float32, order='C', copy=True)

    _check_array(ref, 'ref')

    with_PBC = (box is not None)
    if box is not None:
        boxtype = _box_check(box)
        # Convert [A,B,C,alpha,beta,gamma] to [[A],[B],[C]]
        if (boxtype == 'tri_box'):
            box = triclinic_vectors(box)
        if (boxtype == 'tri_vecs_bad'):
            box = triclinic_vectors(triclinic_box(box[0], box[1], box[2]))

    refnum = ref.shape[0]
    distnum = refnum * (refnum - 1) // 2

    if result is not None:
        _check_results_array(result, (distnum,))
        distances = np.asarray(result)
    else:
        distances = np.zeros((distnum,), np.float64)

    if box is not None:
        if boxtype == 'ortho':
            _run("calc_self_distance_array_ortho",
                   args=(ref, box, distances),
                   backend=backend)
        else:
            _run("calc_self_distance_array_triclinic",
                   args=(ref, box, distances),
                   backend=backend)
    else:
        _run("calc_self_distance_array",
               args=(ref, distances),
               backend=backend)

    return distances


def capped_distance(reference, configuration, max_cutoff, min_cutoff=None,
                    box=None, method=None):
    """Calculates the pairs and distances within a specified distance

    If a *box* is supplied, then a minimum image convention is used
    to evaluate the distances.

    An automatic guessing of optimized method to calculate the distances is
    included in the function. An optional keyword for the method is also
    provided. Users can override the method with this functionality.
    Currently pkdtree and bruteforce are implemented.


    Parameters
    -----------
    reference : array
        reference coordinates array with shape ``reference.shape = (3,)``
        or ``reference.shape = (len(reference), 3)``
    configuration : array
        Configuration coordinate array with shape ``reference.shape = (3,)``
        or ``reference.shape = (len(reference), 3)``
    max_cutoff : float
        Maximum cutoff distance between the reference and configuration
    min_cutoff : (optional) float
        Minimum cutoff distance between reference and configuration [None]
    box : (optional) array or None
        The dimensions, if provided, must be provided in the same
        The unitcell dimesions for this system format as returned
        by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx,ly, lz, alpha, beta, gamma]``. Minimum image convention
        is applied if the box is provided [None]
    method : (optional) 'bruteforce' or 'pkdtree' or 'None'
        Keyword to override the automatic guessing of method built-in
        in the function [None]

    Returns
    -------
    pairs : array
        Pair of indices, one from each reference and configuration such that
        distance between them is  within the ``max_cutoff`` and ``min_cutoff``
        pairs[i,j] contains the indices i from reference coordinates, and
        j from configuration
    distances : array
        Distances corresponding to each pair of indices.
        d[k] corresponding to the pairs[i,j] gives the distance between
        i-th and j-th coordinate in reference and configuration respectively

        .. code-block:: python

            pairs, distances = capped_distances(reference, coordinates, max_cutoff)
            for indx, [a,b] in enumerate(pairs):
                coord1 = reference[a]
                coord2 = configuration[b]
                distance = distances[indx]

    Note
    -----
    Currently only supports brute force and Periodic KDtree

    .. SeeAlso:: :func:'MDAnalysis.lib.distances.distance_array'
    .. SeeAlso:: :func:'MDAnalysis.lib.pkdtree.PeriodicKDTree.search'
    .. SeeAlso:: :class:'MDAnalysis.lib.nsgrid.FastNS.search'
    """
    if box is not None:
        if box.shape[0] != 6:
            raise ValueError('Box Argument is of incompatible type. The dimension'
                         'should be either None or '
                         'of the type [lx, ly, lz, alpha, beta, gamma]')
    method = _determine_method(reference, configuration,
                               max_cutoff, min_cutoff=min_cutoff,
                               box=box, method=method)
    pairs, dist = method(reference, configuration, max_cutoff,
                         min_cutoff=min_cutoff, box=box)

    return np.asarray(pairs), np.asarray(dist)


def _determine_method(reference, configuration, max_cutoff, min_cutoff=None,
                      box=None, method=None):
    """
    Switch between different methods based on the the optimized time.
    All the rules to select the method based on the input can be
    incorporated here.

    Parameters
    ----------
    reference : array
        reference coordinates array with shape ``reference.shape = (3,)``
        or ``reference.shape = (len(reference), 3)``
    configuration : array
        Configuration coordinate array with shape ``reference.shape = (3,)``
        or ``reference.shape = (len(reference), 3)``
    max_cutoff : float
        Maximum cutoff distance between the reference and configuration
    min_cutoff : (optional) float
        Minimum cutoff distance between reference and configuration [None]
    box : (optional) array or None
        The dimensions, if provided, must be provided in the same
        The unitcell dimesions for this system format as returned
        by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx,ly, lz, alpha, beta, gamma]``. Minimum image convention
        is applied if the box is provided [None]
    method : (optional) 'bruteforce' or 'pkdtree' or 'None'
        Keyword to override the automatic guessing of method built-in
        in the function [None]

    Returns
    -------
    Method : Function object
        Returns function object based on the rules and specified method

    Note
    ----
    Currently implemented methods are present in the ``methods`` dictionary
        bruteforce : returns ``_bruteforce_capped``
        PKDtree : return ``_pkdtree_capped`
        NSGrid : return ``_nsgrid_capped`

    """
    methods = {'bruteforce': _bruteforce_capped,
               'pkdtree': _pkdtree_capped,
               'nsgrid': _nsgrid_capped}

    if method is not None:
        return methods[method]

    if len(reference) > 5000 and len(configuration) > 5000:
        if box is None:
            min_dim = np.array([reference.min(axis=0),
                               configuration.min(axis=0)])
            max_dim = np.array([reference.max(axis=0),
                               configuration.max(axis=0)])
            size = max_dim.max(axis=0) - min_dim.min(axis=0)
        elif np.allclose(box[3:], 90):
            size = box[:3]
        else:
            tribox = triclinic_vectors(box)
            size = tribox.max(axis=0) - tribox.min(axis=0)

        if ((np.any(size < 10.0*max_cutoff) and
                   len(reference) > 100000 and
                   len(configuration) > 100000)):
            return methods['bruteforce']
        else:
            return methods['pkdtree']
    return methods['bruteforce']


def _bruteforce_capped(reference, configuration, max_cutoff,
                       min_cutoff=None, box=None):
    """Internal method for bruteforce calculations

    Uses naive distance calulations and returns a list
    containing the indices with one from each
    reference and configuration arrays, such that the distance between
    them is less than the specified cutoff distance

    Returns
    -------
    pairs : list
        List of ``[(i, j)]`` pairs such that atom-index ``i`` is
        from reference and ``j`` from configuration array
    distance: list
         Distance between ``reference[i]`` and ``configuration[j]``
         atom coordinate

    """
    pairs, distance = [], []

    reference = np.asarray(reference, dtype=np.float32)
    configuration = np.asarray(configuration, dtype=np.float32)

    if reference.shape == (3, ):
        reference = reference[None, :]
    if configuration.shape == (3, ):
        configuration = configuration[None, :]

    _check_array(reference, 'reference')
    _check_array(configuration, 'configuration')

    for i, coords in enumerate(reference):
        dist = distance_array(coords[None, :], configuration, box=box)[0]
        if min_cutoff is not None:
            idx = np.where((dist < max_cutoff) & (dist > min_cutoff))[0]
        else:
            idx = np.where((dist < max_cutoff))[0]
        for j in idx:
            pairs.append((i, j))
            distance.append(dist[j])
    return pairs, distance


def _pkdtree_capped(reference, configuration, max_cutoff,
                    min_cutoff=None, box=None):
    """ Capped Distance evaluations using KDtree.

    Uses minimum image convention if *box* is specified

    Returns:
    --------
    pairs : list
        List of atom indices which are within the specified cutoff distance.
        pairs `(i, j)` corresponds to i-th particle in reference and
        j-th particle in configuration
    distance : list
        Distance between two atoms corresponding to the (i, j) indices
        in pairs.

    """
    from .pkdtree import PeriodicKDTree

    pairs, distances = [], []

    reference = np.asarray(reference, dtype=np.float32)
    configuration = np.asarray(configuration, dtype=np.float32)

    if reference.shape == (3, ):
        reference = reference[None, :]
    if configuration.shape == (3, ):
        configuration = configuration[None, :]

    _check_array(reference, 'reference')
    _check_array(configuration, 'configuration')

    kdtree = PeriodicKDTree(box=box)
    cut = max_cutoff if box is not None else None
    kdtree.set_coords(configuration, cutoff=cut)
    # Search for every query point
    for idx, centers in enumerate(reference):
        kdtree.search(centers, max_cutoff)
        indices = kdtree.get_indices()
        dist = distance_array(centers.reshape((1, 3)),
                              configuration[indices], box=box)[0]
        if min_cutoff is not None:
            mask = np.where(dist > min_cutoff)[0]
            dist = dist[mask]
            indices = [indices[mask[i]] for i in range(len(mask))]
        if len(indices) != 0:
            for num, j in enumerate(indices):
                pairs.append((idx, j))
                distances.append(dist[num])
    return pairs, distances


def _nsgrid_capped(reference, configuration, max_cutoff, min_cutoff=None,
                   box=None):
    """Search all the pairs in *reference* and *configuration* within
    a specified distance using Grid Search


    Parameters
    -----------
    reference : array
        reference coordinates array with shape ``reference.shape = (3,)``
        or ``reference.shape = (len(reference), 3)``.
    configuration : array
        Configuration coordinate array with shape ``reference.shape = (3,)``
        or ``reference.shape = (len(reference), 3)``
    max_cutoff : float
        Maximum cutoff distance between the reference and configuration
    min_cutoff : (optional) float
        Minimum cutoff distance between reference and configuration [None]
    box :  array
        The dimensions, if provided, must be provided in the same
        The unitcell dimesions for this system format as returned
        by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx,ly, lz, alpha, beta, gamma]``. Minimum image convention
        is applied if the box is provided

    """
    from .nsgrid import FastNS
    if reference.shape == (3, ):
        reference = reference[None, :]
    if configuration.shape == (3, ):
        configuration = configuration[None, :]

    if box is None:
        # create a pseudobox
        # define the max range
        # and supply the pseudobox
        # along with only one set of coordinates
        pseudobox = np.zeros(6, dtype=np.float32)
        all_coords = np.concatenate([reference, configuration])
        lmax = all_coords.max(axis=0)
        lmin = all_coords.min(axis=0)
        # Using maximum dimension as the box size
        boxsize = (lmax-lmin).max()
        # to avoid failures of very close particles
        # but with larger cutoff
        if boxsize < 2*max_cutoff:
            # just enough box size so that NSGrid doesnot fails
            sizefactor = 2.2*max_cutoff/boxsize
        else:
            sizefactor = 1.2
        pseudobox[:3] = sizefactor*boxsize
        pseudobox[3:] = 90.
        shiftref, shiftconf = reference.copy(), configuration.copy()
        # Extra padding near the origin
        shiftref -= lmin - 0.1*boxsize
        shiftconf -= lmin - 0.1*boxsize
        gridsearch = FastNS(max_cutoff, shiftconf, box=pseudobox, pbc=False)
        results = gridsearch.search(shiftref)
    else:
        gridsearch = FastNS(max_cutoff, configuration, box=box)
        results = gridsearch.search(reference)

    pairs = results.get_pairs()
    pair_distance = results.get_pair_distances()

    if min_cutoff is not None:
        idx = pair_distance > min_cutoff
        pairs, pair_distance = pairs[idx], pair_distance[idx]
    return pairs, pair_distance


def self_capped_distance(reference, max_cutoff, min_cutoff=None,
                         box=None, method=None):
    """Finds all the pairs and respective distances within a specified cutoff
    for a configuration *reference*

    If a *box* is supplied, then a minimum image convention is used
    to evaluate the distances.

    An automatic guessing of optimized method to calculate the distances is
    included in the function. An optional keyword for the method is also
    provided. Users can override the method with this functionality.
    Currently pkdtree and bruteforce are implemented.

    Parameters
    -----------
    reference : array
        reference coordinates array with shape ``reference.shape = (3,)``
        or ``reference.shape = (len(reference), 3)``
    max_cutoff : float
        Maximum cutoff distance to check the neighbors with itself
    min_cutoff : (optional) float
        Minimum cutoff distance [None]
    box : (optional) array or None
        The dimensions, if provided, must be provided in the same
        The unitcell dimesions for this system format as returned
        by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx,ly, lz, alpha, beta, gamma]``. Minimum image convention
        is applied if the box is provided [None]
    method : (optional) 'bruteforce' or 'pkdtree' or 'None'
        Keyword to override the automatic guessing of method built-in
        in the function [None]

    Returns
    -------
    pairs : array
        Pair of indices such that distance between them is
        within the ``max_cutoff`` and ``min_cutoff``
    distances : array
        Distances corresponding to each pair of indices.
        d[k] corresponding to the pairs[i,j] gives the distance between
        i-th and j-th coordinate in reference

        .. code-block:: python

            pairs, distances = self_capped_distances(reference, max_cutoff)
            for indx, [a,b] in enumerate(pairs):
                coord1, coords2 = reference[a], reference[b]
                distance = distances[indx]

    Note
    -----
    Currently only supports brute force, Periodic KDtree and Grid Search

    .. SeeAlso:: :func:'MDAnalysis.lib.distances.self_distance_array'
    .. SeeAlso:: :func:'MDAnalysis.lib.pkdtree.PeriodicKDTree.search'
    .. SeeAlso:: :func:'MDAnalysis.lib.nsgrid.FastNS.self_search'
    """
    if box is not None:
        if box.shape[0] != 6:
            raise ValueError('Box Argument is of incompatible type. The dimension'
                         'should be either None or '
                         'of the type [lx, ly, lz, alpha, beta, gamma]')
    method = _determine_method_self(reference, max_cutoff,
                                    min_cutoff=min_cutoff,
                                    box=box, method=method)
    pairs, dist = method(reference,  max_cutoff,
                         min_cutoff=min_cutoff, box=box)

    return np.asarray(pairs), np.asarray(dist)


def _determine_method_self(reference, max_cutoff, min_cutoff=None,
                           box=None, method=None):
    """
    Switch between different methods based on the the optimized time.
    All the rules to select the method based on the input can be
    incorporated here.

    Parameters
    ----------
    reference : array
        reference coordinates array with shape ``reference.shape = (3,)``
        or ``reference.shape = (len(reference), 3)``
    max_cutoff : float
        Maximum cutoff distance
    min_cutoff : (optional) float
        Minimum cutoff distance [None]
    box : (optional) array or None
        The dimensions, if provided, must be provided in the same
        The unitcell dimesions for this system format as returned
        by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx,ly, lz, alpha, beta, gamma]``. Minimum image convention
        is applied if the box is provided [None]
    method : (optional) 'bruteforce' or 'pkdtree' or 'None'
        Keyword to override the automatic guessing of method built-in
        in the function [None]

    Returns
    -------
    Method : Function object
        Returns function object based on the rules and specified method

    Note
    ----
    Currently implemented methods are present in the ``methods`` dictionary
        bruteforce : returns ``_bruteforce_capped_self``
        PKDtree : return ``_pkdtree_capped_self``
        NSGrid : return ``_nsgrid_capped_self``

    """
    methods = {'bruteforce': _bruteforce_capped_self,
            'pkdtree': _pkdtree_capped_self,
            'nsgrid': _nsgrid_capped_self}

    if method is not None:
        return methods[method]

    if len(reference) > 5000:
        if box is None:
            min_dim = np.array([reference.min(axis=0)])
            max_dim = np.array([reference.max(axis=0)])
            size = max_dim.max(axis=0) - min_dim.min(axis=0)
        elif np.allclose(box[3:], 90):
            size = box[:3]
        else:
            tribox = triclinic_vectors(box)
            size = tribox.max(axis=0) - tribox.min(axis=0)

        if ((np.any(size < 10.0*max_cutoff) and
                   (len(reference) > 100000))):
            return methods['bruteforce']
        else:
            return methods['pkdtree']
    return methods['bruteforce']


def _bruteforce_capped_self(reference, max_cutoff, min_cutoff=None,
                            box=None):
    """Finds all the pairs among the *reference* coordinates within
    a fixed distance using brute force method

    Internal method using brute force method to evaluate all the pairs
    of atoms within a fixed distance.

    Returns
    -------
    pairs : array
        Arrray of ``[i, j]`` pairs such that atom-index ``i``
        and ``j`` from reference array are within the fixed distance
    distance: array
         Distance between ``reference[i]`` and ``reference[j]``
         atom coordinate

    """
    pairs, distance = [], []

    reference = np.asarray(reference, dtype=np.float32)
    if reference.shape == (3, ):
        reference = reference[None, :]
    for i, coords in enumerate(reference):
        # Each pair of atoms needs to be checked only once.
        # Only calculate distance for atomA and atomB
        # if atomidA < atomidB
        dist = distance_array(coords[None, :], reference[i+1:],
                             box=box)[0]

        if min_cutoff is not None:
            idx = np.where((dist < max_cutoff) & (dist > min_cutoff))[0]
        else:
            idx = np.where((dist < max_cutoff))[0]
        for other_idx in idx:
            # Actual atomid for atomB
            # can be direclty obtained in this way
            j = other_idx + 1 + i
            pairs.append((i, j))
            distance.append(dist[other_idx])
    return np.asarray(pairs), np.asarray(distance)


def _pkdtree_capped_self(reference, max_cutoff, min_cutoff=None,
                         box=None):
    """Finds all the pairs among the coordinates within a fixed distance
    using PeriodicKDTree

    Internal method using PeriodicKDTree method to evaluate all the pairs
    of atoms within a fixed distance.

    Returns
    -------
    pairs : array
        Array of ``[(i, j)]`` pairs such that atom-index ``i``
        and ``j`` from reference array are within the fixed distance
    distance: array
         Distance between ``reference[i]`` and ``reference[j]``
         atom coordinate

    """
    from .pkdtree import PeriodicKDTree

    reference = np.asarray(reference, dtype=np.float32)
    if reference.shape == (3, ):
        reference = reference[None, :]

    pairs, distance = [], []
    kdtree = PeriodicKDTree(box=box)
    cut = max_cutoff if box is not None else None
    kdtree.set_coords(reference, cutoff=cut)
    pairs = kdtree.search_pairs(max_cutoff)
    if pairs.size > 0:
        refA, refB = pairs[:, 0], pairs[:, 1]
        distance = calc_bonds(reference[refA], reference[refB], box=box)
        if min_cutoff is not None:
            mask = np.where(distance > min_cutoff)[0]
            pairs, distance = pairs[mask], distance[mask]
    return np.asarray(pairs), np.asarray(distance)


def _nsgrid_capped_self(reference, max_cutoff, min_cutoff=None,
                        box=None):
    """Finds all the pairs among the *reference* coordinates within
    a fixed distance using gridsearch

    Returns
    -------
    pairs : array
        Arrray of ``[i, j]`` pairs such that atom-index ``i``
        and ``j`` from reference array are within a fixed distance
    distance: array
         Distance between ``reference[i]`` and ``reference[j]``
         atom coordinate

    """
    from .nsgrid import FastNS

    reference = np.asarray(reference, dtype=np.float32)
    if reference.shape == (3, ) or len(reference) == 1:
        return [], []

    if box is None:
        # create a pseudobox
        # define the max range
        # and supply the pseudobox
        # along with only one set of coordinates
        pseudobox = np.zeros(6, dtype=np.float32)
        lmax = reference.max(axis=0)
        lmin = reference.min(axis=0)
        # Using maximum dimension as the box size
        boxsize = (lmax-lmin).max()
        # to avoid failures of very close particles
        # but with larger cutoff
        if boxsize < 2*max_cutoff:
            # just enough box size so that NSGrid doesnot fails
            sizefactor = 2.2*max_cutoff/boxsize
        else:
            sizefactor = 1.2
        pseudobox[:3] = sizefactor*boxsize
        pseudobox[3:] = 90.
        shiftref = reference.copy()
        # Extra padding near the origin
        shiftref -= lmin - 0.1*boxsize
        gridsearch = FastNS(max_cutoff, shiftref, box=pseudobox, pbc=False)
        results = gridsearch.self_search()
    else:
        gridsearch = FastNS(max_cutoff, reference, box=box)
        results = gridsearch.self_search()

    pairs = results.get_pairs()[::2, :]
    pair_distance = results.get_pair_distances()[::2]

    if min_cutoff is not None:
        idx = pair_distance > min_cutoff
        pairs, pair_distance = pairs[idx], pair_distance[idx]
    return pairs, pair_distance


def transform_RtoS(inputcoords, box, backend="serial"):
    """Transform an array of coordinates from real space to S space (aka lambda space)

    S space represents fractional space within the unit cell for this system

    Reciprocal operation to :meth:`transform_StoR`

    Parameters
    ----------
    inputcoords : array
        A (3,) coordinate or (n x 3) array of coordinates (``dtype`` is
        arbitrary, will be converted to ``dtype=numpy.float32`` internally)
    box : array
        The unitcell dimesions for this system.
        The dimensions must be provided in the same format as returned
        by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`: ``[lx,
        ly, lz, alpha, beta, gamma]``.
    backend : str, optional
        Select the type of acceleration; "serial" is always available. Other
        possibilities are "OpenMP" (OpenMP).

    Returns
    -------
    outcoords : array
        A n x 3 array of fractional coordiantes.


    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    .. versionchanged:: 0.19.0
       Internal dtype conversion of input coordinates to ``numpy.float32``.
    """
    coords = inputcoords.astype(np.float32, order='C', copy=True)

    is_1d = False  # True if only one vector coordinate
    if len(coords.shape) == 1:
        coords = coords.reshape(1, len(coords))
        is_1d = True

    boxtype = _box_check(box)
    # Convert [A,B,C,alpha,beta,gamma] to [[A],[B],[C]]
    if (boxtype == 'tri_box'):
        box = triclinic_vectors(box)
    if (boxtype == 'tri_vecs_bad'):
        box = triclinic_vectors(triclinic_box(box[0], box[1], box[2]))
    elif (boxtype == 'ortho'):
        box = np.array([[box[0], 0.0, 0.0],
                           [0.0, box[1], 0.0],
                           [0.0, 0.0, box[2]]], dtype=np.float32)

    # Create inverse matrix of box
    # need order C here
    inv = np.array(np.matrix(box).I, dtype=np.float32, order='C')

    _run("coord_transform",
           args=(coords, inv),
           backend=backend)

    if is_1d:
        coords = coords.reshape(len(coords[0]), )
    return coords


def transform_StoR(inputcoords, box, backend="serial"):
    """Transform an array of coordinates from S space into real space.

    S space represents fractional space within the unit cell for this system

    Reciprocal operation to :meth:`transform_RtoS`

    Parameters
    ----------
    inputcoords : array
        A (3,) coordinate or (n x 3) array of coordinates (``dtype`` is
        arbitrary, will be converted to ``dtype=numpy.float32`` internally)
    box : array
        The unitcell dimesions for this system.
        The dimensions must be provided in the same format as returned
        by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`: ``[lx,
        ly, lz, alpha, beta, gamma]``.
    backend : str, optional
        Select the type of acceleration; "serial" is always available. Other
        possibilities are "OpenMP" (OpenMP).

    Returns
    -------
    outcoords : array
        A n x 3 array of fractional coordiantes.


    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    .. versionchanged:: 0.19.0
       Internal dtype conversion of input coordinates to ``numpy.float32``.
    """
    coords = inputcoords.astype(np.float32, order='C', copy=True)

    is_1d = False  # True if only one vector coordinate
    if len(coords.shape) == 1:
        coords = coords.reshape(1, len(coords))
        is_1d = True

    boxtype = _box_check(box)
    # Convert [A,B,C,alpha,beta,gamma] to [[A],[B],[C]]
    if (boxtype == 'tri_box'):
        box = triclinic_vectors(box)
    elif (boxtype == 'ortho'):
        box = np.array([[box[0], 0.0, 0.0],
                        [0.0, box[1], 0.0],
                        [0.0, 0.0, box[2]]], dtype=np.float32)

    _run("coord_transform",
           args=(coords, box),
           backend=backend)

    if is_1d:
        coords = coords.reshape(len(coords[0]),)
    return coords


def calc_bonds(coords1, coords2, box=None, result=None, backend="serial"):
    """
    Calculate all distances between a pair of atoms.  *atom1* and *atom2* are both
    arrays of coordinates, where atom1[i] and atom2[i] represent a bond.

    In comparison to distance_array and self_distance_array which calculate distances
    between all combinations of coordinates, calc_bonds can be used to calculate distance
    between pairs of objects, similar to::

       numpy.linalg.norm(a - b) for a, b in zip(coords1, coords2)

    The optional argument *box* applies minimum image convention if supplied.
    *box* can be either orthogonal or triclinic

    If a 1D numpy array of dtype ``numpy.float64`` with ``len(atom1)`` elements is
    provided in *result* then this preallocated array is filled. This can speed
    up calculations.

    bondlengths = calc_bonds(coords1, coords2 [, box [,result=bondlengths]])

    Parameters
    ----------
    coords1 : array
        An array of coordinates for one half of the bond (``dtype`` is
        arbitrary, will be converted to ``dtype=numpy.float32`` internally)
    coords2 : array
        An array of coordinates for the other half of bond (``dtype`` is
        arbitrary, will be converted to ``dtype=numpy.float32`` internally)
    box : array
        The unitcell dimesions for this system.
        The dimensions must be provided in the same format as returned
        by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`: ``[lx,
        ly, lz, alpha, beta, gamma]``.
    result : array, optional
        Preallocated result array which must be same length as coord
        arrays and ``dtype=numpy.float64``. Avoids creating the
        array which saves time when the function is called repeatedly. [None]
    backend : str
        Select the type of acceleration; "serial" is always available. Other
        possibilities are "OpenMP" (OpenMP).

    Returns
    -------
    bondlengths : array
        The length between each pair in coords1 and coords2


    .. versionadded:: 0.8
    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    .. versionchanged:: 0.19.0
       Internal dtype conversion of input coordinates to ``numpy.float32``.
    """
    atom1 = coords1.astype(np.float32, order='C', copy=True)
    atom2 = coords2.astype(np.float32, order='C', copy=True)

    _check_array(atom1, 'atom1')
    _check_array(atom2, 'atom2')
    _check_lengths_match(atom1, atom2)

    if box is not None:
        boxtype = _box_check(box)
        # Convert [A,B,C,alpha,beta,gamma] to [[A],[B],[C]]
        if (boxtype == 'tri_box'):
            box = triclinic_vectors(box)
        if (boxtype == 'tri_vecs_bad'):
            box = triclinic_vectors(triclinic_box(box[0], box[1], box[2]))

    numatom = atom1.shape[0]

    if result is not None:
        _check_results_array(result, (numatom,))
        distances = np.asarray(result)
    else:
        distances = np.zeros((numatom,), np.float64)

    if box is not None:
        if boxtype == 'ortho':
            _run("calc_bond_distance_ortho",
                   args=(atom1, atom2, box, distances),
                   backend=backend)
        else:
            _run("calc_bond_distance_triclinic",
                   args=(atom1, atom2, box, distances),
                   backend=backend)
    else:
        _run("calc_bond_distance",
               args=(atom1, atom2, distances),
               backend=backend)

    return distances


def calc_angles(coords1, coords2, coords3, box=None, result=None, backend="serial"):
    """
    Calculates the angle formed between three atoms, over a list of coordinates.
    All *atom* inputs are lists of coordinates of equal length, with *atom2*
    representing the apex of the angle.

    If a 1D numpy array of dtype ``numpy.float64`` with ``len(atom1)`` elements is
    provided in *result* then this preallocated array is filled. This can speed
    up calculations.

    The optional argument ``box`` ensures that periodic boundaries are taken into account when
    constructing the connecting vectors between atoms, ie that the vector between atoms 1 & 2
    goes between coordinates in the same image.

    angles = calc_angles(coords1, coords2, coords3, [[box=None],result=angles])

    Parameters
    ----------
    coords1 : numpy.ndarray
        Coordinate array of one side of angles (``dtype`` is arbitrary, will be
        converted to ``dtype=numpy.float32`` internally)
    coords2 : numpy.ndarray
        Coordinate array of apex of angles (``dtype`` is arbitrary, will be
        converted to ``dtype=numpy.float32`` internally)
    coords3 : numpy.ndarray
        Coordinate array of other side of angles (``dtype`` is arbitrary, will be
        converted to ``dtype=numpy.float32`` internally)
    box : numpy.ndarray, optional
        The unitcell dimesions for this system.
        The dimensions must be provided in the same format as returned
        by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`: ``[lx,
        ly, lz, alpha, beta, gamma]``.
    result : numpy.ndarray, optional
        Preallocated result array which must be same length as coord
        arrays and ``dtype=numpy.float64``. Avoids creating the
        array which saves time when the function is called repeatedly. [None]
    backend : str
        Select the type of acceleration; "serial" is always available. Other
        possibilities are "OpenMP" (OpenMP).

    Returns
    -------
    angles : numpy.ndarray
        An array of angles in radians.


    .. versionadded:: 0.8
    .. versionchanged:: 0.9.0
       Added optional box argument to account for periodic boundaries in calculation
    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    .. versionchanged:: 0.19.0
       Internal dtype conversion of input coordinates to ``numpy.float32``.
    """
    atom1 = coords1.astype(np.float32, order='C', copy=True)
    atom2 = coords2.astype(np.float32, order='C', copy=True)
    atom3 = coords3.astype(np.float32, order='C', copy=True)
    numatom = atom1.shape[0]

    _check_array(atom1, 'coords1')
    _check_array(atom2, 'coords2')
    _check_array(atom3, 'coords3')
    _check_lengths_match(atom1, atom2, atom3)

    if box is not None:
        boxtype = _box_check(box)
        # Convert [A,B,C,alpha,beta,gamma] to [[A],[B],[C]]
        if (boxtype == 'tri_box'):
            box = triclinic_vectors(box)
        if (boxtype == 'tri_vecs_bad'):
            box = triclinic_vectors(triclinic_box(box[0], box[1], box[2]))

    if result is not None:
        _check_results_array(result, (numatom,))
        angles = np.asarray(result)
    else:
        angles = np.zeros((numatom,), np.float64)

    if box is not None:
        if boxtype == 'ortho':
            _run("calc_angle_ortho",
                   args=(atom1, atom2, atom3, box, angles),
                   backend=backend)
        else:
            _run("calc_angle_triclinic",
                   args=(atom1, atom2, atom3, box, angles),
                   backend=backend)
    else:
        _run("calc_angle",
               args=(atom1, atom2, atom3, angles),
               backend=backend)

    return angles


def calc_dihedrals(coords1, coords2, coords3, coords4, box=None, result=None,
                   backend="serial"):
    """
    Calculate the dihedral angle formed by four atoms, over a list of coordinates.

    Dihedral angle around axis connecting atoms 1 and 2 (i.e. the angle
    between the planes spanned by atoms (0,1,2) and (1,2,3))::

                  3
                  |
            1-----2
           /
          0

    If a 1D numpy array of dtype ``numpy.float64`` with ``len(atom1)`` elements
    is provided in *result* then this preallocated array is filled. This can
    speed up calculations.

    The optional argument ``box`` ensures that periodic boundaries are taken
    into account when constructing the connecting vectors between atoms, ie
    that the vector between atoms 1 & 2 goes between coordinates in the same
    image::

        angles = calc_dihedrals(coords1, coords2, coords3, coords4 [,box=box, result=angles])

    Parameters
    ----------
    coords1 : array
        Coordinate array of 1st atom in dihedrals (``dtype`` is arbitrary, will
        be converted to ``dtype=numpy.float32`` internally)
    coords2 : array
        Coordinate array of 2nd atom in dihedrals (``dtype`` is arbitrary, will
        be converted to ``dtype=numpy.float32`` internally)
    coords3 : array
        Coordinate array of 3rd atom in dihedrals (``dtype`` is arbitrary, will
        be converted to ``dtype=numpy.float32`` internally)
    coords4 : array
        Coordinate array of 4th atom in dihedrals (``dtype`` is arbitrary, will
        be converted to ``dtype=numpy.float32`` internally)
    box : array
        The unitcell dimesions for this system.
        The dimensions must be provided in the same format as returned
        by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`: ``[lx,
        ly, lz, alpha, beta, gamma]``.
    result : array, optional
        Preallocated result array which must be same length as coord
        arrays and ``dtype=numpy.float64``. Avoids creating the
        array which saves time when the function is called repeatedly. [None]
    backend : str
        Select the type of acceleration; "serial" is always available. Other
        possibilities are "OpenMP" (OpenMP).

    Returns
    -------
    angles : array
        A numpy.array of angles in radians.


    .. versionadded:: 0.8
    .. versionchanged:: 0.9.0
       Added optional box argument to account for periodic boundaries in calculation
    .. versionchanged:: 0.11.0
       Renamed from calc_torsions to calc_dihedrals
    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    .. versionchanged:: 0.19.0
       Internal dtype conversion of input coordinates to ``numpy.float32``.
    """
    atom1 = coords1.astype(np.float32, order='C', copy=True)
    atom2 = coords2.astype(np.float32, order='C', copy=True)
    atom3 = coords3.astype(np.float32, order='C', copy=True)
    atom4 = coords4.astype(np.float32, order='C', copy=True)

    _check_array(atom1, 'atom1')
    _check_array(atom2, 'atom2')
    _check_array(atom3, 'atom3')
    _check_array(atom4, 'atom4')
    _check_lengths_match(atom1, atom2, atom3, atom4)

    numatom = atom1.shape[0]

    if box is not None:
        boxtype = _box_check(box)
        # Convert [A,B,C,alpha,beta,gamma] to [[A],[B],[C]]
        if (boxtype == 'tri_box'):
            box = triclinic_vectors(box)
        if (boxtype == 'tri_vecs_bad'):
            box = triclinic_vectors(triclinic_box(box[0], box[1], box[2]))

    if result is not None:
        _check_results_array(result, (numatom,))
        angles = np.asarray(result)
    else:
        angles = np.zeros((numatom,), np.float64)

    if box is not None:
        if boxtype == 'ortho':
            _run("calc_dihedral_ortho",
                   args=(atom1, atom2, atom3, atom4, box, angles),
                   backend=backend)
        else:
            _run("calc_dihedral_triclinic",
                   args=(atom1, atom2, atom3, atom4, box, angles),
                   backend=backend)
    else:
        _run("calc_dihedral",
               args=(atom1, atom2, atom3, atom4, angles),
               backend=backend)

    return angles


def apply_PBC(incoords, box, backend="serial"):
    """Moves a set of coordinates to all be within the primary unit cell

    newcoords = apply_PBC(coords, box)

    Parameters
    ----------
    incoords : numpy.ndarray
        Coordinate array of shape ``(n, 3)`` (``dtype`` is arbitrary, will be
        converted to ``dtype=numpy.float32`` internally)
    box : array
        The unitcell dimesions for this system; can be either orthogonal or
        triclinic information. The dimensions must be provided in the same
        format as returned by
        :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`: ``[lx, ly, lz,
        alpha, beta, gamma]``.
    backend : str
        Select the type of acceleration; ``"serial"`` is always available. Other
        possibilities are ``"OpenMP"`` (OpenMP).

    Returns
    -------
    newcoords : numpy.ndarray(dtype=numpy.float32)
        Coordinates that are now all within the primary unit cell, as defined
        by box.


    .. versionadded:: 0.8
    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    .. versionchanged:: 0.19.0
       Internal dtype conversion of input coordinates to ``numpy.float32``.
    """
    coords = incoords.astype(np.float32, order='C', copy=True)

    _check_array(coords, 'coords')

    coordnum = coords.shape[0]

    # determine boxtype
    boxtype = _box_check(box)
    # Convert [A,B,C,alpha,beta,gamma] to [[A],[B],[C]]
    if boxtype == 'tri_box':
        box = triclinic_vectors(box)
    if boxtype == 'tri_vecs_bad':
        box = triclinic_vectors(triclinic_box(box[0], box[1], box[2]))

    box_inv = np.zeros((3), dtype=np.float32)
    if boxtype == 'ortho':
        box_inv[0] = 1.0 / box[0]
        box_inv[1] = 1.0 / box[1]
        box_inv[2] = 1.0 / box[2]
        _run("ortho_pbc",
               args=(coords, box, box_inv),
               backend=backend)
    else:
        box_inv[0] = 1.0 / box[0][0]
        box_inv[1] = 1.0 / box[1][1]
        box_inv[2] = 1.0 / box[2][2]
        _run("triclinic_pbc",
               args=(coords, box, box_inv),
               backend=backend)

    return coords


def calc_distance(a, b, box=None):
    """Distance between a and b

    Parameters
    ----------
    a, b : numpy.ndarray
        single coordinate vectors
    box : numpy.ndarray, optional
        simulation box, if given periodic boundary conditions will be applied


    .. versionadded:: 0.18.1
    """
    return calc_bonds(a[None, :], b[None, :], box=box)[0]


def calc_angle(a, b, c, box=None):
    """Angle (in degrees) between a, b and c, where b is apex of angle

    Parameters
    ----------
    a, b, c : numpy.ndarray
        single coordinate vectors
    box : numpy.ndarray
        simulation box if given periodic boundary conditions will be applied to
        the vectors between atoms


    .. versionadded:: 0.18.1
    """
    return np.rad2deg(calc_angles(a[None, :], b[None, :], c[None, :], box=box)[0])


def calc_dihedral(a, b, c, d, box=None):
    """Dihedral angle (in degrees) between planes (a, b, c) and (b, c, d)

    Parameters
    ----------
    a, b, c, d : numpy.ndarray
        single coordinate vectors
    box : numpy.ndarray, optional
        simulation box, if given periodic boundary conditions will be applied


    .. versionadded:: 0.18.1
    """
    return np.rad2deg(
        calc_dihedrals(a[None, :], b[None, :], c[None, :], d[None, :], box)[0])
