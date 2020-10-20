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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
#

"""Fast distance array computation --- :mod:`MDAnalysis.lib.distances`
===================================================================

Fast C-routines to calculate arrays of distances or angles from coordinate
arrays. Many of the functions also exist in parallel versions, which typically
provide higher performance than the serial code.
The boolean attribute `MDAnalysis.lib.distances.USED_OPENMP` can be checked to
see if OpenMP was used in the compilation of MDAnalysis.

Selection of acceleration ("backend")
-------------------------------------

All functions take the optional keyword `backend`, which determines the type of
acceleration. Currently, the following choices are implemented (`backend` is
case-insensitive):

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
.. autofunction:: distance_array
.. autofunction:: self_distance_array
.. autofunction:: capped_distance
.. autofunction:: self_capped_distance
.. autofunction:: calc_bonds
.. autofunction:: calc_angles
.. autofunction:: calc_dihedrals
.. autofunction:: apply_PBC
.. autofunction:: transform_RtoS
.. autofunction:: transform_StoR
.. autofunction:: augment_coordinates(coordinates, box, r)
.. autofunction:: undo_augment(results, translation, nreal)
"""
import numpy as np
from numpy.lib.utils import deprecate

from .util import check_coords, check_box
from .mdamath import triclinic_vectors
from ._augment import augment_coordinates, undo_augment
from .nsgrid import FastNS

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
    """Helper function to select a backend function `funcname`."""
    args = args if args is not None else tuple()
    kwargs = kwargs if kwargs is not None else dict()
    backend = backend.lower()
    try:
        func = getattr(_distances[backend], funcname)
    except KeyError:
        errmsg = (f"Function {funcname} not available with backend {backend} "
                  f"try one of: {_distances.keys()}")
        raise ValueError(errmsg) from None
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


def _check_result_array(result, shape):
    """Check if the result array is ok to use.

    The `result` array must meet the following requirements:
      * Must have a shape equal to `shape`.
      * Its dtype must be ``numpy.float64``.

    Paramaters
    ----------
    result : numpy.ndarray or None
        The result array to check. If `result` is `None``, a newly created
        array of correct shape and dtype ``numpy.float64`` will be returned.
    shape : tuple
        The shape expected for the `result` array.

    Returns
    -------
    result : numpy.ndarray (``dtype=numpy.float64``, ``shape=shape``)
        The input array or a newly created array if the input was ``None``.

    Raises
    ------
    ValueError
        If `result` is of incorrect shape.
    TypeError
        If the dtype of `result` is not ``numpy.float64``.
    """
    if result is None:
        return np.zeros(shape, dtype=np.float64)
    if result.shape != shape:
        raise ValueError("Result array has incorrect shape, should be {0}, got "
                         "{1}.".format(shape, result.shape))
    if result.dtype != np.float64:
        raise TypeError("Result array must be of type numpy.float64, got {}."
                        "".format(result.dtype))
# The following two lines would break a lot of tests. WHY?!
#    if not coords.flags['C_CONTIGUOUS']:
#        raise ValueError("{0} is not C-contiguous.".format(desc))
    return result


@check_coords('reference', 'configuration', reduce_result_if_single=False,
              check_lengths_match=False)
def distance_array(reference, configuration, box=None, result=None,
                   backend="serial"):
    """Calculate all possible distances between a reference set and another
    configuration.

    If there are ``n`` positions in `reference` and ``m`` positions in
    `configuration`, a distance array of shape ``(n, m)`` will be computed.

    If the optional argument `box` is supplied, the minimum image convention is
    applied when calculating distances. Either orthogonal or triclinic boxes are
    supported.

    If a 2D numpy array of dtype ``numpy.float64`` with the shape ``(n, m)``
    is provided in `result`, then this preallocated array is filled. This can
    speed up calculations.

    Parameters
    ----------
    reference : numpy.ndarray
        Reference coordinate array of shape ``(3,)`` or ``(n, 3)`` (dtype is
        arbitrary, will be converted to ``numpy.float32`` internally).
    configuration : numpy.ndarray
        Configuration coordinate array of shape ``(3,)`` or ``(m, 3)`` (dtype is
        arbitrary, will be converted to ``numpy.float32`` internally).
    box : array_like, optional
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx, ly, lz, alpha, beta, gamma]``.
    result : numpy.ndarray, optional
        Preallocated result array which must have the shape ``(n, m)`` and dtype
        ``numpy.float64``.
        Avoids creating the array which saves time when the function
        is called repeatedly.
    backend : {'serial', 'OpenMP'}, optional
        Keyword selecting the type of acceleration.

    Returns
    -------
    d : numpy.ndarray (``dtype=numpy.float64``, ``shape=(n, m)``)
        Array containing the distances ``d[i,j]`` between reference coordinates
        ``i`` and configuration coordinates ``j``.


    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    .. versionchanged:: 0.19.0
       Internal dtype conversion of input coordinates to ``numpy.float32``.
       Now also accepts single coordinates as input.
    """
    confnum = configuration.shape[0]
    refnum = reference.shape[0]

    distances = _check_result_array(result, (refnum, confnum))
    if len(distances) == 0:
        return distances

    if box is not None:
        boxtype, box = check_box(box)
        if boxtype == 'ortho':
            _run("calc_distance_array_ortho",
                 args=(reference, configuration, box, distances),
                 backend=backend)
        else:
            _run("calc_distance_array_triclinic",
                 args=(reference, configuration, box, distances),
                 backend=backend)
    else:
        _run("calc_distance_array",
             args=(reference, configuration, distances),
             backend=backend)

    return distances


@check_coords('reference', reduce_result_if_single=False)
def self_distance_array(reference, box=None, result=None, backend="serial"):
    """Calculate all possible distances within a configuration `reference`.

    If the optional argument `box` is supplied, the minimum image convention is
    applied when calculating distances. Either orthogonal or triclinic boxes are
    supported.

    If a 1D numpy array of dtype ``numpy.float64`` with the shape
    ``(n*(n-1)/2,)`` is provided in `result`, then this preallocated array is
    filled. This can speed up calculations.

    Parameters
    ----------
    reference : numpy.ndarray
        Reference coordinate array of shape ``(3,)`` or ``(n, 3)`` (dtype is
        arbitrary, will be converted to ``numpy.float32`` internally).
    box : array_like, optional
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx, ly, lz, alpha, beta, gamma]``.
    result : numpy.ndarray, optional
        Preallocated result array which must have the shape ``(n*(n-1)/2,)`` and
        dtype ``numpy.float64``. Avoids creating the array which saves time when
        the function is called repeatedly.
    backend : {'serial', 'OpenMP'}, optional
        Keyword selecting the type of acceleration.

    Returns
    -------
    d : numpy.ndarray (``dtype=numpy.float64``, ``shape=(n*(n-1)/2,)``)
        Array containing the distances ``dist[i,j]`` between reference
        coordinates ``i`` and ``j`` at position ``d[k]``. Loop through ``d``:

        .. code-block:: python

            for i in range(n):
                for j in range(i + 1, n):
                    k += 1
                    dist[i, j] = d[k]


    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    .. versionchanged:: 0.19.0
       Internal dtype conversion of input coordinates to ``numpy.float32``.
    """
    refnum = reference.shape[0]
    distnum = refnum * (refnum - 1) // 2

    distances = _check_result_array(result, (distnum,))
    if len(distances) == 0:
        return distances

    if box is not None:
        boxtype, box = check_box(box)
        if boxtype == 'ortho':
            _run("calc_self_distance_array_ortho",
                 args=(reference, box, distances),
                 backend=backend)
        else:
            _run("calc_self_distance_array_triclinic",
                 args=(reference, box, distances),
                 backend=backend)
    else:
        _run("calc_self_distance_array",
             args=(reference, distances),
             backend=backend)

    return distances


def capped_distance(reference, configuration, max_cutoff, min_cutoff=None,
                    box=None, method=None, return_distances=True):
    """Calculates pairs of indices corresponding to entries in the `reference`
    and `configuration` arrays which are separated by a distance lying within
    the specified cutoff(s). Optionally, these distances can be returned as
    well.

    If the optional argument `box` is supplied, the minimum image convention is
    applied when calculating distances. Either orthogonal or triclinic boxes are
    supported.

    An automatic guessing of the optimal method to calculate the distances is
    included in the function. An optional keyword for the method is also
    provided. Users can enforce a particular method with this functionality.
    Currently brute force, grid search, and periodic KDtree methods are
    implemented.

    Parameters
    -----------
    reference : numpy.ndarray
        Reference coordinate array with shape ``(3,)`` or ``(n, 3)``.
    configuration : numpy.ndarray
        Configuration coordinate array with shape ``(3,)`` or ``(m, 3)``.
    max_cutoff : float
        Maximum cutoff distance between the reference and configuration.
    min_cutoff : float, optional
        Minimum cutoff distance between reference and configuration.
    box : array_like, optional
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx, ly, lz, alpha, beta, gamma]``.
    method : {'bruteforce', 'nsgrid', 'pkdtree'}, optional
        Keyword to override the automatic guessing of the employed search
        method.
    return_distances : bool, optional
        If set to ``True``, distances will also be returned.

    Returns
    -------
    pairs : numpy.ndarray (``dtype=numpy.int64``, ``shape=(n_pairs, 2)``)
        Pairs of indices, corresponding to coordinates in the `reference` and
        `configuration` arrays such that the distance between them lies within
        the interval (`min_cutoff`, `max_cutoff`].
        Each row in `pairs` is an index pair ``[i, j]`` corresponding to the
        ``i``-th coordinate in `reference` and the ``j``-th coordinate in
        `configuration`.
    distances : numpy.ndarray (``dtype=numpy.float64``, ``shape=(n_pairs,)``), optional
        Distances corresponding to each pair of indices. Only returned if
        `return_distances` is ``True``. ``distances[k]`` corresponds to the
        ``k``-th pair returned in `pairs` and gives the distance between the
        coordinates ``reference[pairs[k, 0]]`` and
        ``configuration[pairs[k, 1]]``.

        .. code-block:: python

            pairs, distances = capped_distances(reference, configuration,
                                                max_cutoff, return_distances=True)
            for k, [i, j] in enumerate(pairs):
                coord1 = reference[i]
                coord2 = configuration[j]
                distance = distances[k]

    Note
    -----
    Currently supports brute force, grid-based, and periodic KDtree search
    methods.

    See Also
    --------
    distance_array
    MDAnalysis.lib.pkdtree.PeriodicKDTree.search
    MDAnalysis.lib.nsgrid.FastNS.search
    """
    if box is not None:
        box = np.asarray(box, dtype=np.float32)
        if box.shape[0] != 6:
            raise ValueError("Box Argument is of incompatible type. The "
                             "dimension should be either None or of the form "
                             "[lx, ly, lz, alpha, beta, gamma]")
    method = _determine_method(reference, configuration, max_cutoff,
                               min_cutoff=min_cutoff, box=box, method=method)
    return method(reference, configuration, max_cutoff, min_cutoff=min_cutoff,
                  box=box, return_distances=return_distances)


def _determine_method(reference, configuration, max_cutoff, min_cutoff=None,
                      box=None, method=None):
    """Guesses the fastest method for capped distance calculations based on the
    size of the coordinate sets and the relative size of the target volume.

    Parameters
    ----------
    reference : numpy.ndarray
        Reference coordinate array with shape ``(3,)`` or ``(n, 3)``.
    configuration : numpy.ndarray
        Configuration coordinate array with shape ``(3,)`` or ``(m, 3)``.
    max_cutoff : float
        Maximum cutoff distance between `reference` and `configuration`
        coordinates.
    min_cutoff : float, optional
        Minimum cutoff distance between `reference` and `configuration`
        coordinates.
    box : numpy.ndarray
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx, ly, lz, alpha, beta, gamma]``.
    method : {'bruteforce', 'nsgrid', 'pkdtree'}, optional
        Keyword to override the automatic guessing of the employed search
        method.

    Returns
    -------
    function : callable
        The function implementing the guessed (or deliberatly chosen) method.
    """
    methods = {'bruteforce': _bruteforce_capped,
               'pkdtree': _pkdtree_capped,
               'nsgrid': _nsgrid_capped}

    if method is not None:
        return methods[method.lower()]

    if len(reference) < 10 or len(configuration) < 10:
        return methods['bruteforce']
    elif len(reference) * len(configuration) >= 1e8:
        # CAUTION : for large datasets, shouldnt go into 'bruteforce'
        # in any case. Arbitrary number, but can be characterized
        return methods['nsgrid']
    else:
        if box is None:
            min_dim = np.array([reference.min(axis=0),
                                configuration.min(axis=0)])
            max_dim = np.array([reference.max(axis=0),
                                configuration.max(axis=0)])
            size = max_dim.max(axis=0) - min_dim.min(axis=0)
        elif np.all(box[3:] == 90.0):
            size = box[:3]
        else:
            tribox = triclinic_vectors(box)
            size = tribox.max(axis=0) - tribox.min(axis=0)
        if np.any(max_cutoff > 0.3*size):
            return methods['bruteforce']
        else:
            return methods['nsgrid']


@check_coords('reference', 'configuration', enforce_copy=False,
              reduce_result_if_single=False, check_lengths_match=False)
def _bruteforce_capped(reference, configuration, max_cutoff, min_cutoff=None,
                       box=None, return_distances=True):
    """Capped distance evaluations using a brute force method.

    Computes and returns an array containing pairs of indices corresponding to
    entries in the `reference` and `configuration` arrays which are separated by
    a distance lying within the specified cutoff(s). Employs naive distance
    computations (brute force) to find relevant distances.

    Optionally, these distances can be returned as well.

    If the optional argument `box` is supplied, the minimum image convention is
    applied when calculating distances. Either orthogonal or triclinic boxes are
    supported.

    Parameters
    ----------
    reference : numpy.ndarray
        Reference coordinate array with shape ``(3,)`` or ``(n, 3)`` (dtype will
        be converted to ``numpy.float32`` internally).
    configuration : array
        Configuration coordinate array with shape ``(3,)`` or ``(m, 3)`` (dtype
        will be converted to ``numpy.float32`` internally).
    max_cutoff : float
        Maximum cutoff distance between `reference` and `configuration`
        coordinates.
    min_cutoff : float, optional
        Minimum cutoff distance between `reference` and `configuration`
        coordinates.
    box : numpy.ndarray, optional
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx, ly, lz, alpha, beta, gamma]``.
    return_distances : bool, optional
        If set to ``True``, distances will also be returned.

    Returns
    -------
    pairs : numpy.ndarray (``dtype=numpy.int64``, ``shape=(n_pairs, 2)``)
        Pairs of indices, corresponding to coordinates in the `reference` and
        `configuration` arrays such that the distance between them lies within
        the interval (`min_cutoff`, `max_cutoff`].
        Each row in `pairs` is an index pair ``[i, j]`` corresponding to the
        ``i``-th coordinate in `reference` and the ``j``-th coordinate in
        `configuration`.
    distances : numpy.ndarray (``dtype=numpy.float64``, ``shape=(n_pairs,)``), optional
        Distances corresponding to each pair of indices. Only returned if
        `return_distances` is ``True``. ``distances[k]`` corresponds to the
        ``k``-th pair returned in `pairs` and gives the distance between the
        coordinates ``reference[pairs[k, 0]]`` and
        ``configuration[pairs[k, 1]]``.
    """
    # Default return values (will be overwritten only if pairs are found):
    pairs = np.empty((0, 2), dtype=np.intp)
    distances = np.empty((0,), dtype=np.float64)

    if len(reference) > 0 and len(configuration) > 0:
        _distances = distance_array(reference, configuration, box=box)
        if min_cutoff is not None:
            mask = np.where((_distances <= max_cutoff) & \
                            (_distances > min_cutoff))
        else:
            mask = np.where((_distances <= max_cutoff))
        if mask[0].size > 0:
            pairs = np.c_[mask[0], mask[1]]
            if return_distances:
                distances = _distances[mask]

    if return_distances:
        return pairs, distances
    else:
        return pairs


@check_coords('reference', 'configuration', enforce_copy=False,
              reduce_result_if_single=False, check_lengths_match=False)
def _pkdtree_capped(reference, configuration, max_cutoff, min_cutoff=None,
                    box=None, return_distances=True):
    """Capped distance evaluations using a KDtree method.

    Computes and returns an array containing pairs of indices corresponding to
    entries in the `reference` and `configuration` arrays which are separated by
    a distance lying within the specified cutoff(s). Employs a (periodic) KDtree
    algorithm to find relevant distances.

    Optionally, these distances can be returned as well.

    If the optional argument `box` is supplied, the minimum image convention is
    applied when calculating distances. Either orthogonal or triclinic boxes are
    supported.

    Parameters
    ----------
    reference : numpy.ndarray
        Reference coordinate array with shape ``(3,)`` or ``(n, 3)`` (dtype will
        be converted to ``numpy.float32`` internally).
    configuration : array
        Configuration coordinate array with shape ``(3,)`` or ``(m, 3)`` (dtype
        will be converted to ``numpy.float32`` internally).
    max_cutoff : float
        Maximum cutoff distance between `reference` and `configuration`
        coordinates.
    min_cutoff : float, optional
        Minimum cutoff distance between `reference` and `configuration`
        coordinates.
    box : numpy.ndarray, optional
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx, ly, lz, alpha, beta, gamma]``.
    return_distances : bool, optional
        If set to ``True``, distances will also be returned.

    Returns
    -------
    pairs : numpy.ndarray (``dtype=numpy.int64``, ``shape=(n_pairs, 2)``)
        Pairs of indices, corresponding to coordinates in the `reference` and
        `configuration` arrays such that the distance between them lies within
        the interval (`min_cutoff`, `max_cutoff`].
        Each row in `pairs` is an index pair ``[i, j]`` corresponding to the
        ``i``-th coordinate in `reference` and the ``j``-th coordinate in
        `configuration`.
    distances : numpy.ndarray (``dtype=numpy.float64``, ``shape=(n_pairs,)``), optional
        Distances corresponding to each pair of indices. Only returned if
        `return_distances` is ``True``. ``distances[k]`` corresponds to the
        ``k``-th pair returned in `pairs` and gives the distance between the
        coordinates ``reference[pairs[k, 0]]`` and
        ``configuration[pairs[k, 1]]``.
    """
    from .pkdtree import PeriodicKDTree  # must be here to avoid circular import

    # Default return values (will be overwritten only if pairs are found):
    pairs = np.empty((0, 2), dtype=np.intp)
    distances = np.empty((0,), dtype=np.float64)

    if len(reference) > 0 and len(configuration) > 0:
        kdtree = PeriodicKDTree(box=box)
        cut = max_cutoff if box is not None else None
        kdtree.set_coords(configuration, cutoff=cut)
        _pairs = kdtree.search_tree(reference, max_cutoff)
        if _pairs.size > 0:
            pairs = _pairs
            if (return_distances or (min_cutoff is not None)):
                refA, refB = pairs[:, 0], pairs[:, 1]
                distances = calc_bonds(reference[refA], configuration[refB],
                                       box=box)
                if min_cutoff is not None:
                    mask = np.where(distances > min_cutoff)
                    pairs, distances = pairs[mask], distances[mask]

    if return_distances:
        return pairs, distances
    else:
        return pairs


@check_coords('reference', 'configuration', enforce_copy=False,
              reduce_result_if_single=False, check_lengths_match=False)
def _nsgrid_capped(reference, configuration, max_cutoff, min_cutoff=None,
                   box=None, return_distances=True):
    """Capped distance evaluations using a grid-based search method.

    Computes and returns an array containing pairs of indices corresponding to
    entries in the `reference` and `configuration` arrays which are separated by
    a distance lying within the specified cutoff(s). Employs a grid-based search
    algorithm to find relevant distances.

    Optionally, these distances can be returned as well.

    If the optional argument `box` is supplied, the minimum image convention is
    applied when calculating distances. Either orthogonal or triclinic boxes are
    supported.

    Parameters
    ----------
    reference : numpy.ndarray
        Reference coordinate array with shape ``(3,)`` or ``(n, 3)`` (dtype will
        be converted to ``numpy.float32`` internally).
    configuration : array
        Configuration coordinate array with shape ``(3,)`` or ``(m, 3)`` (dtype
        will be converted to ``numpy.float32`` internally).
    max_cutoff : float
        Maximum cutoff distance between `reference` and `configuration`
        coordinates.
    min_cutoff : float, optional
        Minimum cutoff distance between `reference` and `configuration`
        coordinates.
    box : numpy.ndarray (``dtype=numpy.float64``, ``shape=(n_pairs,)``), optional
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx, ly, lz, alpha, beta, gamma]``.
    return_distances : bool, optional
        If set to ``True``, distances will also be returned.

    Returns
    -------
    pairs : numpy.ndarray (``dtype=numpy.int64``, ``shape=(n_pairs, 2)``)
        Pairs of indices, corresponding to coordinates in the `reference` and
        `configuration` arrays such that the distance between them lies within
        the interval (`min_cutoff`, `max_cutoff`].
        Each row in `pairs` is an index pair ``[i, j]`` corresponding to the
        ``i``-th coordinate in `reference` and the ``j``-th coordinate in
        `configuration`.
    distances : numpy.ndarray, optional
        Distances corresponding to each pair of indices. Only returned if
        `return_distances` is ``True``. ``distances[k]`` corresponds to the
        ``k``-th pair returned in `pairs` and gives the distance between the
        coordinates ``reference[pairs[k, 0]]`` and
        ``configuration[pairs[k, 1]]``.
    """
    # Default return values (will be overwritten only if pairs are found):
    pairs = np.empty((0, 2), dtype=np.intp)
    distances = np.empty((0,), dtype=np.float64)

    if len(reference) > 0 and len(configuration) > 0:
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
            # to avoid failures for very close particles but with
            # larger cutoff
            boxsize = np.maximum(boxsize, 2 * max_cutoff)
            pseudobox[:3] = boxsize + 2.2*max_cutoff
            pseudobox[3:] = 90.
            shiftref, shiftconf = reference.copy(), configuration.copy()
            # Extra padding near the origin
            shiftref -= lmin - 0.1*max_cutoff
            shiftconf -= lmin - 0.1*max_cutoff
            gridsearch = FastNS(max_cutoff, shiftconf, box=pseudobox, pbc=False)
            results = gridsearch.search(shiftref)
        else:
            gridsearch = FastNS(max_cutoff, configuration, box=box)
            results = gridsearch.search(reference)

        pairs = results.get_pairs()
        if return_distances or (min_cutoff is not None):
            distances = results.get_pair_distances()
            if min_cutoff is not None:
                idx = distances > min_cutoff
                pairs, distances = pairs[idx], distances[idx]

    if return_distances:
        return pairs, distances
    else:
        return pairs


def self_capped_distance(reference, max_cutoff, min_cutoff=None, box=None,
                         method=None, return_distances=True):
    """Calculates pairs of indices corresponding to entries in the `reference`
    array which are separated by a distance lying within the specified
    cutoff(s). Optionally, these distances can be returned as well.

    If the optional argument `box` is supplied, the minimum image convention is
    applied when calculating distances. Either orthogonal or triclinic boxes are
    supported.

    An automatic guessing of the optimal method to calculate the distances is
    included in the function. An optional keyword for the method is also
    provided. Users can enforce a particular method with this functionality.
    Currently brute force, grid search, and periodic KDtree methods are
    implemented.

    Parameters
    -----------
    reference : numpy.ndarray
        Reference coordinate array with shape ``(3,)`` or ``(n, 3)``.
    max_cutoff : float
        Maximum cutoff distance between `reference` coordinates.
    min_cutoff : float, optional
        Minimum cutoff distance between `reference` coordinates.
    box : array_like, optional
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx, ly, lz, alpha, beta, gamma]``.
    method : {'bruteforce', 'nsgrid', 'pkdtree'}, optional
        Keyword to override the automatic guessing of the employed search
        method.
    return_distances : bool, optional
        If set to ``True``, distances will also be returned.

    Returns
    -------
    pairs : numpy.ndarray (``dtype=numpy.int64``, ``shape=(n_pairs, 2)``)
        Pairs of indices, corresponding to coordinates in the `reference` array
        such that the distance between them lies within the interval
        (`min_cutoff`, `max_cutoff`].
        Each row in `pairs` is an index pair ``[i, j]`` corresponding to the
        ``i``-th and the ``j``-th coordinate in `reference`.
    distances : numpy.ndarray (``dtype=numpy.float64``, ``shape=(n_pairs,)``)
        Distances corresponding to each pair of indices. Only returned if
        `return_distances` is ``True``. ``distances[k]`` corresponds to the
        ``k``-th pair returned in `pairs` and gives the distance between the
        coordinates ``reference[pairs[k, 0]]`` and ``reference[pairs[k, 1]]``.

        .. code-block:: python

            pairs, distances = self_capped_distances(reference, max_cutoff,
                                                     return_distances=True)
            for k, [i, j] in enumerate(pairs):
                coord1 = reference[i]
                coord2 = reference[j]
                distance = distances[k]


    Note
    -----
    Currently supports brute force, grid-based, and periodic KDtree search
    methods.

    See Also
    --------
    self_distance_array
    MDAnalysis.lib.pkdtree.PeriodicKDTree.search
    MDAnalysis.lib.nsgrid.FastNS.self_search

    .. versionchanged:: 0.20.0
       Added `return_distances` keyword.
    """
    if box is not None:
        box = np.asarray(box, dtype=np.float32)
        if box.shape[0] != 6:
            raise ValueError("Box Argument is of incompatible type. The "
                             "dimension should be either None or of the form "
                             "[lx, ly, lz, alpha, beta, gamma]")
    method = _determine_method_self(reference, max_cutoff,
                                    min_cutoff=min_cutoff,
                                    box=box, method=method)
    return method(reference,  max_cutoff, min_cutoff=min_cutoff, box=box,
                  return_distances=return_distances)


def _determine_method_self(reference, max_cutoff, min_cutoff=None, box=None,
                           method=None):
    """Guesses the fastest method for capped distance calculations based on the
    size of the `reference` coordinate set and the relative size of the target
    volume.

    Parameters
    ----------
    reference : numpy.ndarray
        Reference coordinate array with shape ``(3,)`` or ``(n, 3)``.
    max_cutoff : float
        Maximum cutoff distance between `reference` coordinates.
    min_cutoff : float, optional
        Minimum cutoff distance between `reference` coordinates.
    box : numpy.ndarray
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx, ly, lz, alpha, beta, gamma]``.
    method : {'bruteforce', 'nsgrid', 'pkdtree'}, optional
        Keyword to override the automatic guessing of the employed search
        method.

    Returns
    -------
    function : callable
        The function implementing the guessed (or deliberatly chosen) method.
    """
    methods = {'bruteforce': _bruteforce_capped_self,
               'pkdtree': _pkdtree_capped_self,
               'nsgrid': _nsgrid_capped_self}

    if method is not None:
        return methods[method.lower()]

    if len(reference) < 100:
        return methods['bruteforce']

    if box is None:
        min_dim = np.array([reference.min(axis=0)])
        max_dim = np.array([reference.max(axis=0)])
        size = max_dim.max(axis=0) - min_dim.min(axis=0)
    elif np.all(box[3:] == 90.0):
        size = box[:3]
    else:
        tribox = triclinic_vectors(box)
        size = tribox.max(axis=0) - tribox.min(axis=0)

    if max_cutoff < 0.03*size.min():
        return methods['pkdtree']
    else:
        return methods['nsgrid']


@check_coords('reference', enforce_copy=False, reduce_result_if_single=False)
def _bruteforce_capped_self(reference, max_cutoff, min_cutoff=None, box=None,
                            return_distances=True):
    """Capped distance evaluations using a brute force method.

    Computes and returns an array containing pairs of indices corresponding to
    entries in the `reference` array which are separated by a distance lying
    within the specified cutoff(s). Employs naive distance computations (brute
    force) to find relevant distances. Optionally, these distances can be
    returned as well.

    If the optional argument `box` is supplied, the minimum image convention is
    applied when calculating distances. Either orthogonal or triclinic boxes are
    supported.

    Parameters
    ----------
    reference : numpy.ndarray
        Reference coordinate array with shape ``(3,)`` or ``(n, 3)`` (dtype will
        be converted to ``numpy.float32`` internally).
    max_cutoff : float
        Maximum cutoff distance between `reference` coordinates.
    min_cutoff : float, optional
        Minimum cutoff distance between `reference` coordinates.
    box : numpy.ndarray, optional
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx, ly, lz, alpha, beta, gamma]``.
    return_distances : bool, optional
        If set to ``True``, distances will also be returned.

    Returns
    -------
    pairs : numpy.ndarray (``dtype=numpy.int64``, ``shape=(n_pairs, 2)``)
        Pairs of indices, corresponding to coordinates in the `reference` array
        such that the distance between them lies within the interval
        (`min_cutoff`, `max_cutoff`].
        Each row in `pairs` is an index pair ``[i, j]`` corresponding to the
        ``i``-th and the ``j``-th coordinate in `reference`.
    distances : numpy.ndarray (``dtype=numpy.float64``, ``shape=(n_pairs,)``), optional
        Distances corresponding to each pair of indices. Only returned if
        `return_distances` is ``True``. ``distances[k]`` corresponds to the
        ``k``-th pair returned in `pairs` and gives the distance between the
        coordinates ``reference[pairs[k, 0]]`` and
        ``configuration[pairs[k, 1]]``.

    .. versionchanged:: 0.20.0
       Added `return_distances` keyword.
    """
    # Default return values (will be overwritten only if pairs are found):
    pairs = np.empty((0, 2), dtype=np.intp)
    distances = np.empty((0,), dtype=np.float64)

    N = len(reference)
    # We're searching within a single coordinate set, so we need at least two
    # coordinates to find distances between them.
    if N > 1:
        distvec = self_distance_array(reference, box=box)
        dist = np.full((N, N), np.finfo(np.float64).max, dtype=np.float64)
        dist[np.triu_indices(N, 1)] = distvec

        if min_cutoff is not None:
            mask = np.where((dist <= max_cutoff) & (dist > min_cutoff))
        else:
            mask = np.where((dist <= max_cutoff))

        if mask[0].size > 0:
            pairs = np.c_[mask[0], mask[1]]
            distances = dist[mask]
    if return_distances:
        return pairs, distances
    return pairs


@check_coords('reference', enforce_copy=False, reduce_result_if_single=False)
def _pkdtree_capped_self(reference, max_cutoff, min_cutoff=None, box=None,
                         return_distances=True):
    """Capped distance evaluations using a KDtree method.

    Computes and returns an array containing pairs of indices corresponding to
    entries in the `reference` array which are separated by a distance lying
    within the specified cutoff(s). Employs a (periodic) KDtree algorithm to
    find relevant distances. Optionally, these distances can be returned as
    well.

    If the optional argument `box` is supplied, the minimum image convention is
    applied when calculating distances. Either orthogonal or triclinic boxes are
    supported.

    Parameters
    ----------
    reference : numpy.ndarray
        Reference coordinate array with shape ``(3,)`` or ``(n, 3)`` (dtype will
        be converted to ``numpy.float32`` internally).
    max_cutoff : float
        Maximum cutoff distance between `reference` coordinates.
    min_cutoff : float, optional
        Minimum cutoff distance between `reference` coordinates.
    box : numpy.ndarray, optional
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx, ly, lz, alpha, beta, gamma]``.
    return_distances : bool, optional
        If set to ``True``, distances will also be returned.

    Returns
    -------
    pairs : numpy.ndarray (``dtype=numpy.int64``, ``shape=(n_pairs, 2)``)
        Pairs of indices, corresponding to coordinates in the `reference` array
        such that the distance between them lies within the interval
        (`min_cutoff`, `max_cutoff`].
        Each row in `pairs` is an index pair ``[i, j]`` corresponding to the
        ``i``-th and the ``j``-th coordinate in `reference`.
    distances : numpy.ndarray (``dtype=numpy.float64``, ``shape=(n_pairs,)``)
        Distances corresponding to each pair of indices. Only returned if
        `return_distances` is ``True``. ``distances[k]`` corresponds to the
        ``k``-th pair returned in `pairs` and gives the distance between
        the coordinates ``reference[pairs[k, 0]]`` and
        ``reference[pairs[k, 1]]``.

    .. versionchanged:: 0.20.0
       Added `return_distances` keyword.
    """
    from .pkdtree import PeriodicKDTree  # must be here to avoid circular import

    # Default return values (will be overwritten only if pairs are found):
    pairs = np.empty((0, 2), dtype=np.intp)
    distances = np.empty((0,), dtype=np.float64)

    # We're searching within a single coordinate set, so we need at least two
    # coordinates to find distances between them.
    if len(reference) > 1:
        kdtree = PeriodicKDTree(box=box)
        cut = max_cutoff if box is not None else None
        kdtree.set_coords(reference, cutoff=cut)
        _pairs = kdtree.search_pairs(max_cutoff)
        if _pairs.size > 0:
            pairs = _pairs
            if (return_distances or (min_cutoff is not None)):
                refA, refB = pairs[:, 0], pairs[:, 1]
                distances = calc_bonds(reference[refA], reference[refB], box=box)
                if min_cutoff is not None:
                    idx = distances > min_cutoff
                    pairs, distances = pairs[idx], distances[idx]
    if return_distances:
        return pairs, distances
    return pairs


@check_coords('reference', enforce_copy=False, reduce_result_if_single=False)
def _nsgrid_capped_self(reference, max_cutoff, min_cutoff=None, box=None,
                        return_distances=True):
    """Capped distance evaluations using a grid-based search method.

    Computes and returns an array containing pairs of indices corresponding to
    entries in the `reference` array which are separated by a distance lying
    within the specified cutoff(s). Employs a grid-based search algorithm to
    find relevant distances. Optionally, these distances can be returned as
    well.

    If the optional argument `box` is supplied, the minimum image convention is
    applied when calculating distances. Either orthogonal or triclinic boxes are
    supported.

    Parameters
    ----------
    reference : numpy.ndarray
        Reference coordinate array with shape ``(3,)`` or ``(n, 3)`` (dtype will
        be converted to ``numpy.float32`` internally).
    max_cutoff : float
        Maximum cutoff distance between `reference` coordinates.
    min_cutoff : float, optional
        Minimum cutoff distance between `reference` coordinates.
    box : numpy.ndarray, optional
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx, ly, lz, alpha, beta, gamma]``.

    Returns
    -------
    pairs : numpy.ndarray (``dtype=numpy.int64``, ``shape=(n_pairs, 2)``)
        Pairs of indices, corresponding to coordinates in the `reference` array
        such that the distance between them lies within the interval
        (`min_cutoff`, `max_cutoff`].
        Each row in `pairs` is an index pair ``[i, j]`` corresponding to the
        ``i``-th and the ``j``-th coordinate in `reference`.
    distances : numpy.ndarray, optional
        Distances corresponding to each pair of indices. Only returned if
        `return_distances` is ``True``. ``distances[k]`` corresponds to the
        ``k``-th pair returned in `pairs` and gives the distance between the
        coordinates ``reference[pairs[k, 0]]`` and
        ``configuration[pairs[k, 1]]``.

    .. versionchanged:: 0.20.0
       Added `return_distances` keyword.
    """
    # Default return values (will be overwritten only if pairs are found):
    pairs = np.empty((0, 2), dtype=np.intp)
    distances = np.empty((0,), dtype=np.float64)

    # We're searching within a single coordinate set, so we need at least two
    # coordinates to find distances between them.
    if len(reference) > 1:
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
        if return_distances or (min_cutoff is not None):
            distances = results.get_pair_distances()[::2]
            if min_cutoff is not None:
                idx = distances > min_cutoff
                pairs, distances = pairs[idx], distances[idx]

    if return_distances:
        return pairs, distances
    return pairs


@check_coords('coords')
def transform_RtoS(coords, box, backend="serial"):
    """Transform an array of coordinates from real space to S space (a.k.a.
    lambda space)

    S space represents fractional space within the unit cell for this system.

    Reciprocal operation to :meth:`transform_StoR`.

    Parameters
    ----------
    coords : numpy.ndarray
        A ``(3,)`` or ``(n, 3)`` array of coordinates (dtype is arbitrary, will
        be converted to ``numpy.float32`` internally).
    box : numpy.ndarray
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx, ly, lz, alpha, beta, gamma]``.
    backend : {'serial', 'OpenMP'}, optional
        Keyword selecting the type of acceleration.

    Returns
    -------
    newcoords : numpy.ndarray (``dtype=numpy.float32``, ``shape=coords.shape``)
        An array containing fractional coordiantes.


    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    .. versionchanged:: 0.19.0
       Internal dtype conversion of input coordinates to ``numpy.float32``.
       Now also accepts (and, likewise, returns) a single coordinate.
    """
    if len(coords) == 0:
        return coords
    boxtype, box = check_box(box)
    if boxtype == 'ortho':
        box = np.diag(box)
    box = box.astype(np.float64)

    # Create inverse matrix of box
    # need order C here
    inv = np.array(np.linalg.inv(box), order='C')

    _run("coord_transform", args=(coords, inv), backend=backend)

    return coords


@check_coords('coords')
def transform_StoR(coords, box, backend="serial"):
    """Transform an array of coordinates from S space into real space.

    S space represents fractional space within the unit cell for this system.

    Reciprocal operation to :meth:`transform_RtoS`

    Parameters
    ----------
    coords : numpy.ndarray
        A ``(3,)`` or ``(n, 3)`` array of coordinates (dtype is arbitrary, will
        be converted to ``numpy.float32`` internally).
    box : numpy.ndarray
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx, ly, lz, alpha, beta, gamma]``.
    backend : {'serial', 'OpenMP'}, optional
        Keyword selecting the type of acceleration.

    Returns
    -------
    newcoords : numpy.ndarray (``dtype=numpy.float32``, ``shape=coords.shape``)
        An array containing real space coordiantes.


    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    .. versionchanged:: 0.19.0
       Internal dtype conversion of input coordinates to ``numpy.float32``.
       Now also accepts (and, likewise, returns) a single coordinate.
    """
    if len(coords) == 0:
        return coords
    boxtype, box = check_box(box)
    if boxtype == 'ortho':
        box = np.diag(box)
    box = box.astype(np.float64)

    _run("coord_transform", args=(coords, box), backend=backend)
    return coords


@check_coords('coords1', 'coords2')
def calc_bonds(coords1, coords2, box=None, result=None, backend="serial"):
    """Calculates the bond lengths between pairs of atom positions from the two
    coordinate arrays `coords1` and `coords2`, which must contain the same
    number of coordinates. ``coords1[i]`` and ``coords2[i]`` represent the
    positions of atoms connected by the ``i``-th bond. If single coordinates are
    supplied, a single distance will be returned.

    In comparison to :meth:`distance_array` and :meth:`self_distance_array`,
    which calculate distances between all possible combinations of coordinates,
    :meth:`calc_bonds` only calculates distances between pairs of coordinates,
    similar to::

       numpy.linalg.norm(a - b) for a, b in zip(coords1, coords2)

    If the optional argument `box` is supplied, the minimum image convention is
    applied when calculating distances. Either orthogonal or triclinic boxes are
    supported.

    If a numpy array of dtype ``numpy.float64`` with shape ``(n,)`` (for ``n``
    coordinate pairs) is provided in `result`, then this preallocated array is
    filled. This can speed up calculations.

    Parameters
    ----------
    coords1 : numpy.ndarray
        Coordinate array of shape ``(3,)`` or ``(n, 3)`` for one half of a
        single or ``n`` bonds, respectively (dtype is arbitrary, will be
        converted to ``numpy.float32`` internally).
    coords2 : numpy.ndarray
        Coordinate array of shape ``(3,)`` or ``(n, 3)`` for the other half of
        a single or ``n`` bonds, respectively (dtype is arbitrary, will be
        converted to ``numpy.float32`` internally).
    box : numpy.ndarray, optional
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx, ly, lz, alpha, beta, gamma]``.
    result : numpy.ndarray, optional
        Preallocated result array of dtype ``numpy.float64`` and shape ``(n,)``
        (for ``n`` coordinate pairs). Avoids recreating the array in repeated
        function calls.
    backend : {'serial', 'OpenMP'}, optional
        Keyword selecting the type of acceleration.

    Returns
    -------
    bondlengths : numpy.ndarray (``dtype=numpy.float64``, ``shape=(n,)``) or numpy.float64
        Array containing the bond lengths between each pair of coordinates. If
        two single coordinates were supplied, their distance is returned as a
        single number instead of an array.


    .. versionadded:: 0.8
    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    .. versionchanged:: 0.19.0
       Internal dtype conversion of input coordinates to ``numpy.float32``.
       Now also accepts single coordinates as input.
    """
    numatom = coords1.shape[0]
    bondlengths = _check_result_array(result, (numatom,))

    if numatom > 0:
        if box is not None:
            boxtype, box = check_box(box)
            if boxtype == 'ortho':
                _run("calc_bond_distance_ortho",
                     args=(coords1, coords2, box, bondlengths),
                     backend=backend)
            else:
                _run("calc_bond_distance_triclinic",
                     args=(coords1, coords2, box, bondlengths),
                     backend=backend)
        else:
            _run("calc_bond_distance",
                 args=(coords1, coords2, bondlengths),
                 backend=backend)

    return bondlengths


@check_coords('coords1', 'coords2', 'coords3')
def calc_angles(coords1, coords2, coords3, box=None, result=None,
                backend="serial"):
    """Calculates the angles formed between triplets of atom positions from the
    three coordinate arrays `coords1`, `coords2`, and `coords3`. All coordinate
    arrays must contain the same number of coordinates.

    The coordinates in `coords2` represent the apices of the angles::

            2---3
           /
          1

    Configurations where the angle is undefined (e.g., when coordinates 1 or 3
    of a triplet coincide with coordinate 2) result in a value of **zero** for
    that angle.

    If the optional argument `box` is supplied, periodic boundaries are taken
    into account when constructing the connecting vectors between coordinates,
    i.e., the minimum image convention is applied for the vectors forming the
    angles. Either orthogonal or triclinic boxes are supported.

    If a numpy array of dtype ``numpy.float64`` with shape ``(n,)`` (for ``n``
    coordinate triplets) is provided in `result`, then this preallocated array
    is filled. This can speed up calculations.

    Parameters
    ----------
    coords1 : numpy.ndarray
        Array of shape ``(3,)`` or ``(n, 3)`` containing the coordinates of one
        side of a single or ``n`` angles, respectively (dtype is arbitrary, will
        be converted to ``numpy.float32`` internally)
    coords2 : numpy.ndarray
        Array of shape ``(3,)`` or ``(n, 3)`` containing the coordinates of the
        apices of a single or ``n`` angles, respectively (dtype is arbitrary,
        will be converted to ``numpy.float32`` internally)
    coords3 : numpy.ndarray
        Array of shape ``(3,)`` or ``(n, 3)`` containing the coordinates of the
        other side of a single or ``n`` angles, respectively (dtype is
        arbitrary, will be converted to ``numpy.float32`` internally)
    box : numpy.ndarray, optional
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx, ly, lz, alpha, beta, gamma]``.
    result : numpy.ndarray, optional
        Preallocated result array of dtype ``numpy.float64`` and shape ``(n,)``
        (for ``n`` coordinate triplets). Avoids recreating the array in repeated
        function calls.
    backend : {'serial', 'OpenMP'}, optional
        Keyword selecting the type of acceleration.

    Returns
    -------
    angles : numpy.ndarray (``dtype=numpy.float64``, ``shape=(n,)``) or numpy.float64
        Array containing the angles between each triplet of coordinates. Values
        are returned in radians (rad). If three single coordinates were
        supplied, the angle is returned as a single number instead of an array.


    .. versionadded:: 0.8
    .. versionchanged:: 0.9.0
       Added optional box argument to account for periodic boundaries in
       calculation
    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    .. versionchanged:: 0.19.0
       Internal dtype conversion of input coordinates to ``numpy.float32``.
       Now also accepts single coordinates as input.
    """
    numatom = coords1.shape[0]
    angles = _check_result_array(result, (numatom,))

    if numatom > 0:
        if box is not None:
            boxtype, box = check_box(box)
            if boxtype == 'ortho':
                _run("calc_angle_ortho",
                       args=(coords1, coords2, coords3, box, angles),
                       backend=backend)
            else:
                _run("calc_angle_triclinic",
                       args=(coords1, coords2, coords3, box, angles),
                       backend=backend)
        else:
            _run("calc_angle",
                   args=(coords1, coords2, coords3, angles),
                   backend=backend)

    return angles


@check_coords('coords1', 'coords2', 'coords3', 'coords4')
def calc_dihedrals(coords1, coords2, coords3, coords4, box=None, result=None,
                   backend="serial"):
    r"""Calculates the dihedral angles formed between quadruplets of positions
    from the four coordinate arrays `coords1`, `coords2`, `coords3`, and
    `coords4`, which must contain the same number of coordinates.

    The dihedral angle formed by a quadruplet of positions (1,2,3,4) is
    calculated around the axis connecting positions 2 and 3 (i.e., the angle
    between the planes spanned by positions (1,2,3) and (2,3,4))::

                  4
                  |
            2-----3
           /
          1

    If all coordinates lie in the same plane, the cis configuration corresponds
    to a dihedral angle of zero, and the trans configuration to :math:`\pi`
    radians (180 degrees). Configurations where the dihedral angle is undefined
    (e.g., when all coordinates lie on the same straight line) result in a value
    of ``nan`` (not a number) for that dihedral.

    If the optional argument `box` is supplied, periodic boundaries are taken
    into account when constructing the connecting vectors between coordinates,
    i.e., the minimum image convention is applied for the vectors forming the
    dihedral angles. Either orthogonal or triclinic boxes are supported.

    If a numpy array of dtype ``numpy.float64`` with shape ``(n,)`` (for ``n``
    coordinate quadruplets) is provided in `result` then this preallocated array
    is filled. This can speed up calculations.

    Parameters
    ----------
    coords1 : numpy.ndarray
        Coordinate array of shape ``(3,)`` or ``(n, 3)`` containing the 1st
        positions in dihedrals (dtype is arbitrary, will be converted to
        ``numpy.float32`` internally)
    coords2 : numpy.ndarray
        Coordinate array of shape ``(3,)`` or ``(n, 3)`` containing the 2nd
        positions in dihedrals (dtype is arbitrary, will be converted to
        ``numpy.float32`` internally)
    coords3 : numpy.ndarray
        Coordinate array of shape ``(3,)`` or ``(n, 3)`` containing the 3rd
        positions in dihedrals (dtype is arbitrary, will be converted to
        ``numpy.float32`` internally)
    coords4 : numpy.ndarray
        Coordinate array of shape ``(3,)`` or ``(n, 3)`` containing the 4th
        positions in dihedrals (dtype is arbitrary, will be converted to
        ``numpy.float32`` internally)
    box : numpy.ndarray, optional
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx, ly, lz, alpha, beta, gamma]``.
    result : numpy.ndarray, optional
        Preallocated result array of dtype ``numpy.float64`` and shape ``(n,)``
        (for ``n`` coordinate quadruplets). Avoids recreating the array in
        repeated function calls.
    backend : {'serial', 'OpenMP'}, optional
        Keyword selecting the type of acceleration.

    Returns
    -------
    dihedrals : numpy.ndarray (``dtype=numpy.float64``, ``shape=(n,)``) or numpy.float64
        Array containing the dihedral angles formed by each quadruplet of
        coordinates. Values are returned in radians (rad). If four single
        coordinates were supplied, the dihedral angle is returned as a single
        number instead of an array.


    .. versionadded:: 0.8
    .. versionchanged:: 0.9.0
       Added optional box argument to account for periodic boundaries in
       calculation
    .. versionchanged:: 0.11.0
       Renamed from calc_torsions to calc_dihedrals
    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    .. versionchanged:: 0.19.0
       Internal dtype conversion of input coordinates to ``numpy.float32``.
       Now also accepts single coordinates as input.
    """
    numatom = coords1.shape[0]
    dihedrals = _check_result_array(result, (numatom,))

    if numatom > 0:
        if box is not None:
            boxtype, box = check_box(box)
            if boxtype == 'ortho':
                _run("calc_dihedral_ortho",
                     args=(coords1, coords2, coords3, coords4, box, dihedrals),
                     backend=backend)
            else:
                _run("calc_dihedral_triclinic",
                     args=(coords1, coords2, coords3, coords4, box, dihedrals),
                     backend=backend)
        else:
            _run("calc_dihedral",
                 args=(coords1, coords2, coords3, coords4, dihedrals),
                 backend=backend)

    return dihedrals


@check_coords('coords')
def apply_PBC(coords, box, backend="serial"):
    """Moves coordinates into the primary unit cell.

    Parameters
    ----------
    coords : numpy.ndarray
        Coordinate array of shape ``(3,)`` or ``(n, 3)`` (dtype is arbitrary,
        will be converted to ``numpy.float32`` internally).
    box : numpy.ndarray
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
        ``[lx, ly, lz, alpha, beta, gamma]``.
    backend : {'serial', 'OpenMP'}, optional
        Keyword selecting the type of acceleration.

    Returns
    -------
    newcoords : numpy.ndarray  (``dtype=numpy.float32``, ``shape=coords.shape``)
        Array containing coordinates that all lie within the primary unit cell
        as defined by `box`.


    .. versionadded:: 0.8
    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    .. versionchanged:: 0.19.0
       Internal dtype conversion of input coordinates to ``numpy.float32``.
       Now also accepts (and, likewise, returns) single coordinates.
    """
    if len(coords) == 0:
        return coords
    boxtype, box = check_box(box)
    if boxtype == 'ortho':
        _run("ortho_pbc", args=(coords, box), backend=backend)
    else:
        _run("triclinic_pbc", args=(coords, box), backend=backend)

    return coords
