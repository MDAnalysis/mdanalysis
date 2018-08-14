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


def _check_box(box):
    """Take a box input and deduce what type of system it represents based on
    the shape of the array and whether all angles are 90 degrees.

    Parameters
    ----------
    box : numpy.ndarray
        Array with dtype ``numpy.float32`` containing box information of unknown
        format.

    Returns
    -------
    boxtype : str
        * ``'ortho'`` orthogonal box
        * ``'tri_vecs'`` triclinic box vectors

    checked_box : numpy.ndarray
        Array of dtype ``numpy.float32`` containing box information:
        * If `boxtype` is ``'ortho'``, `cecked_box` will have the shape ``(3,)``
          containing the x-, y-, and z-dimensions of the orthogonal box.
        * If  `boxtype` is ``'tri_vecs'``, `cecked_box` will have the shape
          ``(3, 3)`` containing the triclinic box vectors in a lower triangular
          matrix as returned by
          :meth:`~MDAnalysis.lib.mdamath.triclinic_vectors`.

    Raises
    ------
    TypeError
        If the dtype of `box` is not ``numpy.float32``.
    ValueError
        If box type cannot be detected.


    .. seealso: :meth:`~MDAnalysis.lib.mdamath.triclinic_vectors`
    .. versionchanged: 0.19.0
       * Added automatic conversion of input to :class:`numpy.ndarray` with
         dtype ``numpy.float32``.
       * Now also returns the box in the format expected by low-level functions
         in :mod:`~MDAnalysis.lib.c_distances`.
       * Removed obsolete box types ``tri_box`` and ``tri_vecs_bad``.
    """
    if box.dtype != np.float32:
        raise TypeError("Invalid box dtype. Must be numpy.float32, got {0}."
                        "".format(box.dtype))
    if box.shape == (6,):
        if np.all(box[3:] == 90.):
            return 'ortho', box
        else:
            return 'tri_vecs', triclinic_vectors(box)
    elif box.shape == (3,):
        return 'ortho', box
    elif box.shape == (3, 3):
        # Check if triclinic box vectors are properly formatted:
        if np.all([box[0][1] == 0.0, box[0][2] == 0.0, box[1][2] == 0.0]):
            return 'tri_vecs', box
        else:
            return 'tri_vecs', triclinic_vectors(triclinic_box(box[0], box[1],
                                                              box[2]))
    raise ValueError("Box input not recognized, must be an array of box "
                     "dimensions.")


def _check_coord_array(coords, desc):
    """Check if an array is a valid coordinate array.

    The `coords` array must meet the following requirements:
      * Must have a shape of ``(n, 3)``.
      * Its dtype must be ``numpy.float32``.

    Parameters
    ----------
    coords : numpy.ndarray
        The coordinate array to check for validity.
    desc : str
        Name of the coordinate array. Only used in error messages.

    Raises
    ------
    ValueError
        If `coords` has a wrong shape.
    TypeError
        If the dtype of `coords` is not ``numpy.float32``.
    """
    _check_coord_array_shape_2d(coords, desc)
    if coords.dtype != np.float32:
        raise TypeError("{0} must be of type numpy.float32".format(desc))
# The following two lines would break a lot of tests. WHY?!
#    if not coords.flags['C_CONTIGUOUS']:
#        raise ValueError("{0} is not C-contiguous.".format(desc))


def _check_coord_array_shape_1d(coord, desc):
    """Check whether a single coordinate array has the correct shape of
    ``(3,)``.
    """
    if coord.shape != (3,):
        raise ValueError("{0} must have a shape of (3,), got {1}."
                         "".format(desc, coord.shape))


def _check_coord_array_shape_2d(coords, desc):
    """Check whether a coordinate array has the correct shape of ``(n, 3)``."""
    if (coords.ndim != 2 or coords.shape[1] != 3):
        raise ValueError("{0} must have a shape of (n, 3), got {1}."
                         "".format(desc, coords.shape))


def _check_coord_array_shape_1d_or_2d(coords, desc):
    """Check if a coordinate array has the correct shape of either ``(3,)`` or
    ``(n, 3)``.
    """
    if (coords.ndim not in (1, 2) or coords.shape[-1] != 3):
        raise ValueError("{0} must have a shape of either (3,) or (n, 3), "
                         "got {1}.".format(desc, coords.shape))


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
    result : numpy.ndarray
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
        raise TypeError("Result array must be of type numpy.float64, got {}"
                        "".format(result.dtype))
# The following two lines would break a lot of tests. WHY?!
#    if not coords.flags['C_CONTIGUOUS']:
#        raise ValueError("{0} is not C-contiguous.".format(desc))
    return result

def _check_lengths_match(*arrays):
    """Check all arrays are same shape"""
    ref = arrays[0].shape

    if not all( a.shape == ref for a in arrays):
        raise ValueError("Input arrays must all have the same shape, got {0}."
                         "".format([a.shape for a in arrays]))

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
    _check_coord_array_shape_2d(configuration, 'configuration')
    _check_coord_array_shape_2d(reference, 'reference')

    ref = reference.astype(np.float32, order='C', copy=False)
    conf = configuration.astype(np.float32, order='C', copy=False)

    confnum = conf.shape[0]
    refnum = ref.shape[0]

    distances = _check_result_array(result, (refnum, confnum))

    if box is not None:
        boxtype, box = _check_box(box)
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
    _check_coord_array_shape_2d(reference, 'reference')

    ref = reference.astype(np.float32, order='C', copy=False)

    refnum = ref.shape[0]
    distnum = refnum * (refnum - 1) // 2

    distances = _check_result_array(result, (distnum,))

    if box is not None:
        boxtype, box = _check_box(box)
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
                    box=None, method=None, return_distances=True):
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
    return_distances : (optional) 'True'/'False'
        Function returns distances if set to 'True'

    Returns
    -------
    pairs : array
        Pair of indices, one from each reference and configuration such that
        distance between them is  within the ``max_cutoff`` and ``min_cutoff``
        pairs[i,j] contains the indices i from reference coordinates, and
        j from configuration
    distances : (optional) array
        Distances corresponding to each pair of indices.
        d[k] corresponding to the pairs[i,j] gives the distance between
        i-th and j-th coordinate in reference and configuration respectively
        Returns only if ``return_distances`` is set to `True`

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

    if return_distances:
        pairs, dist = method(reference, configuration, max_cutoff,
                         min_cutoff=min_cutoff, box=box,
                         return_distances=return_distances)
        return np.asarray(pairs), np.asarray(dist)
    else:
        pairs = method(reference, configuration, max_cutoff,
                         min_cutoff=min_cutoff, box=box,
                         return_distances=return_distances)

        return np.asarray(pairs)


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

    if len(reference) < 10 or len(configuration) < 10:
        return methods['bruteforce']
    elif len(reference)*len(configuration) >= 1e8:
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
        elif np.allclose(box[3:], 90):
            size = box[:3]
        else:
            tribox = triclinic_vectors(box)
            size = tribox.max(axis=0) - tribox.min(axis=0)
        if np.any(max_cutoff > 0.3*size):
            return methods['bruteforce']
        else:
            return methods['nsgrid']


def _bruteforce_capped(reference, configuration, max_cutoff,
                       min_cutoff=None, box=None, return_distances=True):
    """Internal method for bruteforce calculations

    Uses naive distance calulations and returns a list
    containing the indices with one from each
    reference and configuration arrays, such that the distance between
    them is less than the specified cutoff distance

    Returns
    -------
    pairs : array
        Every ``[i, j]`` pair is such that atom-index ``i`` is
        from reference and ``j`` from configuration array
    distance: (optional) array
         Distance between ``reference[i]`` and ``configuration[j]``
         atom coordinate. Only returns when ``return_distances``
         is set to ``True``

    """
    pairs, distance = [], []
    reference = np.asarray(reference, dtype=np.float32)
    configuration = np.asarray(configuration, dtype=np.float32)

    if reference.shape == (3, ):
        reference = reference[None, :]
    if configuration.shape == (3, ):
        configuration = configuration[None, :]

    _check_coord_array(reference, 'reference')
    _check_coord_array(configuration, 'configuration')

    distance = distance_array(reference, configuration, box=box)
    if min_cutoff is not None:
        mask = np.where((distance <= max_cutoff) & (distance > min_cutoff))
    else:
        mask = np.where((distance <= max_cutoff))

    if mask[0].size > 0:
        pairs = np.c_[mask[0], mask[1]]

    if return_distances:
        distance = distance[mask]
        return pairs, distance
    else:
        return pairs


def _pkdtree_capped(reference, configuration, max_cutoff,
                    min_cutoff=None, box=None, return_distances=True):
    """ Capped Distance evaluations using KDtree.

    Uses minimum image convention if *box* is specified

    Returns:
    --------
    pairs : array
        Array of atom indices which are within the specified cutoff distance.
        Every pairs `(i, j)` corresponds to i-th particle in reference and
        j-th particle in configuration
    distance : (optional) array
        Distance between two atoms corresponding to the (i, j) indices
        in pairs. Only returns when ``return_distances``
         is set to ``True``

    """
    from .pkdtree import PeriodicKDTree

    reference = np.asarray(reference, dtype=np.float32)
    configuration = np.asarray(configuration, dtype=np.float32)

    if reference.shape == (3, ):
        reference = reference[None, :]
    if configuration.shape == (3, ):
        configuration = configuration[None, :]

    _check_coord_array(reference, 'reference')
    _check_coord_array(configuration, 'configuration')
    kdtree = PeriodicKDTree(box=box)
    cut = max_cutoff if box is not None else None
    kdtree.set_coords(configuration, cutoff=cut)
    pairs = kdtree.search_tree(reference, max_cutoff)
    if (return_distances or (min_cutoff is not None)) and pairs.size > 0:
        refA, refB = pairs[:, 0], pairs[:, 1]
        distance = calc_bonds(reference[refA], configuration[refB], box=box)
        if min_cutoff is not None:
            mask = np.where(distance > min_cutoff)
            pairs, distance = pairs[mask], distance[mask]
    else:
        distance = []

    if return_distances:
        return pairs, distance
    else:
        return pairs


def _nsgrid_capped(reference, configuration, max_cutoff, min_cutoff=None,
                   box=None, return_distances=True):
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

    Returns
    -------
    pairs : array
        Array of atom indices which are within the specified cutoff distance.
        Every pairs `(i, j)` corresponds to i-th particle in reference and
        j-th particle in configuration
    distance : (optional) array
        Distance between two atoms corresponding to the (i, j) indices
        in pairs. Only returns when ``return_distances``
         is set to ``True``
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
    if return_distances or (min_cutoff is not None):
        pair_distance = results.get_pair_distances()
        if min_cutoff is not None:
            idx = pair_distance > min_cutoff
            pairs, pair_distance = pairs[idx], pair_distance[idx]

    if return_distances:
        return pairs, pair_distance
    else:
        return pairs


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

    if box is None:
        min_dim = np.array([reference.min(axis=0)])
        max_dim = np.array([reference.max(axis=0)])
        size = max_dim.max(axis=0) - min_dim.min(axis=0)
    elif np.allclose(box[3:], 90):
        size = box[:3]
    else:
        tribox = triclinic_vectors(box)
        size = tribox.max(axis=0) - tribox.min(axis=0)

    if len(reference) < 100:
        return methods['bruteforce']
    elif max_cutoff < 0.03*size.min():
        return methods['pkdtree']
    else:
        return methods['nsgrid']


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

    N = len(reference)
    distvec = np.zeros((N*(N-1)//2), dtype=np.float64)
    self_distance_array(reference, box=box, result=distvec)

    distance = np.ones((N, N), dtype=np.float32)*max_cutoff
    distance[np.triu_indices(N, 1)] = distvec

    if min_cutoff is not None:
        mask = np.where((distance < max_cutoff) & (distance > min_cutoff))
    else:
        mask = np.where((distance < max_cutoff))

    if mask[0].size > 0:
        pairs = np.c_[mask[0], mask[1]]
        distance = distance[mask]
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
        :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:\n
        ``[lx, ly, lz, alpha, beta, gamma]``.
    backend : str, optional
        Select the type of acceleration; ``'serial'`` is always available.
        Another possibility is ``'OpenMP'``.

    Returns
    -------
    newcoords : numpy.ndarray
        An array of dtype ``numpy.float32`` with the same shape as `coords`
        containing fractional coordiantes.


    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    .. versionchanged:: 0.19.0
       Internal dtype conversion of input coordinates to ``numpy.float32``.
    """
    _check_coord_array_shape_1d_or_2d(coords, 'coords')

    # Here, copy=True is necessary to not modify the input coordinates:
    newcoords = coords.astype(np.float32, order='C', copy=True)

    is_1d = False  # True if only one vector coordinate
    if newcoords.ndim == 1:
        newcoords = newcoords[None, :]
        is_1d = True

    boxtype, box = _check_box(box)
    if boxtype == 'ortho':
        box = np.array([[box[0], 0.0, 0.0],
                        [0.0, box[1], 0.0],
                        [0.0, 0.0, box[2]]], dtype=np.float32)

    # Create inverse matrix of box
    # need order C here
    inv = np.array(np.matrix(box).I, dtype=np.float32, order='C')

    _run("coord_transform", args=(newcoords, inv), backend=backend)

    if is_1d:
        return newcoords[0]
    return newcoords


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
        :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:\n
        ``[lx, ly, lz, alpha, beta, gamma]``.
    backend : str, optional
        Select the type of acceleration; ``'serial'`` is always available.
        Another possibility is ``'OpenMP'``.

    Returns
    -------
    newcoords : numpy.ndarray
        An array of dtype ``numpy.float32`` with the same shape as `coords`
        containing real space coordiantes.


    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    .. versionchanged:: 0.19.0
       Internal dtype conversion of input coordinates to ``numpy.float32``.
    """
    _check_coord_array_shape_1d_or_2d(coords, 'coords')

    # Here, copy=True is necessary to not modify the input coordinates:
    newcoords = coords.astype(np.float32, order='C', copy=True)

    is_1d = False  # True if only one vector coordinate
    if len(newcoords.shape) == 1:
        newcoords = newcoords[None, :]
        is_1d = True

    boxtype, box = _check_box(box)
    if boxtype == 'ortho':
        box = np.array([[box[0], 0.0, 0.0],
                        [0.0, box[1], 0.0],
                        [0.0, 0.0, box[2]]], dtype=np.float32)

    _run("coord_transform", args=(newcoords, box), backend=backend)

    if is_1d:
        return newcoords[0]
    return newcoords


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
    _check_coord_array_shape_2d(coords1, 'coords1')
    _check_coord_array_shape_2d(coords2, 'coords2')
    _check_lengths_match(coords1, coords2)

    atom1 = coords1.astype(np.float32, order='C', copy=True)
    atom2 = coords2.astype(np.float32, order='C', copy=True)
    numatom = atom1.shape[0]

    bondlengths = _check_result_array(result, (numatom,))

    if box is not None:
        boxtype, box = _check_box(box)
        if boxtype == 'ortho':
            _run("calc_bond_distance_ortho",
                 args=(atom1, atom2, box, bondlengths),
                 backend=backend)
        else:
            _run("calc_bond_distance_triclinic",
                 args=(atom1, atom2, box, bondlengths),
                 backend=backend)
    else:
        _run("calc_bond_distance",
             args=(atom1, atom2, bondlengths),
             backend=backend)

    return bondlengths


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
    _check_coord_array_shape_2d(coords1, 'coords1')
    _check_coord_array_shape_2d(coords2, 'coords2')
    _check_coord_array_shape_2d(coords3, 'coords3')
    _check_lengths_match(coords1, coords2, coords3)

    atom1 = coords1.astype(np.float32, order='C', copy=True)
    atom2 = coords2.astype(np.float32, order='C', copy=True)
    atom3 = coords3.astype(np.float32, order='C', copy=True)
    numatom = atom1.shape[0]

    angles = _check_result_array(result, (numatom,))

    if box is not None:
        boxtype, box = _check_box(box)
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
    _check_coord_array_shape_2d(coords1, 'coords1')
    _check_coord_array_shape_2d(coords2, 'coords2')
    _check_coord_array_shape_2d(coords3, 'coords3')
    _check_coord_array_shape_2d(coords4, 'coords4')
    _check_lengths_match(coords1, coords2, coords3, coords4)

    atom1 = coords1.astype(np.float32, order='C', copy=True)
    atom2 = coords2.astype(np.float32, order='C', copy=True)
    atom3 = coords3.astype(np.float32, order='C', copy=True)
    atom4 = coords4.astype(np.float32, order='C', copy=True)

    numatom = atom1.shape[0]

    dihedrals = _check_result_array(result, (numatom,))

    if box is not None:
        boxtype, box = _check_box(box)
        if boxtype == 'ortho':
            _run("calc_dihedral_ortho",
                 args=(atom1, atom2, atom3, atom4, box, dihedrals),
                 backend=backend)
        else:
            _run("calc_dihedral_triclinic",
                 args=(atom1, atom2, atom3, atom4, box, dihedrals),
                 backend=backend)
    else:
        _run("calc_dihedral",
             args=(atom1, atom2, atom3, atom4, dihedrals),
             backend=backend)

    return dihedrals


def apply_PBC(coords, box, backend="serial"):
    """Moves coordinates into the primary unit cell.

    Parameters
    ----------
    coords : numpy.ndarray
        Coordinate array of shape ``(n, 3)`` (dtype is arbitrary, will be
        converted to ``numpy.float32`` internally).
    box : numpy.ndarray
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:\n
        ``[lx, ly, lz, alpha, beta, gamma]``.
    backend : str, optional
        Select the type of acceleration; ``'serial'`` is always available.
        Another possibility is ``'OpenMP'``.

    Returns
    -------
    newcoords : numpy.ndarray
        Array of dtype ``numpy.float32`` containing coordinates that all lie
        within the primary unit cell as defined by `box`.


    .. versionadded:: 0.8
    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    .. versionchanged:: 0.19.0
       Internal dtype conversion of input coordinates to ``numpy.float32``.
    """
    _check_coord_array_shape_2d(coords, 'coords')

    newcoords = coords.astype(np.float32, order='C', copy=True)

    boxtype, box = _check_box(box)
    if boxtype == 'ortho':
        box_inv = 1.0 / box
        _run("ortho_pbc", args=(newcoords, box, box_inv), backend=backend)
    else:
        box_inv = 1.0 / np.diagonal(box)
        _run("triclinic_pbc", args=(newcoords, box, box_inv), backend=backend)

    return newcoords


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
    _check_coord_array_shape_1d(a, 'a')
    _check_coord_array_shape_1d(b, 'b')

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
    _check_coord_array_shape_1d(a, 'a')
    _check_coord_array_shape_1d(b, 'b')
    _check_coord_array_shape_1d(c, 'c')

    return np.rad2deg(
        calc_angles(a[None, :], b[None, :], c[None, :], box=box)[0])


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
    _check_coord_array_shape_1d(a, 'a')
    _check_coord_array_shape_1d(b, 'b')
    _check_coord_array_shape_1d(c, 'c')
    _check_coord_array_shape_1d(d, 'd')

    return np.rad2deg(
        calc_dihedrals(a[None, :], b[None, :], c[None, :], d[None, :], box)[0])
