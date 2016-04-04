# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2014 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see AUTHORS for the full list)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

"""Fast distance array computation --- :mod:`MDAnalysis.lib.distances`
===================================================================

Fast C-routines to calculate distance arrays from coordinate
arrays. Many of the functions also exist in parallel versions, that
typically provide higher performance than the serial code.

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
.. autofunction:: transform_RtoS(coordinates, box [, backend])
.. autofunction:: transform_StoR(coordinates, box [,backend])

"""
from six.moves import range

import numpy as np
from numpy.lib.utils import deprecate

from .mdamath import triclinic_vectors, triclinic_box


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


def _box_check(box):
    """Take a box input and deduce what type of system it represents based
    on the shape of the array and whether all angles are 90.

    :Arguments:
      *box*
          box information of unknown format

    :Returns:
      * ``ortho`` orthogonal box
      * ``tri_vecs`` triclinic box vectors
      * ``tri_box`` triclinic box lengths and angles

    :Raises:
      * ``TypeError`` if box is not float32
      * ``ValueError`` if box type not detected
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

    d = distance_array(reference, configuration[,box[,result=d]])

    :Arguments:
        *reference*
             reference coordinate array (must be numpy.float32)
        *configuration*
             configuration coordinate array (must be numpy.float32)
        *box*
             cell dimensions (minimum image convention is applied)
             or None [``None``].
             Cell dimensions must be in an identical to format to those returned
             by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`,
             [lx, ly, lz, alpha, beta, gamma]
        *result*
             optional preallocated result array which must have the
             shape (len(ref), len(conf)) and dtype=numpy.float64.
             Avoids creating the array which saves time when the function
             is called repeatedly. [``None``]
        *backend*
             select the type of acceleration; "serial" is always available. Other
             possibilities are "OpenMP" (OpenMP).
    :Returns:
         *d*
             (len(reference),len(configuration)) numpy array with the distances d[i,j]
             between reference coordinates i and configuration coordinates j

    .. Note:: This method is slower than it could be because internally we need to
              make copies of the ref and conf arrays.

    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    """
    ref = reference.copy('C')
    conf = configuration.copy('C')

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

    :Arguments:
        *ref*
             reference coordinate array with N=len(ref) coordinates
        *box*
             cell dimensions (minimum image convention is applied)
             or None [``None``]
             Cell dimensions must be in an identical to format to those returned
             by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`,
             [lx, ly, lz, alpha, beta, gamma]
        *result*
             optional preallocated result array which must have the shape
             (N*(N-1)/2,) and dtype ``numpy.float64``. Avoids creating
             the array which saves time when the function is called
             repeatedly. [``None``]
        *backend*
             select the type of acceleration; "serial" is always available. Other
             possibilities are "OpenMP" (OpenMP).
    :Returns:
        *d*
             N*(N-1)/2 numpy 1D array with the distances dist[i,j] between ref
             coordinates i and j at position d[k]. Loop through d::

                 for i in range(N):
                     for j in range(i+1, N):
                         k += 1
                         dist[i,j] = d[k]

    .. Note:: This method is slower than it could be because internally we need to
              make copies of the coordinate arrays.

    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    """
    ref = reference.copy('C')

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
    distnum = refnum * (refnum - 1) / 2

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


def transform_RtoS(inputcoords, box, backend="serial"):
    """Transform an array of coordinates from real space to S space (aka lambda space)

    S space represents fractional space within the unit cell for this system

    Reciprocal operation to :meth:`transform_StoR`

    :Arguments:
      *inputcoords*
          An n x 3 array of coordinate data, of type np.float32
      *box*
          The unitcell dimesions for this system
          Cell dimensions must be in an identical to format to those returned
          by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`,
          [lx, ly, lz, alpha, beta, gamma]
      *backend*
          select the type of acceleration; "serial" is always available. Other
          possibilities are "OpenMP" (OpenMP).

    :Returns:
       *outcoords*
          An n x 3 array of fractional coordiantes

    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    """
    coords = inputcoords.copy('C')
    numcoords = coords.shape[0]

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

    return coords


def transform_StoR(inputcoords, box, backend="serial"):
    """Transform an array of coordinates from S space into real space.

    S space represents fractional space within the unit cell for this system

    Reciprocal operation to :meth:`transform_RtoS`

    :Arguments:
      *inputcoords*
           An n x 3 array of coordinate data, of type np.float32
      *box*
           The unitcell dimesions for this system
           Cell dimensions must be in an identical to format to those returned
           by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`,
           [lx, ly, lz, alpha, beta, gamma]
      *backend*
          select the type of acceleration; "serial" is always available. Other
          possibilities are "OpenMP" (OpenMP).

    :Returns:
       *outcoords*
            An n x 3 array of fracional coordiantes

    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    """
    coords = inputcoords.copy('C')
    numcoords = coords.shape[0]

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

    :Arguments:
       *coords1*
          An array of coordinates for one half of the bond
       *coords2*
          An array of coordinates for the other half of bond
       *box*
          Unit cell information if periodic boundary conditions are required [None]
          Cell dimensions must be in an identical to format to those returned
          by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`,
          [lx, ly, lz, alpha, beta, gamma]
       *result*
          optional preallocated result array which must be same length as coord
          arrays and dtype=numpy.float64. Avoids creating the
          array which saves time when the function is called repeatedly. [None]
       *backend*
          select the type of acceleration; "serial" is always available. Other
          possibilities are "OpenMP" (OpenMP).

    :Returns:
       *bondlengths*
          numpy array with the length between each pair in coords1 and coords2

    .. versionadded:: 0.8
    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    """
    atom1 = coords1.copy('C')
    atom2 = coords2.copy('C')

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

    :Arguments:
        *coords1*
            coordinate array of one side of angles
        *coords2*
            coordinate array of apex of angles
        *coords3*
            coordinate array of other side of angles
        *box*
            optional unit cell information.  This ensures that the connecting vectors between
            atoms respect minimum image convention.  This is import when the angle might
            be between atoms in different images.
            Cell dimensions must be in an identical to format to those returned
            by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`,
            [lx, ly, lz, alpha, beta, gamma]
        *result*
            optional preallocated results array which must have same length as coordinate
            array and dtype=numpy.float64.
        *backend*
            select the type of acceleration; "serial" is always available. Other
            possibilities are "OpenMP" (OpenMP).

    :Returns:
        *angles*
            A numpy.array of angles in radians

    .. versionadded:: 0.8
    .. versionchanged:: 0.9.0
       Added optional box argument to account for periodic boundaries in calculation
    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    """
    atom1 = coords1.copy('C')
    atom2 = coords2.copy('C')
    atom3 = coords3.copy('C')
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

    If a 1D numpy array of dtype ``numpy.float64`` with ``len(atom1)`` elements is
    provided in *result* then this preallocated array is filled. This can speed
    up calculations.

    The optional argument ``box`` ensures that periodic boundaries are taken into
    account when constructing the connecting vectors between atoms, ie that the vector
    between atoms 1 & 2 goes between coordinates in the same image::

      angles = calc_dihedrals(coords1, coords2, coords3, coords4 [,box=box, result=angles])

    :Arguments:
        *coords1*
            coordinate array of 1st atom in dihedrals
        *coords2*
            coordinate array of 2nd atom in dihedrals
        *coords3*
            coordinate array of 3rd atom in dihedrals
        *coords4*
            coordinate array of 4th atom in dihedrals
        *box*
            optional unit cell information.  This ensures that the connecting vectors
            between atoms respect minimum image convention.  This is import when the
            angle might be between atoms in different images.
            Cell dimensions must be in an identical to format to those returned
            by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`,
            [lx, ly, lz, alpha, beta, gamma]
        *result*
            optional preallocated results array which must have same length as
            coordinate array and dtype=numpy.float64.
        *backend*
            select the type of acceleration; "serial" is always available. Other
            possibilities are "OpenMP" (OpenMP).

    :Returns:
        *angles*
            A numpy.array of angles in radians

    .. versionadded:: 0.8
    .. versionchanged:: 0.9.0
       Added optional box argument to account for periodic boundaries in calculation
    .. versionchanged:: 0.11.0
       Renamed from calc_torsions to calc_dihedrals
    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    """
    atom1 = coords1.copy('C')
    atom2 = coords2.copy('C')
    atom3 = coords3.copy('C')
    atom4 = coords4.copy('C')

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

calc_torsions = deprecate(calc_dihedrals,
                          old_name='calc_torsions',
                          new_name='calc_dihedrals',
                          message="This will be removed in v0.15.0")


def apply_PBC(incoords, box, backend="serial"):
    """Moves a set of coordinates to all be within the primary unit cell

    newcoords = apply_PBC(coords, box)

    :Arguments:
        *coords*
           coordinate array (of type numpy.float32)
        *box*
           box dimensions, can be either orthogonal or triclinic information
           Cell dimensions must be in an identical to format to those returned
           by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`,
           [lx, ly, lz, alpha, beta, gamma]
        *backend*
           select the type of acceleration; "serial" is always available. Other
           possibilities are "OpenMP" (OpenMP).

    :Returns:
        *newcoords*
           coordinates that are now all within the primary unit cell,
           as defined by box

    .. versionadded:: 0.8
    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    """
    coords = incoords.copy('C')

    _check_array(coords, 'coords')

    coordnum = coords.shape[0]

    # determine boxtype
    boxtype = _box_check(box)
    # Convert [A,B,C,alpha,beta,gamma] to [[A],[B],[C]]
    if (boxtype == 'tri_box'):
        box = triclinic_vectors(box)
    if (boxtype == 'tri_vecs_bad'):
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

applyPBC = deprecate(apply_PBC,
                     old_name='applyPBC',
                     new_name='apply_PBC',
                     message="This will be removed in v0.15.0")
