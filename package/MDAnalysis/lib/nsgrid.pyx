# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
#
# Copyright (C) 2013-2018  Sébastien Buchoux <sebastien.buchoux@gmail.com>
# Copyright (c) 2018 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v3 or any higher version
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


# cython: cdivision=True
# cython: boundscheck=False
# cython: initializedcheck=False
# cython: embedsignature=True

"""
Neighbor search library --- :mod:`MDAnalysis.lib.nsgrid`
========================================================


About the code
--------------

This Neighbor search library is a serialized Cython version greatly
inspired by the NS grid search implemented in
`GROMACS <http://www.gromacs.org/>`_ .

GROMACS 4.x code (more precisely
`nsgrid.c <https://github.com/gromacs/gromacs/commits/master/src/mdlib/nsgrid.c>`_
and `ns.c <https://github.com/gromacs/gromacs/commits/master/src/mdlib/ns.c>`_ )
was used as reference to write this file.

GROMACS 4.x code is released under the GNU Public Licence v2.

About the algorithm
-------------------

The neighbor search implemented here is based on
`cell lists <https://en.wikipedia.org/wiki/Cell_lists>`_ which allow
computation of pairs [#]_ with a cost of :math:`O(N)`, instead
of :math:`O(N^2)`. The basic algorithm is described in
Appendix F,  Page 552 of
``Understanding Molecular Dynamics: From Algorithm to Applications`` by Frenkel and Smit.

In brief, the algorithm divides the domain into smaller subdomains called `cells`
and distributes every particle to these cells based on their positions. Subsequently,
any distance based query first identifies the corresponding cell position in the
domain followed by distance evaluations within the identified cell and
neighboring cells only. Care must be taken to ensure that `cellsize` is
greater than the desired search distance, otherwise all of the neighbours might
not reflect in the results.


.. [#] a pair correspond to two particles that are considered as neighbors .


.. versionadded:: 0.19.0
"""

from MDAnalysis.lib.distances import _check_array
# Used to handle memory allocation
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from libc.math cimport sqrt
import numpy as np
from libcpp.vector cimport vector
cimport numpy as np

# Preprocessor DEFs
DEF DIM = 3
DEF XX = 0
DEF YY = 1
DEF ZZ = 2
DEF EPSILON = 1e-5

ctypedef np.int_t ns_int
ctypedef np.float32_t real
ctypedef real rvec[DIM]
ctypedef ns_int ivec[DIM]
ctypedef real matrix[DIM][DIM]

ctypedef vector[ns_int] intvec
ctypedef vector[real] realvec

# Useful Functions
cdef real rvec_norm2(const rvec a) nogil:
    return a[XX]*a[XX] + a[YY]*a[YY] + a[ZZ]*a[ZZ]

cdef void rvec_clear(rvec a) nogil:
    a[XX] = 0.0
    a[YY] = 0.0
    a[ZZ] = 0.0

###############################
# Utility class to handle PBC #
###############################
cdef struct cPBCBox_t:
    matrix     box
    rvec       fbox_diag
    rvec       hbox_diag
    rvec       mhbox_diag
    real       max_cutoff2


# Class to handle PBC calculations
cdef class PBCBox(object):
    """
    Cython implementation of
    `PBC-related <https://en.wikipedia.org/wiki/Periodic_boundary_conditions>`_
    operations. This class is used by classes :class:`FastNS`
    and :class:`NSGrid` to put all particles inside a brick-shaped box
    and to compute PBC-aware distance. The class can also handle
    non-PBC aware distance evaluations through ``periodic`` argument.

    .. warning::
        This class is not meant to be used by end users.

    .. warning::
        Even if MD triclinic boxes can be handled by this class,
        internal optimization is made based on the assumption that
        particles are inside a brick-shaped box. When this is not
        the case, calculated distances are not
        warranted to be exact.
    """

    cdef cPBCBox_t c_pbcbox
    cdef bint is_triclinic
    cdef bint periodic

    def __init__(self, real[:, ::1] box, bint periodic):
        """
        Parameters
        ----------
        box : numpy.ndarray
            box vectors of shape ``(3, 3)`` or
            as returned by ``MDAnalysis.lib.mdamath.triclinic_vectors``
            ``dtype`` must be ``numpy.float32``
        periodic : boolean
            ``True`` for PBC-aware calculations
            ``False`` for non PBC aware calculations
        """

        self.periodic = periodic
        self.update(box)

    cdef void fast_update(self, real[:, ::1] box) nogil:
        """
        Updates the internal box parameters for
        PBC-aware distance calculations. The internal
        box parameters are used to define the brick-shaped
        box which is eventually used for distance calculations.

        """
        cdef ns_int i, j
        cdef real min_hv2, min_ss, tmp

        # Update matrix
        self.is_triclinic = False
        for i in range(DIM):
            for j in range(DIM):
                self.c_pbcbox.box[i][j] = box[i, j]

                if i != j:
                    if box[i, j] > EPSILON:
                        self.is_triclinic = True

        # Update diagonals
        for i in range(DIM):
            self.c_pbcbox.fbox_diag[i] = box[i, i]
            self.c_pbcbox.hbox_diag[i] = self.c_pbcbox.fbox_diag[i] * 0.5
            self.c_pbcbox.mhbox_diag[i] = - self.c_pbcbox.hbox_diag[i]

        # Update maximum cutoff

        # Physical limitation of the cut-off
        # by half the length of the shortest box vector.
        min_hv2 = min(0.25 * rvec_norm2(&box[XX, XX]), 0.25 * rvec_norm2(&box[YY, XX]))
        min_hv2 = min(min_hv2, 0.25 * rvec_norm2(&box[ZZ, XX]))

        # Limitation to the smallest diagonal element due to optimizations:
        # checking only linear combinations of single box-vectors (2 in x)
        # in the grid search and pbc_dx is a lot faster
        # than checking all possible combinations.
        tmp = box[YY, YY]
        if box[ZZ, YY] < 0:
            tmp -= box[ZZ, YY]
        else:
            tmp += box[ZZ, YY]

        min_ss = min(box[XX, XX], min(tmp, box[ZZ, ZZ]))
        self.c_pbcbox.max_cutoff2 = min(min_hv2, min_ss * min_ss)

    def update(self, real[:, ::1] box):
        """
        Updates internal MD box representation and parameters used for calculations.

        Parameters
        ----------
        box : numpy.ndarray
            Describes the MD box vectors as returned by
            :func:`MDAnalysis.lib.mdamath.triclinic_vectors`.
            `dtype` must be :class:`numpy.float32`

        Note
        ----
        Call to this method is only needed when the MD box is changed
        as it always called when class is instantiated.

        """

        if box.shape[0] != DIM or box.shape[1] != DIM:
            raise ValueError("Box must be a {} x {} matrix. Got: {} x {})".format(
                DIM, DIM, box.shape[0], box.shape[1]))
        if (box[XX, XX] == 0) or (box[YY, YY] == 0) or (box[ZZ, ZZ] == 0):
            raise ValueError("Box does not correspond to PBC=xyz")
        self.fast_update(box)

    cdef void fast_pbc_dx(self, rvec ref, rvec other, rvec dx) nogil:
        """Dislacement between two points for both
        PBC and non-PBC conditions

        Modifies the displacement vector between two points based
        on the minimum image convention for PBC aware calculations.

        For non-PBC aware distance evaluations, calculates the
        displacement vector without any modifications
        """

        cdef ns_int i, j

        for i in range(DIM):
            dx[i] = other[i] - ref[i]

        if self.periodic:
            for i in range(DIM-1, -1, -1):
                while dx[i] > self.c_pbcbox.hbox_diag[i]:
                    for j in range(i, -1, -1):
                        dx[j] -= self.c_pbcbox.box[i][j]

                while dx[i] <= self.c_pbcbox.mhbox_diag[i]:
                    for j in range(i, -1, -1):
                        dx[j] += self.c_pbcbox.box[i][j]

    cdef real fast_distance2(self, rvec a, rvec b) nogil:
        """Distance calculation between two points
        for both PBC and non-PBC aware calculations

        Returns the distance obeying minimum
        image convention if periodic is set to ``True`` while
        instantiating the ``PBCBox`` object.
        """

        cdef rvec dx
        self.fast_pbc_dx(a, b, dx)
        return rvec_norm2(dx)

    cdef real[:, ::1]fast_put_atoms_in_bbox(self, real[:, ::1] coords) nogil:
        """Shifts all ``coords`` to an orthogonal brick shaped box

        All the coordinates are brought into an orthogonal
        box. The box vectors for the brick-shaped box
        are defined in ``fast_update`` method.

        """

        cdef ns_int i, m, d, natoms
        cdef real[:, ::1] bbox_coords

        natoms = coords.shape[0]
        with gil:
            bbox_coords = coords.copy()

        if self.periodic:
            if self.is_triclinic:
                for i in range(natoms):
                    for m in range(DIM - 1, -1, -1):
                        while bbox_coords[i, m] < 0:
                            for d in range(m+1):
                                bbox_coords[i, d] += self.c_pbcbox.box[m][d]
                        while bbox_coords[i, m] >= self.c_pbcbox.box[m][m]:
                            for d in range(m+1):
                                bbox_coords[i, d] -= self.c_pbcbox.box[m][d]
            else:
                for i in range(natoms):
                    for m in range(DIM):
                        while bbox_coords[i, m] < 0:
                            bbox_coords[i, m] += self.c_pbcbox.box[m][m]
                        while bbox_coords[i, m] >= self.c_pbcbox.box[m][m]:
                            bbox_coords[i, m] -= self.c_pbcbox.box[m][m]

        return bbox_coords

#########################
# Neighbor Search Stuff #
#########################

cdef class NSResults(object):
    """Class to store the results

    All the required outputs from :class:`FastNS` is stored in the
    instance of this class. All the methods of :class:`FastNS` returns
    an instance of this class, which can be used to  generate the desired
    results on demand.
    """

    cdef readonly real cutoff
    cdef ns_int npairs

    cdef real[:, ::1] coords  # shape: size, DIM
    cdef real[:, ::1] searchcoords

    cdef vector[intvec] indices_buffer
    cdef vector[realvec] distances_buffer
    cdef vector[ns_int] pairs_buffer
    cdef vector[real] pair_distances_buffer
    cdef vector[real] pair_distances2_buffer

    def __init__(self, real cutoff, real[:, ::1]coords, real[:, ::1]searchcoords):
        """
        Parameters
        ----------
        cutoff : float
            Specified cutoff distance
        coords : numpy.ndarray
            Array with coordinates of atoms of shape ``(N, 3)`` for
            ``N`` particles. ``dtype`` must be ``numpy.float32``
        searchcoords : numpy.ndarray
            Array with query coordinates. Shape must be ``(M, 3)``
            for ``M`` queries. ``dtype`` must be ``numpy.float32``
        """

        self.cutoff = cutoff
        self.coords = coords
        self.searchcoords = searchcoords

        self.npairs = 0

    cdef void add_neighbors(self, ns_int beadid_i, ns_int beadid_j, real distance2) nogil:
        """Internal function to add pairs and distances to buffers

        The buffers populated using this method are used by
        other methods of this class. This is the
        primary function used by :class:`FastNS` to save all
        the pair of atoms,
        which are considered as neighbors.
        """

        self.pairs_buffer.push_back(beadid_i)
        self.pairs_buffer.push_back(beadid_j)
        self.pair_distances2_buffer.push_back(distance2)
        self.npairs += 1

    def get_pairs(self):
        """Returns all the pairs within the desired cutoff distance

        Returns an array of shape ``(N, 2)``, where N is the number of pairs
        between ``reference`` and ``configuration`` within the specified distance.
        For every pair ``(i, j)``, ``reference[i]`` and ``configuration[j]`` are
        atom positions such that ``reference`` is the position of query
        atoms while ``configuration`` coontains the position of group of
        atoms used to search against  the query atoms.

        Returns
        -------
        pairs : numpy.ndarray
            pairs of atom indices of neighbors from query
            and initial atom coordinates of shape ``(N, 2)``
        """

        return np.asarray(self.pairs_buffer).reshape(self.npairs, 2)

    def get_pair_distances(self):
        """Returns all the distances corresponding to each pair of neighbors

        Returns an array of shape ``N`` where N is the number of pairs
        among the query atoms and initial atoms within a specified distance.
        Every element ``[i]`` corresponds to the distance between
        ``pairs[i, 0]`` and ``pairs[i, 1]``, where pairs is the array
        obtained from ``get_pairs()``

        Returns
        -------
        distances : numpy.ndarray
            distances between pairs of query and initial
            atom coordinates of shape ``N``

        See Also
        --------
        MDAnalysis.lib.nsgrid.NSResults.get_pairs

        """

        self.pair_distances_buffer = np.sqrt(self.pair_distances2_buffer)
        return np.asarray(self.pair_distances_buffer)

    cdef void create_buffers(self) nogil:
        """
        Creates buffers to get individual neighbour list and distances
        of the query atoms.

        """

        cdef ns_int i, beadid_i, beadid_j
        cdef ns_int idx, nsearch
        cdef real dist2, dist

        nsearch = len(self.searchcoords)

        self.indices_buffer = vector[intvec]()
        self.distances_buffer = vector[realvec]()

        # initialize rows corresponding to search
        for i in range(nsearch):
            self.indices_buffer.push_back(intvec())
            self.distances_buffer.push_back(realvec())

        for i in range(0, 2*self.npairs, 2):
            beadid_i = self.pairs_buffer[i]
            beadid_j = self.pairs_buffer[i + 1]

            dist2 = self.pair_distances2_buffer[i//2]

            self.indices_buffer[beadid_i].push_back(beadid_j)

            dist = sqrt(dist2)

            self.distances_buffer[beadid_i].push_back(dist)

    def get_indices(self):
        """Individual neighbours of query atom

        For every queried atom ``i``, an array of all its neighbors
        indices can be obtained from ``get_indices()[i]``

        Returns
        -------
        indices : list
            Indices of neighboring atoms.
            Every element i.e. ``indices[i]`` will be a list of
            size ``m`` where m is the number of neighbours of
            query atom[i].

        .. code-block:: python

                results = NSResults()
                indices = results.get_indices()

        ``indices[i]`` will output a list of neighboring
        atoms of ``atom[i]`` from query atoms ``atom``.
        ``indices[i][j]`` will give the atom-id of initial coordinates
        such that ``initial_atom[indices[i][j]]`` is a neighbor of ``atom[i]``

        """

        if self.indices_buffer.empty():
            self.create_buffers()
        return list(self.indices_buffer)

    def get_distances(self):
        """Distance corresponding to individual neighbors of query atom

        For every queried atom ``i``, a list of all the distances
        from its neighboring atoms can be obtained from ``get_distances()[i]``.
        Every ``distance[i][j]`` will correspond
        to the distance between atoms ``atom[i]`` from the query
        atoms and ``atom[indices[j]]`` from the initialized
        set of coordinates, where ``indices`` can be obtained
        by ``get_indices()``

        Returns
        -------
        distances : np.ndarray
            Every element i.e. ``distances[i]`` will be an array of
            shape ``m`` where m is the number of neighbours of
            query atom[i].

        .. code-block:: python

                results = NSResults()
                distances = results.get_distances()


        atoms of ``atom[i]`` and query atoms ``atom``.
        ``indices[i][j]`` will give the atom-id of initial coordinates
        such that ``initial_atom[indices[i][j]]`` is a neighbor of ``atom[i]``

        See Also
        --------
        MDAnalysis.lib.nsgrid.NSResults.get_indices

        """

        if self.distances_buffer.empty():
            self.create_buffers()
        return list(self.distances_buffer)

cdef class NSGrid(object):
    """Constructs a uniform cuboidal grid for a brick-shaped box

    This class uses :class:`PBCBox` to define the brick shaped box
    It is essential to initialize the box with :class:`PBCBox`
    inorder to form the grid.

    The domain is subdivided into number of cells based on the desired search
    radius. Ideally cellsize should be equal to the search radius, but small
    search radius leads to large cell-list data strucutres.
    An optimization of cutoff is imposed to limit the size of data
    structure such that the cellsize is always greater than or
    equal to cutoff distance.

    Note
    ----
    This class assumes that all the coordinates are already
    inside the brick shaped box. Care must be taken to ensure
    all the particles are within the brick shaped box as
    defined by :class:`PBCBox`. This can be ensured by using
    :func:`~MDAnalysis.lib.nsgrid.PBCBox.fast_put_atoms_in_bbox`

    .. warning::
        This class is not meant to be used by end users.

    """

    cdef readonly real cutoff  # cutoff
    cdef ns_int size  # total cells
    cdef ns_int ncoords  # number of coordinates
    cdef ns_int[DIM] ncells  # individual cells in every dimension
    cdef ns_int[DIM] cell_offsets  # Cell Multipliers
    cdef real[DIM] cellsize  # cell size in every dimension
    cdef ns_int nbeads_per_cell  # maximum beads
    cdef ns_int *nbeads  # size (Number of beads in every cell)
    cdef ns_int *beadids  # size * nbeads_per_cell (Beadids in every cell)
    cdef ns_int *cellids  # ncoords (Cell occupation id for every atom)
    cdef bint force  # To negate the effects of optimized cutoff

    def __init__(self, ncoords, cutoff, PBCBox box, max_size, force=False):
        """
        Parameters
        ----------
        ncoords : int
            Number of coordinates to fill inside the brick shaped box
        cutoff : float
            Desired cutoff radius
        box : PBCBox
            Instance of :class:`PBCBox`
        max_size : int
            Maximum total number of cells
        force : boolean
            Optimizes cutoff if set to ``False`` [False]
        """

        cdef ns_int i
        cdef ns_int ncellx, ncelly, ncellz, size, nbeadspercell
        cdef ns_int xi, yi, zi
        cdef real bbox_vol

        self.ncoords = ncoords

        # Calculate best cutoff
        self.cutoff = cutoff
        if not force:
            bbox_vol = box.c_pbcbox.box[XX][XX] * box.c_pbcbox.box[YY][YY] * box.c_pbcbox.box[YY][YY]
            size = bbox_vol/cutoff**3
            nbeadspercell = ncoords/size
            while bbox_vol/self.cutoff**3 > max_size:
                self.cutoff *= 1.2

        for i in range(DIM):
            self.ncells[i] = <ns_int> (box.c_pbcbox.box[i][i] / self.cutoff)
            self.cellsize[i] = box.c_pbcbox.box[i][i] / self.ncells[i]
        self.size = self.ncells[XX] * self.ncells[YY] * self.ncells[ZZ]

        self.cell_offsets[XX] = 0
        self.cell_offsets[YY] = self.ncells[XX]
        self.cell_offsets[ZZ] = self.ncells[XX] * self.ncells[YY]

        # Allocate memory
        # Number of beads in every cell
        self.nbeads = <ns_int *> PyMem_Malloc(sizeof(ns_int) * self.size)
        if not self.nbeads:
            raise MemoryError("Could not allocate memory from NSGrid.nbeads ({} bits requested)".format(sizeof(ns_int) * self.size))
        self.beadids = NULL
        # Cellindex of every bead
        self.cellids = <ns_int *> PyMem_Malloc(sizeof(ns_int) * self.ncoords)
        if not self.cellids:
            raise MemoryError("Could not allocate memory from NSGrid.cellids ({} bits requested)".format(sizeof(ns_int) * self.ncoords))
        self.nbeads_per_cell = 0

        for i in range(self.size):
            self.nbeads[i] = 0

    def __dealloc__(self):
        PyMem_Free(self.nbeads)
        PyMem_Free(self.beadids)
        PyMem_Free(self.cellids)

    cdef ns_int coord2cellid(self, rvec coord) nogil:
        """Finds the cell-id for the given coordinate inside the brick shaped box

        Note
        ----
        Assumes the coordinate is already inside the brick shaped box.
        Return wrong cell-id if this is not the case
        """
        return <ns_int> (coord[ZZ] / self.cellsize[ZZ]) * (self.cell_offsets[ZZ]) +\
               <ns_int> (coord[YY] / self.cellsize[YY]) * self.cell_offsets[YY] + \
               <ns_int> (coord[XX] / self.cellsize[XX])

    cdef bint cellid2cellxyz(self, ns_int cellid, ivec cellxyz) nogil:
        """Finds actual cell position `(x, y, z)` from a cell-id
        """

        if cellid < 0:
            return False
        if cellid >= self.size:
            return False

        cellxyz[ZZ] = <ns_int> (cellid / self.cell_offsets[ZZ])
        cellid -= cellxyz[ZZ] * self.cell_offsets[ZZ]

        cellxyz[YY] = <ns_int> (cellid / self.cell_offsets[YY])
        cellxyz[XX] = cellid - cellxyz[YY] * self.cell_offsets[YY]

        return True

    cdef fill_grid(self, real[:, ::1] coords):
        """Sorts atoms into cells based on their position in the brick shaped box

        Every atom inside the brick shaped box is assigned a
        cell-id based on its position. Another list ``beadids``
        sort the atom-ids in each cell.

        Note
        ----
        The method fails if any coordinate is outside the brick shaped box.

        """

        cdef ns_int i, cellindex = -1
        cdef ns_int ncoords = coords.shape[0]
        cdef ns_int[:] beadcounts = np.empty(self.size, dtype=np.int)

        with nogil:
            # Initialize buffers
            for i in range(self.size):
                beadcounts[i] = 0

            # First loop: find cellindex for each bead
            for i in range(ncoords):
                cellindex = self.coord2cellid(&coords[i, 0])

                self.nbeads[cellindex] += 1
                self.cellids[i] = cellindex

                if self.nbeads[cellindex] > self.nbeads_per_cell:
                    self.nbeads_per_cell = self.nbeads[cellindex]

        # Allocate memory
        self.beadids = <ns_int *> PyMem_Malloc(sizeof(ns_int) * self.size * self.nbeads_per_cell)  # np.empty((self.size, nbeads_max), dtype=np.int)
        if not self.beadids:
            raise MemoryError("Could not allocate memory for NSGrid.beadids ({} bits requested)".format(sizeof(ns_int) * self.size * self.nbeads_per_cell))

        with nogil:
            # Second loop: fill grid
            for i in range(ncoords):

                # Add bead to grid cell
                cellindex = self.cellids[i]
                self.beadids[cellindex * self.nbeads_per_cell + beadcounts[cellindex]] = i
                beadcounts[cellindex] += 1


cdef class FastNS(object):
    """Grid based search between two group of atoms

    Instantiates a class object which uses :class:`PBCBox` and
    :class:`NSGrid` to construct a cuboidal
    grid in an orthogonal brick shaped box.

    Minimum image convention is used for distance evaluations
    if pbc is set to ``True``.
    """
    cdef PBCBox box
    cdef real[:, ::1] coords
    cdef real[:, ::1] coords_bbox
    cdef readonly real cutoff
    cdef NSGrid grid
    cdef ns_int max_gridsize
    cdef bint periodic

    def __init__(self, cutoff, coords, box, max_gridsize=5000, pbc=True):
        """
        Initialize the grid and sort the coordinates in respective
        cells by shifting the coordinates in a brick shaped box.
        The brick shaped box is defined by :class:`PBCBox`
        and cuboidal grid is initialize by :class:`NSGrid`.
        If box is supplied, periodic shifts along box vectors are used
        to contain all the coordinates inside the brick shaped box.
        If box is not supplied, the range of coordinates i.e.
        ``[xmax, ymax, zmax] - [xmin, ymin, zmin]`` should be used
        to construct a pseudo box. Subsequently, the origin should also be
        shifted to ``[xmin, ymin, zmin]``. These arguments must be provided
        to the function.

        Parameters
        ----------
        cutoff : float
            Desired cutoff distance
        coords : numpy.ndarray
            atom coordinates of shape ``(N, 3)`` for ``N`` atoms.
            ``dtype=numpy.float32``. For Non-PBC calculations,
            all the coords must be within the bounding box specified
            by ``box``
        box : numpy.ndarray
            Box dimension of shape (6, ). The dimensions must be
            provided in the same format as returned
            by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`:
            ``[lx, ly, lz, alpha, beta, gamma]``. For non-PBC
            evaluations, provide an orthogonal bounding box
            (dtype = numpy.float32)
        max_gridsize : int
            maximum number of cells in the grid. This parameter
            can be tuned for superior performance.
        pbc : boolean
            Handle to switch periodic boundary conditions on/off [True]

        Note
        ----
        * ``pbc=False`` Only works for orthogonal boxes.
        * Care must be taken such that all particles are inside
          the bounding box as defined by the box argument for non-PBC
          calculations.
        * In case of Non-PBC calculations, a bounding box must be provided
          to encompass all the coordinates as well as the search coordinates.
          The dimension should be similar to ``box`` argument but for
          an orthogonal box. For instance, one valid set of argument
          for ``box`` for the case of no PBC could be
          ``[10, 10, 10, 90, 90, 90]``
        * Following operations are advisable for non-PBC calculations

        ..code-block:: python

            lmax = all_coords.max(axis=0)
            lmin = all_coords.min(axis=0)
            pseudobox[:3] = 1.1*(lmax - lmin)
            pseudobox[3:] = 90.
            shift = all_coords.copy()
            shift -= lmin
            gridsearch = FastNS(max_cutoff, shift, box=pseudobox, pbc=False)

        """

        from MDAnalysis.lib.mdamath import triclinic_vectors


        _check_array(coords, 'coords')

        if np.allclose(box[:3], 0.0):
            raise ValueError("Any of the box dimensions cannot be 0")

        self.periodic = pbc
        self.coords = coords.copy()

        if box.shape != (3, 3):
            box = triclinic_vectors(box)

        self.box = PBCBox(box, self.periodic)

        if cutoff < 0:
            raise ValueError("Cutoff must be positive!")
        if cutoff * cutoff > self.box.c_pbcbox.max_cutoff2:
            raise ValueError("Cutoff greater than maximum cutoff ({:.3f}) given the PBC")

        self.coords_bbox = self.box.fast_put_atoms_in_bbox(self.coords)

        self.cutoff = cutoff
        self.max_gridsize = max_gridsize
        # Note that self.cutoff might be different from self.grid.cutoff
        # due to optimization
        self.grid = NSGrid(self.coords_bbox.shape[0], self.cutoff, self.box, self.max_gridsize)

        self.grid.fill_grid(self.coords_bbox)

    def search(self, search_coords):
        """Search a group of atoms against initialized coordinates

        Creates a new grid with the query atoms and searches
        against the initialized coordinates. The search is exclusive
        i.e. only the pairs ``(i, j)`` such that ``atom[i]`` from query atoms
        and ``atom[j]`` from the initialized set of coordinates is stored as
        neighbors.

        PBC-aware/non PBC-aware calculations are automatically enabled during
        the instantiation of :class:FastNS.

        Parameters
        ----------
        search_coords : numpy.ndarray
            Query coordinates of shape ``(N, 3)`` where
            ``N`` is the number of queries

        Returns
        -------
        results : NSResults object
           The object from :class:NSResults
           contains ``get_indices``, ``get_distances``.
           ``get_pairs``, ``get_pair_distances``

        Note
        ----
        For non-PBC aware calculations, the current implementation doesn't work
        if any of the query coordinates is beyond the range specified in
        ``box`` in :func:`MDAnalysis.lib.nsgrid.FastNS`.

        See Also
        --------
        MDAnalysis.lib.nsgrid.NSResults

        """

        cdef ns_int i, j, size_search
        cdef ns_int d, m
        cdef ns_int current_beadid, bid
        cdef ns_int cellindex, cellindex_probe
        cdef ns_int xi, yi, zi

        cdef NSResults results

        cdef real d2
        cdef rvec probe

        cdef real[:, ::1] searchcoords
        cdef real[:, ::1] searchcoords_bbox
        cdef NSGrid searchgrid
        cdef bint check

        cdef real cutoff2 = self.cutoff * self.cutoff
        cdef ns_int npairs = 0
        _check_array(search_coords, 'search_coords')

        # Generate another grid to search
        searchcoords = np.ascontiguousarray(search_coords, dtype=np.float32)
        searchcoords_bbox = self.box.fast_put_atoms_in_bbox(searchcoords)
        searchgrid = NSGrid(searchcoords_bbox.shape[0], self.grid.cutoff, self.box, self.max_gridsize, force=True)
        searchgrid.fill_grid(searchcoords_bbox)

        size_search = searchcoords.shape[0]

        results = NSResults(self.cutoff, self.coords, searchcoords)

        with nogil:
            for i in range(size_search):
                # Start with first search coordinate
                current_beadid = i
                # find the cellindex of the coordinate
                cellindex = searchgrid.cellids[current_beadid]
                for xi in range(DIM):
                    for yi in range(DIM):
                        for zi in range(DIM):
                            check = True
                            #Probe the search coordinates in a brick shaped box
                            probe[XX] = searchcoords_bbox[current_beadid, XX] + (xi - 1) * searchgrid.cellsize[XX]
                            probe[YY] = searchcoords_bbox[current_beadid, YY] + (yi - 1) * searchgrid.cellsize[YY]
                            probe[ZZ] = searchcoords_bbox[current_beadid, ZZ] + (zi - 1) * searchgrid.cellsize[ZZ]
                            # Make sure the probe coordinates is inside the brick-shaped box
                            if self.periodic:
                                for m in range(DIM - 1, -1, -1):
                                    while probe[m] < 0:
                                        for d in range(m+1):
                                            probe[d] += self.box.c_pbcbox.box[m][d]
                                    while probe[m] >= self.box.c_pbcbox.box[m][m]:
                                        for d in range(m+1):
                                            probe[d] -= self.box.c_pbcbox.box[m][d]
                            else:
                                for m in range(DIM -1, -1, -1):
                                    if probe[m] < 0:
                                        check = False
                                        break
                                    if probe[m] > self.box.c_pbcbox.box[m][m]:
                                        check = False
                                        break
                            if not check:
                                continue
                            # Get the cell index corresponding to the probe
                            cellindex_probe = self.grid.coord2cellid(probe)
                            # for this cellindex search in grid
                            for j in range(self.grid.nbeads[cellindex_probe]):
                                bid = self.grid.beadids[cellindex_probe * self.grid.nbeads_per_cell + j]
                                # find distance between search coords[i] and coords[bid]
                                d2 = self.box.fast_distance2(&searchcoords_bbox[current_beadid, XX], &self.coords_bbox[bid, XX])
                                if d2 < cutoff2:
                                    results.add_neighbors(current_beadid, bid, d2)
                                    npairs += 1
        return results

    def self_search(self):
        """Searches all the pairs within the initialized coordinates

        All the pairs among the initialized coordinates are registered
        in hald the time. Although the algorithm is still the same, but
        the distance checks can be reduced to half in this particular case
        as every pair need not be evaluated twice.

        Returns
        -------
        results : NSResults object
           The object from :class:NSResults
           contains ``get_indices``, ``get_distances``.
           ``get_pairs``, ``get_pair_distances``

        See Also
        --------
        MDAnalysis.lib.nsgrid.NSResults

        """

        cdef ns_int i, j, size_search
        cdef ns_int d, m
        cdef ns_int current_beadid, bid
        cdef ns_int cellindex, cellindex_probe
        cdef ns_int xi, yi, zi

        cdef NSResults results
        cdef real d2
        cdef rvec probe

        cdef real cutoff2 = self.cutoff * self.cutoff
        cdef ns_int npairs = 0
        cdef bint check

        size_search = self.coords.shape[0]

        results = NSResults(self.cutoff, self.coords, self.coords)

        with nogil:
            for i in range(size_search):
                # Start with first search coordinate
                current_beadid = i
                # find the cellindex of the coordinate
                cellindex = self.grid.cellids[current_beadid]
                for xi in range(DIM):
                    for yi in range(DIM):
                        for zi in range(DIM):
                            check = True
                            # Calculate and/or reinitialize shifted coordinates
                            # Probe the search coordinates in a brick shaped box
                            probe[XX] = self.coords_bbox[current_beadid, XX] + (xi - 1) * self.grid.cellsize[XX]
                            probe[YY] = self.coords_bbox[current_beadid, YY] + (yi - 1) * self.grid.cellsize[YY]
                            probe[ZZ] = self.coords_bbox[current_beadid, ZZ] + (zi - 1) * self.grid.cellsize[ZZ]
                            # Make sure the shifted coordinates is inside the brick-shaped box
                            if self.periodic:
                                for m in range(DIM - 1, -1, -1):
                                    while probe[m] < 0:
                                        for d in range(m+1):
                                            probe[d] += self.box.c_pbcbox.box[m][d]
                                    while probe[m] >= self.box.c_pbcbox.box[m][m]:
                                        for d in range(m+1):
                                            probe[d] -= self.box.c_pbcbox.box[m][d]
                            else:
                                for m in range(DIM -1, -1, -1):
                                    if probe[m] < 0:
                                        check = False
                                        break
                                    elif probe[m] >= self.box.c_pbcbox.box[m][m]:
                                        check = False
                                        break
                            if not check:
                                continue
                            # Get the cell index corresponding to the probe
                            cellindex_probe = self.grid.coord2cellid(probe)
                            # for this cellindex search in grid
                            for j in range(self.grid.nbeads[cellindex_probe]):
                                bid = self.grid.beadids[cellindex_probe * self.grid.nbeads_per_cell + j]
                                if bid < current_beadid:
                                    continue
                                # find distance between search coords[i] and coords[bid]
                                d2 = self.box.fast_distance2(&self.coords_bbox[current_beadid, XX], &self.coords_bbox[bid, XX])
                                if d2 < cutoff2 and d2 > EPSILON:
                                    results.add_neighbors(current_beadid, bid, d2)
                                    results.add_neighbors(bid, current_beadid, d2)
                                    npairs += 2
        return results
