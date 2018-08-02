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

This Neighbor search library is a serialized Cython version greatly inspired by the  NS grid search implemented in
`GROMACS <http://www.gromacs.org/>`_ .

GROMACS 4.x code (more precisely `nsgrid.c <https://github.com/gromacs/gromacs/commits/master/src/mdlib/nsgrid.c>`_  and
`ns.c <https://github.com/gromacs/gromacs/commits/master/src/mdlib/ns.c>`_ ) was used as reference to write this file.

GROMACS 4.x code is released under the GNU Public Licence v2.


About the algorithm
-------------------

The neighbor search implemented here is based on `cell lists <https://en.wikipedia.org/wiki/Cell_lists>`_ which allow
computation of pairs [#]_ with a cost of :math:`O(N)`, instead of :math:`O(N^2)`.
The basic algorithm is described in the following references:`

Examples
--------


.. [#] a pair correspond to two particles that are considered as neighbors .

"""


# Preprocessor DEFs
DEF DIM = 3
DEF XX = 0
DEF YY = 1
DEF ZZ = 2

DEF EPSILON = 1e-5

DEF BOX_MARGIN=1.0010
DEF MAX_NTRICVEC=12

DEF OK=0
DEF ERROR=1

# Used to handle memory allocation
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from libc.math cimport sqrt
import numpy as np
cimport numpy as np

ctypedef np.int_t ns_int
ctypedef np.float32_t real
ctypedef real rvec[DIM]
ctypedef ns_int ivec[DIM]
ctypedef ns_int ipair[2]
ctypedef real matrix[DIM][DIM]


# Useful Functions
cdef real rvec_norm2(const rvec a) nogil:
    return a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ]

cdef void rvec_clear(rvec a) nogil:
    a[XX]=0.0
    a[YY]=0.0
    a[ZZ]=0.0

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
    Cython implementation of `PBC-related <https://en.wikipedia.org/wiki/Periodic_boundary_conditions>`_ operations.
    This class is used by classes :class:`FastNS` and :class:`NSGrid` to put all particles inside a brick-shaped box
    and to compute PBC-aware distance.

    .. warning::

       This class is not meant to be used by end users.

    .. warning::

       Even if MD triclinic boxes can be handle by this class, internal optimization is made based on the
       assumption that particles are inside a brick-shaped box. When this is not the case, calculated distances are not
       warranted to be exact.

    """
    cdef cPBCBox_t c_pbcbox
    cdef rvec center
    cdef rvec bbox_center
    cdef bint is_triclinic

    def __init__(self, real[:,::1] box):
        """

        Parameters
        ----------

        box : :class:`numpy.ndarray`
          the MD box vectors as returned by
          :func:`MDAnalysis.lib.mdamath.triclinic_vectors`. dtype must be :class:`numpy.float32`

        """
        self.update(box)

    cdef void fast_update(self, real[:,::1] box) nogil:
        cdef ns_int i, j, k, d, jc, kc, shift
        cdef real d2old, d2new, d2new_c
        cdef rvec trial, pos
        cdef ns_int ii, jj ,kk
        cdef ns_int *order = [0, -1, 1, -2, 2]
        cdef bint use
        cdef real min_hv2, min_ss, tmp

        rvec_clear(self.center)
        # Update matrix
        self.is_triclinic = False
        for i in range(DIM):
            for j in range(DIM):
                self.c_pbcbox.box[i][j] = box[i, j]
                self.center[j] += 0.5 * box[i, j]

                if i != j:
                    if box[i, j] > EPSILON:
                        self.is_triclinic = True
            self.bbox_center[i] = 0.5 *  box[i, i]

        # Update diagonals
        for i in range(DIM):
            self.c_pbcbox.fbox_diag[i]  =  box[i, i]
            self.c_pbcbox.hbox_diag[i]  =  self.c_pbcbox.fbox_diag[i] * 0.5
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

    def update(self, real[:,::1] box):
        """

        Updates internal MD box representation and parameters used for calculations.

        .. note::

            Call to this method is only needed when the MD box is changed as it always called when class is
            instantiated.

        Parameters
        ----------

        box : a :class:`numpy.ndarray` that describes the MD box vectors as returned by
            :func:`MDAnalysis.lib.mdamath.triclinic_vectors`.
            `dtype` must be :class:`numpy.float32`

        """

        if box.shape[0] != DIM or box.shape[1] != DIM:
            raise ValueError("Box must be a {} x {} matrix. Got: {} x {})".format(
                DIM, DIM, box.shape[0], box.shape[1]))
        if (box[XX, XX] == 0) or (box[YY, YY] == 0) or (box[ZZ, ZZ] == 0):
            raise ValueError("Box does not correspond to PBC=xyz")
        self.fast_update(box)


    cdef void fast_pbc_dx(self, rvec ref, rvec other, rvec dx) nogil:
        cdef ns_int i, j
        cdef rvec dx_start, trial

        for i in range(DIM):
            dx[i] = other[i] - ref[i]

        for i in range (DIM-1, -1, -1):
            while dx[i] > self.c_pbcbox.hbox_diag[i]:
                for j in range (i, -1, -1):
                    dx[j] -= self.c_pbcbox.box[i][j]

            while dx[i] <= self.c_pbcbox.mhbox_diag[i]:
                for j in range (i, -1, -1):
                    dx[j] += self.c_pbcbox.box[i][j]

    def dx(self, real[:] a, real[:] b):
        """


        Returns the distance vector :math:`dx` between two coordinates :math:`a` and :math:`b` .
        if the minimum image convention is respected by :math:`a` and :math:`b` then :math:`dx = ab` .
        if not, :math:`dx` is modified so it is always compliant with the minimum image convention.

        Parameters
        ----------
        a : :class:`numpy.ndarray`
          coordinates with a `shape` of 3 and a `dtype` of :class:`numpy.float32`
        b : :class:`numpy.ndarray`
          coordinates with a `shape` of 3 and a `dtype` of :class:`numpy.float32`

        """
        cdef rvec dx
        if a.shape[0] != DIM or b.shape[0] != DIM:
            raise ValueError("Not 3 D coordinates")

        self.fast_pbc_dx(&a[XX], &b[XX], dx)

        return np.array([dx[XX], dx[YY], dx[ZZ]], dtype=np.float32)


    cdef real fast_distance2(self, rvec a, rvec b) nogil:
        cdef rvec dx
        self.fast_pbc_dx(a, b, dx)
        return rvec_norm2(dx)

    def distance2(self, real[:] a, real[:] b):
        """
        Returns the distance vector :math:`dx` between two coordinates :math:`a` and :math:`b` .
        if the minimum image convention is respected by :math:`a` and :math:`b` then :math:`dx = ab` .
        if not, :math:`dx` is modified so it is always compliant with the minimum image convention.

        Parameters
        ----------
        a : :class:`numpy.ndarray`
          coordinates with a `shape` of 3 and a `dtype` of :class:`numpy.float32`
        b : :class:`numpy.ndarray`
          coordinates with a `shape` of 3 and a `dtype` of :class:`numpy.float32`

        """
        if a.shape[0] != DIM or b.shape[0] != DIM:
            raise ValueError("Not 3 D coordinates")
        return self.fast_distance2(&a[XX], &b[XX])

    cdef real fast_distance(self, rvec a, rvec b) nogil:
        return sqrt(self.fast_distance2(a,b))

    def distance(self, real[:] a, real[:] b):
        if a.shape[0] != DIM or b.shape[0] != DIM:
            raise ValueError("Not 3 D coordinates")
        return self.fast_distance(&a[XX], &b[XX])

    cdef real[:, ::1]fast_put_atoms_in_bbox(self, real[:,::1] coords) nogil:
        cdef ns_int i, m, d, natoms, wd = 0
        cdef real[:,::1] bbox_coords

        natoms = coords.shape[0]
        with gil:
            bbox_coords = coords.copy()

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

    def put_atoms_in_bbox(self, real[:,::1] coords):
        return np.asarray(self.fast_put_atoms_in_bbox(coords))


#########################
# Neighbor Search Stuff #
#########################
cdef class NSResults(object):
    cdef readonly real cutoff
    cdef ns_int npairs
    cdef bint debug

    cdef real[:, ::1] coords # shape: size, DIM
    cdef ns_int[:] search_ids

    cdef ns_int allocation_size
    cdef ipair *pairs # shape: pair_allocation
    cdef real *pair_distances2 # shape: pair_allocation

    cdef list indices_buffer
    cdef list coordinates_buffer
    cdef list distances_buffer
    cdef np.ndarray pairs_buffer
    cdef np.ndarray pair_distances_buffer
    cdef np.ndarray pair_coordinates_buffer

    def __init__(self, real cutoff, real[:, ::1]coords, ns_int[:] search_ids, debug=False):
        self.debug = debug
        self.cutoff = cutoff
        self.coords = coords
        self.search_ids = search_ids

        # Preallocate memory
        self.allocation_size = search_ids.shape[0] + 1
        if not self.pairs and not self.pair_distances2:
            self.pairs = <ipair *> PyMem_Malloc(sizeof(ipair) * self.allocation_size)
            if not self.pairs:
                MemoryError("Could not allocate memory for NSResults.pairs "
                            "({} bits requested)".format(sizeof(ipair) * self.allocation_size))
            self.pair_distances2 = <real *> PyMem_Malloc(sizeof(real) * self.allocation_size)
            if not self.pair_distances2:
                raise MemoryError("Could not allocate memory for NSResults.pair_distances2 "
                                  "({} bits requested)".format(sizeof(real) * self.allocation_size))
        else:
            if self.resize(self.allocation_size) != OK:
                raise MemoryError("foo")

        self.npairs = 0

        # Buffer
        self.indices_buffer = None
        self.coordinates_buffer = None
        self.distances_buffer = None
        self.pairs_buffer = None
        self.pair_coordinates_buffer = None

    def __dealloc__(self):
        PyMem_Free(self.pairs)
        PyMem_Free(self.pair_distances2)

    cdef int add_neighbors(self, ns_int beadid_i, ns_int beadid_j, real distance2) nogil:
        # Important: If this function returns ERROR, it means that memory allocation failed

        # Reallocate memory if needed
        if self.npairs >= self.allocation_size:
            # We need to reallocate memory
            if self.resize(self.allocation_size + <ns_int> (self.allocation_size * 0.5 + 1)) != OK:
                return ERROR

        # Actually store pair and distance squared
        if beadid_i < beadid_j:
            self.pairs[self.npairs][0] = beadid_i
            self.pairs[self.npairs][1] = beadid_j
        else:
            self.pairs[self.npairs][1] = beadid_i
            self.pairs[self.npairs][0] = beadid_j
        self.pair_distances2[self.npairs] = distance2
        self.npairs += 1

        return OK

    cdef int resize(self, ns_int new_size) nogil:
        # Important: If this function returns 0, it means that memory allocation failed

        if new_size < self.npairs:
            # Silently ignored the request
            return OK

        if self.allocation_size >= new_size:
            if self.debug:
                with gil:
                    print("NSresults: Reallocation requested but not needed ({} requested but {} already allocated)".format(new_size, self.allocation_size))
            return OK

        self.allocation_size = new_size

        if self.debug:
            with gil:
                print("NSresults: Reallocated to {} pairs".format(self.allocation_size))

        # Allocating memory
        with gil:
            self.pairs = <ipair *> PyMem_Realloc(self.pairs, sizeof(ipair) * self.allocation_size)
            self.pair_distances2 = <real *> PyMem_Realloc(self.pair_distances2, sizeof(real) * self.allocation_size)

        if not self.pairs:
            return ERROR

        if not self.pair_distances2:
            return ERROR

        return OK

    def get_pairs(self):
        cdef ns_int i

        if self.pairs_buffer is None:
            self.pairs_buffer = np.empty((self.npairs, 2), dtype=np.int)
            for i in range(self.npairs):
                self.pairs_buffer[i, 0] = self.pairs[i][0]
                self.pairs_buffer[i, 1] = self.pairs[i][1]
        return self.pairs_buffer

    def get_pair_distances(self):
        cdef ns_int i
        if self.pair_coordinates_buffer is None:
            self.pair_coordinates_buffer = np.empty(self.npairs, dtype=np.float32)
            for i in range(self.npairs):
                self.pair_coordinates_buffer[i] = self.pair_distances2[i]
            self.pair_coordinates_buffer = np.sqrt(self.pair_coordinates_buffer)
        return self.pair_coordinates_buffer

    def get_pair_coordinates(self):
        cdef ns_int i, j, bead_i, bead_j
        if self.pair_coordinates_buffer is None:
            self.pair_coordinates_buffer = np.empty((self.npairs, 2, DIM), dtype=np.float32)
            for i in range(self.npairs):
                bead_i = self.pairs[i][0]
                bead_j = self.pairs[i][1]

                for j in range(DIM):
                    self.pair_coordinates_buffer[i, 0, j] = self.coords[bead_i, j]
                    self.pair_coordinates_buffer[i, 1, j] = self.coords[bead_j, j]
        return self.pair_coordinates_buffer

    cdef create_buffers(self):
        cdef ns_int i, beadid_i, beadid_j
        cdef real dist2
        cdef real[:] coord_i, coord_j
        from collections import defaultdict

        indices_buffer = defaultdict(list)
        coords_buffer = defaultdict(list)
        dists_buffer = defaultdict(list)

        for i in range(self.npairs):
            beadid_i = self.pairs[i][0]
            beadid_j = self.pairs[i][1]

            dist2 = self.pair_distances2[i]
            coord_i = self.coords[beadid_i]
            coord_j = self.coords[beadid_j]

            indices_buffer[beadid_i].append(beadid_j)
            indices_buffer[beadid_j].append(beadid_i)

            coords_buffer[beadid_i].append(coord_j)
            coords_buffer[beadid_j].append((coord_i))

            dists_buffer[beadid_i].append(dist2)
            dists_buffer[beadid_j].append(dist2)

        self.indices_buffer = []
        self.coordinates_buffer = []
        self.distances_buffer = []

        for elm in self.search_ids:
            sorted_indices = np.argsort(indices_buffer[elm])
            self.indices_buffer.append(np.array(indices_buffer[elm])[sorted_indices])
            self.coordinates_buffer.append(np.array(coords_buffer[elm])[sorted_indices])
            self.distances_buffer.append(np.sqrt(dists_buffer[elm])[sorted_indices])

    def get_indices(self):
        if self.indices_buffer is None:
            self.create_buffers()
        return self.indices_buffer

    def get_distances(self):
        if self.distances_buffer is None:
            self.create_buffers()
        return self.distances_buffer

    def get_coordinates(self):
        if self.coordinates_buffer is None:
            self.create_buffers()
        return self.coordinates_buffer



cdef class NSGrid(object):
    cdef bint debug
    cdef readonly real cutoff
    cdef ns_int size
    cdef ns_int ncoords
    cdef ns_int[DIM] ncells
    cdef ns_int[DIM] cell_offsets
    cdef real[DIM] cellsize
    cdef ns_int nbeads_per_cell
    cdef ns_int *nbeads # size
    cdef ns_int *beadids # size * nbeads_per_cell
    cdef ns_int *cellids # ncoords

    def __init__(self, ncoords, cutoff, PBCBox box, max_size, debug=False):
        cdef ns_int i, x, y, z
        cdef ns_int ncellx, ncelly, ncellz, size, nbeadspercell
        cdef real bbox_vol
        self.debug = debug

        self.ncoords = ncoords

        # Calculate best cutoff
        self.cutoff = cutoff
        bbox_vol = box.c_pbcbox.box[XX][XX] * box.c_pbcbox.box[YY][YY] * box.c_pbcbox.box[YY][YY]
        size = bbox_vol/cutoff**3
        nbeadspercell = ncoords/size
        while bbox_vol/self.cutoff**3 > max_size:
            self.cutoff *= 1.2

        for i in range(DIM):
            self.ncells[i] = <ns_int> (box.c_pbcbox.box[i][i] / self.cutoff)
            if self.ncells[i] == 0:
                self.ncells[i] = 1
            self.cellsize[i] = box.c_pbcbox.box[i][i] / self.ncells[i]
        self.size = self.ncells[XX] * self.ncells[YY] * self.ncells[ZZ]

        if self.debug:
            print("NSGrid: Requested cutoff: {:.3f} (Ncells={}, Avg # of beads per cell={}), Optimized cutoff= {:.3f} (Ncells={}, Avg # of beads per cell={})".format(
                cutoff, size, nbeadspercell,
                self.cutoff, self.size, <ns_int> (ncoords / self.size)
            ))
            print("NSGrid: Size={}x{}x{}={}".format(self.ncells[XX], self.ncells[YY], self.ncells[ZZ], self.size))

        self.cell_offsets[XX] = 0
        self.cell_offsets[YY] = self.ncells[XX]
        self.cell_offsets[ZZ] = self.ncells[XX] * self.ncells[YY]

        # Allocate memory
        self.nbeads = <ns_int *> PyMem_Malloc(sizeof(ns_int) * self.size)
        if not self.nbeads:
            raise MemoryError("Could not allocate memory from NSGrid.nbeads ({} bits requested)".format(sizeof(ns_int) * self.size))
        self.beadids = NULL
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
        return <ns_int> (coord[ZZ] / self.cellsize[ZZ]) * (self.ncells[XX] * self.ncells[YY]) +\
               <ns_int> (coord[YY] / self.cellsize[YY]) * self.ncells[XX] + \
               <ns_int> (coord[XX] / self.cellsize[XX])

    cdef bint cellid2cellxyz(self, ns_int cellid, ivec cellxyz) nogil:
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
        self.beadids = <ns_int *> PyMem_Malloc(sizeof(ns_int) * self.size * self.nbeads_per_cell) #np.empty((self.size, nbeads_max), dtype=np.int)
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
    cdef bint debug
    cdef PBCBox box
    cdef real[:, ::1] coords
    cdef real[:, ::1] coords_bbox
    cdef readonly real cutoff
    cdef bint prepared
    cdef NSGrid grid

    def __init__(self, u, cutoff, coords=None, prepare=True, debug=False, max_gridsize=5000):
        import MDAnalysis as mda
        from MDAnalysis.lib.mdamath import triclinic_vectors

        self.debug = debug

        if not isinstance(u, mda.Universe):
            raise TypeError("FastNS class must be initialized with a valid MDAnalysis.Universe instance")
        box = triclinic_vectors(u.dimensions)

        self.box = PBCBox(box)


        if coords is None:
            coords = u.atoms.positions

        self.coords = coords.copy()

        self.coords_bbox =  self.box.fast_put_atoms_in_bbox(coords)

        if cutoff < 0:
            raise ValueError("Cutoff must be positive!")
        if cutoff * cutoff > self.box.c_pbcbox.max_cutoff2:
            raise ValueError("Cutoff greater than maximum cutoff ({:.3f}) given the PBC")
        self.cutoff = cutoff

        self.grid = NSGrid(self.coords_bbox.shape[0], cutoff, self.box, max_gridsize, debug=debug)
        self.prepared = False
        if prepare:
            self.prepare()


    def prepare(self, force=False):
        if self.prepared and not force:
            return

        self.grid.fill_grid(self.coords_bbox)

        self.prepared = True

    def search(self, search_ids=None):
        cdef ns_int i, j, size_search
        cdef ns_int d, m
        cdef NSResults results
        cdef ns_int size = self.coords_bbox.shape[0]

        cdef ns_int current_beadid, bid
        cdef rvec current_coords

        cdef ns_int cellindex, cellindex_adjacent, cellindex_probe
        cdef ivec cellxyz, debug_cellxyz

        cdef real[:, ::1] search_coords
        cdef ns_int[:] search_ids_view

        cdef ns_int xi, yi, zi
        cdef real d2
        cdef rvec shifted_coord, probe, dx

        cdef ns_int nchecked = 0
        cdef ns_int[:] checked = np.zeros(size, dtype=np.int)

        cdef real cutoff2 = self.cutoff * self.cutoff
        cdef ns_int npairs = 0

        if not self.prepared:
            self.prepare()

        if search_ids is None:
            search_ids=np.arange(size)
        elif type(search_ids) == np.int:
            search_ids = np.array([search_ids,], dtype=np.int)
        elif type(search_ids) != np.ndarray:
            search_ids = np.array(search_ids, dtype=np.int)

        search_ids_view = search_ids
        size_search = search_ids.shape[0]

        results = NSResults(self.cutoff, self.coords, search_ids, self.debug)

        cdef bint memory_error = False

        with nogil:
            for i in range(size_search):
                if memory_error:
                    break
                current_beadid = search_ids_view[i]
                cellindex = self.grid.cellids[current_beadid]
                self.grid.cellid2cellxyz(cellindex, cellxyz)
                for xi in range(DIM):
                    if memory_error:
                        break
                    for yi in range(DIM):
                        if memory_error:
                            break
                        for zi in range(DIM):
                            if memory_error:
                                break
                            # Calculate and/or reinitialize shifted coordinates
                            shifted_coord[XX] = self.coords[current_beadid, XX] + (xi - 1) * self.grid.cellsize[XX]
                            shifted_coord[YY] = self.coords[current_beadid, YY] + (yi - 1) * self.grid.cellsize[YY]
                            shifted_coord[ZZ] = self.coords[current_beadid, ZZ] + (zi - 1) * self.grid.cellsize[ZZ]
                            probe[XX] = self.coords[current_beadid, XX] + (xi - 1) * self.cutoff
                            probe[YY] = self.coords[current_beadid, YY] + (yi - 1) * self.cutoff
                            probe[ZZ] = self.coords[current_beadid, ZZ] + (zi - 1) * self.cutoff
                            # Make sure the shifted coordinates is inside the brick-shaped box
                            for m in range(DIM - 1, -1, -1):
                                while shifted_coord[m] < 0:
                                    for d in range(m+1):
                                        shifted_coord[d] += self.box.c_pbcbox.box[m][d]
                                while shifted_coord[m] >= self.box.c_pbcbox.box[m][m]:
                                    for d in range(m+1):
                                        shifted_coord[d] -= self.box.c_pbcbox.box[m][d]
                                while probe[m] < 0:
                                    for d in range(m+1):
                                        probe[d] += self.box.c_pbcbox.box[m][d]
                                while probe[m] >= self.box.c_pbcbox.box[m][m]:
                                    for d in range(m+1):
                                        probe[d] -= self.box.c_pbcbox.box[m][d]

                            # Get the cell index corresponding to the coord
                            cellindex_adjacent = self.grid.coord2cellid(shifted_coord)
                            cellindex_probe = self.grid.coord2cellid(probe)

                            for j in range(self.grid.nbeads[cellindex_adjacent]):
                                bid = self.grid.beadids[cellindex_adjacent * self.grid.nbeads_per_cell + j]
                                if checked[bid] != 0:
                                    continue
                                d2 = self.box.fast_distance2(&self.coords_bbox[current_beadid, XX], &self.coords_bbox[bid, XX])
                                if d2 < cutoff2:
                                    if d2 < EPSILON:
                                        continue
                                    elif results.add_neighbors(current_beadid, bid, d2) != OK:
                                        memory_error = True
                                        break
                                    npairs += 1
                checked[current_beadid] = 1
        if memory_error:
            raise MemoryError("Could not allocate memory to store NS results")
        return results
