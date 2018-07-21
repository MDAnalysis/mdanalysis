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

"""
Neighbor search library --- :mod:`MDAnalysis.lib.grid`
======================================================

This Neighbor search library is a serialized Cython port of the NS grid search implemented in GROMACS.
"""


# Preprocessor DEFs
DEF DIM = 3
DEF XX = 0
DEF YY = 1
DEF ZZ = 2

DEF EPSILON = 1e-5

DEF BOX_MARGIN=1.0010
DEF MAX_NTRICVEC=12

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
    cdef cPBCBox_t c_pbcbox
    cdef rvec center
    cdef rvec bbox_center
    cdef bint is_triclinic

    def __init__(self, real[:,::1] box):
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
            if self.resize(self.allocation_size + <ns_int> (self.allocation_size * 0.5 + 1)) == 0:
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
        # Important: If this function returns 0, it means that memory allocation failed

        # Reallocate memory if needed
        if self.npairs >= self.allocation_size:
            # We need to reallocate memory
            if self.resize(self.allocation_size + <ns_int> (self.allocation_size * 0.5 + 1)) == 0:
                return 0

        # Actually store pair and distance squared
        if beadid_i < beadid_j:
            self.pairs[self.npairs][0] = beadid_i
            self.pairs[self.npairs][1] = beadid_j
        else:
            self.pairs[self.npairs][1] = beadid_i
            self.pairs[self.npairs][0] = beadid_j
        self.pair_distances2[self.npairs] = distance2
        self.npairs += 1

        return self.npairs

    cdef int resize(self, ns_int new_size) nogil:
        # Important: If this function returns 0, it means that memory allocation failed

        if new_size < self.npairs:
            # Silently ignored the request
            return 1

        if self.allocation_size >= new_size:
            if self.debug:
                with gil:
                    print("NSresults: Reallocation requested but not needed ({} requested but {} already allocated)".format(new_size, self.allocation_size))
            return 1

        self.allocation_size = new_size

        if self.debug:
            with gil:
                print("NSresults: Reallocated to {} pairs".format(self.allocation_size))

        # Allocating memory
        with gil:
            self.pairs = <ipair *> PyMem_Realloc(self.pairs, sizeof(ipair) * self.allocation_size)
            self.pair_distances2 = <real *> PyMem_Realloc(self.pair_distances2, sizeof(real) * self.allocation_size)

        if not self.pairs:
            return 0

        if not self.pair_distances2:
            return 0

        return 1

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
        cdef ns_int *beadcounts = NULL

        # Allocate memory
        beadcounts = <ns_int *> PyMem_Malloc(sizeof(ns_int) * self.size)
        if not beadcounts:
            raise MemoryError("Could not allocate memory for bead count buffer ({} bits requested)".format(sizeof(ns_int) * self.size))

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
            with gil:
                self.beadids = <ns_int *> PyMem_Malloc(sizeof(ns_int) * self.size * self.nbeads_per_cell) #np.empty((self.size, nbeads_max), dtype=np.int)
                if not self.beadids:
                    raise MemoryError("Could not allocate memory for NSGrid.beadids ({} bits requested)".format(sizeof(ns_int) * self.size * self.nbeads_per_cell))

            # Second loop: fill grid
            for i in range(ncoords):

                # Add bead to grid cell
                cellindex = self.cellids[i]
                self.beadids[cellindex * self.nbeads_per_cell + beadcounts[cellindex]] = i
                beadcounts[cellindex] += 1

        # Now we can free the allocation buffer
        PyMem_Free(beadcounts)


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

        cdef real cutoff2 = self.cutoff * self.cutoff
        cdef ns_int[:] checked
        cdef ns_int npairs = 0

        #cdef bint debug=False


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

        checked = np.zeros(size, dtype=np.int)

        results = NSResults(self.cutoff, self.coords, search_ids, self.debug)

        cdef bint memory_error = False

        # if self.debug and debug:
        #     print("FastNS: Debug flag is set to True for FastNS.search()")

        with nogil:
            for i in range(size_search):
                if memory_error:
                    break
                current_beadid = search_ids_view[i]

                cellindex = self.grid.cellids[current_beadid]
                self.grid.cellid2cellxyz(cellindex, cellxyz)

                # if self.debug and debug:
                #     with gil:
                #         print("FastNS: Checking neighbors for bead #{} ({:.3f},{:.3f},{:.3f})->rect({:.3f},{:.3f},{:.3f}) - cell[{},{},{}]:" .format(
                #             current_beadid,
                #             self.coords[current_beadid, XX], self.coords[current_beadid, YY], self.coords[current_beadid, ZZ],
                #             self.coords_bbox[current_beadid, XX], self.coords_bbox[current_beadid, YY], self.coords_bbox[current_beadid, ZZ],
                #             cellxyz[XX], cellxyz[YY], cellxyz[ZZ]))


                for xi in range(DIM):
                    if memory_error:
                        break

                    if not self.box.is_triclinic:
                        # If box is not triclinic (ie rect), when can already check if the shift can be skipped (ie cutoff inside the cell)
                        if xi == 0:
                            if self.coords_bbox[current_beadid, XX] - self.cutoff > self.grid.cellsize[XX] * cellxyz[XX]:
                                # if self.debug and debug:
                                #     with gil:
                                #         print("FastNS:   -> Bead X={:.3f}, Cell X={:.3f}, cutoff={:.3f} -> -X shift ignored".format(
                                #             self.coords_bbox[current_beadid, XX],
                                #             self.grid.cellsize[XX] * cellxyz[XX],
                                #             self.cutoff
                                #         ))
                                continue

                        if xi == 2:
                            if self.coords_bbox[current_beadid, XX] + self.cutoff < self.grid.cellsize[XX] * (cellxyz[XX] + 1):
                                # if self.debug and debug:
                                #     with gil:
                                #         print(
                                #             "FastNS:   -> Bead X={:.3f}, Next cell X={:.3f}, cutoff={:.3f} -> +X shift ignored".format(
                                #                 self.coords_bbox[current_beadid, XX],
                                #                 self.grid.cellsize[XX] * (cellxyz[XX] + 1),
                                #                 self.cutoff
                                #             ))
                                continue

                    for yi in range(DIM):
                        if memory_error:
                            break

                        if not self.box.is_triclinic:
                            if yi == 0:
                                if self.coords_bbox[current_beadid, YY] - self.cutoff > self.grid.cellsize[YY] * cellxyz[YY]:
                                    # if self.debug and debug:
                                    #     with gil:
                                    #         print("FastNS:   -> Bead Y={:.3f}, Cell Y={:.3f}, cutoff={:.3f} -> -Y shift is ignored".format(
                                    #             self.coords_bbox[current_beadid, YY],
                                    #             self.grid.cellsize[YY] * cellxyz[YY],
                                    #             self.cutoff,
                                    #         ))
                                    continue

                            if yi == 2:
                                if self.coords_bbox[current_beadid, YY] + self.cutoff < self.grid.cellsize[YY] * (cellxyz[YY] + 1):
                                    # if self.debug and debug:
                                    #     with gil:
                                    #         print("FastNS:   -> Bead Y={:.3f}, Next cell Y={:.3f}, cutoff={:.3f} -> +Y shift is ignored".format(
                                    #             self.coords_bbox[current_beadid, YY],
                                    #             self.grid.cellsize[YY] * (cellxyz[YY] +1),
                                    #             self.cutoff,
                                    #         ))
                                    continue


                        for zi in range(DIM):
                            if not self.box.is_triclinic:
                                if zi == 0:
                                    if self.coords_bbox[current_beadid, ZZ] - self.cutoff > self.grid.cellsize[ZZ] * cellxyz[ZZ]:
                                        if self.coords_bbox[current_beadid, ZZ] - self.cutoff > 0:
                                            # if self.debug and debug:
                                            #     with gil:
                                            #         print("FastNS:   -> Bead Z={:.3f}, Cell Z={:.3f}, cutoff={:.3f} -> -Z shift is ignored".format(self.coords_bbox[current_beadid, ZZ], self.grid.cellsize[ZZ] * cellxyz[ZZ], self.cutoff))
                                            continue

                                if zi == 2:
                                    if self.coords_bbox[current_beadid, ZZ] + self.cutoff < self.grid.cellsize[ZZ] * (cellxyz[ZZ] + 1):
                                        if self.coords_bbox[current_beadid, ZZ] + self.cutoff < self.box.c_pbcbox.box[ZZ][ZZ]:
                                            # if self.debug and debug:
                                            #     with gil:
                                            #         print("FastNS:   -> Bead Z={:.3f}, Next cell Z={:.3f}, cutoff={:.3f} -> +Z shift is ignored".format(self.coords_bbox[current_beadid, ZZ], self.grid.cellsize[XX] * (cellxyz[ZZ] + 1), self.cutoff))
                                            continue

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

                            if cellindex == cellindex_probe and xi != 1 and yi != 1 and zi != 1:
                                # if self.debug and debug:
                                #     with gil:
                                #         print("FastNS:   Grid shift [{}][{}][{}]: Cutoff is inside current cell -> This shift is ignored".format(
                                #             xi - 1,
                                #             yi -1,
                                #             zi -1
                                #         ))
                                continue

                            # if self.debug and debug:
                            #     self.grid.cellid2cellxyz(cellindex_adjacent, debug_cellxyz)
                            #     with gil:
                            #         dist_shift = self.box.fast_distance(&self.coords[current_beadid, XX], shifted_coord)
                            #         grid_shift = np.array([(xi - 1) * self.grid.cellsize[XX],
                            #                                (yi - 1) * self.grid.cellsize[YY],
                            #                                (zi - 1) * self.grid.cellsize[ZZ]])
                            #         print("FastNS:   -> Checking cell#{} ({},{},{}) for neighbors (dshift={:.3f}, grid_shift=({:.3f},{:.3f},{:.3f}->{:.3f})".format(
                            #             cellindex,
                            #             debug_cellxyz[XX], debug_cellxyz[YY], debug_cellxyz[ZZ],
                            #             dist_shift,
                            #             grid_shift[XX], grid_shift[YY], grid_shift[ZZ],
                            #             np.sqrt(np.sum(grid_shift**2))
                            #         ))


                            for j in range(self.grid.nbeads[cellindex_adjacent]):
                                bid = self.grid.beadids[cellindex_adjacent * self.grid.nbeads_per_cell + j]

                                if checked[bid] != 0:
                                    continue

                                d2 = self.box.fast_distance2(&self.coords_bbox[current_beadid, XX], &self.coords_bbox[bid, XX])

                                # if self.debug:
                                #     self.grid.cellid2cellxyz(cellindex, debug_cellxyz)
                                #     with gil:
                                #         print(
                                #             "Beads #{} (cell[{},{},{}]-coords[{:.3f},{:.3f},{:.3f}]) and #{} (cell[{},{},{}]-coords[{:.3f},{:.3f},{:.3f}]) are tested (d2={:.3f})".format(
                                #                 current_beadid,
                                #                 cellxyz[XX], cellxyz[YY], cellxyz[ZZ],
                                #                 self.coords_bbox[current_beadid, XX],
                                #                 self.coords_bbox[current_beadid, YY],
                                #                 self.coords_bbox[current_beadid, ZZ],
                                #                 bid,
                                #                 debug_cellxyz[XX], debug_cellxyz[YY], debug_cellxyz[ZZ],
                                #                 self.coords_bbox[bid, XX], self.coords_bbox[bid, YY],
                                #                 self.coords_bbox[bid, ZZ],
                                #                 d2))

                                if d2 < cutoff2:

                                    if d2 < EPSILON:
                                        continue

                                    # if self.debug and debug:
                                    #     self.grid.cellid2cellxyz(cellindex, debug_cellxyz)
                                    #     with gil:
                                    #         self.box.fast_pbc_dx(&self.coords[current_beadid, XX], &self.coords[bid, XX], dx)
                                    #         dx_py = np.array([dx[XX], dx[YY], dx[ZZ]])
                                    #         print("FastNS:       \_ Neighbor found:  bead#{} (cell[{},{},{}]) -> dx={} -> d={:.3f}".format(bid, debug_cellxyz[XX], debug_cellxyz[YY], debug_cellxyz[ZZ], dx_py, np.sqrt(d2)))
                                    if results.add_neighbors(current_beadid, bid, d2) == 0:
                                        memory_error = True
                                        break
                                    npairs += 1
                checked[current_beadid] = 1

        if memory_error:
            raise MemoryError("Could not allocate memory to store NS results")


        if self.debug:
            print("Total number of pairs={}".format(npairs))

            # ref_bead = 13937
            # beads = np.array([4398, 4401, 13939, 13940, 13941, 17987, 23518, 23519, 23521, 23734, 47451]) - 1
            # for bid in beads:
            #     self.box.fast_pbc_dx(&self.coords[ref_bead, XX], &self.coords[bid, XX], dx)
            #     dx_py = np.array([dx[XX], dx[YY], dx[ZZ]])
            #     self.box.fast_pbc_dx(&self.coords_bbox[ref_bead, XX], &self.coords_bbox[bid, XX], dx)
            #     rect_dx_py = np.array([dx[XX], dx[YY], dx[ZZ]])
            #     self.grid.cellid2cellxyz(self.grid.coord2cellid(&self.coords_bbox[bid, XX]), cellxyz)
            #     print("Bead #{} ({:.3f},{:.3f},{:.3f})->rect({:.3f},{:.3f},{:.3f}) - cell[{},{},{}]: dx=[{:.3f},{:.3f},{:.3f}] -> dist: {:.3f} ({})".format(
            #         bid,
            #         self.coords[bid, XX], self.coords[bid, YY], self.coords[bid,ZZ],
            #         self.coords_bbox[bid, XX], self.coords_bbox[bid, YY], self.coords_bbox[bid,ZZ],
            #         cellxyz[XX], cellxyz[YY], cellxyz[ZZ],
            #         dx[XX], dx[YY], dx[ZZ],
            #         np.sqrt(np.sum(dx_py**2)),
            #         self.box.fast_distance(&self.coords[ref_bead, XX], &self.coords[bid, XX]) <= self.cutoff,
            #      ))

        return results