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

#cython: cdivision=True
#cython: boundscheck=False
#cython: initializedcheck=False

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

DEF NEIGHBORHOOD_ALLOCATION_INCREMENT = 50

DEF BOX_MARGIN=1.0010
DEF MAX_NTRICVEC=12

from libc.stdlib cimport malloc, realloc, free
from libc.math cimport sqrt

import numpy as np
cimport numpy as np

ctypedef np.int_t ns_int
ctypedef np.float32_t real
ctypedef real rvec[DIM]
ctypedef real matrix[DIM][DIM]

cdef struct ns_grid:
    ns_int size
    ns_int[DIM] ncells
    real[DIM] cellsize
    ns_int *nbeads
    ns_int **beadids

cdef struct ns_neighborhood:
    real cutoff
    ns_int allocated_size
    ns_int size
    ns_int *beadids

cdef struct ns_neighborhood_holder:
    ns_int size
    ns_neighborhood **neighborhoods

# Useful stuff

cdef real rvec_norm2(const rvec a) nogil:
    return a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ]

cdef void rvec_clear(rvec a) nogil:
    a[XX]=0.0
    a[YY]=0.0
    a[ZZ]=0.0


cdef struct cPBCBox_t:
    matrix     box
    rvec       fbox_diag
    rvec       hbox_diag
    rvec       mhbox_diag
    real       max_cutoff2
    ns_int        ntric_vec
    ns_int[DIM]   tric_shift[MAX_NTRICVEC]
    real[DIM]  tric_vec[MAX_NTRICVEC]


# Class to handle PBC calculations
cdef class PBCBox(object):
    cdef cPBCBox_t c_pbcbox
    cdef rvec center
    cdef rvec bbox_center

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
        for i in range(DIM):
            for j in range(DIM):
                self.c_pbcbox.box[i][j] = box[i, j]
                self.center[j] += 0.5 * box[i, j]
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

        # Update shift vectors
        self.c_pbcbox.ntric_vec = 0

        # We will only use single shifts
        for kk in range(3):
            k = order[kk]

            for jj in range(3):
                j = order[jj]

                for ii in range(3):
                    i = order[ii]

                    # A shift is only useful when it is trilinic
                    if j != 0 or k != 0:
                        d2old = 0
                        d2new = 0

                        for d in range(DIM):
                            trial[d] = i*box[XX, d] + j*box[YY, d] + k*box[ZZ, d]

                            # Choose the vector within the brick around 0,0,0 that
                            # will become the shortest due to shift try.

                            if d == DIM:
                                trial[d] = 0
                                pos[d] = 0
                            else:
                                if trial[d] < 0:
                                    pos[d] = min(self.c_pbcbox.hbox_diag[d], -trial[d])
                                else:
                                    pos[d] = max(-self.c_pbcbox.hbox_diag[d], -trial[d])

                            d2old += pos[d]**2
                            d2new += (pos[d] + trial[d])**2

                        if BOX_MARGIN*d2new < d2old:
                            if  not (j < -1 or j > 1 or k < -1 or k > 1):
                                use = True

                                for dd in range(DIM):
                                    if dd == 0:
                                        shift = i
                                    elif dd == 1:
                                        shift = j
                                    else:
                                        shift = k

                                    if shift:
                                        d2new_c = 0

                                        for d in range(DIM):
                                            d2new_c += (pos[d] + trial[d] - shift*box[dd, d])**2

                                        if d2new_c <= BOX_MARGIN*d2new:
                                            use = False

                                if use: # Accept this shift vector.
                                    if self.c_pbcbox.ntric_vec >= MAX_NTRICVEC:
                                        with gil:
                                            print("\nWARNING: Found more than %d triclinic "
                                                  "correction vectors, ignoring some."
                                                  % MAX_NTRICVEC)
                                            print("  There is probably something wrong with "
                                                  "your box.")
                                            print(np.array(box))

                                            for i in range(self.c_pbcbox.ntric_vec):
                                                print(" -> shift #{}: [{}, {}, {}]".format(i+1,
                                                                                           self.c_pbcbox.tric_shift[i][XX],
                                                                                           self.c_pbcbox.tric_shift[i][YY],
                                                                                           self.c_pbcbox.tric_shift[i][ZZ]))
                                    else:
                                        for d in range(DIM):
                                            self.c_pbcbox.tric_vec[self.c_pbcbox.ntric_vec][d] = \
                                                trial[d]
                                        self.c_pbcbox.tric_shift[self.c_pbcbox.ntric_vec][XX] = i
                                        self.c_pbcbox.tric_shift[self.c_pbcbox.ntric_vec][YY] = j
                                        self.c_pbcbox.tric_shift[self.c_pbcbox.ntric_vec][ZZ] = k
                                        self.c_pbcbox.ntric_vec += 1


    def update(self, real[:,::1] box):
        if box.shape[0] != DIM or box.shape[1] != DIM:
            raise ValueError("Box must be a %i x %i matrix. (shape: %i x %i)" %
                             (DIM, DIM, box.shape[0], box.shape[1]))
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

    cdef real fast_distance2(self, rvec a, rvec b) nogil:
        cdef rvec dx
        self.fast_pbc_dx(a, b, dx)
        return rvec_norm2(dx)

    cdef real fast_distance(self, rvec a, rvec b) nogil:
        return sqrt(self.fast_distance2(a,b))

    cdef real[:, ::1]fast_put_atoms_in_bbox(self, real[:,::1] coords) nogil:
        cdef ns_int i, m, d, natoms, wd = 0
        cdef real[:,::1] bbox_coords

        natoms = coords.shape[0]
        with gil:
            bbox_coords = coords.copy()

        for i in range(natoms):
            for m in range(DIM - 1, -1, -1):
                while bbox_coords[i, m] < 0:
                    for d in range(m+1):
                        bbox_coords[i, d] += self.c_pbcbox.box[m][d]
                while bbox_coords[i, m] >= self.c_pbcbox.box[m][m]:
                    for d in range(m+1):
                        bbox_coords[i, d] -= self.c_pbcbox.box[m][d]
        return bbox_coords

    def put_atoms_in_bbox(self, real[:,::1] coords):
        if coords.shape[0] == 0:
            return np.zeros((0, DIM), dtype=np.float32)
        return np.asarray(self.fast_put_atoms_in_bbox(coords))

########################################################################################################################
#
# Neighbor Search Stuff
#
########################################################################################################################
cdef struct ns_grid:
    ns_int size
    ns_int[DIM] ncells
    real[DIM] cellsize
    ns_int *nbeads
    ns_int **beadids

cdef ns_neighborhood_holder *create_neighborhood_holder() nogil:
    cdef ns_neighborhood_holder *holder

    holder = <ns_neighborhood_holder *> malloc(sizeof(ns_neighborhood_holder))
    holder.size = 0
    holder.neighborhoods = NULL

    return holder

cdef void free_neighborhood_holder(ns_neighborhood_holder *holder) nogil:
    cdef ns_int i

    if holder == NULL:
        return

    for i in range(holder.size):
        if holder.neighborhoods[i].beadids != NULL:
            free(holder.neighborhoods[i].beadids)
        free(holder.neighborhoods[i])

    if holder.neighborhoods != NULL:
        free(holder.neighborhoods)
    free(holder)

cdef ns_neighborhood *retrieve_neighborhood(rvec current_coords, real[:, ::1]neighborcoords, ns_grid *grid, PBCBox box, real cutoff2) nogil:
    cdef ns_int d, m
    cdef ns_int xi, yi, zi, bid
    cdef real d2
    cdef rvec shifted_coord, dx, neighbor_coord, corrected_coords

    cdef bint already_checked[27]
    cdef bint skip
    cdef ns_int nchecked = 0, icheck
    cdef ns_int cell_index

    cdef ns_neighborhood *neighborhood = <ns_neighborhood *> malloc(sizeof(ns_neighborhood)) 
    if neighborhood == NULL:
        return NULL

    neighborhood.size = 0
    neighborhood.allocated_size = NEIGHBORHOOD_ALLOCATION_INCREMENT
    neighborhood.beadids = <ns_int *> malloc(NEIGHBORHOOD_ALLOCATION_INCREMENT * sizeof(ns_int))

    if neighborhood.beadids == NULL:
        free(neighborhood)
        return NULL

    for zi in range(3):
        for yi in range(3):
            for xi in range(3):
                # Calculate and/or reinitialize shifted coordinates
                shifted_coord[XX] = current_coords[XX] + (xi - 1) * grid.cellsize[XX]
                shifted_coord[YY] = current_coords[YY] + (yi - 1) * grid.cellsize[YY]
                shifted_coord[ZZ] = current_coords[ZZ] + (zi - 1) * grid.cellsize[ZZ]

                # Make sure the shifted coordinates is inside the brick-shaped box
                for m in range(DIM - 1, -1, -1):

                    while shifted_coord[m] < 0:
                        for d in range(m+1):
                            shifted_coord[d] += box.c_pbcbox.box[m][d]


                    while shifted_coord[m] >= box.c_pbcbox.box[m][m]:
                        for d in range(m+1):
                            shifted_coord[d] -= box.c_pbcbox.box[m][d]

                # Get the cell index corresponding to the coord
                cell_index = <ns_int> (shifted_coord[ZZ] / grid.cellsize[ZZ]) * (grid.ncells[XX] * grid.ncells[YY]) +\
                             <ns_int> (shifted_coord[YY] / grid.cellsize[YY]) * grid.ncells[XX] + \
                             <ns_int> (shifted_coord[XX] / grid.cellsize[XX])

                # Just a safeguard
                if cell_index >= grid.size:
                    continue

                # Check the cell index was not already selected
                skip = False
                for icheck in range(nchecked):
                    if already_checked[icheck] == cell_index:
                        skip = True
                        break
                if skip:
                    continue

                # Search for neighbors inside this cell
                for i_bead in range(grid.nbeads[cell_index]):
                    bid = grid.beadids[cell_index][i_bead]

                    box.fast_pbc_dx(current_coords, &neighborcoords[bid, XX], dx)

                    d2 = rvec_norm2(dx)

                    if d2 < cutoff2:
                        if d2 < EPSILON: # Don't add the current bead as its own neighbor!
                            continue

                        # Update neighbor lists
                        neighborhood.beadids[neighborhood.size] = bid
                        neighborhood.size += 1
                        
                        if neighborhood.size >= neighborhood.allocated_size:
                            neighborhood.allocated_size += NEIGHBORHOOD_ALLOCATION_INCREMENT
                            neighborhood.beadids = <ns_int *> realloc(<void*> neighborhood.beadids, neighborhood.allocated_size * sizeof(ns_int))

                            if neighborhood.beadids == NULL:
                                free(neighborhood)
                                return NULL

                # Register the cell as checked
                already_checked[nchecked] = cell_index
                nchecked += 1

    return neighborhood


cdef ns_neighborhood_holder *ns_core(real[:, ::1] refcoords,
                                     real[:, ::1] neighborcoords,
                                     ns_grid *grid,
                                     PBCBox box,
                                     real cutoff) nogil:
    cdef ns_int coordid, i, j
    cdef ns_int ncoords = refcoords.shape[0]
    cdef ns_int ncoords_neighbors = neighborcoords.shape[0]
    cdef real cutoff2 = cutoff * cutoff
    cdef ns_neighborhood_holder *holder

    cdef ns_int *neighbor_buf
    cdef ns_int buf_size, ibuf

    holder = create_neighborhood_holder()
    if holder == NULL:
        return NULL

    holder.neighborhoods = <ns_neighborhood **> malloc(sizeof(ns_neighborhood *) * ncoords)
    if holder.neighborhoods == NULL:
        free_neighborhood_holder(holder)
        return NULL

    # Here starts the real core and the iteration over coordinates
    for coordid in range(ncoords):
        holder.neighborhoods[coordid] = retrieve_neighborhood(&refcoords[coordid, XX],
                                                              neighborcoords,
                                                              grid,
                                                              box,
                                                              cutoff2)
        if holder.neighborhoods[coordid] == NULL:
            free_neighborhood_holder(holder)
            return NULL

        holder.neighborhoods[coordid].cutoff = cutoff
        holder.size += 1

    return holder

cdef class NSResults(object):
    """
    Class used to store results returned by `MDAnalysis.lib.grid.FastNS.search`
    """
    cdef PBCBox box
    cdef readonly real cutoff
    cdef np.ndarray grid_coords
    cdef np.ndarray ref_coords
    cdef ns_int **nids
    cdef ns_int *nsizes
    cdef ns_int size
    cdef list indices
    cdef list coordinates
    cdef list distances

    def __dealloc__(self):
        if self.nids != NULL:
            for i in range(self.size):
                if self.nids[i] != NULL:
                    free(self.nids[i])
            free(self.nids)

        if self.nsizes != NULL:
            free(self.nsizes)

    def __init__(self, PBCBox box, real cutoff):
        self.box = box
        self.cutoff = cutoff

        self.size = 0
        self.nids = NULL
        self.nsizes = NULL

        self.grid_coords = None
        self.ref_coords = None

        self.indices = None
        self.coordinates = None
        self.distances = None


    cdef populate(self, ns_neighborhood_holder *holder, grid_coords, ref_coords):
        cdef ns_int nid, i
        cdef ns_neighborhood *neighborhood

        self.grid_coords = np.asarray(grid_coords)
        self.ref_coords = np.asarray(ref_coords)

        # Allocate memory
        self.nsizes = <ns_int *> malloc(sizeof(ns_int) * holder.size)
        if self.nsizes == NULL:
            raise MemoryError("Could not allocate memory for NSResults")

        self.nids = <ns_int **> malloc(sizeof(ns_int *) * holder.size)
        if self.nids == NULL:
            raise MemoryError("Could not allocate memory for NSResults")

        for nid in range(holder.size):
            neighborhood = holder.neighborhoods[nid]

            self.nsizes[nid] = neighborhood.size

            self.nids[nid] = <ns_int *> malloc(sizeof(ns_int *) * neighborhood.size)
            if self.nids[nid] == NULL:
                raise MemoryError("Could not allocate memory for NSResults")

        with nogil:
            for nid in range(holder.size):
                neighborhood = holder.neighborhoods[nid]

                for i in range(neighborhood.size):
                    self.nids[nid][i] = neighborhood.beadids[i]

        self.size = holder.size

    def get_indices(self):
        """
        Return Neighbors indices.

        :return: list of indices
        """
        cdef ns_int i, nid, size

        if self.indices is None:
            indices = []

            for nid in range(self.size):
                size = self.nsizes[nid]

                tmp_incides = np.empty((size), dtype=np.int)

                for i in range(size):
                    tmp_incides[i] = self.nids[nid][i]

                indices.append(tmp_incides)

            self.indices = indices

        return self.indices


    def get_coordinates(self):
        """
        Return coordinates of neighbors.

        :return: list of coordinates
        """
        cdef ns_int i, nid, size, beadid

        if self.coordinates is None:
            coordinates = []

            for nid in range(self.size):
                size = self.nsizes[nid]

                tmp_values = np.empty((size, DIM), dtype=np.float32)

                for i in range(size):
                    beadid = self.nids[nid][i]
                    tmp_values[i] = self.grid_coords[beadid]

                coordinates.append(tmp_values)

            self.coordinates = coordinates

        return self.coordinates


    def get_distances(self):
        """
        Return coordinates of neighbors.

        :return: list of distances
        """
        cdef ns_int i, nid, size, j, beadid
        cdef rvec ref, other, dx
        cdef real dist
        cdef real[:, ::1] ref_coords = self.ref_coords
        cdef real[:, ::1] grid_coords = self.grid_coords

        if self.distances is None:
            distances = []

            for nid in range(self.size):
                size = self.nsizes[nid]

                tmp_values = np.empty((size), dtype=np.float32)
                ref = <rvec> &ref_coords[nid, 0]

                for i in range(size):
                    beadid = self.nids[nid][i]
                    other = <rvec> &grid_coords[beadid, 0]

                    tmp_values[i] = self.box.fast_distance(ref, other)

                distances.append(tmp_values)

            self.distances = distances

        return self.distances


# Python interface
cdef class FastNS(object):
    cdef PBCBox box
    cdef readonly real[:, ::1] coords
    cdef real[:, ::1] coords_bbox
    cdef readonly real cutoff
    cdef bint prepared
    cdef ns_grid grid

    def __init__(self, u):
        import MDAnalysis as mda
        from MDAnalysis.lib.mdamath import triclinic_vectors

        if not isinstance(u, mda.Universe):
            raise TypeError("FastNS class must be initialized with a valid MDAnalysis.Universe instance")
        box = triclinic_vectors(u.dimensions)

        if box.shape != (3, 3):
            raise ValueError("Box must be provided as triclinic_dimensions (a 3x3 numpy.ndarray of unit cell vectors")

        self.box = PBCBox(box)

        self.coords = None
        self.coords_bbox = None

        self.cutoff = -1

        self.prepared = False

        self.grid.size = 0


    def __dealloc__(self):
        cdef ns_int i
        # Deallocate NS grid
        if self.grid.nbeads != NULL:
            free(self.grid.nbeads)

        for i in range(self.grid.size):
            if self.grid.beadids[i] != NULL:
                free(self.grid.beadids[i])
        free(self.grid.beadids)

        self.grid.size = 0


    def set_coords(self, real[:, ::1] coords):
        self.coords = coords

        # Make sure atoms are inside the brick-shaped box
        self.coords_bbox = self.box.fast_put_atoms_in_bbox(coords)

        self.prepared = False


    def set_cutoff(self, real cutoff):
        self.cutoff = cutoff

        self.prepared = False


    def prepare(self, force=False):
        cdef ns_int i, cellindex = -1
        cdef ns_int *allocated_size = NULL
        cdef ns_int ncoords = self.coords.shape[0]
        cdef ns_int allocation_guess
        cdef rvec *coords = <rvec *>  &self.coords_bbox[0, 0]

        if self.prepared and not force:
            print("NS already prepared, nothing to do!")

        if self.coords is None:
            raise ValueError("Coordinates must be set before NS preparation!")

        if self.cutoff < 0:
            raise ValueError("Cutoff must be set before NS preparation!")

        with nogil:

            # Initializing grid
            for i in range(DIM):
                self.grid.ncells[i] = <ns_int> (self.box.c_pbcbox.box[i][i] / self.cutoff)
                if self.grid.ncells[i] == 0:
                    self.grid.ncells[i] = 1
                self.grid.cellsize[i] = self.box.c_pbcbox.box[i][i] / self.grid.ncells[i]
            self.grid.size = self.grid.ncells[XX] * self.grid.ncells[YY] * self.grid.ncells[ZZ]

            # This is just a guess on how much memory we might for each grid cell:
            # we just assume an average bead density and we take four times this density just to be safe
            allocation_guess = <ns_int> (4 * (ncoords / self.grid.size + 1))

            # Allocate memory for the grid
            self.grid.nbeads = <ns_int *> malloc(sizeof(ns_int) * self.grid.size)
            if self.grid.nbeads == NULL:
                with gil:
                    raise MemoryError("Could not allocate memory for NS grid")

            # Allocate memory from temporary allocation counter
            allocated_size = <ns_int *> malloc(sizeof(ns_int) * self.grid.size)
            if allocated_size == NULL:
                # No need to free grid.nbeads as it will be freed by destroy_nsgrid called by  __dealloc___
                with gil:
                    raise MemoryError("Could not allocate memory for allocation buffer")

            # Pre-allocate some memory for grid cells
            for i in range(self.grid.size):
                self.grid.nbeads[i] = 0
                allocated_size[i] = allocation_guess

            self.grid.beadids = <ns_int **> malloc(sizeof(ns_int *) * self.grid.size)
            if self.grid.beadids == NULL:
                with gil:
                    raise MemoryError("Could not allocate memory for grid cells")

            for i in range(self.grid.size):
                self.grid.beadids[i] = <ns_int *> malloc(sizeof(ns_int) * allocated_size[i])
                if self.grid.beadids[i] == NULL:
                    with gil:
                        raise MemoryError("Could not allocate memory for grid cell")

            # Populate grid cells using the coordinates (ie do the heavy work)
            for i in range(ncoords):
                cellindex = <ns_int> (coords[i][ZZ] / self.grid.cellsize[ZZ]) * (self.grid.ncells[XX] * self.grid.ncells[YY]) +\
                            <ns_int> (coords[i][YY] / self.grid.cellsize[YY]) * self.grid.ncells[XX] + \
                            <ns_int> (coords[i][XX] / self.grid.cellsize[XX])

                self.grid.beadids[cellindex][self.grid.nbeads[cellindex]] = i
                self.grid.nbeads[cellindex] += 1

                # We need to allocate more memory (simply double the amount of memory as
                # 1. it should barely be needed
                # 2. the size should stay fairly reasonable
                if self.grid.nbeads[cellindex] >= allocated_size[cellindex]:
                    allocated_size[cellindex] *= 2
                    self.grid.beadids[cellindex] = <ns_int *> realloc(<void *> self.grid.beadids[cellindex], sizeof(ns_int) * allocated_size[cellindex])

            # Now we can free the allocation buffer
            free(allocated_size)
        self.prepared = True


    def search(self, search_coords):
        cdef real[:, ::1] search_coords_bbox
        cdef real[:, ::1] search_coords_view
        cdef ns_int nid, i, j
        cdef ns_neighborhood_holder *holder
        cdef ns_neighborhood *neighborhood

        if not self.prepared:
            self.prepare()

        # Check the shape of search_coords as a array of 3D coords if needed
        shape = search_coords.shape
        if len(shape) == 1:
            if not shape[0] == 3:
                raise ValueError("Coordinates must be 3D")
            else:
                search_coords_view = np.array([search_coords,], dtype=np.float32)
        else:
            search_coords_view = search_coords

        # Make sure atoms are inside the brick-shaped box
        search_coords_bbox = self.box.fast_put_atoms_in_bbox(search_coords_view)

        with nogil:
            holder = ns_core(search_coords_bbox, self.coords_bbox, &self.grid, self.box, self.cutoff)

        if holder == NULL:
            raise MemoryError("Could not allocate memory to run NS core")

        results = NSResults(self.box, self.cutoff)
        results.populate(holder, self.coords, search_coords_view)

        # Free memory allocated to holder
        free_neighborhood_holder(holder)

        return results
