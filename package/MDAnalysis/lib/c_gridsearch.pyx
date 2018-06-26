# -*- coding: utf-8; Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
#  Copyright (C) 2013-2016  Sébastien Buchoux <sebastien.buchoux@gmail.com>
#
#    This file is part of FATSLiM.
#
#    FATSLiM is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    FATSLiM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with FATSLiM.  If not, see <http://www.gnu.org/licenses/>.
#cython: cdivision=True
#cython: boundscheck=False

# Preprocessor DEFs
DEF DIM = 3
DEF XX = 0
DEF YY = 1
DEF ZZ = 2
DEF RET_OK = 1
DEF RET_ERROR = 0
DEF EPSILON = 1e-5
DEF NEIGHBORHOOD_ALLOCATION_INCREMENT = 50
DEF GRID_ALLOCATION_INCREMENT = 50

DEF BOX_MARGIN=1.0010
DEF MAX_NTRICVEC=12


# Cython C imports (no Python here!)
from cython.parallel cimport prange
from libc.stdlib cimport malloc, realloc, free, abort
from libc.stdio cimport fprintf, stderr
from libc.math cimport sqrt
from libc.math cimport abs as real_abs

cimport openmp


# Python imports
import numpy as np
cimport numpy as np

# Ctypes
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
    real *beaddist
    ###
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

# noinspection PyNoneFunctionAssignment
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
        # We will only use single shifts, but we will check a few
        # more shifts to see if there is a limiting distance
        # above which we can not be sure of the correct distance.
        for kk in range(5):
            k = order[kk]

            for jj in range(5):
                j = order[jj]

                for ii in range(5):
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

                            d2old += sqrt(pos[d])
                            d2new += sqrt(pos[d] + trial[d])

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
                                            d2new_c += sqrt(pos[d] + trial[d] - shift*box[dd, d])

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
                                            print(box)
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

    cdef real[:, ::1]fast_put_atoms_in_bbox(self, real[:,::1] coords) nogil:
        cdef ns_int i, m, d, natoms, wd = 0
        cdef real[:,::1] bbox_coords

        natoms = coords.shape[0]
        with gil:
            if natoms == 0:
                bbox_coords = np.empty((0, DIM))
            else:
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

cdef ns_grid initialize_nsgrid(matrix box,
                               float cutoff) nogil:
    cdef ns_grid grid
    cdef ns_int i

    for i in range(DIM):
        grid.ncells[i] = <ns_int> (box[i][i] / cutoff)
        if grid.ncells[i] == 0:
            grid.ncells[i] = 1
        grid.cellsize[i] = box[i][i] / grid.ncells[i]

    grid.size = grid.ncells[XX] * grid.ncells[YY] * grid.ncells[ZZ]
    return grid

cdef ns_int populate_grid(ns_grid *grid,
                       real[:,::1] coords) nogil:
    cdef ns_int ncoords = coords.shape[0]
    cdef bint ret_val

    ret_val = populate_grid_array(grid,
                                  <rvec *> &coords[0, 0],
                                  ncoords)

    return ret_val

cdef ns_int populate_grid_array(ns_grid *grid,
                             rvec *coords,
                             ns_int ncoords) nogil:
    cdef ns_int i, cellindex = -1
    cdef ns_int grid_size = grid.size
    cdef ns_int *allocated_size = NULL

    if grid_size != grid.ncells[XX] * grid.ncells[YY] * grid.ncells[ZZ]: # Grid not initialized
        return RET_ERROR


    # Allocate memory
    grid.nbeads = <ns_int *> malloc(sizeof(ns_int) * grid_size)
    if grid.nbeads == NULL:
        fprintf(stderr,"FATAL: Could not allocate memory for NS grid.nbeads (requested: %i bytes)\n",
                sizeof(ns_int) * grid_size)
        abort()

    allocated_size = <ns_int *> malloc(sizeof(ns_int) * grid_size)
    if allocated_size == NULL:
        fprintf(stderr,"FATAL: Could not allocate memory for NS allocated_size (requested: %i bytes)\n",
                sizeof(ns_int) * grid_size)
        abort()

    for i in range(grid_size):
        grid.nbeads[i] = 0
        allocated_size[i] = GRID_ALLOCATION_INCREMENT

    grid.beadids = <ns_int **> malloc(sizeof(ns_int *) * grid_size)
    if grid.beadids == NULL:
        fprintf(stderr,"FATAL: Could not allocate memory for NS grid.beadids (requested: %i bytes)\n",
                sizeof(ns_int *) * grid_size)
        abort()

    for i in range(grid_size):
        grid.beadids[i] = <ns_int *> malloc(sizeof(ns_int) * allocated_size[i])
        if grid.beadids[i] == NULL:
            fprintf(stderr,"FATAL: Could not allocate memory for NS grid.beadids[i] (requested: %i bytes)\n",
                sizeof(ns_int) * allocated_size[i])
            abort()

    # Get cell indices for coords
    for i in range(ncoords):
        cellindex = <ns_int> (coords[i][ZZ] / grid.cellsize[ZZ]) * (grid.ncells[XX] * grid.ncells[YY]) +\
                    <ns_int> (coords[i][YY] / grid.cellsize[YY]) * grid.ncells[XX] + \
                    <ns_int> (coords[i][XX] / grid.cellsize[XX])

        grid.beadids[cellindex][grid.nbeads[cellindex]] = i
        grid.nbeads[cellindex] += 1

        if grid.nbeads[cellindex] >= allocated_size[cellindex]:
            allocated_size[cellindex] += GRID_ALLOCATION_INCREMENT
            grid.beadids[cellindex] = <ns_int *> realloc(<void *> grid.beadids[cellindex], sizeof(ns_int) * allocated_size[cellindex])
    free(allocated_size)
    return RET_OK

cdef void destroy_nsgrid(ns_grid *grid) nogil:
    cdef ns_int i
    if grid.nbeads != NULL:
        free(grid.nbeads)

    for i in range(grid.size):
        if grid.beadids[i] != NULL:
            free(grid.beadids[i])
    free(grid.beadids)


cdef ns_neighborhood_holder *create_neighborhood_holder() nogil:
    cdef ns_neighborhood_holder *holder

    holder = <ns_neighborhood_holder *> malloc(sizeof(ns_neighborhood_holder))

    return holder

cdef void free_neighborhood_holder(ns_neighborhood_holder *holder) nogil:
    cdef ns_int i

    if holder == NULL:
        return

    for i in range(holder.size):
        if holder.neighborhoods[i].beadids != NULL:
            free(holder.neighborhoods[i].beadids)
        free(holder.neighborhoods[i])
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
        abort()

    neighborhood.size = 0
    neighborhood.allocated_size = NEIGHBORHOOD_ALLOCATION_INCREMENT
    neighborhood.beadids = <ns_int *> malloc(NEIGHBORHOOD_ALLOCATION_INCREMENT * sizeof(ns_int))
    ###Modified here
    neighborhood.beaddist = <real *> malloc(NEIGHBORHOOD_ALLOCATION_INCREMENT * sizeof(real))
    ### 
    if neighborhood.beadids == NULL:
        abort()

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
                        ### Modified here
                        neighborhood.beaddist[neighborhood.size] = d2
                        ###
                        neighborhood.size += 1
                        
                        if neighborhood.size >= neighborhood.allocated_size:
                            neighborhood.allocated_size += NEIGHBORHOOD_ALLOCATION_INCREMENT
                            neighborhood.beadids = <ns_int *> realloc(<void*> neighborhood.beadids, neighborhood.allocated_size * sizeof(ns_int))
                            ###Modified here
                            neighborhood.beaddist = <real *> realloc(<void*> neighborhood.beaddist, neighborhood.allocated_size * sizeof(real))
                            ###
                            if neighborhood.beadids == NULL:
                                abort()
                            ###Modified
                            if neighborhood.beaddist == NULL:
                                abort()
                            ###
                # Register the cell as checked
                already_checked[nchecked] = cell_index
                nchecked += 1

    return neighborhood


cdef ns_neighborhood_holder *ns_core_parallel(real[:, ::1] refcoords,
                                     real[:, ::1] neighborcoords,
                                     ns_grid *grid,
                                     PBCBox box,
                                     real cutoff,
                                     int nthreads=-1) nogil:
    cdef ns_int coordid, i, j
    cdef ns_int ncoords = refcoords.shape[0]
    cdef ns_int ncoords_neighbors = neighborcoords.shape[0]
    cdef real cutoff2 = cutoff * cutoff
    cdef ns_neighborhood_holder *holder

    cdef ns_int *neighbor_buf
    cdef ns_int buf_size, ibuf

    if nthreads < 0:
        nthreads = openmp.omp_get_num_threads()

    holder = create_neighborhood_holder()
    if holder == NULL:
        fprintf(stderr,"FATAL: Could not allocate memory for NS holder\n",
                sizeof(ns_int) * ncoords)
        abort()

    holder.size = ncoords
    holder.neighborhoods = <ns_neighborhood **> malloc(sizeof(ns_neighborhood *) * ncoords)
    if holder.neighborhoods == NULL:
        fprintf(stderr,"FATAL: Could not allocate memory for NS holder.neighborhoods (requested: %i bytes)\n",
                sizeof(ns_neighborhood) * ncoords)
        abort()

    # Here starts the real core and the iteration over coordinates
    for coordid in prange(ncoords, schedule='dynamic', num_threads=nthreads):
        holder.neighborhoods[coordid] = retrieve_neighborhood(&refcoords[coordid, XX],
                                                              neighborcoords,
                                                              grid,
                                                              box,
                                                              cutoff2)
        holder.neighborhoods[coordid].cutoff = cutoff

    return holder

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
        fprintf(stderr,"FATAL: Could not allocate memory for NS holder\n",
                sizeof(ns_int) * ncoords)
        abort()

    holder.size = ncoords
    holder.neighborhoods = <ns_neighborhood **> malloc(sizeof(ns_neighborhood *) * ncoords)
    if holder.neighborhoods == NULL:
        fprintf(stderr,"FATAL: Could not allocate memory for NS holder.neighborhoods (requested: %i bytes)\n",
                sizeof(ns_neighborhood) * ncoords)
        abort()

    # Here starts the real core and the iteration over coordinates
    for coordid in range(ncoords):
        holder.neighborhoods[coordid] = retrieve_neighborhood(&refcoords[coordid, XX],
                                                              neighborcoords,
                                                              grid,
                                                              box,
                                                              cutoff2)
        holder.neighborhoods[coordid].cutoff = cutoff

    return holder


# Python interface
cdef class FastNS(object):
    cdef PBCBox box
    cdef readonly int nthreads
    cdef readonly real[:, ::1] coords
    cdef real[:, ::1] coords_bbox
    cdef readonly real cutoff
    cdef bint prepared
    cdef ns_grid *grid


    def __init__(self, box):
        if box.shape != (3, 3):
            raise ValueError("Box must be provided as triclinic_dimensions (a 3x3 numpy.ndarray of unit cell vectors")

        self.box = PBCBox(box)

        self.nthreads = 1

        self.coords = None
        self.coords_bbox = None

        self.cutoff = -1

        self.prepared = False

        self.grid = <ns_grid *> malloc(sizeof(ns_grid))


    def __dealloc__(self):
        #destroy_nsgrid(self.grid)
        self.grid.size = 0

        #free(self.grid)

    def set_nthreads(self, nthreads, silent=False):
        import multiprocessing

        if nthreads > multiprocessing.cpu_count():
            print("Warning: the number of threads requested if greater than the number of cores available. Performances may not be optimal!")

        if not silent:
            print("Number of threads for NS adjusted to {}.".format(nthreads))

        self.nthreads = nthreads


    def set_coords(self, real[:, ::1] coords):
        self.coords = coords

        # Make sure atoms are inside the brick-shaped box
        self.coords_bbox = self.box.fast_put_atoms_in_bbox(coords)

        self.prepared = False


    def set_cutoff(self, real cutoff):
        self.cutoff = cutoff

        self.prepared = False


    def prepare(self, force=False):
        cdef ns_int i
        cdef bint initialization_ok

        if self.prepared and not force:
            print("NS already prepared, nothing to do!")

        if self.coords is None:
            raise ValueError("Coordinates must be set before NS preparation!")

        if self.cutoff < 0:
            raise ValueError("Cutoff must be set before NS preparation!")

        with nogil:
            initialization_ok = False

            # Initializing grid
            for i in range(DIM):
                self.grid.ncells[i] = <ns_int> (self.box.c_pbcbox.box[i][i] / self.cutoff)
                if self.grid.ncells[i] == 0:
                    self.grid.ncells[i] = 1
                self.grid.cellsize[i] = self.box.c_pbcbox.box[i][i] / self.grid.ncells[i]

            self.grid.size = self.grid.ncells[XX] * self.grid.ncells[YY] * self.grid.ncells[ZZ]

            # Populating grid
            if populate_grid(self.grid, self.coords_bbox) == RET_OK:
                initialization_ok = True


        if initialization_ok:
            self.prepared = True
        else:
            raise RuntimeError("Could not initialize NS grid")


    def search(self, real[:, ::1]search_coords, return_ids=False):
        cdef real[:, ::1] search_coords_bbox
        cdef ns_int nid, i, j
        cdef ns_neighborhood_holder *holder
        cdef ns_neighborhood *neighborhood

        if not self.prepared:
            self.prepare()


        # Make sure atoms are inside the brick-shaped box
        search_coords_bbox = self.box.fast_put_atoms_in_bbox(search_coords)


        with nogil:
            # Retrieve neighbors from grid
            if self.nthreads == 1:
                holder = ns_core(search_coords_bbox, self.coords_bbox, self.grid, self.box, self.cutoff)
            else:
                holder = ns_core_parallel(search_coords_bbox, self.coords_bbox, self.grid, self.box, self.cutoff, self.nthreads)



        neighbors = []
        ###Modify for distance
        sqdist = []
        indx = []
        ###
        for nid in range(holder.size):
            neighborhood = holder.neighborhoods[nid]

            if return_ids:
                neighborhood_py = np.empty(neighborhood.size, dtype=np.int64)
                ###Modify for distance
                neighborhood_dis = np.empty(neighborhood.size, dtype=np.float32)
                neighborhood_indx = np.empty(neighborhood.size, dtype=np.int64)
                ###
                for i in range(neighborhood.size):
                    neighborhood_py[i] = neighborhood.beadids[i]
                    ###Modify for distance
                    neighborhood_dis[i] = neighborhood.beaddist[i]
                    ###
            else:
                neighborhood_py = np.empty((neighborhood.size, DIM), dtype=np.float32)
                ###Modify for distance
                neighborhood_dis = np.empty((neighborhood.size), dtype=np.float32)
                neighborhood_indx = np.empty(neighborhood.size, dtype=np.int64)
                ###
                for i in range(neighborhood.size):
                    ###Modify for distance
                    neighborhood_dis[i] = neighborhood.beaddist[i]
                    neighborhood_indx[i] = neighborhood.beadids[i]
                    ###

                    for j in range(DIM):
                        neighborhood_py[i,j] = self.coords[neighborhood.beadids[i], j]
                    
            neighbors.append(neighborhood_py)
            sqdist.append(neighborhood_dis)
            indx.append(neighborhood_indx)

        # Free Memory
        free_neighborhood_holder(holder)

        return neighbors, sqdist, indx



__version__ = "26"