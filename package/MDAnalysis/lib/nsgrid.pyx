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


# Used to handle memory allocation
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from libc.math cimport sqrt
import numpy as np
cimport numpy as np
from libcpp.vector cimport vector


ctypedef np.int_t ns_int
ctypedef np.float32_t real
ctypedef real rvec[DIM]
ctypedef ns_int ivec[DIM]
ctypedef real matrix[DIM][DIM]

ctypedef vector[ns_int] intvec
ctypedef vector[real] realvec

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
    cdef bint is_triclinic
    cdef bint periodic

    def __init__(self, real[:,::1] box, bint periodic):
        self.periodic = periodic
        self.update(box)


    cdef void fast_update(self, real[:,::1] box) nogil:

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

        for i in range(DIM):
            dx[i] = other[i] - ref[i]

        if self.periodic:
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
        cdef ns_int i, m, d, natoms
        cdef real[:,::1] bbox_coords

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

    def put_atoms_in_bbox(self, real[:,::1] coords):
        return np.asarray(self.fast_put_atoms_in_bbox(coords))

#########################
# Neighbor Search Stuff #
#########################

cdef class NSResults(object):
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
        self.cutoff = cutoff
        self.coords = coords
        self.searchcoords = searchcoords

        self.npairs = 0
     
    cdef void add_neighbors(self, ns_int beadid_i, ns_int beadid_j, real distance2) nogil:
        # Important: If this function returns ERROR, it means that memory allocation failed

        self.pairs_buffer.push_back(beadid_i)
        self.pairs_buffer.push_back(beadid_j)
        self.pair_distances2_buffer.push_back(distance2)
        self.npairs += 1
    
    def get_pairs(self):
        return np.asarray(self.pairs_buffer).reshape(self.npairs, 2)

    def get_pair_distances(self):
        self.pair_distances_buffer = np.sqrt(self.pair_distances2_buffer)
        return np.asarray(self.pair_distances_buffer)
    
    cdef void create_buffers(self) nogil:
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
        if self.indices_buffer.empty():
            self.create_buffers()
        return np.ascontiguousarray(self.indices_buffer)

    def get_distances(self):
        if self.distances_buffer.empty():
            self.create_buffers()
        return np.ascontiguousarray(self.distances_buffer)
    
cdef class NSGrid(object):
    cdef readonly real cutoff  # cutoff
    cdef ns_int size  # total cells
    cdef ns_int ncoords  # number of coordinates
    cdef ns_int[DIM] ncells # individual cells in every dimension
    cdef ns_int[DIM] cell_offsets # Cell Multipliers
    cdef real[DIM] cellsize # cell size in every dimension
    cdef ns_int nbeads_per_cell # maximum beads
    cdef ns_int *nbeads # size (Number of beads in every cell)
    cdef ns_int *beadids # size * nbeads_per_cell (Beadids in every cell)
    cdef ns_int *cellids # ncoords (Cell occupation id for every atom)
    cdef bint force  # To negate the effects of optimized cutoff

    def __init__(self, ncoords, cutoff, PBCBox box, max_size, force=False):
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
        return <ns_int> (coord[ZZ] / self.cellsize[ZZ]) * (self.cell_offsets[ZZ]) +\
               <ns_int> (coord[YY] / self.cellsize[YY]) * self.cell_offsets[YY] + \
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
   
    cdef PBCBox box
    cdef real[:, ::1] coords
    cdef real[:, ::1] coords_bbox
    cdef readonly real cutoff
    cdef NSGrid grid
    cdef ns_int max_gridsize
    cdef bint periodic
    

    def __init__(self, cutoff, coords, box=None, max_gridsize=5000):
        import MDAnalysis as mda
        from MDAnalysis.lib.mdamath import triclinic_vectors

        cdef real[:] pseudobox = np.zeros(6, dtype=np.float32)
        cdef real[DIM] bmax, bmin
        cdef ns_int i
        self.periodic = True

        if (box is None) or np.allclose(box[:3], 0.):
            bmax = np.max(coords, axis=0)
            bmin = np.min(coords, axis=0)
            for i in range(DIM):
                pseudobox[i] = 1.1*(bmax - bmin)
                pseudobox[DIM + i] = 90.
            box = pseudobox
            # shift the origin
            coords -= bmin
            self.periodic = False



        if box.shape != (3,3):
            box = triclinic_vectors(box)

        self.box = PBCBox(box, self.periodic)

        if cutoff < 0:
            raise ValueError("Cutoff must be positive!")
        if cutoff * cutoff > self.box.c_pbcbox.max_cutoff2:
            raise ValueError("Cutoff greater than maximum cutoff ({:.3f}) given the PBC")
        
        self.coords = coords.copy()

        self.coords_bbox =  self.box.fast_put_atoms_in_bbox(coords)

        self.cutoff = cutoff
        self.max_gridsize = max_gridsize
        # Note that self.cutoff might be different from self.grid.cutoff
        self.grid = NSGrid(self.coords_bbox.shape[0], self.cutoff, self.box, self.max_gridsize)
        
        self.grid.fill_grid(self.coords_bbox)

        

    def search(self, search_coords):
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

        # Generate another grid to search
        searchcoords = np.asarray(search_coords, dtype=np.float32)
        searchcoords_bbox =  self.box.fast_put_atoms_in_bbox(searchcoords)
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
                                for m in range(DIM, -1, -1, -1):
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
                            #for this cellindex search in grid
                            for j in range(self.grid.nbeads[cellindex_probe]):
                                bid = self.grid.beadids[cellindex_probe * self.grid.nbeads_per_cell + j]
                                #find distance between search coords[i] and coords[bid]
                                d2 = self.box.fast_distance2(&searchcoords_bbox[current_beadid, XX], &self.coords_bbox[bid, XX])
                                if d2 < cutoff2 and d2 > EPSILON:
                                    results.add_neighbors(current_beadid, bid, d2)
                                    npairs += 1
        return results


    def self_search(self):
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
                            #Probe the search coordinates in a brick shaped box
                            probe[XX] = self.coords_bbox[current_beadid, XX] + (xi - 1) * self.grid.cellsize[XX]
                            probe[YY] = self.coords_bbox[current_beadid, YY] + (yi - 1) * self.grid.cellsize[XX]
                            probe[ZZ] = self.coords_bbox[current_beadid, ZZ] + (zi - 1) * self.grid.cellsize[XX]
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
                                for m in range(DIM, -1, -1, -1):
                                    if probe[m] < 0:
                                        check = False
                                        break
                                    elif probe[m] > self.box.c_pbcbox.box[m][m]:
                                        check = False
                                        break
                            if not check:
                                continue
                            # Get the cell index corresponding to the probe
                            cellindex_probe = self.grid.coord2cellid(probe)
                            #for this cellindex search in grid
                            for j in range(self.grid.nbeads[cellindex_probe]):
                                bid = self.grid.beadids[cellindex_probe * self.grid.nbeads_per_cell + j]
                                if bid < current_beadid:
                                    continue
                                #find distance between search coords[i] and coords[bid]
                                d2 = self.box.fast_distance2(&self.coords_bbox[current_beadid, XX], &self.coords_bbox[bid, XX])
                                if d2 < cutoff2 and d2 > EPSILON:
                                    results.add_neighbors(current_beadid, bid, d2)
                                    results.add_neighbors(bid, current_beadid, d2)
                                    npairs += 1
        return results
