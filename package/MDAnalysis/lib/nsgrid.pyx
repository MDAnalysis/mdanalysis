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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
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
.. versionchanged:: 1.0.2
   Rewrote module
.. versionchanged:: 2.1.0
   Capped max grid dimension to 1290, which when cubed is the max value of
   a 32 bit signed integer.

Classes
-------
"""
import numpy as np

from libcpp.vector cimport vector
from libc cimport math

cdef int END = -1

cdef int XX = 0
cdef int XY = 3
cdef int YY = 4
cdef int XZ = 6
cdef int YZ = 7
cdef int ZZ = 8

# Cube root of the maximum size of a 32 bit signed integer. If the system is divided into more
# grids than this, integer overflow will occur.
cdef int  MAX_GRID_DIM = 1290

ctypedef float coordinate[3]

cdef extern from "calc_distances.h" nogil:
    void _minimum_image_ortho_lazy(double* x, float* box, float* half_box)
    void minimum_image_triclinic(double* dx, float* box)
    void _ortho_pbc(coordinate* coords, int numcoords, float* box)
    void _triclinic_pbc(coordinate* coords, int numcoords, float* box)


cdef inline float fmax(float a, float b):
    if a > b:
        return a
    else:
        return b


cdef inline float fmin(float a, float b):
    if a < b:
        return a
    else:
        return b

cdef inline float degsin(float deg):
    # sin in degrees
    deg *= math.M_PI / 180.
    return math.sin(deg)


cdef class NSResults(object):
    """Class to store the results

    All outputs from :class:`FastNS` are stored in an instance of this class.
    All methods of :class:`FastNS` return an instance of this class, which can
    be used to  generate the desired results on demand.
    """

    cdef vector[int] pairs
    cdef vector[double] distances2

    cdef void add_neighbors(self, int beadid_i, int beadid_j, double distance2) nogil:
        """Internal function to add pairs and distances to buffers

        The buffers populated using this method are used by
        other methods of this class. This is the
        primary function used by :class:`FastNS` to save all
        the pair of atoms,
        which are considered as neighbors.
        """

        self.pairs.push_back(beadid_i)
        self.pairs.push_back(beadid_j)
        self.distances2.push_back(distance2)

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
        return np.asarray(self.pairs, dtype=np.intp).reshape(-1, 2)

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
        :meth:`~NSResults.get_pairs`

        """
        dist2 = np.asarray(self.distances2)

        return np.sqrt(dist2)


cdef class FastNS(object):
    """Grid based search between positions

    Minimum image convention is used for distance evaluations
    if pbc is set to ``True``.

    .. versionchanged::  1.0.2
       Rewrote to fix bugs with triclinic boxes
    """
    cdef readonly double cutoff
    cdef float[:, ::1] coords_bbox

    cdef int[3] ncells  # individual cells in every dimension
    cdef int[3] cell_offsets  # Cell Multipliers
    # cellsize MUST be double precision, otherwise coord2cellid() may fail for
    # coordinates very close to the upper box boundaries! See Issue #2132
    # diagonal stores the cell width, off diagonal elements are the "tilt"
    # i.e. cellsize[3] is the dxdy tilt (box[XY] / box[YY])
    cdef double[9] cellsize

    cdef int[::1] head_id  # first coord id for a given cell
    cdef int[::1] next_id  # next coord id after a given cell
    cdef bint triclinic
    cdef float[6] dimensions
    cdef float[3] half_dimensions
    cdef float[3] inverse_dimensions
    cdef float[9] triclinic_dimensions
    # are we periodic in the X, Y and Z dimension?
    cdef bint pbc  # periodic at all?
    cdef bint periodic[3]

    def __init__(self, cutoff, coords, box, pbc=True):
        """
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

        .. code-block:: python

            lmax = all_coords.max(axis=0)
            lmin = all_coords.min(axis=0)
            pseudobox[:3] = 1.1*(lmax - lmin)
            pseudobox[3:] = 90.
            shift = all_coords.copy()
            shift -= lmin
            gridsearch = FastNS(max_cutoff, shift, box=pseudobox, pbc=False)

        """
        if (coords.ndim != 2 or coords.shape[1] != 3):
            raise ValueError("coords must have a shape of (n, 3), got {}."
                             "".format(coords.shape))
        if box.shape != (6,):
            raise ValueError("Box must be a numpy array of [lx, ly, lz, alpha, beta, gamma], got {}"
                             "".format(box))
        if (box[:3] == 0.0).any():
            raise ValueError("Any of the box dimensions cannot be 0")
        if cutoff < 0:
            raise ValueError("Cutoff must be positive")
        self.cutoff = cutoff
        max_cutoff = self._prepare_box(box, pbc)
        if cutoff > max_cutoff:
            raise ValueError("Cutoff {} too large for box (max {})".format(cutoff, max_cutoff))
        self._pack_grid(coords)

    cdef float _prepare_box(self, box, bint pbc):
        """
        Parameters
        ----------
        box : numpy ndarray shape=(6,)
            Box info, [lx, ly, lz, alpha, beta, gamma]
        pbc : bool
            is this NSGrid periodic at all?

        Returns
        -------
        max_cutoff : float
           the maximum allowable cutoff given the box shape and size
        """
        cdef float cutoff, min_cellsize, max_cutoff, new_cellsize
        cdef int i

        from MDAnalysis.lib.mdamath import triclinic_vectors

        for i in range(3):
            self.dimensions[i] = box[i]
            self.half_dimensions[i] = 0.5 * box[i]
            self.inverse_dimensions[i] = 1.0 / box[i]
        self.dimensions[3] = box[3]
        self.dimensions[4] = box[4]
        self.dimensions[5] = box[5]

        self.triclinic_dimensions = triclinic_vectors(box).reshape((9,))
        self.triclinic = (self.triclinic_dimensions[XY] != 0 or
                          self.triclinic_dimensions[XZ] != 0 or
                          self.triclinic_dimensions[YZ] != 0)
        cutoff = max(self.cutoff, 1.0)  # TODO: Figure out max ncells and stick to that
        # Find maximum allowable cutoff
        # For orthogonal cell this is just half the shortest boxlength (sine values all 1.0)
        # For triclinic the box tilt creates a shorter path across the box we must account for
        # alpha
        max_cutoff = self.triclinic_dimensions[YY] * degsin(self.dimensions[3])
        max_cutoff = fmin(max_cutoff, self.triclinic_dimensions[ZZ] * degsin(self.dimensions[3]))
        # beta
        max_cutoff = fmin(max_cutoff, self.triclinic_dimensions[XX] * degsin(self.dimensions[4]))
        max_cutoff = fmin(max_cutoff, self.triclinic_dimensions[ZZ] * degsin(self.dimensions[4]))
        # gamma
        max_cutoff = fmin(max_cutoff, self.triclinic_dimensions[XX] * degsin(self.dimensions[5]))
        max_cutoff = fmin(max_cutoff, self.triclinic_dimensions[YY] * degsin(self.dimensions[5]))
        max_cutoff /= 2

        # for triclinic cells, we need to worry about the shortest path across the cells
        min_cellsize = cutoff
        if self.triclinic:
            for i in range(3, 6):
                # cutoff/sin(theta) to elongate the XX/YY/ZZ dimension to make smallest diagonal large enough
                new_cellsize = cutoff / degsin(self.dimensions[i])
                min_cellsize = fmax(new_cellsize, min_cellsize)

        # add 0.001 here to avoid floating point errors
        # will make cells slightly too large as a result, ah well
        min_cellsize += 0.001
        # If the cell size is too small, indexing overflow will occur. Limit the number
        # of cells in any dimension to the cube root of the maximum of 32 bit integer values.
        self.ncells[0] = <int> min(math.floor(self.triclinic_dimensions[XX] / min_cellsize), MAX_GRID_DIM)
        self.ncells[1] = <int> min(math.floor(self.triclinic_dimensions[YY] / min_cellsize), MAX_GRID_DIM)
        self.ncells[2] = <int> min(math.floor(self.triclinic_dimensions[ZZ] / min_cellsize), MAX_GRID_DIM)

        self.pbc = pbc
        # If there aren't enough cells in a given dimension it's equivalent to one
        # this prevents double counting of results if a cell would have the same
        # neighbour above and below
        if pbc:
            for i in range(3):
                if self.ncells[i] <= 3:
                    self.ncells[i] = 1
                    self.periodic[i] = False
                else:
                    self.periodic[i] = True
        else:
            for i in range(3):
                if self.ncells[i] <= 2:
                    self.ncells[i] = 1
                self.periodic[i] = False
        # Off diagonal cellsizes are actually the tilt, i.e. dy/dz or similar
        self.cellsize[XX] = self.triclinic_dimensions[XX] / <double> self.ncells[0]
        # [YX] and [ZX] are 0
        self.cellsize[XY] = self.triclinic_dimensions[XY] / self.triclinic_dimensions[YY]
        self.cellsize[YY] = self.triclinic_dimensions[YY] / <double> self.ncells[1]
        # [ZY] is zero
        self.cellsize[XZ] = self.triclinic_dimensions[XZ] / self.triclinic_dimensions[ZZ]
        self.cellsize[YZ] = self.triclinic_dimensions[YZ] / self.triclinic_dimensions[ZZ]
        self.cellsize[ZZ] = self.triclinic_dimensions[ZZ] / <double> self.ncells[2]

        self.cell_offsets[0] = 0
        self.cell_offsets[1] = self.ncells[0]
        self.cell_offsets[2] = self.ncells[0] * self.ncells[1]

        return max_cutoff

    cdef void _pack_grid(self, float[:, :] coords):
        """Assigns coordinates into cells

        Parameters
        ----------
        coords : np.ndarray float of shape (ncoords, 3)
            Coordinates to populate the box
        """
        cdef int i, j

        # Linked list for each cell
        # Starting coordinate index for each cell (END if empty cell)
        self.head_id = np.full(self.cell_offsets[2] * self.ncells[2], END, dtype=np.int32, order='C')
        # Next coordinate index in cell for each coordinate (END if end of sequence)
        self.next_id = np.full(coords.shape[0], END, dtype=np.int32, order='C')

        self.coords_bbox = coords.copy()
        with nogil:
            if self.triclinic:
                _triclinic_pbc(<coordinate*>&self.coords_bbox[0][0],
                               self.coords_bbox.shape[0],
                               &self.triclinic_dimensions[0])
            else:
                _ortho_pbc(<coordinate*>&self.coords_bbox[0][0],
                           self.coords_bbox.shape[0],
                           &self.dimensions[0])

            for i in range(self.coords_bbox.shape[0]):
                j = self.coord2cellid(&self.coords_bbox[i][0])
                self.next_id[i] = self.head_id[j]
                self.head_id[j] = i

    cdef int coord2cellid(self, const float* coord) nogil:
        """Finds the cell-id for the given coordinate

        Note
        ----
        Assumes the coordinate is already inside the primary unit cell.
        Return wrong cell-id if this is not the case
        """
        cdef int xyz[3]

        self.coord2cellxyz(coord, xyz)

        return xyz[0] + xyz[1] * self.cell_offsets[1] + xyz[2] * self.cell_offsets[2]

    cdef void coord2cellxyz(self, const float* coord, int* xyz) nogil:
        """Calculate cell coordinate for coord"""
        # This assumes coordinate is inside the primary unit cell
        xyz[2] = <int> (coord[2] / self.cellsize[ZZ])
        xyz[1] = <int> ((coord[1] - coord[2] * self.cellsize[YZ]) / self.cellsize[YY])
        xyz[0] = <int> ((coord[0] - coord[1] * self.cellsize[XY]
                         - coord[2] * self.cellsize[XZ]) / self.cellsize[XX])
        # Make sure cell coordinate indices are within the primary unit cell
        # (better safe than sorry):
        xyz[0] %= self.ncells[0]
        xyz[1] %= self.ncells[1]
        xyz[2] %= self.ncells[2]

    cdef int cellxyz2cellid(self, int cx, int cy, int cz) nogil:
        """Convert cell coordinate to cell id, END for out of bounds"""
        if cx < 0:
            if self.periodic[0]:
                cx = self.ncells[0] - 1
            else:
                return END
        elif cx == self.ncells[0]:
            if self.periodic[0]:
                cx = 0
            else:
                return END
        if cy < 0:
            if self.periodic[1]:
                cy = self.ncells[1] - 1
            else:
                return END
        elif cy == self.ncells[1]:
            if self.periodic[1]:
                cy = 0
            else:
                return END
        if cz < 0:
            if self.periodic[2]:
                cz = self.ncells[2] - 1
            else:
                return END
        elif cz == self.ncells[2]:
            if self.periodic[2]:
                cz = 0
            else:
                return END

        return cx + cy * self.cell_offsets[1] + cz * self.cell_offsets[2]

    cdef double calc_distsq(self, const float* a, const float* b) nogil:
        cdef double dx[3]

        dx[0] = a[0] - b[0]
        dx[1] = a[1] - b[1]
        dx[2] = a[2] - b[2]

        if self.pbc:
            if self.triclinic:
                minimum_image_triclinic(dx, &self.triclinic_dimensions[0])
            else:
                _minimum_image_ortho_lazy(dx, &self.dimensions[0],
                                          &self.half_dimensions[0])

        return dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]

    def search(self, float[:, :] search_coords):
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
        results : NSResults
           An :class:`NSResults` object holding neighbor search results, which
           can be accessed by its methods :meth:`~NSResults.get_pairs` and
           :meth:`~NSResults.get_pair_distances`.

        Note
        ----
        For non-PBC aware calculations, the current implementation doesn't work
        if any of the query coordinates lies outside the `box` supplied to
        :class:`~MDAnalysis.lib.nsgrid.FastNS`.
        """
        cdef int i, j, size_search
        cdef int cx, cy, cz
        cdef int cellid
        cdef int xi, yi, zi
        cdef int cellcoord[3]
        cdef float tmpcoord[3]

        cdef NSResults results = NSResults()
        cdef double d2, cutoff2

        cutoff2 = self.cutoff * self.cutoff

        if (search_coords.ndim != 2 or search_coords.shape[1] != 3):
            raise ValueError("search_coords must have a shape of (n, 3), got "
                             "{}.".format(search_coords.shape))

        with nogil:
            size_search = search_coords.shape[0]
            for i in range(size_search):
                tmpcoord[0] = search_coords[i][0]
                tmpcoord[1] = search_coords[i][1]
                tmpcoord[2] = search_coords[i][2]
                if self.triclinic:
                    _triclinic_pbc(<coordinate*>&tmpcoord[0], 1,
                                   &self.triclinic_dimensions[0])
                else:
                    _ortho_pbc(<coordinate*>&tmpcoord[0], 1,
                               &self.dimensions[0])
                # which cell is atom *i* in
                self.coord2cellxyz(&tmpcoord[0], cellcoord)
                # loop over all 27 neighbouring cells
                for xi in range(3):
                    for yi in range(3):
                        for zi in range(3):
                            cx = cellcoord[0] - 1 + xi
                            cy = cellcoord[1] - 1 + yi
                            cz = cellcoord[2] - 1 + zi
                            cellid = self.cellxyz2cellid(cx, cy, cz)

                            if cellid == END:  # out of bounds
                                continue
                            # for loop over atoms in searchcoord
                            j = self.head_id[cellid]
                            while (j != END):
                                d2 = self.calc_distsq(&tmpcoord[0],
                                                      &self.coords_bbox[j][0])
                                if d2 <= cutoff2:
                                    # place search_coords then self.bbox_coords
                                    results.add_neighbors(i, j, d2)
                                j = self.next_id[j]
        return results

    def self_search(self):
        """Searches all the pairs within the initialized coordinates

        All the pairs among the initialized coordinates are registered
        in hald the time. Although the algorithm is still the same, but
        the distance checks can be reduced to half in this particular case
        as every pair need not be evaluated twice.

        Returns
        -------
        results : NSResults
           An :class:`NSResults` object holding neighbor search results, which
           can be accessed by its methods :meth:`~NSResults.get_pairs` and
           :meth:`~NSResults.get_pair_distances`.
        """
        cdef int cx, cy, cz, ox, oy, oz
        cdef int ci, cj, i, j, nj
        cdef NSResults results = NSResults()
        cdef double d2
        cdef double cutoff2 = self.cutoff * self.cutoff
        # route over 13 neighbouring cells
        cdef int[13][3] route = [[1, 0, 0], [1, 1, 0], [0, 1, 0], [-1, 1, 0],
                                 [1, 0, -1], [1, 1, -1], [0, 1, -1], [-1, 1, -1],
                                 [1, 0, 1], [1, 1, 1], [0, 1, 1], [-1, 1, 1], [0, 0, 1]]

        for cx in range(self.ncells[0]):
            for cy in range(self.ncells[1]):
                for cz in range(self.ncells[2]):
                    ci = self.cellxyz2cellid(cx, cy, cz)

                    i = self.head_id[ci]
                    while (i != END):
                        # pairwise within this cell
                        j = self.next_id[i]
                        while (j != END):
                            d2 = self.calc_distsq(&self.coords_bbox[i][0],
                                                  &self.coords_bbox[j][0])
                            if d2 <= cutoff2:
                                results.add_neighbors(i, j, d2)
                            j = self.next_id[j]

                        # loop over 13 neighbouring cells
                        for nj in range(13):
                            ox = cx + route[nj][0]
                            oy = cy + route[nj][1]
                            oz = cz + route[nj][2]

                            cj = self.cellxyz2cellid(ox, oy, oz)
                            if cj == END:
                                continue

                            j = self.head_id[cj]
                            while (j != END):
                                d2 = self.calc_distsq(&self.coords_bbox[i][0],
                                                      &self.coords_bbox[j][0])
                                if d2 <= cutoff2:
                                    results.add_neighbors(i, j, d2)
                                j = self.next_id[j]

                        # move to next position in cell *ci*
                        i = self.next_id[i]

        return results
