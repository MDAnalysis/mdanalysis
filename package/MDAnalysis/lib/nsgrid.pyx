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

Classes
-------
"""

# Used to handle memory allocation
import numpy as np
from libcpp.vector cimport vector
from libc.math cimport floor

DEF END = -1

DEF XX = 0
DEF XY = 3
DEF YY = 4
DEF XZ = 6
DEF YZ = 7
DEF ZZ = 8


cdef extern from "calc_distances.h" nogil:
    void minimum_image(double* x, float* box, float* inverse_box)
    void minimum_image_triclinic(double* dx, float* box)


cdef class NSResults(object):
    """Class to store the results

    All outputs from :class:`FastNS` are stored in an instance of this class.
    All methods of :class:`FastNS` return an instance of this class, which can
    be used to  generate the desired results on demand.
    """

    cdef vector[int] pairs
    cdef vector[double] distances2

    def __init__(self):
        pass

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


cdef class _NSGrid(object):
    """Constructs a uniform cuboidal grid for a brick-shaped box

    The domain is subdivided into number of cells based on the desired search
    radius.  An optimization of cutoff is imposed to limit the size of data
    structure such that the cellsize is always greater than or
    equal to cutoff distance.

    .. warning::
        This class is not meant to be used by end users.

    """
    cdef float[:, ::1] coords_bbox

    cdef int size  # total cells
    cdef int[3] ncells  # individual cells in every dimension
    cdef int[3] cell_offsets  # Cell Multipliers
    # cellsize MUST be double precision, otherwise coord2cellid() may fail for
    # coordinates very close to the upper box boundaries! See Issue #2132
    cdef double[9] cellsize  # cell size in every dimension
    cdef double max_cutoff  # maximum allowable cutoff
    
    cdef int[::1] head_id  # first coord id for a given cell
    cdef int[::1] next_id  # next coord id after a given cell
    cdef bint triclinic
    cdef float[6] dimensions
    cdef float[3] inverse_dimensions
    cdef float[9] triclinic_dimensions
    # are we periodic in the X, Y and Z dimension?
    cdef bint periodic[3]
    
    def __init__(self, float[:, ::1] coords, double cutoff, box, bint pbc):
        """
        Parameters
        ----------
        coords : np.ndarray float of shape (ncoords, 3)
            Coordinates to fill inside the brick shaped box
        cutoff : float
            Minimum desired cutoff radius
        box : numpy ndarray shape=(6,)
            Box info, [lx, ly, lz, alpha, beta, gamma]
        pbc : bool
            is this NSGrid periodic at all?
        """
        cdef int i, j
        cdef int xyz[3]
        cdef float coordcopy[3]

        if box.shape != (6,):
            raise ValueError("Box must be a numpy array of [lx, ly, lz, alpha, beta, gamma]")

        box = box.astype(np.float32)
        
        from MDAnalysis.lib.mdamath import triclinic_vectors

        for i in range(3):
            self.dimensions[i] = box[i]
            self.inverse_dimensions[i] = 1.0 / box[i]
        self.dimensions[3] = box[3]
        self.dimensions[4] = box[4]
        self.dimensions[5] = box[5]

        self.triclinic_dimensions = triclinic_vectors(box).reshape((9,))
        self.triclinic = (self.triclinic_dimensions[XY] > 0 or
                          self.triclinic_dimensions[XZ] > 0 or
                          self.triclinic_dimensions[YZ] > 0)
        if cutoff < 0:
            raise ValueError("Cutoff must be positive")
        cutoff = max(cutoff, 0.01)

        self.ncells[0] = <int> floor(self.triclinic_dimensions[XX] / (cutoff + 0.001))
        self.ncells[1] = <int> floor(self.triclinic_dimensions[YY] / (cutoff + 0.001))
        self.ncells[2] = <int> floor(self.triclinic_dimensions[ZZ] / (cutoff + 0.001))
        # If there aren't enough cells in a given dimension it's equivalent to one
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

        self.cellsize[XX] = self.triclinic_dimensions[XX] / <double> self.ncells[0]
        # [YX] and [ZX] are 0
        self.cellsize[XY] = self.triclinic_dimensions[XY] / <double> self.ncells[1]
        self.cellsize[YY] = self.triclinic_dimensions[YY] / <double> self.ncells[1]
        # [ZY] is zero
        self.cellsize[XZ] = self.triclinic_dimensions[XZ] / <double> self.ncells[2]
        self.cellsize[YZ] = self.triclinic_dimensions[YZ] / <double> self.ncells[2]
        self.cellsize[ZZ] = self.triclinic_dimensions[ZZ] / <double> self.ncells[2]

        self.max_cutoff = self.cellsize[ZZ]
        if self.cellsize[YY] < self.max_cutoff:
            self.max_cutoff = self.cellsize[YY]
        if self.cellsize[XX] < self.max_cutoff:
            self.max_cutoff = self.cellsize[XX]

        self.size = self.ncells[0] * self.ncells[1] * self.ncells[2]
        
        self.cell_offsets[0] = 0
        self.cell_offsets[1] = self.ncells[0]
        self.cell_offsets[2] = self.ncells[0] * self.ncells[1]

        # Assign coordinates into cells
        # Linked list for each cell
        self.head_id = np.full(self.size, END, dtype=np.int32)
        self.next_id = np.full(coords.shape[0], END, dtype=np.int32)

        self.coords_bbox = np.empty_like(coords)
        for i in range(coords.shape[0]):
            self.coords_bbox[i][0] = coords[i][0]
            self.coords_bbox[i][1] = coords[i][1]
            self.coords_bbox[i][2] = coords[i][2]

            self.coordintoprimarycell(&self.coords_bbox[i][0])

            j = self.coord2cellid(&self.coords_bbox[i][0])

            self.next_id[i] = self.head_id[j]
            self.head_id[j] = i

    cdef int coord2cellid(self, const float* coord) nogil:
        """Finds the cell-id for the given coordinate inside the brick shaped box

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
        cdef float dydz, dxdz, dxdy

        dydz = self.triclinic_dimensions[YZ] / self.triclinic_dimensions[ZZ]
        dxdz = self.triclinic_dimensions[XZ] / self.triclinic_dimensions[ZZ]
        dxdy = self.triclinic_dimensions[XY] / self.triclinic_dimensions[YY]

        # This assumes coordinate is inside the primary unit cell
        xyz[2] = <int> (coord[2] / self.cellsize[ZZ])
        xyz[1] = <int> ((coord[1] - coord[2] * dydz) / self.cellsize[YY])
        xyz[0] = <int> ((coord[0] - coord[1] * dxdy - coord[2] * dydz) / self.cellsize[XX])

    cdef void coordintoprimarycell(self, float* coord) nogil:
        cdef float dydz, dxdz, dxdy
        cdef float adj_cy, adj_cx

        dydz = self.triclinic_dimensions[YZ] / self.triclinic_dimensions[ZZ]
        dxdz = self.triclinic_dimensions[XZ] / self.triclinic_dimensions[ZZ]
        dxdy = self.triclinic_dimensions[XY] / self.triclinic_dimensions[YY]

        while (coord[2] >= self.triclinic_dimensions[ZZ]):
            coord[2] -= self.triclinic_dimensions[ZZ]
            coord[1] -= self.triclinic_dimensions[YZ]
            coord[0] -= self.triclinic_dimensions[XZ]
        while (coord[2] < 0):
            coord[2] += self.triclinic_dimensions[ZZ]
            coord[1] += self.triclinic_dimensions[YZ]
            coord[0] += self.triclinic_dimensions[XZ]

        adj_cy = coord[2] * dydz
        while ((coord[1] - adj_cy) >= self.triclinic_dimensions[YY]):
            coord[1] -= self.triclinic_dimensions[YY]
            coord[0] -= self.triclinic_dimensions[XY]
        while ((coord[1] - adj_cy) < 0):
            coord[1] += self.triclinic_dimensions[YY]
            coord[0] += self.triclinic_dimensions[XY]

        adj_cx = coord[2] * dxdz + coord[1] * dxdy
        while ((coord[0] - adj_cx) >= self.triclinic_dimensions[XX]):
            coord[0] -= self.triclinic_dimensions[XX]
        while ((coord[0] - adj_cx) < 0):
            coord[0] += self.triclinic_dimensions[XX]

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

    cdef double calc_distsq(self, const float* a, const float* b):
        cdef double dx[3]

        dx[0] = a[0] - b[0]
        dx[1] = a[1] - b[1]
        dx[2] = a[2] - b[2]

        if self.triclinic:
            minimum_image_triclinic(dx, &self.triclinic_dimensions[0])
        else:
            minimum_image(dx, &self.dimensions[0], &self.inverse_dimensions[0])

        return dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]


cdef class FastNS(object):
    """Grid based search between two group of atoms

    Minimum image convention is used for distance evaluations
    if pbc is set to ``True``.
    """
    cdef readonly double cutoff

    cdef _NSGrid grid

    def __init__(self, cutoff, coords, box, pbc=True):
        """
        Initialize the grid and sort the coordinates in respective
        cells by shifting the coordinates in a brick shaped box.
        The brick shaped box is defined by :class:`_PBCBox`
        and cuboidal grid is initialize by :class:`_NSGrid`.
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
        box : numpy.ndarray or None
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

        if np.allclose(box[:3], 0.0):
            raise ValueError("Any of the box dimensions cannot be 0")

        self.cutoff = cutoff
        # Note that self.cutoff might be different from self.grid.cutoff
        # due to optimization
        self.grid = _NSGrid(coords, self.cutoff, box, pbc)
        
        if self.cutoff > self.grid.max_cutoff:
            raise ValueError("Cutoff greater than maximum cutoff ({:.3f}) given the PBC"
                             "".format(self.grid.max_cutoff))

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
           can be accessed by its methods :meth:`~NSResults.get_indices`,
           :meth:`~NSResults.get_distances`, :meth:`~NSResults.get_pairs`, and
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
        cdef int cellcoord[3], searchcoord[3]
        cdef float tmpcoord[3]
        
        cdef NSResults results
        cdef double d2, cutoff2

        cutoff2 = self.cutoff * self.cutoff

        if (search_coords.ndim != 2 or search_coords.shape[1] != 3):
            raise ValueError("search_coords must have a shape of (n, 3), got "
                             "{}.".format(search_coords.shape))

        size_search = search_coords.shape[0]
        results = NSResults()

        for i in range(size_search):
            for j in range(3):
                tmpcoord[j] = search_coords[i][j]
            self.grid.coordintoprimarycell(tmpcoord)
            # which cell is atom *i* in
            self.grid.coord2cellxyz(&tmpcoord[0], &cellcoord[0])
            # loop over all 27 neighbouring cells
            for xi in range(3):
                for yi in range(3):
                    for zi in range(3):
                        cx = cellcoord[0] - 1 + xi
                        cy = cellcoord[1] - 1 + yi
                        cz = cellcoord[2] - 1 + zi
                        cellid = self.grid.cellxyz2cellid(cx, cy, cz)

                        if cellid == END:  # out of bounds
                            continue
                        # for loop over atoms in searchcoord
                        j = self.grid.head_id[cellid]
                        while (j != END):
                            d2 = self.grid.calc_distsq(&tmpcoord[0], &self.grid.coords_bbox[j][0])
                            if d2 <= cutoff2:
                                # place self.coords index first then search_coords
                                results.add_neighbors(j, i, d2)
                            j = self.grid.next_id[j]

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
           can be accessed by its methods :meth:`~NSResults.get_indices`,
           :meth:`~NSResults.get_distances`, :meth:`~NSResults.get_pairs`, and
           :meth:`~NSResults.get_pair_distances`.
        """
        cdef int cx, cy, cz, ox, oy, oz
        cdef int ci, cj, i, j, nj
        cdef int cellindex, cellindex_probe
        cdef int xi, yi, zi
        cdef int[:, ::1] route
        cdef NSResults results = NSResults()
        cdef double d2
        cdef double cutoff2 = self.cutoff * self.cutoff

        # route over 13 neighbouring cells
        route = np.array([[1, 0, 0], [1, 1, 0], [0, 1, 0], [-1, 1, 0],
                          [1, 0, -1], [1, 1, -1], [0, 1, -1], [-1, 1, -1],
                          [1, 0, 1], [1, 1, 1], [0, 1, 1], [-1, 1, 1], [0, 0, 1]], dtype=np.int32)

        for cx in range(self.grid.ncells[0]):
            for cy in range(self.grid.ncells[1]):
                for cz in range(self.grid.ncells[2]):
                    ci = self.grid.cellxyz2cellid(cx, cy, cz)

                    i = self.grid.head_id[ci]
                    while (i != END):
                        # pairwise within this cell
                        j = self.grid.next_id[i]
                        while (j != END):
                            d2 = self.grid.calc_distsq(&self.grid.coords_bbox[i][0],
                                                       &self.grid.coords_bbox[j][0])
                            if d2 <= cutoff2:
                                results.add_neighbors(i, j, d2)
                            j = self.grid.next_id[j]

                        # loop over 13 neighbouring cells
                        for nj in range(13):
                            ox = cx + route[nj][0]
                            oy = cy + route[nj][1]
                            oz = cz + route[nj][2]

                            cj = self.grid.cellxyz2cellid(ox, oy, oz)
                            if cj == END:
                                continue

                            j = self.grid.head_id[cj]
                            while (j != END):
                                d2 = self.grid.calc_distsq(&self.grid.coords_bbox[i][0],
                                                           &self.grid.coords_bbox[j][0])
                                if d2 <= cutoff2:
                                    results.add_neighbors(i, j, d2)
                                j = self.grid.next_id[j]

                        # move to next position in cell *ci*
                        i = self.grid.next_id[i]


        return results
