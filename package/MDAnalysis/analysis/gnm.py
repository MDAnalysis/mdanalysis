# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

#Analyse a trajectory using elastic network models, following the approach of Hall et al (JACS 2007)
#Ben Hall (benjamin.a.hall@ucl.ac.uk) is to blame
#Copyright 2011; Consider under GPL v2 or later
r"""
Elastic network analysis of MD trajectories --- :mod:`MDAnalysis.analysis.gnm`
==============================================================================

:Author: Benjamin Hall <benjamin.a.hall@ucl.ac.uk>
:Year: 2011
:Copyright: GNU Public License v2 or later


Analyse a trajectory using elastic network models, following the approach of
:footcite:p:`Hall2007`.

An example is provided in the MDAnalysis Cookbook_, listed as GNMExample_.

.. _GNMExample: https://github.com/MDAnalysis/MDAnalysisCookbook/blob/master/examples/GNMExample.py
.. _Cookbook: https://github.com/MDAnalysis/MDAnalysisCookbook

The basic approach is to pass a trajectory to :class:`GNMAnalysis` and then run
the analysis:

.. code-block:: python

    u = MDAnalysis.Universe(PSF, DCD)
    C = MDAnalysis.analysis.gnm.GNMAnalysis(u, ReportVector="output.txt")

    C.run()
    output = zip(*C.results)

    with open("eigenvalues.dat", "w") as outputfile:
        for item in output[1]:
            outputfile.write(item + "\n")


The results are found in :attr:`GNMAnalysis.results`, which can be
used for further processing (see :footcite:p:`Hall2007`).

.. rubric:: References

.. footbibliography::


Analysis tasks
--------------

.. autoclass:: GNMAnalysis
   :members:
.. autoclass:: closeContactGNMAnalysis
   :members:

Utility functions
-----------------

The following functions are used internally and are typically not
directly needed to perform the analysis.

.. autofunction:: generate_grid
.. autofunction:: order_list

.. versionchanged:: 0.16.0
   removed unused function :func:`backup_file`

"""
import itertools
import logging
import warnings

import numpy as np

from .base import AnalysisBase, ResultsGroup


from MDAnalysis.analysis.base import Results

logger = logging.getLogger('MDAnalysis.analysis.GNM')


def _dsq(a, b):
    diff = (a - b)
    return np.dot(diff, diff)


def generate_grid(positions, cutoff):
    """Simple grid search.

    An alternative to searching the entire list of each atom; divide the
    structure into `cutoff` sized boxes This way, for each particle you only need
    to search the neighbouring boxes to find the particles within the `cutoff`.

    Observed a 6x speed up for a smallish protein with ~300 residues; this
    should get better with bigger systems.

    Parameters
    ----------
    positions : array
        coordinates of the atoms
    cutoff : float
        find particles with distance less than `cutoff` from each other; the
        grid will consist of boxes with sides of at least length `cutoff`

    """
    positions = np.asarray(positions)

    x, y, z = positions.T
    high_x = x.max()
    high_y = y.max()
    high_z = z.max()
    low_x = x.min()
    low_y = y.min()
    low_z = z.min()
    #Ok now generate a list with 3 dimensions representing boxes in x, y and z
    grid = [[[[] for i in range(int((high_z - low_z) / cutoff) + 1)]
             for j in range(int((high_y - low_y) / cutoff) + 1)]
            for k in range(int((high_x - low_x) / cutoff) + 1)]
    for i, pos in enumerate(positions):
        x_pos = int((pos[0] - low_x) / cutoff)
        y_pos = int((pos[1] - low_y) / cutoff)
        z_pos = int((pos[2] - low_z) / cutoff)
        grid[x_pos][y_pos][z_pos].append(i)
    return grid


def neighbour_generator(positions, cutoff):
    """
    return atom pairs that are in neighboring regions of space from a verlet-grid

    Parameters
    ----------
    positions : ndarray
        atom positions
    cutoff : float
        size of grid box

    Yields
    ------
    i_atom, j_atom
        indices of close atom pairs
    """
    grid = generate_grid(positions, cutoff)
    n_x = len(grid)
    n_y = len(grid[0])
    n_z = len(grid[0][0])
    for cell_x, cell_y, cell_z in itertools.product(
            range(n_x), range(n_y), range(n_z)):
        atoms = grid[cell_x][cell_y][cell_z]
        # collect all atoms in own cell and neighboring cell
        all_atoms = []
        nei_cells = (-1, 0, 1)
        for x, y, z in itertools.product(nei_cells, nei_cells, nei_cells):
            gx = cell_x + x
            gy = cell_y + y
            gz = cell_z + z
            if 0 <= gx < n_x and 0 <= gy < n_y and 0 <= gz < n_z:
                all_atoms += grid[gx][gy][gz]
        # return all possible atom pairs in current cell
        for i_atom in atoms:
            for j_atom in all_atoms:
                yield i_atom, j_atom


def order_list(w):
    """Returns a dictionary showing the order of eigenvalues (which are reported scrambled normally)"""
    ordered = list(w)
    unordered = list(w)
    ordered.sort()
    list_map = {}
    for i in range(len(w)):
        list_map[i] = unordered.index(ordered[i])
    return list_map


class GNMAnalysis(AnalysisBase):
    """Basic tool for GNM analysis.

    Each frame is treated as a novel structure and the GNM
    calculated.  By default, this stores the dominant eigenvector
    and its associated eigenvalue; either can be used to monitor
    conformational change in a simulation.

    Parameters
    ----------
    universe : Universe
          Analyze the full trajectory in the universe.
    select : str (optional)
          MDAnalysis selection string
    cutoff : float (optional)
          Consider selected atoms within the cutoff as neighbors for the
          Gaussian network model.
    ReportVector : str (optional)
          filename to write eigenvectors to, by default no output is written
    Bonus_groups : tuple
          This is a tuple of selection strings that identify additional groups
          (such as ligands). The center of mass of each group will be added as
          a single point in the ENM (it is a popular way of treating small
          ligands such as drugs). You need to ensure that none of the atoms in
          `Bonus_groups` is contained in `selection` as this could lead to
          double counting. No checks are applied.

    Attributes
    ----------
    results.times : numpy.ndarray
            simulation times used in analysis
    results.eigenvalues : numpy.ndarray
            calculated eigenvalues
    results.eigenvectors : numpy.ndarray
            calculated eigenvectors

    See Also
    --------
    :class:`closeContactGNMAnalysis`


    .. versionchanged:: 0.16.0
       Made :meth:`generate_output` a private method :meth:`_generate_output`.

    .. versionchanged:: 1.0.0
       Changed `selection` keyword to `select`

    .. versionchanged:: 2.0.0
       Use :class:`~MDAnalysis.analysis.AnalysisBase` as parent class and
       store results as attributes ``times``, ``eigenvalues`` and
       ``eigenvectors`` of the ``results`` attribute.

    .. versionchanged:: 2.8.0
       introduced a :meth:`get_supported_backends` allowing for execution on with
       ``multiprocessing`` and ``dask`` backends.
    """

    _analysis_algorithm_is_parallelizable = True

    @classmethod
    def get_supported_backends(cls):
        return ("serial", "multiprocessing", "dask")
    
    def __init__(self,
                 universe,
                 select='protein and name CA',
                 cutoff=7.0,
                 ReportVector=None,
                 Bonus_groups=None):
        super(GNMAnalysis, self).__init__(universe.trajectory)
        self.u = universe
        self.select = select
        self.cutoff = cutoff
        self.results = Results()
        self.results.eigenvalues = []
        self.results.eigenvectors = []
        self._timesteps = None  # time for each frame
        self.ReportVector = ReportVector
        self.Bonus_groups = [self.u.select_atoms(item) for item in Bonus_groups] \
                            if Bonus_groups else []
        self.ca = self.u.select_atoms(self.select)

    def _generate_output(self, w, v, outputobject,
                         ReportVector=None, counter=0):
        """Appends time, eigenvalues and eigenvectors to results.

        This generates the output by adding eigenvalue and
        eigenvector data to an appendable object and optionally
        printing some of the results to file. This is the function
        to replace if you want to generate a more complex set of
        outputs
        """
        list_map = order_list(w)
        if ReportVector:
            with open(ReportVector, "a") as oup:
                for item in enumerate(v[list_map[1]]):
                    print(
                        "",
                        counter,
                        item[0] + 1,
                        w[list_map[1]],
                        item[1],
                        file=oup)

        outputobject.eigenvalues.append(w[list_map[1]])
        outputobject.eigenvectors.append(v[list_map[1]])

    def generate_kirchoff(self):
        """Generate the Kirchhoff matrix of contacts.

        This generates the neighbour matrix by generating a grid of
        near-neighbours and then calculating which are are within
        the cutoff.

        Returns
        -------
        array
                the resulting Kirchhoff matrix
        """
        positions = self.ca.positions

        #add the com from each bonus group to the ca_positions list
        for item in self.Bonus_groups:
            #bonus = self.u.select_atoms(item)
            positions = np.vstack((positions, item.center_of_mass()))

        natoms = len(positions)
        matrix = np.zeros((natoms, natoms), np.float64)

        cutoffsq = self.cutoff**2

        for i_atom, j_atom in neighbour_generator(positions, self.cutoff):
            if j_atom > i_atom and _dsq(positions[i_atom],
                                        positions[j_atom]) < cutoffsq:
                matrix[i_atom][j_atom] = -1.0
                matrix[j_atom][i_atom] = -1.0
                matrix[i_atom][i_atom] = matrix[i_atom][i_atom] + 1
                matrix[j_atom][j_atom] = matrix[j_atom][j_atom] + 1

        return matrix

    def _single_frame(self):
        matrix = self.generate_kirchoff()
        try:
            _, w, v = np.linalg.svd(matrix)
        except np.linalg.LinAlgError:
            msg = f"SVD with cutoff {self.cutoff} failed to converge. "
            msg += f"Skip frame at {self._ts.time}."
            warnings.warn(msg)
            logger.warning(msg)
            return
        # Save the results somewhere useful in some useful format. Usefully.
        self._generate_output(
            w,
            v,
            self.results,
            ReportVector=self.ReportVector,
            counter=self._ts.frame)

    def _conclude(self):
        self.results.times = self.times
        self.results.eigenvalues = np.asarray(self.results.eigenvalues)
        self.results.eigenvectors = np.asarray(self.results.eigenvectors)

    def _get_aggregator(self):
        return ResultsGroup(
            lookup={
                "eigenvectors": ResultsGroup.ndarray_hstack,
                "eigenvalues": ResultsGroup.ndarray_hstack,
                "times": ResultsGroup.ndarray_hstack,
            }
        )


class closeContactGNMAnalysis(GNMAnalysis):
    r"""GNMAnalysis only using close contacts.

    This is a version of the GNM where the Kirchoff matrix is
    constructed from the close contacts between individual atoms
    in different residues.

    Parameters
    ----------
    universe : Universe
          Analyze the full trajectory in the universe.
    select : str (optional)
          MDAnalysis selection string
    cutoff : float (optional)
          Consider selected atoms within the cutoff as neighbors for the
          Gaussian network model.
    ReportVector : str (optional)
          filename to write eigenvectors to, by default no output is written
    weights : {"size", None} (optional)
          If set to "size" (the default) then weight the contact by
          :math:`1/\sqrt{N_i N_j}` where :math:`N_i` and :math:`N_j` are the
          number of atoms in the residues :math:`i` and :math:`j` that contain
          the atoms that form a contact.

    Attributes
    ----------
    results.times : numpy.ndarray
            simulation times used in analysis
    results.eigenvalues : numpy.ndarray
            calculated eigenvalues
    results.eigenvectors : numpy.ndarray
            calculated eigenvectors

    Notes
    -----
    The `MassWeight` option has now been removed.

    See Also
    --------
    :class:`GNMAnalysis`


    .. versionchanged:: 0.16.0
       Made :meth:`generate_output` a private method :meth:`_generate_output`.

    .. deprecated:: 0.16.0
       Instead of ``MassWeight=True`` use ``weights="size"``.

    .. versionchanged:: 1.0.0
       MassWeight option (see above deprecation entry).
       Changed `selection` keyword to `select`

    .. versionchanged:: 2.0.0
       Use :class:`~MDAnalysis.analysis.AnalysisBase` as parent class and
       store results as attributes ``times``, ``eigenvalues`` and
       ``eigenvectors`` of the `results` attribute.
    """

    def __init__(self,
                 universe,
                 select='protein',
                 cutoff=4.5,
                 ReportVector=None,
                 weights="size"):
        super(closeContactGNMAnalysis, self).__init__(universe,
                                                      select,
                                                      cutoff,
                                                      ReportVector)
        self.weights = weights

    def generate_kirchoff(self):
        nresidues = self.ca.n_residues
        positions = self.ca.positions
        residue_index_map = [
            resnum
            for [resnum, residue] in enumerate(self.ca.residues)
            for atom in residue.atoms
        ]
        matrix = np.zeros((nresidues, nresidues), dtype=np.float64)
        cutoffsq = self.cutoff**2

        # cache sqrt of residue sizes (slow) so that sr[i]*sr[j] == sqrt(r[i]*r[j])
        inv_sqrt_res_sizes = np.ones(len(self.ca.residues))
        if self.weights == 'size':
            inv_sqrt_res_sizes = 1 / np.sqrt(
                [r.atoms.n_atoms for r in self.ca.residues])

        for i_atom, j_atom in neighbour_generator(positions, self.cutoff):
            if j_atom > i_atom and _dsq(positions[i_atom],
                                        positions[j_atom]) < cutoffsq:
                iresidue = residue_index_map[i_atom]
                jresidue = residue_index_map[j_atom]
                contact = (inv_sqrt_res_sizes[iresidue] *
                           inv_sqrt_res_sizes[jresidue])
                matrix[iresidue][jresidue] -= contact
                matrix[jresidue][iresidue] -= contact
                matrix[iresidue][iresidue] += contact
                matrix[jresidue][jresidue] += contact

        return matrix
