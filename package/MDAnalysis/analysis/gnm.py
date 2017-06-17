# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
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


Analyse a trajectory using elastic network models, following the approach of [Hall2007]_.

An example is provided in the MDAnalysis Cookbook_, listed as GNMExample_.

.. _GNMExample: https://github.com/MDAnalysis/MDAnalysisCookbook/blob/master/examples/GNMExample.py
.. _Cookbook: https://github.com/MDAnalysis/MDAnalysisCookbook

The basic approach is to pass a trajectory to :class:`GNMAnalysis` and then run
the analysis::

    u = MDAnalysis.Universe(PSF, DCD)
    C = MDAnalysis.analysis.gnm.GNMAnalysis(u, ReportVector="output.txt")

    C.run()
    output = zip(*C.results)

    with open("eigenvalues.dat", "w") as outputfile:
        for item in output[1]:
            outputfile.write(item + "\n")


The results are found in :attr:`GNMAnalysis.results`, which can be
used for further processing (see [Hall2007]_).

.. rubric:: References

.. [Hall2007]  Benjamin A. Hall, Samantha L. Kaye, Andy Pang, Rafael Perera, and
               Philip C. Biggin. Characterization of Protein Conformational
               States by Normal-Mode Frequencies. *JACS* 129 (2007), 11394--11401.


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
   removed un-unsed function :func:`backup_file`

"""

from __future__ import print_function, division, absolute_import
from six.moves import range

import numpy as np

import warnings
import logging

logger = logging.getLogger('MDAnalysis.analysis.GNM')


def _dsq(a, b):
    return ((a - b)**2).sum()

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
    natoms = len(positions)
    #Ok now generate a list with 3 dimensions representing boxes in x, y and z
    grid = [[[
        [] for i in range(int((high_z - low_z) / cutoff) + 1)] for j in range(int((high_y - low_y) / cutoff) + 1)]
        for k in range(int((high_x - low_x) / cutoff) + 1)]
    res_positions = []
    for i in range(natoms):
        x_pos = int((positions[i][0] - low_x) / cutoff)
        y_pos = int((positions[i][1] - low_y) / cutoff)
        z_pos = int((positions[i][2] - low_z) / cutoff)
        grid[x_pos][y_pos][z_pos].append(i)
        res_positions.append([x_pos, y_pos, z_pos])
    return (res_positions, grid, low_x, low_y, low_z)


def order_list(w):
    """Returns a dictionary showing the order of eigenvalues (which are reported scrambled normally)"""
    ordered = list(w)
    unordered = list(w)
    ordered.sort()
    list_map = {}
    for i in range(len(w)):
        list_map[i] = unordered.index(ordered[i])
    return list_map


class GNMAnalysis(object):
    """Basic tool for GNM analysis.

    Each frame is treated as a novel structure and the GNM
    calculated.  By default, this stores the dominant eigenvector
    and its associated eigenvalue; either can be used to monitor
    conformational change in a simulation.

    Parameters
    ----------
    universe : Universe
          Analyze the full trajectory in the universe.
    selection : str (optional)
          MDAnalysis selection string, default "protein and name CA"
    cutoff : float (optional)
          Consider selected atoms within the cutoff as neighbors for the
          Gaussian network model.
    ReportVector : str (optional)
          filename to write eigenvectors to, by default no output is written
          (``None``)
    Bonus_groups : tuple
          This is a tuple of selection strings that identify additional groups
          (such as ligands). The center of mass of each group will be added as
          a single point in the ENM (it is a popular way of treating small
          ligands such as drugs). You need to ensure that none of the atoms in
          `Bonus_groups` is contained in `selection` as this could lead to
          double counting. No checks are applied. Default is ``None``.

    See Also
    --------
    :class:`closeContactGNMAnalysis`


    .. versionchanged:: 0.16.0
       Made :meth:`generate_output` a private method :meth:`_generate_output`.

    """

    def __init__(self, universe, selection='protein and name CA', cutoff=7.0,
                 ReportVector=None, Bonus_groups=None):
        self.u = universe
        self.selection = selection
        self.cutoff = cutoff
        self.results = []  # final result
        self._timesteps = None  # time for each frame
        self.ReportVector = ReportVector
        self.Bonus_groups = [self.u.select_atoms(item) for item in Bonus_groups] \
                            if Bonus_groups else []
        self.ca = self.u.select_atoms(self.selection)

    def _generate_output(self, w, v, outputobject, time, matrix, nmodes=2,
                         ReportVector=None, counter=0):
        """Appends eigenvalues and eigenvectors to results.

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
                    print("", counter, time, item[0] + 1,
                          w[list_map[1]], item[1], file=oup)
        outputobject.append((time, w[list_map[1]], v[list_map[1]]))
        #outputobject.append((time, [ w[list_map[i]] for i in range(nmodes) ], [ v[list_map[i]] for i in range(
        # nmodes) ] ))

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

        natoms = len(positions)

        #add the com from each bonus group to the ca_positions list
        for item in self.Bonus_groups:
            #bonus = self.u.select_atoms(item)
            positions = np.vstack((positions, item.center_of_mass()))
            natoms += 1

        matrix = np.zeros((natoms, natoms), "float")
        res_positions, grid, low_x, low_y, low_z = generate_grid(positions, self.cutoff)
        icounter = 0

        cutoffsq = self.cutoff ** 2

        for icounter in range(natoms):
            #find neighbours from the grid
            neighbour_atoms = []
            for x in (-1, 0, 1):
                #print icounter, natoms, len(positions), len(res_positions)
                gx = res_positions[icounter][0] + x
                if 0 <= gx < len(grid):
                    for y in (-1, 0, 1):
                        gy = res_positions[icounter][1] + y
                        if 0 <= gy < len(grid[0]):
                            for z in (-1, 0, 1):
                                gz = res_positions[icounter][2] + z
                                if 0 <= gz < len(grid[0][0]):
                                    neighbour_atoms += grid[gx][gy][gz]

            #for jcounter in range(icounter+1,natoms):
            for jcounter in neighbour_atoms:
                if jcounter > icounter and _dsq(positions[icounter], positions[jcounter]) < cutoffsq:
                    matrix[icounter][jcounter] = -1.0
                    matrix[jcounter][icounter] = -1.0
                    matrix[icounter][icounter] = matrix[icounter][icounter] + 1
                    matrix[jcounter][jcounter] = matrix[jcounter][jcounter] + 1
        return matrix

    def run(self, start=None, stop=None, step=None):
        """Analyze trajectory and produce timeseries.

        Parameters
        ----------
        start : int (optional)
        stop : int (optional)
        step : int (optional)

        Returns
        -------
        results : list
            GNM results per frame::

                results = [(time,eigenvalues[1],eigenvectors[1]),(time,eigenvalues[1],eigenvectors[1])... ]

        .. versionchanged:: 0.16.0
           use start, stop, step instead of skip
        """
        logger.info("GNM analysis: starting")

        self.timeseries = []
        self._timesteps = []

        for ts in self.u.trajectory[start:stop:step]:
            self._timesteps.append(ts.time)

            matrix = self.generate_kirchoff()
            try:
                u, w, v = np.linalg.svd(matrix)
            except np.linalg.LinAlgError:
                print("\nFrame skip at", ts.time,
                      "(SVD failed to converge). Cutoff", self.cutoff)
                continue
            #Save the results somewhere useful in some useful format. Usefully.
            self._generate_output(w, v, self.results, ts.time, matrix, ReportVector=self.ReportVector,
                                  counter=ts.frame)


class closeContactGNMAnalysis(GNMAnalysis):
    """GNMAnalysis only using close contacts.

    This is a version of the GNM where the Kirchoff matrix is
    constructed from the close contacts between individual atoms
    in different residues.

    Parameters
    ----------
    universe : Universe
          Analyze the full trajectory in the universe.
    selection : str (optional)
          MDAnalysis selection string, default "protein"
    cutoff : float (optional)
          Consider selected atoms within the cutoff as neighbors for the
          Gaussian network model [4.5 Å].
    ReportVector : str (optional)
          filename to write eigenvectors to, by default no output is written
          (``None``)
    weights : {"size", None} (optional)
          If set to "size" (the default) then weight the contact by
          :math:`1/\sqrt{N_i N_j}` where :math:`N_i` and :math:`N_j` are the
          number of atoms in the residues :math:`i` and :math:`j` that contain
          the atoms that form a contact.
    MassWeight : bool (deprecated, optional)
          if set to ``True`` equivalent to `weights` set to "size".

    Notes
    -----
    The `MassWeight` option does not perform a true mass weighting but
    weighting by the number of atoms in each residue; the name of the parameter
    exists for historical reasons and will be removed in 0.17.0. Until then,
    setting `MassWeight` to anything but ``None`` will override `weights`.

    See Also
    --------
    :class:`GNMAnalysis`


    .. versionchanged:: 0.16.0
       Made :meth:`generate_output` a private method :meth:`_generate_output`.

    .. deprecated:: 0.16.0
       Instead of ``MassWeight=True`` use ``weights="size"``.
    """

    def __init__(self, universe, selection='protein', cutoff=4.5, ReportVector=None,
                 weights="size", MassWeight=None):
        self.u = universe
        self.selection = selection
        self.cutoff = cutoff
        self.results = []  # final result
        self._timesteps = None  # time for each frame
        self.ReportVector = ReportVector
        self.ca = self.u.select_atoms(self.selection)

        self.weights = weights
        # remove MassWeight in 0.17.0
        if MassWeight is not None:
            warnings.warn("MassWeight=True|False is deprecated in favor of weights='size'|None "
                          "and will be removed in 0.17.0",
                          category=DeprecationWarning)
            self.weights = "size" if MassWeight else None


    def generate_kirchoff(self):
        natoms = self.ca.n_atoms
        nresidues = self.ca.n_residues
        positions = self.ca.positions
        res_positions, grid, low_x, low_y, low_z = generate_grid(positions, self.cutoff)
        residue_index_map = [resnum for [resnum, residue] in enumerate(self.ca.residues) for atom in residue.atoms]
        matrix = np.zeros((nresidues, nresidues), dtype=np.float_)
        cutoffsq = self.cutoff ** 2

        # cache sqrt of residue sizes (slow) so that sr[i]*sr[j] == sqrt(r[i]*r[j])
        sqrt_res_sizes = np.sqrt([r.atoms.n_atoms for r in  self.ca.residues]) \
                         if self.weights == "size" else None

        for icounter in range(natoms):
            neighbour_atoms = []
            for x in (-1, 0, 1):
                gx = res_positions[icounter][0] + x
                if 0 <= gx < len(grid):
                    for y in (-1, 0, 1):
                        gy = res_positions[icounter][1] + y
                        if 0 <= gy < len(grid[0]):
                            for z in (-1, 0, 1):
                                gz = res_positions[icounter][2] + z
                                if 0 <= gz < len(grid[0][0]):
                                    neighbour_atoms += grid[gx][gy][gz]
            for jcounter in neighbour_atoms:
                if jcounter > icounter and _dsq(positions[icounter], positions[jcounter]) <= cutoffsq:
                    iresidue, jresidue = residue_index_map[icounter], residue_index_map[jcounter]
                    contact = 1.0 / (sqrt_res_sizes[iresidue] * sqrt_res_sizes[jresidue]) \
                              if self.weights == "size" else 1.0
                    matrix[iresidue][jresidue] -= contact
                    matrix[jresidue][iresidue] -= contact
                    matrix[iresidue][iresidue] += contact
                    matrix[jresidue][jresidue] += contact
        return matrix
