# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

#Analyse a trajectory using elastic network models, following the approach of Hall et al (JACS 2007)
#Ben Hall (benjamin.a.hall@ucl.ac.uk) is to blame
#Copyright 2011; Consider under GPL v2 or later

"""
Elastic network analysis of MD trajectories --- :mod:`MDAnalysis.analysis.gnm`
==============================================================================

:Author: Benjamin Hall <benjamin.a.hall@ucl.ac.uk>
:Year: 2011
:Copyright: GNU Public License v2 or later


Analyse a trajectory using elastic network models, following the approach of [Hall2007]_.

An example is provided in :file:`examples/GNMExample.py`. The basic
approach is to pass a trajectory to :class:`GNMAnalysis` and then run
the analysis::

   u = MDAnalysis.Universe(PSF,DCD)
   C = MDAnalysis.analysis.gnm.GNMAnalysis(u,ReportVector="output.txt")

   C.run()
   output = zip(*C.results)

   outputfile = open("eigenvalues.dat","w")
   for item in output[1]:
      print >> outputfile, item
   outputfile.close()

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

.. autofunction:: backup_file
.. autofunction:: generate_grid
.. autofunction:: order_list

"""

# import copy #unused

import numpy as np
from numpy import linalg

import os

#import warnings #unused
import logging

logger = logging.getLogger('MDAnalysis.analysis.GNM')


def backup_file(filename):
    '''
    This function helps prevent overwriting default named files
    '''
    if os.path.exists(filename):
        target_name = "#" + filename
        failure = True
        if not os.path.exists(target_name):
            os.rename(filename, target_name)
            failure = False
        else:
            for i in range(20):
                alt_target_name = target_name + "." + str(i)
                if os.path.exists(alt_target_name):
                    continue
                else:
                    os.rename(filename, alt_target_name)
                    failure = False
                    break
        if failure:
            print "Too many backups. Clean up and try again"
            exit()


def generate_grid(positions, cutoff):
    '''
    An alternative to searching the entire list of each atom; divide the structure into cutoff sized boxes
    This way, for each particle you only need to search the neighbouring boxes to find the particles within the cutoff
    Observed a 6x speed up for a smallish protein with ~300 residues; this should get better with bigger systems.
    '''
    [x, y, z] = zip(*positions)
    high_x = max(x)
    high_y = max(y)
    high_z = max(z)
    low_x = min(x)
    low_y = min(y)
    low_z = min(z)
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
    '''
    Returns a dictionary showing the order of eigenvalues (which are reported scrambled normally)
    '''
    ordered = list(w)
    unordered = list(w)
    ordered.sort()
    list_map = {}
    for i in range(len(w)):
        list_map[i] = unordered.index(ordered[i])
    return list_map


class GNMAnalysis(object):
    '''Basic tool for GNM analysis.

    Each frame is treated as a novel structure and the GNM
    calculated.  By default, this stores the dominant eigenvector
    and its associated eigenvalue; either can be used to monitor
    conformational change in a simulation.
    '''

    def __init__(self, universe, selection='protein and name CA', cutoff=7.0, ReportVector=None, Bonus_groups=()):
        self.u = universe
        self.selection = selection
        self.cutoff = cutoff
        self.results = []  # final result
        self._timesteps = None  # time for each frame
        self.ReportVector = ReportVector
        # this is a tuple of selection groups, the com of each which will be added as a single point in the ENM
        # (its a popular way of treating small ligands eg drugs)
        # be careful that it doesn't have any atoms which can be caught by the global selection as this could lead to
        #  double counting
        self.Bonus_groups = [self.u.select_atoms(item) for item in Bonus_groups]
        self.ca = self.u.select_atoms(self.selection)

    def generate_output(self, w, v, outputobject, time, matrix, nmodes=2, ReportVector=None, counter=0):
        '''Appends eigenvalues and eigenvectors to results.

        This generates the output by adding eigenvalue and
        eigenvector data to an appendable object and optionally
        printing some of the results to file. This is the function
        to replace if you want to generate a more complex set of
        outputs
        '''
        list_map = order_list(w)
        #print round(time), w[list_map[1]]
        if ReportVector:
            with open(ReportVector, "a") as oup:
                for item in enumerate(v[list_map[1]]):
                    print >> oup, "", counter, time, item[0] + 1, w[list_map[1]], item[1]
        outputobject.append((time, w[list_map[1]], v[list_map[1]]))
        #outputobject.append((time, [ w[list_map[i]] for i in range(nmodes) ], [ v[list_map[i]] for i in range(
        # nmodes) ] ))

    def generate_kirchoff(self):
        '''Generate the Kirchhoff matrix of contacts.

        This generates the neighbour matrix by generating a grid of
        near-neighbours and then calculating which are are within
        the cutoff. Returns the resulting matrix
        '''
        #ca = self.u.select_atoms(self.selection)
        positions = self.ca.coordinates()

        natoms = len(positions)

        #add the com from each bonus group to the ca_positions list
        for item in self.Bonus_groups:
            #bonus = self.u.select_atoms(item)
            positions = np.vstack((positions, item.center_of_mass()))
            natoms += 1

        matrix = np.zeros((natoms, natoms), "float")
        [res_positions, grid, low_x, low_y, low_z] = generate_grid(positions, self.cutoff)
        icounter = 0
        for icounter in range(natoms):
            #find neighbours from the grid
            neighbour_atoms = []
            for x in (-1, 0, 1):
                #print icounter, natoms, len(positions), len(res_positions)
                if (res_positions[icounter][0] + x) >= 0 and (res_positions[icounter][0] + x) < len(grid):
                    for y in (-1, 0, 1):
                        if (res_positions[icounter][1] + y) >= 0 and (res_positions[icounter][1] + y) < len(grid[0]):
                            for z in (-1, 0, 1):
                                if (res_positions[icounter][2] + z) >= 0 and (res_positions[icounter][2] + z) < len(
                                        grid[0][0]):
                                    neighbour_atoms += grid[
                                        res_positions[icounter][0] + x][res_positions[icounter][1]
                                                                        + y][res_positions[icounter][2] + z]
            #for jcounter in range(icounter+1,natoms):
            for jcounter in neighbour_atoms:
                if jcounter > icounter and ((positions[
                    icounter][0] - positions[jcounter][0]) ** 2 +
                        (positions[icounter][1] - positions[jcounter][1]) ** 2 +
                        (positions[icounter][2] - positions[jcounter][2]) ** 2) <= self.cutoff ** 2:
                    matrix[icounter][jcounter] = -1.0
                    matrix[jcounter][icounter] = -1.0
                    matrix[icounter][icounter] = matrix[icounter][icounter] + 1
                    matrix[jcounter][jcounter] = matrix[jcounter][jcounter] + 1
        return (matrix)

    def run(self, skip=1):
        '''Analyze trajectory and produce timeseries.

        Returns GNM results per frame::

          results = [(time,eigenvalues[1],eigenvectors[1]),(time,eigenvalues[1],eigenvectors[1])... ]

        '''
        logger.info("GNM analysis: starting")
        counter = 0

        self.timeseries = []
        self._timesteps = []

        try:
            self.u.trajectory.time

            def _get_timestep():
                return self.u.trajectory.time

            logger.debug("GNM analysis is recording time step")
        except NotImplementedError:
            # chained reader or xyz(?) cannot do time yet
            def _get_timestep():
                return self.u.trajectory.frame

            logger.warn("GNM analysis is recording frame number instead of time step")

        for ts in self.u.trajectory:
            if counter % skip != 0:
                counter += 1
                continue
            counter += 1
            frame = ts.frame
            timestep = _get_timestep()
            self._timesteps.append(timestep)

            matrix = self.generate_kirchoff()
            try:
                [u, w, v] = linalg.svd(matrix)
            except:
                print "\nFrame skip at", timestep, "(SVD failed to converge). Cutoff", self.cutoff
                continue
            #Save the results somewhere useful in some useful format. Usefully.
            self.generate_output(w, v, self.results, timestep, matrix, ReportVector=self.ReportVector, counter=counter)


class closeContactGNMAnalysis(GNMAnalysis):
    """GNMAnalysis only using close contacts.

    This is a version of the GNM where the Kirchoff matrix is
    constructed from the close contacts between individual atoms
    in different residues
    """

    def __init__(self, universe, selection='protein', cutoff=4.5, ReportVector=None, MassWeight=True):
        self.u = universe
        self.selection = selection
        self.cutoff = cutoff
        self.results = []  # final result
        self._timesteps = None  # time for each frame
        self.ReportVector = ReportVector
        # no bonus groups in this version of the GNM analysis tool; this is because this version doesn't use CA atoms
        #  or centroids, and so doesn't need them
        self.ca = self.u.select_atoms(self.selection)
        self.MassWeight = MassWeight

    def generate_kirchoff(self):
        natoms = len(self.ca.atoms)
        nresidues = len(self.ca.residues)
        positions = self.ca.coordinates()
        [res_positions, grid, low_x, low_y, low_z] = generate_grid(positions, self.cutoff)
        residue_index_map = [resnum for [resnum, residue] in enumerate(self.ca.residues) for atom in residue]
        matrix = np.zeros((nresidues, nresidues), "float")
        for icounter in range(natoms):
            neighbour_atoms = []
            for x in (-1, 0, 1):
                if (res_positions[icounter][0] + x) >= 0 and (res_positions[icounter][0] + x) < len(grid):
                    for y in (-1, 0, 1):
                        if (res_positions[icounter][1] + y) >= 0 and (res_positions[icounter][1] + y) < len(grid[0]):
                            for z in (-1, 0, 1):
                                if (res_positions[icounter][2] + z) >= 0 and (res_positions[icounter][2] + z) < len(
                                        grid[0][0]):
                                    neighbour_atoms += grid[res_positions[icounter][0] + x][
                                        res_positions[icounter][1] + y][res_positions[icounter][2] + z]
            for jcounter in neighbour_atoms:
                if jcounter > icounter and ((positions[icounter][
                    0] - positions[jcounter][0]) ** 2 +
                        (positions[icounter][1] - positions[jcounter][1]) ** 2 +
                        (positions[icounter][2] - positions[jcounter][2]) ** 2) <= self.cutoff ** 2:
                    [iresidue, jresidue] = [residue_index_map[icounter], residue_index_map[jcounter]]
                    if self.MassWeight:
                        contact = 1.0 / ((len(self.ca.residues[iresidue])) * (len(self.ca.residues[jresidue]))) ** 0.5
                    else:
                        contact = 1.0
                    matrix[iresidue][jresidue] = matrix[iresidue][jresidue] - contact
                    matrix[jresidue][iresidue] = matrix[jresidue][iresidue] - contact
                    matrix[iresidue][iresidue] = matrix[iresidue][iresidue] + contact
                    matrix[jresidue][jresidue] = matrix[jresidue][jresidue] + contact
        return (matrix)
