# Ensemble.py --- Representation of a protein ensemble
# Copyright (C) 2014 Wouter Boomsma, Matteo Tiberti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Ensemble representation --- :mod:`MDAnalysis.analysis.ensemble.ensemble`
=====================================================================

This module contains the Ensemble class allowing for easy reading in
and alignment of the ensemble contained in one or more trajectory files.
Trajectory files can be specified in several formats, including the popular
xtc and dcd, as well as experimental multiple-conformation pdb files, i.e.
those coming from NMR structure resoltion experiments.

.. autoclass:: Ensemble
:Author: Matteo Tiberti, Wouter Boomsma, Tone Bengtsen
:Year: 2015--2016
:Copyright: GNU Public License v3
:Maintainer: Wouter Boomsma <wb@bio.ku.dk>, wouterboomsma on github

.. versionadded:: 0.14.0

"""

import numpy as np

import MDAnalysis
import MDAnalysis.analysis
import MDAnalysis.analysis.align
from MDAnalysis.coordinates.memory import MemoryReader


class Ensemble(MDAnalysis.Universe):
    """
    A wrapper class around Universe providing functionality for aligning
    all frames in a trajectory, and providing easy access to the underlying
    array of coordinates. This class makes use of the MemoryReader
    trajectory reader to store the entire trajectory in a numpy array, in
    which coordinates can be manipulated upon alignment. The frame_interval
    option makes it possible to read in a lower number of frames (e.g. with
    frame-interval=2 only every second frame will be loaded).

    The align method takes an atom selection string, using the MDAnalysis
    syntax for selections
    (see http://mdanalysis.googlecode.com/git/package/doc/html/ \
    documentation_pages/selections.html for details). By default all the
    alpha carbons ("CA") are considered. Frames in an Ensemble object can be
    superimposed to a reference conformation using the reference argument.


    Examples
    --------

    The examples show how to use ENCORE to initiate an Ensemble object.
    The topology- and trajectory files are obtained from the MDAnalysis
    test suite for a simulation of the protein AdK. To run the
    example some imports first need to be executed: ::

        >>> import MDAnalysis.analysis.encore as encore
        >>> from MDAnalysis.tests.datafiles import PDB_small, DCD
        >>> ens = encore.Ensemble(topology=PDB_small,trajectory=DCD)

    In addition, to decrease the computations the :class:`Ensemble` object
    can be initialized by only loading every nth frame from the trajectory
    using the parameter `frame_interval`: ::

        >>> ens = encore.Ensemble(topology=PDB_small, trajectory=DCD,
                                  frame_interval=3)


    """

    def __init__(self,
                 topology=None,
                 trajectory=None,
                 frame_interval=1,
                 **kwargs):
        """
        Constructor for the Ensemble class. See the module description for more
        details.

        Parameters
        ----------

        topology : str
            Topology file name

        trajectory : iterable or str
            One or more Trajectory file name(s)

        frame_interval : int
            Interval at which frames should be included

        """

        # Chained trajectories cannot use TimeSeries functionality
        # and the analysis is therefore slower - we therefore use a
        # single trajectory value when possible
        if len(trajectory) == 1:
            trajectory = trajectory[0]
        MDAnalysis.Universe.__init__(self, topology, trajectory,
                                     **kwargs)

        if kwargs.get('format') != MemoryReader:

            # Try to extract coordinates using Timeseries object
            # This is significantly faster, but only implemented for certain
            # trajectory file formats
            try:
                # frame_interval already takes into account
                coordinates = self.universe.trajectory.timeseries(
                    self.atoms, format='afc', skip=frame_interval)

            # if the Timeseries extraction fails,
            # fall back to a slower approach
            except AttributeError:
                coordinates = np.zeros(
                    tuple([self.universe.trajectory.n_frames/frame_interval]) +
                    self.atoms.coordinates().shape)

                k = 0
                for i, time_step in enumerate(self.universe.trajectory):
                    if (i+1) % frame_interval == 0:
                        coordinates[k] = self.atoms.coordinates(time_step)
                        k += 1
                coordinates = np.swapaxes(coordinates, 0, 1)

            # Overwrite trajectory in universe with an MemoryReader
            # object, to provide fast access and allow coordinates
            # to be manipulated
            self.trajectory = MemoryReader(coordinates)

    def get_coordinates(self, selection="", format='afc'):
        """
        Convenience method for extracting array of coordinates. If no
        selection is provided, this version is slightly faster than accessing
        the coordinates through the timeseries interface (which always takes
        a copy of the array).

        Parameters
        ----------

        selection : str
            Atom selection string in the MDAnalysis format.

        *format*
           the order/shape of the return data array, corresponding
           to (a)tom, (f)rame, (c)oordinates all six combinations
           of 'a', 'f', 'c' are allowed ie "fac" - return array
           where the shape is (frame, number of atoms,
           coordinates)

        """
        # if selection == "":
        #     # If no selection is applied, return raw array
        #     return self.trajectory.get_array(format=format)
        # else:
        return self.trajectory.timeseries(self.select_atoms(selection),
                                          format=format)

    def align(self, selection="name CA", reference=None, weighted=True):
        """
        Least-square superimposition of the Ensemble coordinates to a reference
         structure.

        Parameters
        ----------

        selection : str
            Atom selection string in the MDAnalysis format. Default is
            "name CA"

        reference : None or MDAnalysis.Universe
            Reference structure on which those belonging to the Ensemble will
            be fitted upon.  It must have the same topology as the Ensemble
            topology. If reference is None, the structure in the first frame of
            the ensemble will be used as reference.

        weighted : bool
            Whether to perform weighted superimposition or not

        """

        coordinates = self.trajectory.timeseries(format='fac')

        alignment_subset_selection = self.select_atoms(selection)
        alignment_subset_coordinates = \
            self.trajectory.timeseries(alignment_subset_selection,
                                       format='fac')

        if weighted:
            alignment_subset_masses = alignment_subset_selection.masses
        else:
            alignment_subset_masses = np.ones(
                alignment_subset_selection.masses.shape[0])

        # Find center of mass of alignment subset for all frames
        alignment_subset_coordinates_center_of_mass = np.average(
            alignment_subset_coordinates,
            axis=1,
            weights=alignment_subset_masses)

        # Move both subset atoms and the other atoms to the center of mass of
        # subset atoms
        coordinates -= alignment_subset_coordinates_center_of_mass[:,
                                                                   np.newaxis]

        # if reference: no offset
        if reference:
            offset = 0
            # Select the same atoms in reference structure
            reference_atom_selection = reference.select_atoms(selection)
            reference_coordinates = reference_atom_selection.atoms.coordinates()
        else:
            reference_atom_selection = self.select_atoms(selection)
            reference_coordinates = alignment_subset_coordinates[0]

            # Skip the first frame, which is used as reference
            offset = 1

        if weighted:
            reference_masses = reference_atom_selection.masses
        else:
            reference_masses = np.ones(
                reference_atom_selection.masses.shape[0])

        # Reference center of mass
        reference_center_of_mass = np.average(reference_coordinates, axis=0,
                                              weights=reference_masses)
        # Move reference structure to its center of mass
        reference_coordinates -= reference_center_of_mass

        import logging
        logging.info("reference_coordinates: " + str(reference_coordinates))

        # Apply optimal rotations for each frame
        for i in range(offset, len(coordinates)):
            # Find rotation matrix on alignment subset
            rotation_matrix = MDAnalysis.analysis.align.rotation_matrix(
                alignment_subset_coordinates[i],
                reference_coordinates,
                alignment_subset_masses)[0]

            # Apply rotation matrix
            coordinates[i][:] = np.transpose(np.dot(rotation_matrix,
                                                    np.transpose(
                                                        coordinates[i][:])))
