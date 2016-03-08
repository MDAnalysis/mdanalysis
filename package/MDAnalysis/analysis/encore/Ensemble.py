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
:Maintainer: Matteo Tiberti <matteo.tiberti@gmail.com>, mtiberti on github

.. versionadded:: 0.14.0

"""

import MDAnalysis
import MDAnalysis.analysis
import MDAnalysis.analysis.align
import numpy
import logging
import numpy as np
import errno
from MDAnalysis.coordinates.array import ArrayReader


class Ensemble(MDAnalysis.Universe):
    """
    A wrapper class around Universe providing
    Ensemble class designed to easily manage more than one trajectory files.
    Users can provide either a topology/trajectory(es) combination or a
    MDAnalysis.Universe object. Topology and trajectory files must have the
    same number of atoms, and order is of course important.

    While creating a new Ensemble object it is possible to load from a
    trajectory a selected subset of atoms, using the MDAnalysis syntax for
    selections
    (see http://mdanalysis.googlecode.com/git/package/doc/html/ \
    documentation_pages/selections.html for details)
    and the atom_selection_string argument. By default all the alpha carbons
    ("CA") are considered. It is also possible to load a lower number of frames
     for each trajectory, by selecting only one frame every frame_interval
     (e.g. with frame-interval=2 only every second frame will be loaded).

    Frames in an Ensemble object can be superimposed to a reference
    conformation (see method align). By default the rotation matrix for this
    superimposition is calculated on all the atoms of the system, as defined
    by the atom_selection_string. However, if the
    superimposition_selection_string is provided, that subset will be used to
    calculate the rotation matrix, which will be applied on the whole
    atom_selection_string. Notice that the set defined by
     superimposition_selection_string is completely independent from the
     atom_selection_string atoms, as it can be a subset or superset of that,
     although it must refer to the same topology.

    Attributes
    ----------

    topology_filename : str
        Topology file name.

    trajectory_filename : str
        Trajectory file name. If more then one are specified, it is a list of
        comma-separated names (e.g. "traj1.xtc,traj2.xtc")

    frame_interval : int
        Keep only one frame every frame_interval (see the package or module
        description)

    selection : str
        Atom selection string in the MDAnalysis format
         (see http://mdanalysis.googlecode.com/git/package/doc/html/documentation_pages/selections.html)

    atom_selection : MDAnalysis.core.AtomGroup
        MDAnalysis atom selection, which corresponds to the selection
        defined by atom_selection_string on universe

    coordinates : (x,N,3) numpy.array
        Array of coordinate which will be used in the calculations, where x is
        the number of frames and N is the number of atoms. Notice that these coordinates may be different from those of universe, because of the atom_selection and frame_interval.

    superimposition_selection_string : str
	Analogous to atom_selection_string, but related to the subset of atoms that
	will be used for 3D superimposition.

    superimposition_selection : MDAnalysis.core.AtomGroup
	Analogous to atom_selection, but related to the subset of atoms that will
	 be used for 3D superimposition.

    superimposition_coordinates : (x,N,3) numpy.array
        Analogous to coordinates, but related to the subset of atoms that will
         be used for 3D superimposition.


    Examples
    --------

	The examples show how to use ENCORE to initiate an Ensemble object.
	The topology- and trajectory files are obtained from the MDAnalysis
	test suite for a simulation of the protein AdK. To run the
	example some imports first need to be executed: ::

	    >>> from MDAnalysis import *
	    >>> from MDAnalysis.analysis.encore.similarity import *
	    >>> from MDAnalysis.tests.datafiles import PDB_small, DCD
	    >>> ens = Ensemble(topology=PDB_small,trajectory=DCD)

	In addition, to decrease the computations the :class:`Ensemble` object
	can be initialized by only loading every nth frame from the trajectory
	using the parameter `frame_interval`: ::

	    >>> ens = Ensemble(topology=PDB_small, trajectory=DCD, frame_interval=3)


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

        selection : str
            Atom selection string in the MDAnalysis format
            (see http://mdanalysis.googlecode.com/git/package/doc/html/documentation_pages/selections.html)

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


        if kwargs.get('format', None) != ArrayReader:

            # Try to extract coordinates using Timeseries object
            # This is significantly faster, but only implemented for certain
            # trajectory file formats
            try:
                # frame_interval already takes into account
                coordinates = self.universe.trajectory.timeseries(
                    self.atoms, format='afc', skip=frame_interval)

            # if the Timeseries extraction fails, fall back to a slower approach
            except AttributeError:
                coordinates = numpy.zeros(
                    tuple([self.universe.trajectory.n_frames]) +
                    self.atoms.coordinates().shape)

                k = 0
                for i, time_step in enumerate(self.universe.trajectory):
                    if i%frame_interval == 0:
                        coordinates[k] = self.atoms.coordinates(time_step)
                        k+=1
                coordinates = np.swapaxes(coordinates,0,1)

            # Overwrite trajectory in universe with an ArrayReader
            # object, to provide fast access and allow coordinates
            # to be manipulated
            self.trajectory = ArrayReader(coordinates)
                # self._get_coordinates(frame_interval=frame_interval),
                # format='afc')

        # # Overwrite atoms selection from Universe
        # self.atoms_selection = self.select_atoms(self.atom_selection_string)
        #
        # # Set the attributes for the atom set on which fitting will be
        # # performed. Fitting and calculation may be performed on two
        # # non-overlapping sets. This is optional.
        # if superimposition_selection_string:
        #     self.superimposition_selection_string \
        #         = superimposition_selection_string
        #     self.superimposition_selection = self.select_atoms(
        #         superimposition_selection_string)
        #     self.superimposition_coordinates = self.get_coordinates(
        #         subset_selection_string=self.superimposition_selection_string)
        # else:
        #     self.superimposition_selection_string = self.atom_selection_string
        #     self.superimposition_selection = self.atoms_selection
        #     self.superimposition_coordinates = numpy.copy(self.trajectory.get_array())

        # # Save trajectories filename for future reference
        # if type(trajectory) == str:
        #     self.trajectory_filename = trajectory
        # else:
        #     self.trajectory_filename = ", ".join(trajectory)
        #
        #     # Save topology filename for future reference
        self.topology_filename = topology

    def get_coordinates(self, selection, format):
        if selection == "":
            # If no selection is applied, return raw array
            return self.trajectory.get_array(format=format)
        else:
            return self.trajectory.timeseries(self.select_atoms(selection),
                                              format=format)

    # def _get_coordinates(self, selection="", frame_interval=1):
    #     """
    #     Get a set of coordinates from Universe.
    #
    #     Parameters
    #     ----------
    #
    #     subset_selection_string : None or str
    #         Selection string that selects the universe atoms whose coordinates
    #         have to be returned. The frame_interval will be automatically
    #         applied. If the argument is None,  the atoms defined in the
    #         atom_selection_string will be considered.
    #
    #     Returns
    #     -------
    #
    #     coordinates : (x,N,3) numpy array
    #         The requested array of coordinates.
    #
    #     """
    #
    #     if selection == "":
    #         atomgroup = self.atoms
    #     else:
    #         atomgroup = self.select_atoms(selection)
    #
    #     # if not subset_selection_string:
    #     #     subset_selection_string = self.atom_selection_string
    #     # subset_selection = self.universe.select_atoms(subset_selection_string)
    #
    #     if len(atomgroup) == 0:
    #         logging.error(
    #             "ERROR: selection \'%s\' not found in topology."
    #             % subset_selection_string)
    #         exit(1)
    #
    #     # Try to extract coordinates using Timeseries object
    #     # This is significantly faster, but only implemented for certain
    #     # trajectory file formats
    #     try:
    #         # frame_interval already takes into account
    #         coordinates = self.universe.trajectory.timeseries(
    #             atomgroup, format='afc', skip=frame_interval)
    #
    #     # if the Timeseries extraction fails, fall back to a slower approach
    #     except:
    #         coordinates = numpy.zeros(
    #             tuple([self.universe.trajectory.n_frames]) +
    #             atomgroup.coordinates().shape)
    #
    #         k = 0
    #         for i, time_step in enumerate(self.universe.trajectory):
    #             if i%frame_interval == 0:
    #                 coordinates[k] = atomgroup.coordinates(time_step)
    #                 k+=1
    #         coordinates = np.swapaxes(coordinates,0,1)
    #     return coordinates

    def align(self, selection="name *", reference=None, weighted=True):
        """
        Least-square superimposition of the Ensemble coordinates to a reference
         structure.

        Parameters
        ----------

        reference : None or MDAnalysis.Universe
            Reference structure on which those belonging to the Ensemble will
            be fitted upon.  It must have the same topology as the Ensemble
            topology. If reference is None, the structure in the first frame of
            the ensemble will be used as reference.

        weighted : bool
            Whether to perform weighted superimposition or not

        """

        coordinates = self.trajectory.get_array(format='fac')

        alignment_subset_selection = self.select_atoms(selection)
        alignment_subset_coordinates = \
            self.trajectory.timeseries(alignment_subset_selection,
                                       format='fac')

        # alignment_subset_atom_selection = self.superimposition_selection
        # alignment_subset_coordinates = self.superimposition_coordinates

        if weighted:
            alignment_subset_masses = alignment_subset_selection.masses
        else:
            alignment_subset_masses = np.ones(
                alignment_subset_selection.masses.shape[0])

        # Find center of mass of alignment subset for all frames
        alignment_subset_coordinates_center_of_mass = numpy.average(
            alignment_subset_coordinates,
            axis=1,
            weights=alignment_subset_masses)

        # Move both subset atoms and the other atoms to the center of mass of
        # subset atoms
        # alignment_subset_coordinates -= \
        #     alignment_subset_coordinates_center_of_mass[ :, numpy.newaxis]
        # print alignment_subset_coordinates[0]
        coordinates -= alignment_subset_coordinates_center_of_mass[:,
                       numpy.newaxis]
        # print coordinates.shape

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
        reference_center_of_mass = numpy.average(reference_coordinates, axis=0,
                                                 weights=reference_masses)
        # Move reference structure to its center of mass
        reference_coordinates -= reference_center_of_mass

        # Apply optimal rotations for each frame
        for i in range(offset, len(coordinates)):
            # Find rotation matrix on alignment subset
            rotation_matrix = MDAnalysis.analysis.align.rotation_matrix(
                alignment_subset_coordinates[i],
                reference_coordinates,
                alignment_subset_masses)[0]

            # Apply rotation matrix
            coordinates[i][:] = numpy.transpose(numpy.dot(rotation_matrix,
                                                          numpy.transpose(
                                                          coordinates[i][:])))
        # self.trajectory.set_array(self.coordinates)
        # for k, ts in enumerate(self.trajectory[:1]):
        #     print k, self.atoms.positions, id(self.trajectory.ts.positions)
        #     # self.trajectory[i].positions = self.coordinates[i]
        #     # self.atoms.set_positions(self.coordinates[i][:])
        #     self.trajectory.ts.positions[:] = self.coordinates[i][:]
        #     print k, "***", self.atoms.positions, id(self.trajectory.ts.positions)
        # for k, ts in enumerate(self.trajectory[:1]):
        #     print k, self.atoms.positions, id(self.trajectory.ts.positions)
        # print "&&&&&&&&&&&&&", self.trajectory.ts[0]
