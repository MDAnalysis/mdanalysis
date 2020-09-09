# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
r""" Average Aggregation number ---
======================================================

This module is able to calculate the number average :math: `g_{n}`, weight
average :math: `g_{w}`, and z average :math: `g_{z}` aggregation number
based on the three different categories including the distance between
closest atom, the distance between the centres of mass, and the distance
between the two certain atoms of two molecules. Also, average gyration radius
of aggregates can be calculated (optional). The monomers are eliminated as a
default since they have not been involve in aggregation, but the user have an
option to consider them by inserting (no_monomer = False) in the
:func:`AggregationSize`. The formula of each type of average is:
..math::
    g_{n} = \frac {{sum_{i=2}^number of molecules} \langle
    number of cluster*number of molecules in clusters
    \rangle}{{sum_{i=2}^number of molecules} \langle
    number of molecules in clusters \rangle}
    g_{w} = \frac {{sum_{i=2}^number of molecules} \langle
    number of cluster*number of molecules in clusters^2
    \rangle}{{sum_{i=2}^number of molecules} \langle
    number of cluster*number of molecules in clusters
    \rangle}
    g_{z} = \frac {{sum_{i=2}^number of molecules} \langle
    number of cluster*number of molecules in clusters^3
    \rangle}{{sum_{i=2}^number of molecules} \langle
    number of cluster*number of molecules in clusters^2
    \rangle}
---------------------
Algorithm
---------------------
1- Import the needed module
2- Calculate the distance between molecules based on the intended
criteria
3- Find the aggregates by checking only one time each using the
"aggregate" subfunction (to be efficient)
4- Check and merge two aggregates as one aggregate if the distance
between any element of them is less than cut_off
5- Eliminate monomers (optional)
6- measure average gyration radius (optional)
7- measure the number of aggregates with different size
---------------------
P.S.
---------------------
1- The box should be unwrapped before running the code. We
recommend doing -pbc whole and -pbc cluster using trjconv keyword available
in GROMCAS or any similar procedure.
2- If you have more than one type of molecule and you want to measure the
aggregation without considering the type of molecule the code will works
perfectly for "Closest" and "COM" Types while for the "Specific" the results
will not be correct. The reasone is that you cannot tell atom_1 and atom_2 is
related to which residue!
3- If you have more than one type of molecule and you want to measure
aggregation for each type of molecule separately, you can keep only the
target molecule in the box and run the code multiple times for different
molecule type.
4- You will get an error if there is no aggregation in the box or if the
aggregation distance choose too short.
"""
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib import distances
from MDAnalysis.analysis.base import AnalysisBase


class AggregationSize(AnalysisBase):
    """Average aggregates size function
        ag_size = AggregationSize(atomgroup, number_of_mol,
        agregation_criteria, cut_off, no_monomer = True,
        ave_gyr_calc = True, atom_1 = "C1", atom_2 = "C1")

    Arguments
    ---------
    atomgroup : selected residue to calculate ave aggregate size
    number_of_mol :  molecule number of selected residue
    agregation_criteria :  the criterial for considering molecules
                        as an aggregate. It can be "Closest", "COM",
                        and "Specific". default "Closest"
    cut_off :  it is the distance bellow which two molecule will consider
            as one aggregate. For "Closest", it is 3.5 Angstrom while for
            other condition it is 8.5 Angstrom
    no_monomer : consider monomer in averaging or not (optional)[True]
    ave_gyr_calc : calculate the average gyration radius of aggregate
                 in each frame or not (optional) [False]
    atom_1 & atom_2 : they are only nessesary to enter for "Specific"
                    agregation_criteria.[None]

    Example
    -------
    Creat :class:`AggregationSize` objects, by supplying an atomgroup
    and set the number of mlecules which is mandatory and other optional
    parameters. Then use the :meth:`run` method ::

        ag_size = AggregationSize(atomgroup, number_of_mol = 30, cut_off = 10
        aggregation_criteria = "COM", no_monomer = False, ave_gyr_calc = True
        ).run()

    Results are available through the :attr:`results` in five coulmns
    including frame number, number average, weight average, z average, and
    average gyration radius::

        ag_size.results
    """
    def __init__(self, atomgroup, number_of_mol,
                 aggregation_criteria, no_monomer=True,
                 ave_gyr_calc=False, cut_off=False,
                 atom_1=None, atom_2=None, **kwargs):
        """Set up the initial analysis parameters.
        """
        trajectory = atomgroup.universe.trajectory
        super(AggregationSize, self).__init__(trajectory)
        self.atomgroup = atomgroup
        self._molnum = number_of_mol
        self._aggregation_criteria = aggregation_criteria
        if (self._aggregation_criteria != "Closest"
            and self._aggregation_criteria != "COM"
                and self._aggregation_criteria != "Specific"):
            raise ValueError("Enter valid aggregation_criteria")
        if (self._aggregation_criteria == "Specific"
                and (not(atom_1 and atom_2))):
            raise ValueError(
                "Enter valid atom_1 & atom_2: atom name and number(e.g. C1)")
        self._cut_off = cut_off
        self._monomer = no_monomer
        self.atom_1 = atom_1
        self.atom_2 = atom_2
        self._ave_gyr_calc = ave_gyr_calc
        self._ave_GYR = []

    def _prepare(self):
        """Create array of zeroes as a placeholder for parameters including results
        _distance_of_mol, size_num, Connected_NLevel.
        """
        self.results = np.zeros((self.n_frames, 5))
        self._distance_of_mol = [0] * self._molnum
        i = 0
        for it in reversed(range(1, self._molnum+1)):
            self._distance_of_mol[i] = np.zeros(it * 1).reshape(it, 1)
            i += 1
        self.size_num = [0] * self.n_frames
        self.Connected_NLevel = []

    def closest_atom(self):
        """
        Calculate the distance between each two molecules based on the closest
        atom of molecules. The default cut_off distance is 3.5 Angstrom but can
        be changed by the user in the :func:`AggregationSize`.
        *cut_off unit is Angstrom

        Example
        -------
            ag_size = AggregationSize(atomgroup, number_of_mol = 30,
            cut_off = 10).run()

        Returns
        -------
            self._distance_of_mol : list
        """
        for i in range(1, self._molnum+1):
            ii = 0
            for j in range(i+1, self._molnum+1):
                residue1 = self.atomgroup.select_atoms(f"resid {i}").positions
                residue2 = self.atomgroup.select_atoms(f"resid {j}").positions
                distance_of_every_atoms = distances.distance_array(
                    residue1, residue2, box=self.atomgroup.universe.dimensions)
                self._distance_of_mol[i][ii] = np.amin(distance_of_every_atoms)
                ii += 1
        return self._distance_of_mol

    def center_of_mass(self):
        """Calculate the distance between the center of mass of two molecules.
        The default cut_off distance is 8.5 Angstrom but can be changed by the
        user in the :func:`AggregationSize`.
        *cut_off unit is Angstrom

        Example
        -------
            ag_size = AggregationSize(atomgroup, number_of_mol = 30,
            cut_off = 10).run()

        Returns
        -------
            self._distance_of_mol : list
        """
        cm = {}
        for x in range((self._molnum)):
            cm[x+1] = (
                self.atomgroup.select_atoms(f"resid {x+1}").center_of_mass())
        cm_coordintes = np.array(list(cm.values()))
        distance_of_coms = distances.distance_array(
                cm_coordintes, cm_coordintes,
                box=self.atomgroup.universe.dimensions)
        for j in range(self._molnum):
            ii = 0
            for i in range(j+1, self._molnum):
                self._distance_of_mol[j+1][ii] = distance_of_coms[i][j]
                ii += 1
        return self._distance_of_mol

    def specific_atoms(self):
        """Calculate the distance between the selected atom of two molecules.
        Besides the cut_off distance, the two atoms' name (atom_1 and atom_2)
        must be defined. The default value for cut_off distance is 8.5 Angstrom
        while there are no default values for atom name.
        *cut_off unit is Angstrom

        Example
        -------
            ag_size = AggregationSize(atomgroup, number_of_mol = 30,
            aggregation_criteria = "Specific", no_monomer = True, atom_1 = "C1"
            , atom_2 = "C1").run()

        Returns
        -------
            self._distance_of_mol : list
        """
        self._cut_off = self._cut_off or 8.5
        residue1 = self.atomgroup.select_atoms(f"name {self.atom_1}").positions
        residue2 = self.atomgroup.select_atoms(f"name {self.atom_2}").positions
        distance_of_atoms = distances.distance_array(
            residue1, residue2, box=self.atomgroup.universe.dimensions)
        for j in range(self._molnum):
            ii = 0
            for i in range(j+1, self._molnum):
                self._distance_of_mol[j+1][ii] = distance_of_atoms[i][j]
                ii += 1
        return self._distance_of_mol

    def merge_aggregates(self, aggregate):
        """Check if two aggregates has any element with distance less than
        cut-off. If they are, they will be merged as a one aggregate.

        Returns
        -------
            aggregate : list
        """
        restart = True
        while restart:
            restart = False
            for i in range(0, len(aggregate)-1):
                aggregate1 = set(aggregate[i])
                for j in range(i+1, len(aggregate)-1):
                    aggregate2 = set(aggregate[j])
                    if len(aggregate1 & aggregate2) > 0:
                        aggregate[i] = aggregate[i] + aggregate[j]
                        del aggregate[j]
                        aggregate[i] = (list(dict.fromkeys(aggregate[i])))
                        check = True
                        break
                    else:
                        check = False
                if check is True:
                    restart = True
                    break
        return aggregate

    def connecting_loop(self, connected, counter, p_level, aggregate):
        """This function with the help of :func:`aggregation` will find the
        connected molecule to the target molecules
        Returns
        -------
            connected : list
            aggregate : list
        """
        while len(connected) > 0:
            adjust = 0
            n_level = np.array([])
            for k in connected:
                k = int(k)
                # "if" will bypass function for last molecule because
                # there is no molecule after that.
                if k != self._molnum:
                    self.aggregation(k, self._frame_index, counter, adjust,
                                     self._cut_off, connected)
                    aggregate[counter-1].extend(self.Connected_NLevel)
                else:
                    aggregate[counter-1].extend([k])
                # n_level find molecule attach to main aggregate for
                # levels after level1
                n_level = np.append(n_level, self.Connected_NLevel)
                adjust += 1
            p_level = np.append(p_level, connected)
            # Delete duplicated number in n_Level
            n_level_no_repetition = list(dict.fromkeys(list(n_level)))
            # Check if the molecule for the next level
            # has been invstigated in pervious level
            for i in set(n_level_no_repetition):
                if i in (list(p_level)):
                    n_level_no_repetition.remove(i)
            connected = np.array(n_level_no_repetition)
            return(connected, aggregate)

    def average_gyration_radius_calculation(self, aggregate):
        """Calculate the average gyration radius of aggregates in each frame
        Returns
        -------
            self._ave_GYR
        """
        gyr_of_each_cluster = []
        for x in range(len(aggregate)-1):
            aggregate[x] = [int(elem) for elem in aggregate[x]]
            if len(aggregate[x]) != 1:
                molecules_in_cluster = ([f'resid {i}' for i in aggregate[x]])
                atomgroup_in_cluster = self.atomgroup.select_atoms(
                    molecules_in_cluster[0])
                for x in molecules_in_cluster[1:]:
                    atomgroup_in_cluster = (
                        atomgroup_in_cluster + self.atomgroup.select_atoms(x))
                gyr_of_each_cluster.append(
                    atomgroup_in_cluster.radius_of_gyration())
            elif self._monomer is False:
                molecules_in_cluster = ([f'resid {i}' for i in aggregate[x]])
                atomgroup_in_cluster = self.atomgroup.select_atoms(
                    molecules_in_cluster[0])
                gyr_of_each_cluster.append(
                    atomgroup_in_cluster.radius_of_gyration())
        gyr_of_each_cluster = np.array(gyr_of_each_cluster)
        return self._ave_GYR.append(np.average(gyr_of_each_cluster))

    def aggregation(self, k, tss, counter, adjust, cut_off, connected):
        """
        This is a subfunction that enables the code to find all connected
        molecules in a cluster.
        k is a connected molecule to the main molecule in the first level.
        self._frame_index is a time frame.
        counter is a step.
        adjust provides a suitable off_set to match the number of molecules.
        Limit is the cut_off provided by the user.
        """
        boolean_connected_mol = self._distance_of_mol[k] < cut_off
        # This "if" will check if there is any uninvestigated molecule left
        if len(boolean_connected_mol) != 0:
            index_connected_mol = np.where(boolean_connected_mol)
            self.Connected_NLevel = np.array(
                [j for j in index_connected_mol[0]] + connected[adjust] + 1)
        else:
            self.Connected_NLevel = []

    def _single_frame(self):
        aggregate = [0] * self._molnum
        if self._aggregation_criteria == "Closest":
            self._cut_off = self._cut_off or 3.5
            self.closest_atom()
        elif self._aggregation_criteria == "COM":
            self._cut_off = self._cut_off or 8.5
            self.center_of_mass()
        elif self._aggregation_criteria == "Specific":
            self._cut_off = self._cut_off or 8.5
            self.specific_atoms()
        # p_level collact all the investigated nodes in past levels
        p_level = []
        # len(Cluster) will be the number of clusters
        clusters = []
        # this for will go through all molecules to
        # check their connection with others
        for counter in range(1, self._molnum):
            # aggergate is a list of list which includes the molecules in
            # each cluster
            aggregate[counter-1] = [counter]
            # this "if" will restrict repeating the calculation of the same
            # molecule to increase efficincy of code
            if counter in p_level:
                continue
            # recognize the molecules that are attached to a molecule in first
            # level
            index_connected_mol = np.where(
                self._distance_of_mol[counter] < self._cut_off)
            # match array number with molecule number: For example if molecule
            # 2 is under investigation and it is connected to molecule 3,
            # the index_connected_mol for the distance between molecule 2 and 3
            # will be 0, while it should be 3 for further calculation
            connected = (
                np.array([i for i in index_connected_mol[0]]) + counter + 1)
            connected_list = connected.tolist()
            # update "aggregate" parameter by adding attached
            # moleculs in the first level
            aggregate[counter-1].extend(connected_list)
            # "while" will allow code investigates attached molecule in
            # level 2, 3, ... as long as there is any molecule which has
            # not been investigated before
            self.connecting_loop(connected, counter, p_level, aggregate)
            aggregate[counter-1] = (
                list(dict.fromkeys(aggregate[counter-1])))
            clusters.append(counter)
        # put all aggregate variable in a single list called collect
        collect = []
        for counter in clusters:
            collect = collect + aggregate[counter-1]
        if self._molnum not in collect:
            aggregate[self._molnum-1] = [self._molnum]
            clusters = np.append(clusters, self._molnum)
        # merge aggregates if they have any in commen molecules
        self.merge_aggregates(aggregate)
        # make a dictionary with the key of number of molecules in clusters
        # and value of number of cluster
        size_agg = [0] * (len(aggregate)-1)
        for counter in range(len(aggregate)-1):
            size_agg[counter] = len(aggregate[counter])
        self.size_num[self._frame_index] = {
            j: size_agg.count(j) for j in size_agg}
        if 0 in self.size_num[self._frame_index].keys():
            del self.size_num[self._frame_index][0]
        # Give option to not considering the monomer in calculation.
        # the default is not considering monomers [True]
        if self._monomer:
            if 1 in self.size_num[self._frame_index].keys():
                del self.size_num[self._frame_index][1]
        # Give option to calculate average cluster gyration radius thoughout
        # the simulation. The default is [False]
        if self._ave_gyr_calc:
            self.average_gyration_radius_calculation(aggregate)
        self.results[self._frame_index, 0] = self._trajectory.time

    def _conclude(self):
        if self._ave_gyr_calc:
            self.results[:, 4] = np.array(self._ave_GYR)
        number_average = []
        weight_average = []
        z_average = []
        for counter in range(0, len(self.size_num)):
            size = list(self.size_num[counter].keys())
            num = list(self.size_num[counter].values())
            # number_average
            number_average.append(np.sum(np.array(num)*(np.array(size)))/
                                  np.sum(np.array(num)))
            # weight_average
            weight_average.append(np.sum(np.array(num)*(np.array(size)**2))/
                                  np.sum(np.array(num)*(np.array(size))))
            # z_average
            z_average.append(np.sum(np.array(num)*(np.array(size)**3))/
                             np.sum(np.array(num)*(np.array(size)**2)))
        self.results[:, 1] = np.array(number_average)
        self.results[:, 2] = np.array(weight_average)
        self.results[:, 3] = np.array(z_average)