r""" Average Aggregation number ---
======================================================
This module is able to calculate the number average :math: 'g_{n}', weight
average math: 'g_{w}', and z average :math: 'g_{z}' aggregation number
based on the three different categories including the distance between
closest atom, the distance between the centres of mass, and the distance
between the two certain atoms of two molecules. In this calculation, the
monomers are not considered as a default, but the user can change it
by inserting (monomer = "Yes") as an argument in the fuction. The formula
of each aggregation type is:
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
--------------------
The 'Closest' Type
--------------------
In this method, the distance between the closest atom of two molecules
will be considered as a criterion for considering two molecules as one
aggregate. The default cut_off distance is 3.5 Angstrom but can be
changed by the user with an argument in the function.
Example: cut_off = 10 => desired value in Angstrom
--------------------
The 'COM' Type
--------------------
In this method, the distance between the centres of mass of two molecules
will be considered as a criterion for considering two molecules as one
aggregate. The default cut_off distance is 8.5 Angstrom but can be
changed by the user with an argument in the function.
Example: cut_off = 10 => desired value in Angstrom
--------------------
The 'Specific' Type
--------------------
In this method, the distance between the two selected atoms of molecules
will be considered as a criterion for considering two molecules as a one
aggregate. Besides the cut_off distance, the two atoms' name (atom_1 and
atom_2) should be defined. The default value for cut_off distance is 8.5
Angstrom while there is no default value for atom names.
---------------------
P.S.
---------------------
1- This code does not consider the residue of molecules. So as a prerequisite,
you may need to eliminate the non-target molecules after simulation and before
runnig the code.
2- The periodic boundary should be eliminated before running the code. We
recommend doing -pbc whole and -pbc cluster using trjconv keyword (available)
in GROMCAS or any similar procedure.
3- If you have more than one type of molecule and you want to measure the
aggregation without considering the type of molecule the code will works
perfectly for "Closest" and "COM" Types while for the "Specific" the results
will not be correct. The reasone is that you cannot tell atom_1 and atom_2 is
related to which residue!
4- If you have more than one type of molecule and you want to measure
aggregation for each type of molecule separately, you can keep only the
target molecule in the box and run the code multiple times for different
molecule type.
5- You will get an error if there is no aggregation in the box or if the
aggregation distance choose too short.
"""
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib import distances
from MDAnalysis.analysis.base import (AnalysisBase)


class Aggregation_size(AnalysisBase):  # subclass AnalysisBase
    def __init__(self, atomgroup, number_of_molecules, Type="Closest",
                 no_monomer=True, Gyr_calc=False, cut_off=False,
                 atom_1=None, atom_2=None, **kwargs):
        """
        Set up the initial analysis parameters.
        """
        # must first run AnalysisBase.__init__ and pass the trajectory
        trajectory = atomgroup.universe.trajectory
        super(Aggregation_size, self).__init__(trajectory)
        # set atomgroup as a property for access in other methods
        self.atomgroup = atomgroup
        self._molnum = number_of_molecules
        self._type = Type
        self._cut_off = cut_off
        self._monomer = no_monomer
        self._Gyr_calc = Gyr_calc
        self._GYR = []

    def closest_atom(self, d):
        """
        This function calculates the distance between the closest atom
        of two molecules.
        """
        self._cut_off = self._cut_off or 3.5
        for i in range(1, self._molnum+1):
            ii = 0
            for j in range(i+1, self._molnum+1):
                residue1 = u.select_atoms(f"resid {i}").positions
                residue2 = u.select_atoms(f"resid {j}").positions
                dists = distances.distance_array(
                    residue1, residue2, box=u.universe.dimensions,
                    backend='OpenMP')
                d[i][ii] = np.amin(dists)
                ii += 1

    def center_of_mass(self, d):
        """
        This function calculates the distance between the COM
        of two molecules.
        """
        self._cut_off = self._cut_off or 8.5
        cm = {}
        for x in range((self._molnum)):
            cm[x+1] = u.select_atoms(f"resid {x+1}").center_of_mass()
        cm_coordintes = np.array(list(cm.values()))
        dists = distances.distance_array(
            cm_coordintes, cm_coordintes, box=u.universe.dimensions)
        for j in range(self._molnum):
            ii = 0
            for i in range(j+1, self._molnum):
                d[j+1][ii] = dists[i][j]
                ii += 1

    def specific_atoms(self, d):
        """
        This function calculates the distance between the selected atom
        of two molecules.
        """
        if (not(atom_1 or atom_2)):
            raise ValueError("enter valid atom_1 and atom_2")
        self._cut_off = self._cut_off or 8.5
        residue1 = u.select_atoms(f"name {atom_1}").positions
        residue2 = u.select_atoms(f"name {atom_2}").positions
        dists = distances.distance_array(
            residue1, residue2, box=u.universe.dimensions)
        for j in range(self._molnum):
            ii = 0
            for i in range(j+1, self._molnum):
                d[j+1][ii] = dists[i][j]
                ii += 1

    def _prepare(self):
        """
        Create array of zeroes as a placeholder for differnet parameters.
        """
        self.results = np.zeros((self.n_frames, 5))
        self.d = [0] * self._molnum
        i = 0
        for it in reversed(range(1, self._molnum+1)):
            self.d[i] = np.zeros(it * 1).reshape(it, 1)
            i += 1
        self.size_num = [0] * self.n_frames

    def aggregation(self, k, tss, counter, m, cut_off, connected):
        """
        This is a subfunction that enables the code to find all
        connected molecules in a cluster.
        k is a connected molecule to the main molecule in the first level.
        self._frame_index is a time frame.
        counter is a step.
        m provides a suitable off_set to match the number of molecules.
        Limit is the cut_off provided by the user.
        """
        bool_sector = self.d[k] < cut_off
        Coordinates = self.d[k][bool_sector]
        # This "if" will check if there is any uninvestigated molecule left
        if len(Coordinates) != 0:
            index1 = np.where(bool_sector)
            listOfCoordinates = list(zip(index1[0], index1[1]))
            connect = ([j[0] for j in listOfCoordinates]
                       + connected[m] + 1)
            Connected_NLevel = np.array(connect)
        else:
            Connected_NLevel = []
        return Connected_NLevel

    def _single_frame(self):
        """
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
        """
        aggregate = [0] * self._molnum
        if self._type == "Closest":
            self.closest_atom(self.d)
        elif self._type == "COM":
            self.center_of_mass(self.d)
        elif self._type == "Specific":
            self.specific_atoms(self.d)
        else:
            raise ValueError("enter valid Type")
        # P_Level collact all the investigated nodes in past levels
        P_Level = []
        # Group will collect the molecule first investigating molecule
        # len(Group) will be the number of clusters
        Group = []
        # this for will go through all molecules to
        # check their connection with others
        for counter in range(1, self._molnum):
            # aggergate parameter include the molecules in each cluster
            aggregate[counter-1] = [counter]
            # this "if" will restrict double calculation of the same molecule
            # to increase efficincy of code
            if counter not in P_Level:
                # define the ones that are attached
                # to a molecule in first level
                bool_sector = self.d[counter] < self._cut_off
                Coordinates = self.d[counter][bool_sector]
                index1 = np.where(bool_sector)
                # list the connected molecules to a molecule in first level
                listOfCoordinates = list(zip(index1[0], index1[1]))
                # match array number with molecule number: 2-3[0] => 2-3 [3]
                connected = (
                    np.array([i[0] for i in listOfCoordinates]) + counter + 1)
                connected_list = connected.tolist()
                # update "aggregate" parameter by adding attached
                # moleculs in the first level
                aggregate[counter-1].extend(connected_list)
                # "while" will allow code investigates attached molecule in
                # level 2, 3, ... as long as there is any molecule which has
                # not been investigated before
                while len(connected) > 0:
                    m = 0
                    N_level = np.array([])
                    for k in connected:
                        k = int(k)
                        # "aggregate" is a function to find connceted molecules
                        # in levels after level1
                        # "if" will neglect function for last molecules because
                        # there is no molecule after that.
                        if k != self._molnum:
                            Connected_NLevel = self.aggregation(
                                k, self._frame_index, counter, m,
                                self._cut_off, connected)
                            aggregate[counter-1].extend(Connected_NLevel)
                        else:
                            aggregate[counter-1].extend([k])
                        # N_level find molecule attach to main aggregate for
                        # levels after level1
                        N_level = np.append(N_level, Connected_NLevel)
                        m += 1
                    P_Level = np.append(P_Level,
                                        connected)
                    # Delete duplicated number in N_Level
                    mylist = list(dict.fromkeys(list(N_level)))
                    # Check if the molecule in the next level
                    # has been invstigated in pervious level
                    for i in set(mylist):
                        if i in (list(P_Level)):
                            mylist.remove(i)
                    # substitute new "conncected" with previous one
                    connected = np.array(mylist)
                aggregate[counter-1] = (
                    list(dict.fromkeys(aggregate[counter-1])))
                Group = np.append(Group, counter)
            else:
                continue
        Group = list(map(int, Group))
        # put all aggregate variable in a single list
        Collect = []
        for counter in Group:
            Collect = Collect + aggregate[counter-1]
        if self._molnum not in Collect:
            aggregate[self._molnum-1] = [self._molnum]
            Group = np.append(Group, self._molnum)
            Group = list(map(int, Group))
        # Check if two aggregates has any element with distance less than
        # cut-off to be considered as a one aggregate
        restart = True
        while restart:
            restart = False
            for j in range(0, len(aggregate)-1):
                aggregate1 = set(aggregate[j])
                for k in range(j+1, len(aggregate)-1):
                    aggregate2 = set(aggregate[k])
                    if len(aggregate1 & aggregate2) > 0:
                        aggregate[j] = aggregate[j] + aggregate[k]
                        del aggregate[k]
                        aggregate[j] = (list(dict.fromkeys(aggregate[j])))
                        check1 = 1
                        break
                    else:
                        check1 = 0
                if check1 == 1:
                    restart = True
                    break
        # make a dictionary with the key of number of molecules in clusters
        # and value of number of cluster
        size_agg = [0] * (len(aggregate)-1)
        for counter in range(len(aggregate)-1):
            size_agg[counter] = len(aggregate[counter])
        self.size_num[self._frame_index] = {
            j:size_agg.count(j) for j in size_agg}
        if 0 in self.size_num[self._frame_index].keys():
            del self.size_num[self._frame_index][0]
        # Give option to not considering the monomer in calculation.
        # the default is not considering monomers
        if self._monomer:
            if 1 in self.size_num[self._frame_index].keys():
                del self.size_num[self._frame_index][1]
        # Give option to calculate average cluster gyration radius thoughout
        # the simulation. The default is No
        print("first", self._Gyr_calc)
        if self._Gyr_calc:
            Gyr = []
            for x in range(len(aggregate)-1):
                aggregate[x] = [
                    int(elem) for elem in aggregate[x]]
                if len(aggregate[x]) != 1:
                    values2 = ([f'resid {i}' for i in aggregate[x]])
                    iio = u.select_atoms(values2[0])
                    for x in values2[1:]:
                        io = (u.select_atoms(x))
                        iio = iio + io
                    Gyr.append(iio.radius_of_gyration())
            Gyr = np.array(Gyr)
            GyR = np.average(Gyr)
            self._GYR.append(GyR)
        self.results[self._frame_index, 0] = self._trajectory.time

    def _conclude(self):
        if self._Gyr_calc:
            self.results[:, 4] = np.array(self._GYR)
        number_average = []
        weight_average = []
        z_average = []
        for counter in range(0, len(self.size_num)):
            size = list(self.size_num[counter].keys())
            num = list(self.size_num[counter].values())
            # number_average
            number_average.append(np.sum(np.array(num)*(np.array(size))) /
                                  np.sum(np.array(num)))
            # weight_average
            weight_average.append(np.sum(np.array(num)*(np.array(size)**2)) /
                                  np.sum(np.array(num)*(np.array(size))))
            # z_average
            z_average.append(np.sum(np.array(num)*(np.array(size)**3)) /
                             np.sum(np.array(num)*(np.array(size)**2)))
        self.results[:, 1] = np.array(number_average)
        self.results[:, 2] = np.array(weight_average)
        self.results[:, 3] = np.array(z_average)
        columns = ['Frame', 'number_average',
                   'weight_average', 'z_average',
                   'Radius of Gyration']
        self.df = pd.DataFrame(self.results, columns=columns)



import os
os.chdir(r"path")
gro_file = "mdrun.gro"
xtc_file = "mdrun.xtc"
u = mda.Universe(gro_file, xtc_file)
Mol = u.select_atoms('resname AsphC')
rog_base = Aggregation_size(Mol, number_of_molecules = 50, Type = "COM", 
no_monomer = True, Gyr_calc=True).run()

print(rog_base.results)
rog_base.df





