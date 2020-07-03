#!/usr/bin/env python
# coding: utf-8

# In[44]:


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


def aggregate(k, tss, counter, m, Limit):
    """
    This is a subfunction that enables the code to find all 
    connected molecules in a cluster.
    k is a connected molecule to the main molecule in the first level.
    tss is a time frame.
    counter is a step.
    m provides a suitable off_set to match the number of molecules.
    Limit is the cut_off provided by the user.
    """
    import numpy as np
    bool_sector1 = globals()["d" + str(k)][:, [tss]] < Limit
    Coordinates1 = globals()["d" + str(k)][:, [tss]][bool_sector1]
    # This "if" will check if there is any uninvestigated molecule left
    if len(Coordinates1) != 0:
        index1 = np.where(bool_sector1)
        listOfCoordinates1 = list(zip(index1[0], index1[1]))
        connect = ([j[0] for j in listOfCoordinates1]
                   + globals()["connected" + str(counter)][m] + 1)
        rt = np.array(connect)
    else:
        rt = []
    return rt


def Aggregation_size(gro_file, xtc_file, number_of_molecules, Type, monomer="No", Gyr_calc="No", cut_off=False,
                        atom_1=None, atom_2=None):
    """
    1- Import the needed module
    2- Calculate the distance between molecules based on the intended criteria
    3- Find the aggregates by checking only one time each using the
    "aggregate" subfunction (to be efficient)
    4- Check and merge two aggregates as one aggregate if the distance
    between any element of them is less than cut_off
    5- Eliminate monomers (optional)
    6- measure average gyration radius (optional)
    7- measure the number of aggregates with different size
    """
    import warnings
    import numpy as np
    import MDAnalysis as mda
    from MDAnalysis.lib import distances
    u = mda.Universe(gro_file, xtc_file)
    molnum = number_of_molecules
    Limit = cut_off
    i = 0
    for it in reversed(range(1, molnum)):
        i += 1
        globals()["d" + str(i)] = np.zeros(
            it * len(u.trajectory)).reshape(it, len(u.trajectory))
    if Type == "Closest":
        cut_off = cut_off or 3.5
        j_time = 0
        for ts in u.trajectory:
            for i in range(1, molnum+1):
                ii = 0
                for j in range(i+1, molnum+1):
                    residue1 = u.select_atoms(f"resid {i}").positions
                    residue2 = u.select_atoms(f"resid {j}").positions
                    dists = distances.distance_array(
                        residue1, residue2, box=u.universe.dimensions,
                        backend='OpenMP')
                    globals()["d" + str(i)][ii][j_time] = np.amin(dists)
                    ii += 1
            j_time += 1
    elif Type == "COM":
        cut_off = cut_off or 8.5
        j_time = 0
        for ts in u.trajectory:
            cm = {}
            for x in range((molnum)):
                cm.update(
                    {x+1:u.select_atoms(f"resid {x+1}").center_of_mass()})
            cm_coordintes = np.array(list(cm.values()))
            dists = distances.distance_array(
                cm_coordintes, cm_coordintes, box=u.universe.dimensions)
            for j in range(molnum):
                ii = 0
                for i in range(j+1, molnum):
                    globals()["d" + str(j+1)][ii][j_time] = dists[i][j]
                    ii += 1
            j_time += 1
    elif Type == "Specific":
        if (not(atom_1 and atom_2)):
            warnings.warn("enter valid atom_1 and atom_2")
        cut_off = cut_off or 8.5
        j_time = 0
        for ts in u.trajectory:
            residue1 = u.select_atoms(f"name {atom_1}").positions
            residue2 = u.select_atoms(f"name {atom_2}").positions
            dists = distances.distance_array(
                residue1, residue2, box=u.universe.dimensions)
            for j in range(molnum):
                ii = 0
                for i in range(j+1, molnum):
                    globals()["d" + str(j+1)][ii][j_time] = dists[i][j]
                    ii += 1
            j_time += 1
    else:
        warnings.warn("enter valid Type")
    GYR = []
    ij = 0
    for ts in u.trajectory:
        u.trajectory[ij]
        ij += 1
        tss = ts.frame
        # P_Level collact all the investigated nodes in past levels
        P_Level = []
        # Group will colloct the molecule first investigating molecule
        # len(Group) will be the number of clusters
        Group = []
        # this for will go through all molecules to
        # check their connection with others
        for counter in range(1, molnum):
            # aggergate parameter include the molecules in each cluster
            globals()["aggregate" + str(counter)] = [counter]
            # this "if" will restrict double calculation of the same molecule
            # to increase efficincy of code
            if counter not in P_Level:
                # define the ones that are attached
                # to a molecule in first level
                bool_sector1 = globals()["d" + str(counter)][:, [tss]] < Limit
                Coordinates1 = (globals()["d" + str(counter)]
                                [:, [tss]][bool_sector1])
                index1 = np.where(bool_sector1)
                # list the connected molecules to a molecule in first level
                listOfCoordinates1 = list(zip(index1[0], index1[1]))
                # match array number with molecule number: 2-3[0] => 2-3 [3]
                globals()["connected" + str(counter)] = (
                    np.array([i[0] for i in listOfCoordinates1]) + counter + 1)
                # update "aggregate" parameter by adding attached
                # moleculs in the first level
                globals()["aggregate" + str(counter)] = (
                    np.append(globals()["aggregate" + str(counter)],
                              globals()["connected" + str(counter)]))
                # "while" will allow code investigates attached molecule in
                # level 2, 3, ... as long as there is any molecule which has
                # not been investigated before
                while len(globals()["connected" + str(counter)]) > 0:
                    m = 0
                    N_level = np.array([])
                    for k in globals()["connected" + str(counter)]:
                        k = int(k)
                        # "aggregate" is a function to find connceted molecules
                        # in levels after level1
                        # "if" will neglect function for last molecules because
                        # there is no molecule after that.
                        if k != molnum:
                            rt = aggregate(k, tss, counter, m, Limit)
                            globals()["aggregate" + str(counter)] = (
                                np.append(
                                    globals()["aggregate" + str(counter)], rt))
                        else:
                            globals()["aggregate" + str(counter)] = (
                                np.append(
                                    globals()["aggregate" + str(counter)], k))
                        # N_level find molecule attach to main aggregate for
                        # levels after level1
                        N_level = np.append(N_level, rt)
                        m += 1
                    P_Level = np.append(P_Level,
                                        globals()["connected" + str(counter)])
                    # Delete duplicated number in N_Level
                    mylist = list(dict.fromkeys(list(N_level)))
                    # Check if the molecule in the next level
                    # has been invstigated in pervious level
                    for i in set(mylist):
                        if i in (list(P_Level)):
                            mylist.remove(i)
                    # substitute new "conncected" with previous one
                    globals()["connected" + str(counter)] = np.array(mylist)
                globals()["aggregate" + str(counter)] = (
                    list(dict.fromkeys(list(globals()["aggregate" +
                                                      str(counter)]))))
                Group = np.append(Group, counter)
            else:
                continue
        Group = list(map(int, Group))
        # put all aggregate variable in a single list
        Colloct = []
        for i in Group:
            globals()["aggregate" + str(i)] = (
                np.array(globals()["aggregate" + str(i)]))
            Colloct = np.append(Colloct, globals()["aggregate" + str(i)])
        if molnum not in Colloct:
            globals()["aggregate" + str(molnum)] = [molnum]
            Group = np.append(Group, molnum)
            Group = list(map(int, Group))
        # Check if two aggregates has any element with distance less than
        # cut-off to be considered as a one aggregate
        restart = True
        while restart:
            restart = False
            for j in range(0, len(Group)-1):
                L1 = set(globals()["aggregate" + str(Group[j])])
                for k in range(j+1, len(Group)):
                    L2 = set(globals()["aggregate" + str(Group[k])])
                    if len(L1 & L2) > 0:
                        globals()["aggregate" + str(Group[j])] = (
                            np.append(globals()["aggregate" + str(Group[j])],
                                      globals()["aggregate" + str(Group[k])]))
                        del globals()["aggregate" + str(Group[k])]
                        Group.remove(Group[k])
                        globals()["aggregate" + str(Group[j])] = (
                            list(dict.fromkeys(list(
                                globals()["aggregate" + str(Group[j])]))))
                        check1 = 1
                        break
                    else:
                        check1 = 0
                if check1 == 1:
                    restart = True
                    break
        # make a dictionary with the key of number of molecules in clusters
        # and value of number of cluster
        size_agg = [0] * len(Group)
        for i in range(len(Group)):
            size_agg[i] = len(globals()["aggregate" + str(Group[i])])
        globals()["size_num" + str(ij)] = {
            j:size_agg.count(j) for j in size_agg}
        # Give option to not considering the monomer in calculation.
        # the default is not considering monomers
        if monomer=="No":
            if 1 in globals()["size_num" + str(ij)].keys():
                del globals()["size_num" + str(ij)][1]
        # Give option to calculate average cluster gyration radius thoughout
        # the simulation. The default is No
        if Gyr_calc == "Yes":
            Gyr = []
            for x in range(len(Group)):
                globals()["aggregate" + str(Group[x])] = [
                    int(elem) for elem in globals()["aggregate" + str(Group[x])]]
                if len(globals()["aggregate" + str(Group[x])]) != 1:
                    values2 = ([f'resid {i}' for i in
                                globals()["aggregate" + str(Group[x])]])
                    iio = u.select_atoms(values2[0])
                    for x in values2[1:]:
                        io = (u.select_atoms(x))
                        iio = iio + io
                    Gyr.append(iio.radius_of_gyration())
            Gyr = np.array(Gyr)
            GyR = np.average(Gyr)
            GYR.append(GyR)
    if Gyr_calc == "Yes":
        shape_data = np.zeros(1*len(u.trajectory))
        shape_data = shape_data.reshape(1, len(u.trajectory))
        GYR = np.array(GYR)
        shape_data[0][:] = GYR
        np.savetxt("shape_data.txt", shape_data, fmt='%.4f')
    number_average = []
    weight_average = []
    z_average = []
    ij = 1
    for ij in range(1, len(u.trajectory)):
        size = list(globals()["size_num" + str(ij)].keys())
        num = list(globals()["size_num" + str(ij)].values())
        # number_average
        number_average.append(np.sum(np.array(num)*(np.array(size))) /
                              np.sum(np.array(num)))
        # weight_average
        weight_average.append(np.sum(np.array(num)*(np.array(size)**2)) /
                              np.sum(np.array(num)*(np.array(size))))
        # z_average
        z_average.append(np.sum(np.array(num)*(np.array(size)**3)) /
                         np.sum(np.array(num)*(np.array(size)**2)))
    np.savetxt("number_average.txt", number_average, fmt="%.4f")
    np.savetxt("weight_average.txt", weight_average, fmt="%.4f")
    np.savetxt("z_average.txt", z_average, fmt="%.4f")


# In[50]:


# import os
# os.chdir(r"path")
# gro_file = "mdrun.gro"
# xtc_file = "mdrun.xtc"

# Aggregation_size(gro_file, xtc_file, number_of_molecules = 50, Type = "Specific", Gyr_calc="Yes", cut_off= 10)


# In[ ]:




