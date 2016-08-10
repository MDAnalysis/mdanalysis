# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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

"""
MSD --- :mod:'MDAnalysis.analysis.MSD'
==============================================================

:Author: Andrew Biedermann
:Year: 2016
:Copyright: GNU Public License v2

This tool uses all available data to calculate the mean squared displacement
(MSD) as a function of time, averaging over all atoms in the atom selection
and all specified restarts. The MSD of any atom selection can be calculated,
as any atoms which remain within the atom group over an entire time interval
are included in the calculation, allowing for highly specific, position-
dependent MSD calculations. The tool is robust to zero-length atom selections,
and allows for the calculation of multiple atom selections simultaneously.

Important Note: Calculations are run on the trajectory as given to the tool,
so some preprocessing may be required to properly account for periodic
boundary conditions. For example, if calculating the MSD of a Gromacs
simulation, the following preprocessing is necessary:
    >>> gmx trjconv -f mysim.xtc -s mysim.tpr -o processed_sim.xtc -pbc nojump

MSD Tutorial
------------

The example uses files provided as part of the MDAnalysis test suite to show
how to use the MSD class.

First, import MDAnalysis, the MSD class, and the datafiles, then create a
universe object:

    >>> import MDAnalysis
    >>> from MDAnalysis.analysis.MSD import MSD
    >>> from MDAnalysisTests.datafiles import waterPSF, waterDCD
    >>>
    >>> my_universe = MDAnalysis.Universe(waterPSF, waterDCD)

Next, create an MSD object and run the calculation:

    >>> my_MSD = MSD(my_universe, ['type OT','type HT'],0,9,5)
    >>> my_MSD.run()

This line initializes an MSD object by passing it the universe object, a list
of selection groups to calculate, the first and last frame numbers to analyze,
and the number of frames between each restart point in the MSD calculation.
Note that this script is capable of calculating the MSD of any atom selection,
regardless of complexity, so the MSD of something like
 'type OT around 4 type HT' can be calculated.

By default, the MSD object will write output in binary to a pickle file named
'msd' which will contain MSD calculated for tau less than 250 frames or the
length of the trajectory, which ever is smaller. Data in this file is in the
form: [msd,n_samples], where msd is a vector of length=2 in this example,
containing vectors of MSD values with indexes corresponding to tau,
and n_samples is a similar array containing the number of distances used
to calculate the msd at the associated index.

To change these defaults, you can use the optional parameters 'out_name' and
'len_msd':

    >>> my_MSD = MSD(my_universe, ['type OT','type HT'],0,9,5,
                     out_name='my_new_name.p',len_msd=5)

Setting write_output=False will cause my_MSD.run() to return the data it would
otherwise write to the output file. Finally, setting max_data=False will allow
the calculation to use data past the specified last frame in the msd
calculation. This can be useful for splitting large jobs over multiple batch
scripts, as it ensures that no data is lost in the process.

"""


from __future__ import absolute_import, print_function, division
from six import range

import sys
import getopt
import warnings
from collections import deque

import numpy as np
import pickle

import MDAnalysis
from MDAnalysis.lib.distances import squared_distance_vector


class MSD(object):
    """Object for calculating mean square displacement (MSD)

    Attributes
    ----------
    pos : list
        contains position information, organized:
        selection X frame X atom X (x, y, z)
    atomids_per_sel_at_t : list
        contains set of atoms which satisfy each selection at each
        frame, organized:
        selection X frame X set of atom ids
    map_atomids_to_pos_at_t : list
        contains dictionaries which map atom ids to corresponding
        pos index, organized:
        selection X frame X dict
    msd : np.array(dtype=float64)
        Contains the sum of squared displacements data, organized:
        selection X tau
        (MSD is calculated as msd/n_samples)
    n_samples : np.array(dtype=int)
        Contains the number of samples corresponding to each msd
        entry, organized:
        selection X tau

    Methods
    -------
    _select(self)
        Generates list of atom groups that satisfy select_list strings
    _init_pos_deque(self):
        Initializes lists necessary for MSD calculations at the first restart.
    _update_deques(self, this_restart):
        Updates lists necessary for MSD calculations
    _process_pos_data(self):
        Runs MSD calculation
    run(self):
        Calls _process_pos_data and saves final output

    """

    def __init__(self, universe, select_list, start_frame, stop_frame,
                 dt_restart, out_name='msd', len_msd=250, write_output=True,
                 max_data=False):
        """
        Parameters
        ----------
        universe : object
            Universe object containing trajectory data
        select_list : list
            List of selection strings
        start_frame : int
            Time when analysis begins (in frames)
        stop_frame : int
            Time when analysis ends (in frames)
        dt_restart : int
            Spacing of restarts (in frames)
        out_name: str
            Name of pickle output file: (out_name).p
        len_msd : int
            Specifies the maximum tau to be calculated (frames)
            (saves memory in array allocation)
        write_output : bool
            Specifies whether to write pickle file
        max_data : bool
            If true, include tau values after stop_frame in calculation (this
            is useful for avoiding data loss when breaking large calculations
            across multiple batch jobs)
        """

        self.universe = universe
        self.select_list = select_list
        if isinstance(select_list, str):
            self.select_list = select_list.split(',')
        print(self.select_list)
        self.start_frame = int(start_frame)
        self.stop_frame = int(stop_frame)
        self.dt_restart = int(dt_restart)
        self.out_name = out_name
        # handling short simulation cases
        n_frames = int(stop_frame) - int(start_frame) + 1
        if int(len_msd) > int(n_frames):
            len_msd = int(n_frames)
        self.len_msd = int(len_msd)
        self.write_output = write_output
        self.max_data = max_data

        # initializing object attributes
        self.pos, self.atomids_per_sel_at_t, self.map_atomids_to_pos_at_t, = \
            self._init_pos_deque()

        n_sel = len(self.pos)
        self.msd = np.zeros((n_sel, self.len_msd), dtype='float64')
        self.n_samples = np.zeros((n_sel, self.len_msd), dtype='int')

        # check that trajectory length is correct
        try:
            assert len(self.universe.trajectory) >= \
                    (self.stop_frame-self.start_frame+1)
        except RuntimeError:
            print("Sample interval exceeds trajectory length. This may result"
                  "from choice of start_frame/stop_frame or insufficient RAM"
                  "when loading the trajectory.")

    def _select(self):
        """Generates list of atom groups that satisfy select_list strings

        Returns
        -------
        selections : list
            list of atom groups
        """
        selections = []
        for i, sel in enumerate(self.select_list):
            selections.append(self.universe.select_atoms(sel))
            if selections[i].n_atoms == 0:
                selections.pop(i)
                selections.append(None)
        
        return selections

    def _init_pos_deque(self):
        """Initializes lists necessary for MSD calculations at the
           first restart.

        Deques with maxlen=len_msd are used to conserve memory.

        Returns:
        --------
        pos : list
            contains position information, organized:
            selection X frame X atom X (x, y, z)
        atomids_per_sel_at_t : list
            contains set of atoms which satisfy each selection at each
            frame, organized:
            selection X frame X set of atom ids
        map_atomids_to_pos_at_t : list
            contains dictionaries which map atom ids to corresponding
            pos index, organized:
            selection X frame X dict
        """

        n_sel = len(self.select_list)
        
        # dictionaries with key: (id) and value: (array index)
        # selection X frame X dictionary
        map_atomids_to_pos_at_t = [[dict() for j in range(self.len_msd)]
                     for i in range(n_sel)]
        # atom ids which satisfy the selection at each frame
        atomids_per_sel_at_t = [[set() for j in range(self.len_msd)]
                    for i in range(n_sel)]

        # pre-allocate position array
        pos = [[np.array([]) for j in range(self.len_msd)]
               for i in range(n_sel)]

        print("Populating deque...")
        for frame, ts in self.universe.trajectory[
                         self.start_frame:self.start_frame + self.len_msd]:
            selections = self._select()

            for i, sel in enumerate(selections):
                if sel is not None:  # if there is data
                    pos[i][frame-self.start_frame] = sel.positions.copy()
                    
                    temp_list = sel.atoms.ids
                    # store set of atom ids which satisfy selection
                    atomids_per_sel_at_t[i][frame-self.start_frame] = \
                        set(temp_list)
                    # link atom id to array index
                    for j in range(len(temp_list)):
                        map_atomids_to_pos_at_t[i][frame-self.start_frame][temp_list[j]] = j
                    
            if frame % (self.len_msd/10) == 0:
                print("Frame: "+str(frame))
        # converting lists to deques
        for i in range(n_sel):
            map_atomids_to_pos_at_t[i] = deque(map_atomids_to_pos_at_t[i],
                                               maxlen=self.len_msd)
            atomids_per_sel_at_t[i] = deque(atomids_per_sel_at_t[i],
                                            maxlen=self.len_msd)
            pos[i] = deque(pos[i], maxlen=self.len_msd)
        print("Complete!")

        return pos, atomids_per_sel_at_t, map_atomids_to_pos_at_t

    def _update_deques(self, this_restart):
        """Updates lists necessary for MSD calculations

        Parameters:
        -----------
        this_restart : int
            frame number of the current restart
        """
        # stop updating when there are no more frames to analyze
        top_frame = this_restart + self.len_msd
        calc_cutoff = self.stop_frame
        if self.max_data:
            calc_cutoff = len(self.universe.trajectory) - 1

        if top_frame+self.dt_restart > calc_cutoff:
            # feed in zeros
            for i in range(len(self.select_list)):
                pos[i].append(np.array([]))
                self.atomids_per_sel_at_t[i].append(set())
                self.map_atomids_to_pos_at_t[i].append(dict())
            return
        
        for ts in self.universe.trajectory[top_frame:top_frame
                                           + self.dt_restart]:
            selections = self._select()

            for i in range(len(selections)):
                if selections[i] is not None:  # if there is data
                    pos[i].append(selections[i].positions)

                    temp_list = selections[i].atoms.ids
                    # store set of atom ids which satisfy selection
                    self.atomids_per_sel_at_t[i].append(set(temp_list))
                    # link atom id to array index
                    temp_dict = dict()
                    for j in range(len(temp_list)):
                        temp_dict[temp_list[j]] = j
                    self.map_atomids_to_pos_at_t[i].append(temp_dict)
                else:
                    pos[i].append(np.array([]))
                    self.atomids_per_sel_at_t[i].append(set())
                    self.map_atomids_to_pos_at_t[i].append(dict())

        return

    def _process_pos_data(self):
        """Runs MSD calculation"""

        n_sel = len(pos)
        
        calc_cutoff = self.stop_frame
        if self.max_data:
            calc_cutoff = len(self.universe.trajectory) - 1

        # for each restart point
        for j in range(self.start_frame, self.stop_frame+1, self.dt_restart):
            for i in range(n_sel):  # for each selection

                # storing initial set at restart
                atoms_in_sel = self.atomids_per_sel_at_t[i][0]

                # for each frame after restart
                for ts in range(j, calc_cutoff+1):
                    # avoid computing irrelevantly long taus
                    if ts-j >= self.len_msd:
                        break

                    # updating restart set to exclude any atoms which have
                    # left the selection since the restart point
                    atoms_in_sel = atoms_in_sel.intersection(
                        self.atomids_per_sel_at_t[i][ts-j])
                    
                    # find mutual atoms at times 0 and ts-j
                    shared0 = [self.map_atomids_to_pos_at_t[i][0][k]
                               for k in atoms_in_sel]
                    shared = [self.map_atomids_to_pos_at_t[i][ts-j][k]
                              for k in atoms_in_sel]
                    
                    # move to next restart if there's nothing to evaluate
                    if len(shared) == 0:
                        break  # skip to next restart
                    
                    self.msd[i][ts-j] = (self.msd[i][ts-j] +
                                         squared_distance_vector(
                                             pos[i][0][shared0],
                                             pos[i][ts-j][shared]).sum(axis=0))
                    self.n_samples[i][ts-j] += len(shared)
            
            if j % 100 == 0:
                print("Frame: "+str(j))
                if self.write_output:
                    pickle.dump([self.msd/self.n_samples, self.n_samples],
                                open(self.out_name, 'wb'))

            # Update deques
            self._update_deques(j)

        return

    def run(self):
        """Calls _process_pos_data and saves final output"""

        self._process_pos_data()

        if self.write_output:
            pickle.dump([self.msd/self.n_samples, self.n_samples],
                        open(self.out_name, 'wb'))
        # for testing purposes
        else:
            return [self.msd/self.n_samples, self.n_samples]
