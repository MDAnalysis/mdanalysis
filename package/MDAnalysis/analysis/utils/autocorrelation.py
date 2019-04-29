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

from __future__ import absolute_import, division
import numpy as np
from copy import deepcopy


def autocorrelation(list_of_sets, tau_max, window_step=1):
    r"""The descrete implementation of the autocorrelation function.

    Here is a random equation that shows something.

    .. math::
       C_{HB}^c(\tau) = \frac{\sum_{ij}h_{ij}(t_0)h'_{ij}(t_0+\tau)}{\sum_{ij}h_{ij}(t_0)}


    fixme -  list_of_sets: Modifies in place!
    Parameters
    ----------
    list_of_sets : list
      List of sets,
    tau_max : int
      The last tau (inclusive) for which to carry out autocorrelation.
    window_step : int, optional
      The step for the t0 to perform autocorrelation (without the overlap). Default is 1.
    fixme - move intermittency to the right place
    intermittency : int, optional
      Helps with the graps in the data. If we want to remove some of the fluctuations and focus on the patterns,
      intermittency removes the gaps. The default intermittency=0 which means that if the datapoint is missing at any
      frame, it is not counted. With the value of 2 the datapoint can be missing for 2 consecutive frames but it will
      be counted as present.

    Returns
    --------
    tau_timeseries : list of int
        the tau for which the autocorrelation was calculated
    timeseries : list of int
        the autocorelation values for each of the tau values
    timeseries_data : list of list of int
        the raw data from which the autocorrelation is computed. The time dependant evolution can be investigated.

    .. versionadded:: 0.19.2
    """
    # fixme - add checks if this is a list of sets
    # fixme - check dimensions
    # fixme - check parameters?

    tau_timeseries = list(range(1, tau_max + 1))
    timeseries_data = [[] for _ in range(tau_max)]

    # calculate autocorrelation
    for t in range(0, len(list_of_sets), window_step):
        Nt = len(list_of_sets[t])

        if Nt == 0:
            continue

        for tau in tau_timeseries:
            if t + tau >= len(list_of_sets):
                break

            # continuous: IDs that survive from t to t + tau and at every frame in between
            Ntau = len(set.intersection(*list_of_sets[t:t + tau + 1]))
            timeseries_data[tau - 1].append(Ntau / float(Nt))

    timeseries = [np.mean(x) for x in timeseries_data]

    # at time 0 the value has to be one
    tau_timeseries.insert(0, 0)
    timeseries.insert(0, 1)

    return tau_timeseries, timeseries, timeseries_data

# fixme - is intermittency consitent with list of sets of sets? (hydrogen bonds)
def correct_intermittency(list_of_sets, intermittency):
    """
    Pre-process Consecutive Intermittency with a single pass over the data.
    If an atom is absent for a number of frames equal or smaller
    than the parameter intermittency, then correct the data and remove the absence.
    ie 7,A,A,7 with intermittency=2 will be replaced by 7,7,7,7, where A=absence

    Parameters
    ----------
    id_list: list of sets
        returns a new list with the IDs with added IDs which disappeared for <= :param intermittency
    intermittency: int
        the max gap allowed and to be corrected
    """
    if intermittency == 0:
        return list_of_sets

    list_of_sets = deepcopy(list_of_sets)

    for i, ids in enumerate(list_of_sets):
        # initially update each frame as seen 0 ago (now)
        seen_frames_ago = {i: 0 for i in ids}
        for j in range(1, intermittency + 2):
            for atomid in seen_frames_ago.keys():
                # no more frames
                if i + j >= len(list_of_sets):
                    continue

                # if the atom is absent now
                if not atomid in list_of_sets[i + j]:
                    # increase its absence counter
                    seen_frames_ago[atomid] += 1
                    continue

                # the atom is found
                if seen_frames_ago[atomid] == 0:
                    # the atom was present in the last frame
                    continue

                # it was absent more times than allowed
                if seen_frames_ago[atomid] > intermittency:
                    continue

                # the atom was absent but returned (within <= intermittency_value)
                # add it to the frames where it was absent.
                # Introduce the corrections.
                for k in range(seen_frames_ago[atomid], 0, -1):
                    list_of_sets[i + j - k].add(atomid)

                seen_frames_ago[atomid] = 0
    return list_of_sets