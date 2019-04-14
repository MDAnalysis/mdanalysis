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

"""
A flexible implementation of an autocorrelation function.
"""

import numpy as np

# TODO add documentation
def autocorrelation(list_of_sets, tau_max, window_jump=1, intermittency=0):
    """

    :param list_of_sets: Modifies in place!
    :param tau_max:
    :param window_jump:
    :param intermittency:
    :return:
    """
    # FIXME - make sure that it is not just a list of numbers

    # correct the dataset for gaps (intermittency)
    _correct_intermittency(intermittency, list_of_sets)

    # calculate Survival Probability
    tau_timeseries = list(range(1, tau_max + 1))
    sp_timeseries_data = [[] for _ in range(tau_max)]

    for t in range(0, len(list_of_sets), window_jump):
        Nt = len(list_of_sets[t])

        if Nt == 0:
            continue

        for tau in tau_timeseries:
            if t + tau >= len(list_of_sets):
                break

            # continuous: IDs that survive from t to t + tau and at every frame in between
            Ntau = len(set.intersection(*list_of_sets[t:t + tau + 1]))
            sp_timeseries_data[tau - 1].append(Ntau / float(Nt))

    sp_timeseries = [np.mean(sp) for sp in sp_timeseries_data]

    # at time 0 the value has to be one
    tau_timeseries.insert(0, 0)
    sp_timeseries.insert(0, 1)

    return tau_timeseries, sp_timeseries, sp_timeseries_data


def _correct_intermittency(intermittency, id_list, verbose=False):
    """
    Pre-process Consecutive Intermittency with a single pass over the data.
    If an atom is absent for a number of frames equal or smaller
    than the parameter intermittency, then correct the data and remove the absence.
    ie 7,A,A,7 with intermittency=2 will be replaced by 7,7,7,7, where A=absence

    :param intermittency: the max gap allowed and to be corrected
    :param id_list: modifies the selecteded IDs in place by adding atoms which left for <= :param intermittency
    """
    if intermittency == 0:
        return

    if verbose:
        print('Correcting the selected IDs for intermittancy (gaps). ')

    for i, ids in enumerate(id_list):
        # initially update each frame as seen 0 ago (now)
        seen_frames_ago = {i: 0 for i in ids}
        for j in range(1, intermittency + 2):
            for atomid in seen_frames_ago.keys():
                # no more frames
                if i + j >= len(id_list):
                    continue

                # if the atom is absent now
                if not atomid in id_list[i + j]:
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

                # the atom was absent but returned (within <= intermittency)
                # add it to the frames where it was absent.
                # Introduce the corrections.
                for k in range(seen_frames_ago[atomid], 0, -1):
                    id_list[i + j - k].add(atomid)

                seen_frames_ago[atomid] = 0