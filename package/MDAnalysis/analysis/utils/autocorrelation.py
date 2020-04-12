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
    r"""Implementation of a discrete autocorrelation function.

    The autocorrelation of a property $x$ from a time $t=t_0$ to $t=t_0 + \tau$
    is given by:
    .. math::
        C(\tau) = \langle \frac{ x(t_0)x(t_0 +\tau) }{ x(t_0)x(t_0) } \rangle

    where $x$ may represent any property of a particle, such as velocity or
    potential energy.

    The survival probability, $S(\tau)$, is a special case of the time
    autocorrelation function in which the property under consideration can
    be encoded with indicator variables, $0$ and $1$, to represent the binary
    state of said property. For instance, in calculating the survival probability
    of water molecules within $5 \rm \AA$, each water molecule will either be
    within this cutoff range ($1$) or not ($0$). The total number of water
    molecules within the cutoff at time $t_0$ will be given by $N(t_0)$.

    The survival probability of a property of a set of particles is
    given by:

    .. math::
        S(\tau) =  \langle \frac{ N(t_0, t_0 + \tau )} { N(t_0) }\rangle

    where $N(t0)$ is the number of particles at time $t_0$ for which the feature
    is observed, and $N(t0, t_0 + \tau)$ is the number of particles for which
    this feature is present at every frame from $t_0$ to $N(t0, t_0 + \tau)$.
    The angular brackets represent an average over all time origins, $t_0$.

    See Araya-Secchi et al., 2014, (https://doi.org/10.1016/j.bpj.2014.05.037)
    for a description survival probability.

    Parameters
    ----------
    list_of_sets : list
      List of sets. Each set corresponds to data from a single frame. Each element in a set
      may be, for example, an atom id or a tuple of atoms ids. In the case of calculating the
      survival probability of water around a protein, these atom ids in a given set will be
      those of the atoms which are within a cutoff distance of the protein at a given frame.
    tau_max : int
      The last tau (lag time, inclusive) for which to calculate the autocorrelation. e.g if tau_max = 20,
      the survival probability will be calculated over 20 frames.
    window_step : int, optional
      The step size for t0 to perform autocorrelation. Ideally, window_step will be larger than
       tau_max to ensure independence of each window for which the calculation is performed.
       Default is 1.

    Returns
    --------
    tau_timeseries : list of int
        the values of tau for which the autocorrelation was calculated
    timeseries : list of int
        the autocorelation values for each of the tau values
    timeseries_data : list of list of int
        the raw data from which the autocorrelation is computed, i.e $S(\tau)$ at each window.
        This allows the time dependant evolution of $S(\tau)$ to be investigated.

    .. versionadded:: 0.19.2
    """

    # check types
    if (type(list_of_sets) != list and len(list_of_sets) != 0) or type(list_of_sets[0]) != set:
        raise TypeError("list_of_sets must be a one-dimensional list of sets")  # pragma: no cover

    # Check dimensions of parameters
    if len(list_of_sets) < tau_max:
        raise ValueError("tau_max cannot be greater than the length of list_of_sets") # pragma: no cover

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


def correct_intermittency(list_of_sets, intermittency):
    """Preprocess data to allow intermittent behaviour prior to calling `autocorrelation`.

    Process consecutive intermittency with a single pass over the data. The returned data can be used as input to
    the function `autocorrelation` in order to calculate the survival probability with a given
    intermittency.

    For example, if an atom is absent for a number of frames equal or smaller than the parameter `intermittency`,
    then correct the data and remove the absence.
    e.g 7,A,A,7 with `intermittency=2` will be replaced by 7,7,7,7, where A=absence.

    See Gowers and Carbonne, 2015, (DOI:10.1063.1.4922445) for a description of
    intermittency in the calculation of hydrogen bond lifetimes.

    # TODO - is intermittency consitent with list of sets of sets? (hydrogen bonds)

    Parameters
    ----------
    list_of_sets: list
        In the simple case of e.g survival probability, a list of sets of atom ids present at each frame, where a
        single set contains atom ids at a given frame, e.g [{0, 1}, {0}, {0}, {0, 1}]
    intermittency : int
      The maximum gap allowed. The default intermittency=0 means that if the datapoint is missing at any frame, no
      changes are made to the data. With the value of `intermittency=2`, all datapoints missing for up to two
       consecutive frames will be instead be considered present.

    Returns
    -------
        returns a new list with the IDs with added IDs which disappeared for <= :param intermittency.
        e.g If [{0, 1}, {0}, {0}, {0, 1}] is a list of sets of atom ids present at each frame and `intermittency=2`,
        both atoms will be considered present throughout and thus the returned list of sets will be
        [{0, 1}, {0, 1}, {0, 1}, {0, 1}].

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
