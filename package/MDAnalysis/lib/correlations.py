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

"""Correlations utilities --- :mod:`MDAnalysis.lib.correlations`
=================================================================================


:Authors: Paul Smith & Mateusz Bieniek
:Year: 2020
:Copyright: GNU Public License v2

.. versionadded:: 1.0.0

This module is primarily for internal use by other analysis modules. It
provides functionality for calculating the time autocorrelation function
of a binary variable (i.e one that is either true or false at each
frame for a given atom/molecule/set of molecules). This module includes
functions for calculating both the time continuous autocorrelation and
the intermittent autocorrelation. The function :func:`autocorrelation`
calculates the continuous autocorrelation only. The data may be
pre-processed using the function :func:`intermittency` in order to
acount for intermittency before passing the results to
:func:`autocorrelation`.

This module is inspired by seemingly disparate analyses that rely on the same
underlying calculation, including the survival probability of water around
proteins :footcite:p:`ArayaSecchi2014`, hydrogen bond lifetimes
:footcite:p:`Gowers2015,ArayaSecchi2014`, and the rate of cholesterol
flip-flop in lipid bilayers :footcite:p:`Gu2019`.

.. seeAlso::

    Analysis tools that make use of modules:

        * :class:`MDAnalysis.analysis.waterdynamics.SurvivalProbability`
            Calculates the continuous or intermittent survival probability
            of an atom group in a region of interest.

        * :class:`MDAnalysis.analysis.hbonds.hbond_analysis`
            Calculates the continuous or intermittent hydrogen bond
            lifetime.

.. rubric:: References

.. footbibliography::

"""

import numpy as np
from copy import deepcopy


def autocorrelation(list_of_sets, tau_max, window_step=1):
    r"""Implementation of a discrete autocorrelation function.

    The autocorrelation of a property :math:`x` from a time :math:`t=t_0` to :math:`t=t_0 + \tau`
    is given by:

    .. math::
        C(\tau) = \langle \frac{ x(t_0)x(t_0 +\tau) }{ x(t_0)x(t_0) } \rangle

    where :math:`x` may represent any property of a particle, such as velocity or
    potential energy.

    This function is an implementation of a special case of the time
    autocorrelation function in which the property under consideration can
    be encoded with indicator variables, :math:`0` and :math:`1`, to represent the binary
    state of said property. This special case is often referred to as the
    survival probability (:math:`S(\tau)`). As an example, in calculating the survival
    probability of water molecules within 5 Å of a protein, each water
    molecule will either be within this cutoff range (:math:`1`) or not (:math:`0`). The
    total number of water molecules within the cutoff at time :math:`t_0` will be
    given by :math:`N(t_0)`. Other cases include the Hydrogen Bond Lifetime as
    well as the translocation rate of cholesterol across a bilayer.

    The survival probability of a property of a set of particles is
    given by:

    .. math::
        S(\tau) =  \langle \frac{ N(t_0, t_0 + \tau )} { N(t_0) }\rangle

    where :math:`N(t0)` is the number of particles at time :math:`t_0` for which the feature
    is observed, and :math:`N(t0, t_0 + \tau)` is the number of particles for which
    this feature is present at every frame from :math:`t_0` to :math:`N(t0, t_0 + \tau)`.
    The angular brackets represent an average over all time origins, :math:`t_0`.

    See :footcite:`ArayaSecchi2014` for a description survival probability.

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
        the raw data from which the autocorrelation is computed, i.e :math:`S(\tau)` at each window.
        This allows the time dependant evolution of :math:`S(\tau)` to be investigated.

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
    r"""Preprocess data to allow intermittent behaviour prior to calling :func:`autocorrelation`.

    Survival probabilty may be calculated either with a strict continuous requirement or
    a less strict intermittency. If calculating the survival probability water around a
    protein for example, in the former case the water must be within a cutoff distance
    of the protein at every frame from :math:`t_0` to :math:`t_0 + \tau` in order for it to be considered
    present at :math:`t_0 + \tau`. In the intermittent case, the water molecule is allowed to
    leave the region of interest for up to a specified consecutive number of frames whilst still
    being considered present at :math:`t_0 + \tau`.

    This function pre-processes data, such as the atom ids of water molecules within a cutoff
    distance of a protein at each frame, in order to allow for intermittent behaviour, with a
    single pass over the data.

    For example, if an atom is absent for a number of frames equal or smaller than the parameter
    :attr:`intermittency`, then this absence will be removed and thus the atom is considered to have
    not left.
    e.g 7,A,A,7 with `intermittency=2` will be replaced by 7,7,7,7, where A=absence.

    The returned data can be used as input to the function :func:`autocorrelation` in order
    to calculate the survival probability with a given intermittency.

    See :footcite:p:`Gowers2015` for a description of intermittency in the
    calculation of hydrogen bond lifetimes.

    # TODO - is intermittency consitent with list of sets of sets? (hydrogen bonds)

    Parameters
    ----------
    list_of_sets: list
        In the simple case of e.g survival probability, a list of sets of atom ids present at each frame, where a
        single set contains atom ids at a given frame, e.g [{0, 1}, {0}, {0}, {0, 1}]
    intermittency : int
        The maximum gap allowed. The default `intermittency=0` means that if the datapoint is missing at any frame, no
        changes are made to the data. With the value of `intermittency=2`, all datapoints missing for up to two
        consecutive frames will be instead be considered present.

    Returns
    -------
    list_of_sets: list
        returns a new list with the IDs with added IDs which disappeared for <= :attr:`intermittency`.
        e.g If [{0, 1}, {0}, {0}, {0, 1}] is a list of sets of atom ids present at each frame and `intermittency=2`,
        both atoms will be considered present throughout and thus the returned list of sets will be
        [{0, 1}, {0, 1}, {0, 1}, {0, 1}].

    """
    if intermittency == 0:
        return list_of_sets

    list_of_sets = deepcopy(list_of_sets)

    # an element (a superset) takes the form of:
    # - an atom pair when computing hydrogen bonds lifetime,
    # - atom ID in the case of water survival probability,
    for i, elements in enumerate(list_of_sets):
        # initially update each frame as seen 0 ago (now)
        seen_frames_ago = {i: 0 for i in elements}
        for j in range(1, intermittency + 2):
            for element in seen_frames_ago.keys():
                # no more frames
                if i + j >= len(list_of_sets):
                    continue

                # if the element is absent now
                if element not in list_of_sets[i + j]:
                    # increase its absence counter
                    seen_frames_ago[element] += 1
                    continue

                # the element is present
                if seen_frames_ago[element] == 0:
                    # the element was present in the last frame
                    continue

                # element was absent more times than allowed
                if seen_frames_ago[element] > intermittency:
                    continue

                # The element was absent but returned, i.e.
                # within <= intermittency_value.
                # Add it to the frames where it was absent.
                # Introduce the corrections.
                for k in range(seen_frames_ago[element], 0, -1):
                    list_of_sets[i + j - k].add(element)

                seen_frames_ago[element] = 0
    return list_of_sets

