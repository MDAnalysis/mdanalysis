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

"""Water dynamics analysis --- :mod:`MDAnalysis.analysis.waterdynamics`
=======================================================================

:Author: Alejandro Bernardin
:Year: 2014-2015
:Copyright: GNU Public License v3

.. versionadded:: 0.11.0

.. deprecated:: 2.8.0
  This module is deprecated in favour of the mdakit
  `waterdynamics <https://www.mdanalysis.org/waterdynamics/>`_ and
  will be removed in MDAnalysis 3.0.0.

See Also
--------
:mod:`waterdynamics.waterdynamics`

"""
from MDAnalysis.lib.correlations import autocorrelation, correct_intermittency
import MDAnalysis.analysis.hbonds
from itertools import zip_longest
import logging
import warnings
import numpy as np


logger = logging.getLogger('MDAnalysis.analysis.waterdynamics')
from MDAnalysis.lib.log import ProgressBar


class WaterOrientationalRelaxation(object):
    r"""Water orientation relaxation analysis

    Function to evaluate the Water Orientational Relaxation proposed by Yu-ling
    Yeh and Chung-Yuan Mou :footcite:p:`Yeh1999`. WaterOrientationalRelaxation
    indicates "how fast" water molecules are rotating or changing direction.
    This is a time correlation function given by:

    .. math::
        C_{\hat u}(\tau)=\langle \mathit{P}_2[\mathbf{\hat{u}}(t_0)\cdot\mathbf{\hat{u}}(t_0+\tau)]\rangle

    where :math:`P_2=(3x^2-1)/2` is the second-order Legendre polynomial and :math:`\hat{u}` is
    a unit vector along HH, OH or dipole vector.


    Parameters
    ----------
    universe : Universe
      Universe object
    selection : str
      Selection string for water [‘byres name OH2’].
    t0 : int
      frame  where analysis begins
    tf : int
      frame where analysis ends
    dtmax : int
      Maximum dt size, `dtmax` < `tf` or it will crash.


    .. versionadded:: 0.11.0

    .. versionchanged:: 1.0.0
       Changed `selection` keyword to `select`
    """

    def __init__(self, universe, select, t0, tf, dtmax, nproc=1):
        self.universe = universe
        self.selection = select
        self.t0 = t0
        self.tf = tf
        self.dtmax = dtmax
        self.nproc = nproc
        self.timeseries = None

    def _repeatedIndex(self, selection, dt, totalFrames):
        """
        Indicates the comparation between all the t+dt.
        The results is a list of list with all the repeated index per frame
        (or time).
        Ex: dt=1, so compare frames (1,2),(2,3),(3,4)...
        Ex: dt=2, so compare frames (1,3),(3,5),(5,7)...
        Ex: dt=3, so compare frames (1,4),(4,7),(7,10)...
        """
        rep = []
        for i in range(int(round((totalFrames - 1) / float(dt)))):
            if (dt * i + dt < totalFrames):
                rep.append(self._sameMolecTandDT(
                    selection, dt * i, (dt * i) + dt))
        return rep

    def _getOneDeltaPoint(self, universe, repInd, i, t0, dt):
        """
        Gives one point to calculate the mean and gets one point of the plot
        C_vect vs t.
        Ex: t0=1 and tau=1 so calculate the t0-tau=1-2 intervale.
        Ex: t0=5 and tau=3 so calcultate the t0-tau=5-8 intervale.
        i = come from getMeanOnePoint (named j) (int)
        """
        valOH = 0
        valHH = 0
        valdip = 0
        n = 0
        for j in range(len(repInd[i]) // 3):
            begj = 3 * j
            universe.trajectory[t0]
            Ot0 = repInd[i][begj]
            H1t0 = repInd[i][begj + 1]
            H2t0 = repInd[i][begj + 2]
            OHVector0 = H1t0.position - Ot0.position
            HHVector0 = H1t0.position - H2t0.position
            dipVector0 = ((H1t0.position + H2t0.position) * 0.5) - Ot0.position

            universe.trajectory[t0 + dt]
            Otp = repInd[i][begj]
            H1tp = repInd[i][begj + 1]
            H2tp = repInd[i][begj + 2]

            OHVectorp = H1tp.position - Otp.position
            HHVectorp = H1tp.position - H2tp.position
            dipVectorp = ((H1tp.position + H2tp.position) * 0.5) - Otp.position

            normOHVector0 = np.linalg.norm(OHVector0)
            normOHVectorp = np.linalg.norm(OHVectorp)
            normHHVector0 = np.linalg.norm(HHVector0)
            normHHVectorp = np.linalg.norm(HHVectorp)
            normdipVector0 = np.linalg.norm(dipVector0)
            normdipVectorp = np.linalg.norm(dipVectorp)

            unitOHVector0 = [OHVector0[0] / normOHVector0,
                             OHVector0[1] / normOHVector0,
                             OHVector0[2] / normOHVector0]
            unitOHVectorp = [OHVectorp[0] / normOHVectorp,
                             OHVectorp[1] / normOHVectorp,
                             OHVectorp[2] / normOHVectorp]
            unitHHVector0 = [HHVector0[0] / normHHVector0,
                             HHVector0[1] / normHHVector0,
                             HHVector0[2] / normHHVector0]
            unitHHVectorp = [HHVectorp[0] / normHHVectorp,
                             HHVectorp[1] / normHHVectorp,
                             HHVectorp[2] / normHHVectorp]
            unitdipVector0 = [dipVector0[0] / normdipVector0,
                              dipVector0[1] / normdipVector0,
                              dipVector0[2] / normdipVector0]
            unitdipVectorp = [dipVectorp[0] / normdipVectorp,
                              dipVectorp[1] / normdipVectorp,
                              dipVectorp[2] / normdipVectorp]

            valOH += self.lg2(np.dot(unitOHVector0, unitOHVectorp))
            valHH += self.lg2(np.dot(unitHHVector0, unitHHVectorp))
            valdip += self.lg2(np.dot(unitdipVector0, unitdipVectorp))
            n += 1
        return  (valOH/n, valHH/n, valdip/n) if n > 0 else (0, 0, 0)


    def _getMeanOnePoint(self, universe, selection1, selection_str, dt,
                         totalFrames):
        """
        This function gets one point of the plot C_vec vs t. It uses the
        _getOneDeltaPoint() function to calculate the average.

        """
        repInd = self._repeatedIndex(selection1, dt, totalFrames)
        sumsdt = 0
        n = 0.0
        sumDeltaOH = 0.0
        sumDeltaHH = 0.0
        sumDeltadip = 0.0

        for j in range(totalFrames // dt - 1):
            a = self._getOneDeltaPoint(universe, repInd, j, sumsdt, dt)
            sumDeltaOH += a[0]
            sumDeltaHH += a[1]
            sumDeltadip += a[2]
            sumsdt += dt
            n += 1

        # if no water molecules remain in selection, there is nothing to get
        # the mean, so n = 0.
        return (sumDeltaOH / n, sumDeltaHH / n, sumDeltadip / n) if n > 0 else (0, 0, 0)

    def _sameMolecTandDT(self, selection, t0d, tf):
        """
        Compare the molecules in the t0d selection and the t0d+dt selection and
        select only the particles that are repeated in both frame. This is to
        consider only the molecules that remains in the selection after the dt
        time has elapsed.
        The result is a list with the indices of the atoms.
        """
        a = set(selection[t0d])
        b = set(selection[tf])
        sort = sorted(list(a.intersection(b)))
        return sort

    def _selection_serial(self, universe, selection_str):
        selection = []
        for ts in ProgressBar(universe.trajectory, verbose=True,
                              total=universe.trajectory.n_frames):
            selection.append(universe.select_atoms(selection_str))
        return selection

    @staticmethod
    def lg2(x):
        """Second Legendre polynomial"""
        return (3*x*x - 1)/2

    def run(self, **kwargs):
        """Analyze trajectory and produce timeseries"""

        # All the selection to an array, this way is faster than selecting
        # later.
        if self.nproc == 1:
            selection_out = self._selection_serial(
                self.universe, self.selection)
        else:
            # selection_out = self._selection_parallel(self.universe,
            # self.selection, self.nproc)
            # parallel selection to be implemented
            selection_out = self._selection_serial(
                self.universe, self.selection)
        self.timeseries = []
        for dt in list(range(1, self.dtmax + 1)):
            output = self._getMeanOnePoint(
                self.universe, selection_out, self.selection, dt, self.tf)
            self.timeseries.append(output)


class AngularDistribution(object):
    r"""Angular distribution function analysis

    The angular distribution function (AD) is defined as the distribution
    probability of the cosine of the :math:`\theta` angle formed by the OH
    vector, HH vector or dipolar vector of water molecules and a vector
    :math:`\hat n` parallel to chosen axis (z is the default value). The cosine
    is define as :math:`\cos \theta = \hat u \cdot \hat n`, where :math:`\hat
    u` is OH, HH or dipole vector.  It creates a histogram and returns a list
    of lists, see Output_. The AD is also know as Angular Probability (AP).


    Parameters
    ----------
    universe : Universe
        Universe object
    select : str
        Selection string to evaluate its angular distribution ['byres name OH2']
    bins : int (optional)
        Number of bins to create the histogram by means of :func:`numpy.histogram`
    axis : {'x', 'y', 'z'} (optional)
        Axis to create angle with the vector (HH, OH or dipole) and calculate
        cosine theta ['z'].


    .. versionadded:: 0.11.0

    .. versionchanged:: 1.0.0
       Changed `selection` keyword to `select`
    """

    def __init__(self, universe, select, bins=40, nproc=1, axis="z"):
        self.universe = universe
        self.selection_str = select
        self.bins = bins
        self.nproc = nproc
        self.axis = axis
        self.graph = None

    def _getCosTheta(self, universe, selection, axis):
        valOH = []
        valHH = []
        valdip = []

        i = 0
        while i <= (len(selection) - 1):
            universe.trajectory[i]
            line = selection[i].positions

            Ot0 = line[::3]
            H1t0 = line[1::3]
            H2t0 = line[2::3]

            OHVector0 = H1t0 - Ot0
            HHVector0 = H1t0 - H2t0
            dipVector0 = (H1t0 + H2t0) * 0.5 - Ot0

            unitOHVector0 = OHVector0 / \
                np.linalg.norm(OHVector0, axis=1)[:, None]
            unitHHVector0 = HHVector0 / \
                np.linalg.norm(HHVector0, axis=1)[:, None]
            unitdipVector0 = dipVector0 / \
                np.linalg.norm(dipVector0, axis=1)[:, None]

            j = 0
            while j < len(line) / 3:
                if axis == "z":
                    valOH.append(unitOHVector0[j][2])
                    valHH.append(unitHHVector0[j][2])
                    valdip.append(unitdipVector0[j][2])

                elif axis == "x":
                    valOH.append(unitOHVector0[j][0])
                    valHH.append(unitHHVector0[j][0])
                    valdip.append(unitdipVector0[j][0])

                elif axis == "y":
                    valOH.append(unitOHVector0[j][1])
                    valHH.append(unitHHVector0[j][1])
                    valdip.append(unitdipVector0[j][1])

                j += 1
            i += 1
        return (valOH, valHH, valdip)

    def _getHistogram(self, universe, selection, bins, axis):
        """
        This function gets a normalized histogram of the cos(theta) values. It
        return a list of list.
        """
        a = self._getCosTheta(universe, selection, axis)
        cosThetaOH = a[0]
        cosThetaHH = a[1]
        cosThetadip = a[2]
        lencosThetaOH = len(cosThetaOH)
        lencosThetaHH = len(cosThetaHH)
        lencosThetadip = len(cosThetadip)
        histInterval = bins
        histcosThetaOH = np.histogram(cosThetaOH, histInterval, density=True)
        histcosThetaHH = np.histogram(cosThetaHH, histInterval, density=True)
        histcosThetadip = np.histogram(cosThetadip, histInterval, density=True)

        return (histcosThetaOH, histcosThetaHH, histcosThetadip)

    def _hist2column(self, aList):
        """
        This function transform from the histogram format
        to a column format.
        """
        a = []
        for x in zip_longest(*aList, fillvalue="."):
            a.append(" ".join(str(i) for i in x))
        return a

    def run(self, **kwargs):
        """Function to evaluate the angular distribution of cos(theta)"""

        if self.nproc == 1:
            selection = self._selection_serial(
                self.universe, self.selection_str)
        else:
            # not implemented yet
            # selection = self._selection_parallel(self.universe,
            # self.selection_str,self.nproc)
            selection = self._selection_serial(
                self.universe, self.selection_str)

        self.graph = []
        output = self._getHistogram(
            self.universe, selection, self.bins, self.axis)
        # this is to format the exit of the file
        # maybe this output could be improved
        listOH = [list(output[0][1]), list(output[0][0])]
        listHH = [list(output[1][1]), list(output[1][0])]
        listdip = [list(output[2][1]), list(output[2][0])]

        self.graph.append(self._hist2column(listOH))
        self.graph.append(self._hist2column(listHH))
        self.graph.append(self._hist2column(listdip))

    def _selection_serial(self, universe, selection_str):
        selection = []
        for ts in ProgressBar(universe.trajectory, verbose=True,
                              total=universe.trajectory.n_frames):
            selection.append(universe.select_atoms(selection_str))
        return selection


class MeanSquareDisplacement(object):
    r"""Mean square displacement analysis

    Function to evaluate the Mean Square Displacement (MSD_). The MSD gives the
    average distance that particles travels. The MSD is given by:

    .. math::
        \langle\Delta r(t)^2\rangle = 2nDt

    where :math:`r(t)` is the position of particle in time :math:`t`,
    :math:`\Delta r(t)` is the displacement after time lag :math:`t`,
    :math:`n` is the dimensionality, in this case :math:`n=3`,
    :math:`D` is the diffusion coefficient and :math:`t` is the time.

    .. _MSD: http://en.wikipedia.org/wiki/Mean_squared_displacement


    Parameters
    ----------
    universe : Universe
      Universe object
    select : str
      Selection string for water [‘byres name OH2’].
    t0 : int
      frame  where analysis begins
    tf : int
      frame where analysis ends
    dtmax : int
      Maximum dt size, `dtmax` < `tf` or it will crash.


    .. versionadded:: 0.11.0

    .. versionchanged:: 1.0.0
       Changed `selection` keyword to `select`
    """

    def __init__(self, universe, select, t0, tf, dtmax, nproc=1):
        self.universe = universe
        self.selection = select
        self.t0 = t0
        self.tf = tf
        self.dtmax = dtmax
        self.nproc = nproc
        self.timeseries = None

    def _repeatedIndex(self, selection, dt, totalFrames):
        """
        Indicate the comparation between all the t+dt.
        The results is a list of list with all the repeated index per frame
        (or time).

        - Ex: dt=1, so compare frames (1,2),(2,3),(3,4)...
        - Ex: dt=2, so compare frames (1,3),(3,5),(5,7)...
        - Ex: dt=3, so compare frames (1,4),(4,7),(7,10)...
        """
        rep = []
        for i in range(int(round((totalFrames - 1) / float(dt)))):
            if (dt * i + dt < totalFrames):
                rep.append(self._sameMolecTandDT(
                    selection, dt * i, (dt * i) + dt))
        return rep

    def _getOneDeltaPoint(self, universe, repInd, i, t0, dt):
        """
        Gives one point to calculate the mean and gets one point of the plot
        C_vect vs t.

        - Ex: t0=1 and dt=1 so calculate the t0-dt=1-2 interval.
        - Ex: t0=5 and dt=3 so calcultate the t0-dt=5-8 interval

        i = come from getMeanOnePoint (named j) (int)
        """
        valO = 0
        n = 0
        for j in range(len(repInd[i]) // 3):
            begj = 3 * j
            universe.trajectory[t0]
            # Plus zero is to avoid 0to be equal to 0tp
            Ot0 = repInd[i][begj].position + 0

            universe.trajectory[t0 + dt]
            # Plus zero is to avoid 0to be equal to 0tp
            Otp = repInd[i][begj].position + 0

            # position oxygen
            OVector = Ot0 - Otp
            # here it is the difference with
            # waterdynamics.WaterOrientationalRelaxation
            valO += np.dot(OVector, OVector)
            n += 1

        # if no water molecules remain in selection, there is nothing to get
        # the mean, so n = 0.
        return valO/n if n > 0 else 0

    def _getMeanOnePoint(self, universe, selection1, selection_str, dt,
                         totalFrames):
        """
        This function gets one point of the plot C_vec vs t. It's uses the
        _getOneDeltaPoint() function to calculate the average.

        """
        repInd = self._repeatedIndex(selection1, dt, totalFrames)
        sumsdt = 0
        n = 0.0
        sumDeltaO = 0.0
        valOList = []

        for j in range(totalFrames // dt - 1):
            a = self._getOneDeltaPoint(universe, repInd, j, sumsdt, dt)
            sumDeltaO += a
            valOList.append(a)
            sumsdt += dt
            n += 1

        # if no water molecules remain in selection, there is nothing to get
        # the mean, so n = 0.
        return sumDeltaO/n if n > 0 else 0

    def _sameMolecTandDT(self, selection, t0d, tf):
        """
        Compare the molecules in the t0d selection and the t0d+dt selection and
        select only the particles that are repeated in both frame. This is to
        consider only the molecules that remains in the selection after the dt
        time has elapsed. The result is a list with the indexs of the atoms.
        """
        a = set(selection[t0d])
        b = set(selection[tf])
        sort = sorted(list(a.intersection(b)))
        return sort

    def _selection_serial(self, universe, selection_str):
        selection = []
        for ts in ProgressBar(universe.trajectory, verbose=True,
                              total=universe.trajectory.n_frames):
            selection.append(universe.select_atoms(selection_str))
        return selection

    def run(self, **kwargs):
        """Analyze trajectory and produce timeseries"""

        # All the selection to an array, this way is faster than selecting
        # later.
        if self.nproc == 1:
            selection_out = self._selection_serial(
                self.universe, self.selection)
        else:
            # parallel not yet implemented
            # selection = selection_parallel(universe, selection_str, nproc)
            selection_out = self._selection_serial(
                self.universe, self.selection)
        self.timeseries = []
        for dt in list(range(1, self.dtmax + 1)):
            output = self._getMeanOnePoint(
                self.universe, selection_out, self.selection, dt, self.tf)
            self.timeseries.append(output)


class SurvivalProbability(object):
    r"""
    Survival Probability (SP) gives the probability for a group of particles to remain
    in a certain region. The SP is given by:

    .. math::
        P(\tau) = \langle \frac{ N(t, t + \tau )} { N(t) }\rangle

    where :math:`\tau` is the timestep, :math:`N(t)` the number of particles at time 
    :math:`t`, and :math:`N(t, t+\tau)` is the number of particles at every frame from 
    :math:`t` to :math:`t + \tau`. The angular brackets represent an average over all time 
    origins, :math:`t`. See :func:`MDAnalysis.lib.correlations.autocorrelation` for
    technical details.


    Parameters
    ----------
    universe : Universe
      Universe object
    select : str
      Selection string; any selection is allowed. With this selection you
      define the region/zone where to analyze, e.g.: "resname SOL and around 5 
      (resid 10)". See `SP-examples`_.
    verbose : Boolean, optional
      When True, prints progress and comments to the console.


    Notes
    -----
    Currently :class:`SurvivalProbability` is the only on in
    :mod:`MDAnalysis.analysis.waterdynamics` to support an `exclusive`
    behaviour (i.e. similar to the current behaviour of :class:`AnalysisBase`
    to the `stop` keyword passed to :meth:`SurvivalProbability.run`. Unlike
    other :mod:`MDAnalysis.analysis.waterdynamics` final frame definitions
    which are `inclusive`.


    .. versionadded:: 0.11.0
    .. versionchanged:: 1.0.0
       Using the MDAnalysis.lib.correlations.py to carry out the intermittency
       and autocorrelation calculations.
       Changed `selection` keyword to `select`.
       Removed support for the deprecated `t0`, `tf`, and `dtmax` keywords.
       These should instead be passed to :meth:`SurvivalProbability.run` as
       the `start`, `stop`, and `tau_max` keywords respectively.
       The `stop` keyword as passed to :meth:`SurvivalProbability.run` has now
       changed behaviour and will act in an `exclusive` manner (instead of it's
       previous `inclusive` behaviour),
    .. versionchanged:: 2.7.0
       Updated docs to align with discrete autocorrelation function.
    """

    def __init__(self, universe, select, verbose=False):
        self.universe = universe
        self.selection = select
        self.verbose = verbose

    def run(self, tau_max=20, start=None, stop=None, step=None, residues=False,
            intermittency=0, verbose=False):
        """
        Computes and returns the Survival Probability (SP) timeseries

        Parameters
        ----------
        start : int, optional
            Zero-based index of the first frame to be analysed, Default: None
            (first frame).
        stop : int, optional
            Zero-based index of the last frame to be analysed (exclusive),
            Default: None (last frame).
        step : int, optional
            Jump every `step`-th frame. This is compatible but independant of
            the taus used, and it is good to consider using the  `step` equal
            to `tau_max` to remove the overlap. Note that `step` and `tau_max`
            work consistently with intermittency. Default: None
            (use every frame).
        tau_max : int, optional
            Survival probability is calculated for the range
            1 <= `tau` <= `tau_max`.
        residues : Boolean, optional
            If true, the analysis will be carried out on the residues
            (.resids) rather than on atom (.ids). A single atom is sufficient
            to classify the residue as within the distance.
        intermittency : int, optional
            The maximum number of consecutive frames for which an atom can
            leave but be counted as present if it returns at the next frame.
            An intermittency of `0` is equivalent to a continuous survival
            probability, which does not allow for the leaving and returning of
            atoms. For example, for `intermittency=2`, any given atom may leave
            a region of interest for up to two consecutive frames yet be
            treated as being present at all frames. The default is continuous
            (0).
        verbose : Boolean, optional
            Print the progress to the console.

        Returns
        -------
        tau_timeseries : list
            tau from 1 to `tau_max`. Saved in the field tau_timeseries.
        sp_timeseries : list
            survival probability for each value of `tau`. Saved in the field
            sp_timeseries.
        sp_timeseries_data: list
            raw datapoints from which the average is taken (sp_timeseries).
            Time dependancy and distribution can be extracted.


        .. versionchanged:: 1.0.0
           To math other analysis methods, the `stop` keyword is now exclusive
           rather than inclusive.
        """

        start, stop, step = self.universe.trajectory.check_slice_indices(
            start,
            stop,
            step
        )

        if tau_max > (stop - start):
            raise ValueError("Too few frames selected for given tau_max.")

        # preload the frames (atom IDs) to a list of sets
        self._selected_ids = []

        # fixme - to parallise: the section should be rewritten so that this loop only creates a list of indices,
        # on which the parallel _single_frame can be applied.

        # skip frames that will not be used in order to improve performance
        # because AtomGroup.select_atoms is the most expensive part of this calculation
        # Example: step 5 and tau 2: LLLSS LLLSS, ... where L = Load, and S = Skip
        # Intermittency means that we have to load the extra frames to know if the atom is actually missing.
        # Say step=5 and tau=1, intermittency=0: LLSSS LLSSS
        # Say step=5 and tau=1, intermittency=1: LLLSL LLLSL
        frame_loaded_counter = 0
        # only for the first window (frames before t are not used)
        frames_per_window = tau_max + 1 + intermittency
        # This number will apply after the first windows was loaded
        frames_per_window_subsequent = (tau_max + 1) + (2 * intermittency)
        num_frames_to_skip = max(step - frames_per_window_subsequent, 0)

        frame_no = start
        while frame_no < stop:      # we have already added 1 to stop, therefore <
            if num_frames_to_skip != 0 and frame_loaded_counter == frames_per_window:
                logger.info("Skipping the next %d frames:", num_frames_to_skip)
                frame_no += num_frames_to_skip
                frame_loaded_counter = 0
                # Correct the number of frames to be loaded after the first window (which starts at t=0, and
                # intermittency does not apply to the frames before)
                frames_per_window = frames_per_window_subsequent
                continue

            # update the frame number
            self.universe.trajectory[frame_no]

            logger.info("Loading frame: %d", self.universe.trajectory.frame)
            atoms = self.universe.select_atoms(self.selection)

            # SP of residues or of atoms
            ids = atoms.residues.resids if residues else atoms.ids
            self._selected_ids.append(set(ids))

            frame_no += 1
            frame_loaded_counter += 1

        # adjust for the frames that were not loaded (step>tau_max + 1),
        # and for extra frames that were loaded (intermittency)
        window_jump = step - num_frames_to_skip

        self._intermittent_selected_ids = correct_intermittency(self._selected_ids, intermittency=intermittency)
        tau_timeseries, sp_timeseries, sp_timeseries_data = autocorrelation(self._intermittent_selected_ids,
                                                                            tau_max, window_jump)

        # warn the user if the NaN are found
        if all(np.isnan(sp_timeseries[1:])):
            logger.warning('NaN Error: Most likely data was not found. Check your atom selections. ')

        # user can investigate the distribution and sample size
        self.sp_timeseries_data = sp_timeseries_data

        self.tau_timeseries = tau_timeseries
        self.sp_timeseries = sp_timeseries
        return self
