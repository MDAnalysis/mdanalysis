# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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
"""\
Hydrogen bond autocorrelation --- :mod:`MDAnalysis.analysis.hydrogenbonds.hbond_autocorrel`
===========================================================================================

:Author: Richard J. Gowers
:Year: 2014
:Copyright: Lesser GNU Public License v2.1+

.. versionadded:: 0.9.0

.. versionchanged:: 2.0.0

   Module moved from :mod:`MDAnalysis.analysis.hbonds.hbond_autocorrel` to
   :mod:`MDAnalysis.analysis.hydrogenbonds.hbond_autocorrel`.


Description
-----------

Calculates the time autocorrelation function, :math:`C_x(t)`, for the hydrogen
bonds in the selections passed to it.  The population of hydrogen bonds at a
given startpoint, :math:`t_0`, is evaluated based on geometric criteria and
then the lifetime of these bonds is monitored over time.  Multiple passes
through the trajectory are used to build an average of the behaviour.

.. math::
   C_x(t) = \\left \\langle \\frac{h_{ij}(t_0) h_{ij}(t_0 + t)}{h_{ij}(t_0)^2} \\right\\rangle

The subscript :math:`x` refers to the definition of lifetime being used, either
continuous or intermittent.  The continuous definition measures the time that
a particular hydrogen bond remains continuously attached, whilst the
intermittent definition allows a bond to break and then subsequently reform and
be counted again.  The relevent lifetime, :math:`\\tau_x`, can then be found
via integration of this function

.. math::
   \\tau_x = \\int_0^\\infty C_x(t) dt`

For this, the observed behaviour is fitted to a multi exponential function,
using 2 exponents for the continuous lifetime and 3 for the intermittent
lifetime.

.. math::
    C_x(t) = A_1 \\exp( - t / \\tau_1)
    + A_2 \\exp( - t / \\tau_2)
    [+ A_3 \\exp( - t / \\tau_3)]

Where the final pre expoential factor :math:`A_n` is subject to the condition:

.. math::
    A_n = 1 - \\sum\\limits_{i=1}^{n-1} A_i

For further details see :footcite:p:`Gowers2015`.

Input
-----

Three AtomGroup selections representing the **hydrogens**, **donors** and
**acceptors** that you wish to analyse.  Note that the **hydrogens** and
**donors** selections must be aligned, that is **hydrogens[0]** and
**donors[0]** must represent a bonded pair.  For systems such as water,
this will mean that each oxygen appears twice in the **donors** AtomGroup.
The function :func:`find_hydrogen_donors` can be used to construct the donor
AtomGroup
::

  import MDAnalysis as mda
  from MDAnalysis.analysis import hydrogenbonds
  from MDAnalysis.tests.datafiles import waterPSF, waterDCD
  u = mda.Universe(waterPSF, waterDCD)
  hydrogens = u.select_atoms('name H*')
  donors = hydrogenbonds.find_hydrogen_donors(hydrogens)


Note that this requires the Universe to have bond information.  If this isn't
present in the topology file, the
:meth:`MDAnalysis.core.groups.AtomGroup.guess_bonds` method can be used
as so
::

  import MDAnalysis as mda
  from MDAnalysis.analysis import hydrogenbonds
  from MDAnalysis.tests.datafiles import GRO
  # we could load the Universe with guess_bonds=True
  # but this would guess **all** bonds
  u = mda.Universe(GRO)
  water = u.select_atoms('resname SOL and not type DUMMY')
  # guess bonds only within our water atoms
  # this adds the bond information directly to the Universe
  water.guess_bonds()
  hydrogens = water.select_atoms('type H')
  # this is now possible as we guessed the bonds
  donors = hydrogenbonds.find_hydrogen_donors(hydrogens)


The keyword **exclusions** allows a tuple of array addresses to be provided,
(Hidx, Aidx),these pairs of hydrogen-acceptor are then not permitted to be
counted as part of the analysis.  This could be used to exclude the
consideration of hydrogen bonds within the same functional group, or to perform
analysis on strictly intermolecular hydrogen bonding.

Hydrogen bonds are defined on the basis of geometric criteria; a
Hydrogen-Acceptor distance of less then **dist_crit** and a
Donor-Hydrogen-Acceptor angle of greater than **angle_crit**.

The length of trajectory to analyse in ps, **sample_time**, is used to choose
what length to analyse.

Multiple passes, controlled by the keyword **nruns**, through the trajectory
are performed and an average calculated.  For each pass, **nsamples** number
of points along the run are calculated.


Output
------

All results of the analysis are available through the *solution* attribute.
This is a dictionary with the following keys

- *results*:  The raw results of the time autocorrelation function.
- *time*:     Time axis, in ps, for the results.
- *fit*:      Results of the exponential curve fitting procedure. For the
             *continuous* lifetime these are (A1, tau1, tau2), for the
             *intermittent* lifetime these are (A1, A2, tau1, tau2, tau3).
- *tau*:      Calculated time constant from the fit.
- *estimate*: Estimated values generated by the calculated fit.

The *results* and *time* values are only filled after the :meth:`run` method,
*fit*, *tau* and *estimate* are filled after the :meth:`solve` method has been
used.


Worked Example for Polyamide
----------------------------

This example finds the continuous hydrogen bond lifetime between N-H..O in a
polyamide system.  This will use the default geometric definition for hydrogen
bonds of length 3.0 Å and angle of 130 degrees.
It will observe a window of 2.0 ps (`sample_time`) and try to gather 1000
sample point within this time window (this relies upon the trajectory being
sampled frequently enough).  This process is repeated for 20 different start
points to build a better average.

::

  import MDAnalysis as mda
  from MDAnalysis.analysis import hydrogenbonds
  from MDAnalysis.tests.datafiles import TRZ_psf, TRZ
  import matplotlib.pyplot as plt
  # load system
  u = mda.Universe(TRZ_psf, TRZ)
  # select atoms of interest into AtomGroups
  H = u.select_atoms('name Hn')
  N = u.select_atoms('name N')
  O = u.select_atoms('name O')
  # create analysis object
  hb_ac = hydrogenbonds.HydrogenBondAutoCorrel(u,
              acceptors=O, hydrogens=H, donors=N,
              bond_type='continuous',
              sample_time=2.0, nsamples=1000, nruns=20)
  # call run to gather results
  hb_ac.run()
  # attempt to fit results to exponential equation
  hb_ac.solve()
  # grab results from inside object
  tau = hb_ac.solution['tau']
  time = hb_ac.solution['time']
  results = hb_ac.solution['results']
  estimate = hb_ac.solution['estimate']
  # plot to check!
  plt.plot(time, results, 'ro')
  plt.plot(time, estimate)
  plt.show()


Functions and Classes
---------------------

.. autofunction:: find_hydrogen_donors

.. autoclass:: HydrogenBondAutoCorrel
   :members:

"""
import numpy as np
import scipy.optimize
import warnings

from MDAnalysis.lib._cutil import _in2d
from MDAnalysis.lib.log import ProgressBar
from MDAnalysis.lib.distances import capped_distance, calc_angles, calc_bonds
from MDAnalysis.core.groups import requires

from MDAnalysis.due import due, Doi

due.cite(
    Doi("10.1063/1.4922445"),
    description="Hydrogen bonding autocorrelation time",
    path="MDAnalysis.analysis.hydrogenbonds.hbond_autocorrel",
)
del Doi


@requires("bonds")
def find_hydrogen_donors(hydrogens):
    """Returns the donor atom for each hydrogen

    Parameters
    ----------
    hydrogens : AtomGroup
      the hydrogens that will form hydrogen bonds

    Returns
    -------
    donors : AtomGroup
      the donor atom for each hydrogen, found via bond information


    .. versionadded:: 0.20.0
    """
    return sum(h.bonded_atoms[0] for h in hydrogens)


class HydrogenBondAutoCorrel(object):
    """Perform a time autocorrelation of the hydrogen bonds in the system.

    Parameters
    ----------
    universe : Universe
        MDAnalysis Universe that all selections belong to
    hydrogens : AtomGroup
        AtomGroup of Hydrogens which can form hydrogen bonds
    acceptors : AtomGroup
        AtomGroup of all Acceptor atoms
    donors : AtomGroup
        The atoms which are connected to the hydrogens.  This group
        must be identical in length to the hydrogen group and matched,
        ie hydrogens[0] is bonded to donors[0].
        For water, this will mean a donor appears twice in this
        group, once for each hydrogen.
    bond_type : str
        Which definition of hydrogen bond lifetime to consider, either
        'continuous' or 'intermittent'.
    exclusions : ndarray, optional
        Indices of Hydrogen-Acceptor pairs to be excluded.
        With nH and nA Hydrogens and Acceptors, a (nH x nA) array of distances
        is calculated, *exclusions* is used as a mask on this array to exclude
        some pairs.
    angle_crit : float, optional
        The angle (in degrees) which all bonds must be greater than [130.0]
    dist_crit : float, optional
        The maximum distance (in Angstroms) for a hydrogen bond [3.0]
    sample_time : float, optional
        The amount of time, in ps, that you wish to observe hydrogen
        bonds for [100]
    nruns : int, optional
        The number of different start points within the trajectory
        to use [1]
    nsamples : int, optional
        Within each run, the number of frames to analyse [50]
    pbc : bool, optional
        Whether to consider periodic boundaries in calculations [``True``]

    ..versionchanged: 1.0.0
      ``save_results()`` method was removed. You can instead use ``np.savez()``
      on :attr:`HydrogenBondAutoCorrel.solution['time']` and
      :attr:`HydrogenBondAutoCorrel.solution['results']` instead.
    """

    def __init__(
        self,
        universe,
        hydrogens=None,
        acceptors=None,
        donors=None,
        bond_type=None,
        exclusions=None,
        angle_crit=130.0,
        dist_crit=3.0,  # geometric criteria
        sample_time=100,  # expected length of the decay in ps
        time_cut=None,  # cutoff time for intermittent hbonds
        nruns=1,  # number of times to iterate through the trajectory
        nsamples=50,  # number of different points to sample in a run
        pbc=True,
    ):

        # warnings.warn("This class is deprecated, use analysis.hbonds.HydrogenBondAnalysis "
        #              "which has .autocorrelation function",
        #              category=DeprecationWarning)

        self.u = universe
        # check that slicing is possible
        try:
            self.u.trajectory[0]
        except Exception:
            raise ValueError("Trajectory must support slicing") from None

        self.h = hydrogens
        self.a = acceptors
        self.d = donors
        if not len(self.h) == len(self.d):
            raise ValueError(
                "Donors and Hydrogen groups must be identical "
                "length.  Try using `find_hydrogen_donors`."
            )

        if exclusions is not None:
            if len(exclusions[0]) != len(exclusions[1]):
                raise ValueError(
                    "'exclusion' must be two arrays of identical length"
                )
            self.exclusions = np.column_stack(
                (exclusions[0], exclusions[1])
            ).astype(np.intp)
        else:
            self.exclusions = None

        self.bond_type = bond_type
        if self.bond_type not in ["continuous", "intermittent"]:
            raise ValueError(
                "bond_type must be either 'continuous' or 'intermittent'"
            )

        self.a_crit = np.deg2rad(angle_crit)
        self.d_crit = dist_crit
        self.pbc = pbc
        self.sample_time = sample_time
        self.nruns = nruns
        self.nsamples = nsamples
        self._slice_traj(sample_time)
        self.time_cut = time_cut

        self.solution = {
            "results": None,  # Raw results
            "time": None,  # Time axis of raw results
            "fit": None,  # coefficients for fit
            "tau": None,  # integral of exponential fit
            "estimate": None,  # y values of fit against time
        }

    def _slice_traj(self, sample_time):
        """Set up start and end points in the trajectory for the
        different passes
        """
        dt = self.u.trajectory.dt  # frame step size in time
        req_frames = int(sample_time / dt)  # the number of frames required

        n_frames = len(self.u.trajectory)
        if req_frames > n_frames:
            warnings.warn(
                "Number of required frames ({}) greater than the"
                " number of frames in trajectory ({})".format(
                    req_frames, n_frames
                ),
                RuntimeWarning,
            )

        numruns = self.nruns
        if numruns > n_frames:
            numruns = n_frames
            warnings.warn(
                "Number of runs ({}) greater than the number of"
                " frames in trajectory ({})".format(self.nruns, n_frames),
                RuntimeWarning,
            )

        self._starts = np.arange(0, n_frames, n_frames / numruns, dtype=int)
        # limit stop points using clip
        self._stops = np.clip(self._starts + req_frames, 0, n_frames)

        self._skip = req_frames // self.nsamples
        if self._skip == 0:  # If nsamples > req_frames
            warnings.warn(
                "Desired number of sample points too high, using {0}".format(
                    req_frames
                ),
                RuntimeWarning,
            )
            self._skip = 1

    def run(self, force=False):
        """Run all the required passes

        Parameters
        ----------
        force : bool, optional
            Will overwrite previous results if they exist
        """
        # if results exist, don't waste any time
        if self.solution["results"] is not None and not force:
            return

        main_results = np.zeros_like(
            np.arange(self._starts[0], self._stops[0], self._skip),
            dtype=np.float32,
        )
        # for normalising later
        counter = np.zeros_like(main_results, dtype=np.float32)

        for i, (start, stop) in ProgressBar(
            enumerate(zip(self._starts, self._stops)),
            total=self.nruns,
            desc="Performing run",
        ):

            # needed else trj seek thinks a np.int64 isn't an int?
            results = self._single_run(int(start), int(stop))

            nresults = len(results)
            if nresults == len(main_results):
                main_results += results
                counter += 1.0
            else:
                main_results[:nresults] += results
                counter[:nresults] += 1.0

        main_results /= counter

        self.solution["time"] = (
            np.arange(len(main_results), dtype=np.float32)
            * self.u.trajectory.dt
            * self._skip
        )
        self.solution["results"] = main_results

    def _single_run(self, start, stop):
        """Perform a single pass of the trajectory"""
        self.u.trajectory[start]

        # Calculate partners at t=0
        box = self.u.dimensions if self.pbc else None

        # 2d array of all distances
        pair = capped_distance(
            self.h.positions,
            self.a.positions,
            max_cutoff=self.d_crit,
            box=box,
            return_distances=False,
        )
        if self.exclusions is not None:
            pair = pair[~_in2d(pair, self.exclusions)]

        hidx, aidx = np.transpose(pair)

        a = calc_angles(
            self.d.positions[hidx],
            self.h.positions[hidx],
            self.a.positions[aidx],
            box=box,
        )
        # from amongst those, who also satisfiess angle crit
        idx2 = np.where(a > self.a_crit)
        hidx = hidx[idx2]
        aidx = aidx[idx2]

        nbonds = len(hidx)  # number of hbonds at t=0
        results = np.zeros_like(
            np.arange(start, stop, self._skip), dtype=np.float32
        )

        if self.time_cut:
            # counter for time criteria
            count = np.zeros(nbonds, dtype=np.float64)

        for i, ts in enumerate(self.u.trajectory[start : stop : self._skip]):
            box = self.u.dimensions if self.pbc else None

            d = calc_bonds(
                self.h.positions[hidx], self.a.positions[aidx], box=box
            )
            a = calc_angles(
                self.d.positions[hidx],
                self.h.positions[hidx],
                self.a.positions[aidx],
                box=box,
            )

            winners = (d < self.d_crit) & (a > self.a_crit)
            results[i] = winners.sum()

            if self.bond_type == "continuous":
                # Remove losers for continuous definition
                hidx = hidx[np.where(winners)]
                aidx = aidx[np.where(winners)]
            elif self.bond_type == "intermittent":
                if self.time_cut:
                    # Add to counter of where losers are
                    count[~winners] += self._skip * self.u.trajectory.dt
                    count[winners] = 0  # Reset timer for winners

                    # Remove if you've lost too many times
                    # New arrays contain everything but removals
                    hidx = hidx[count < self.time_cut]
                    aidx = aidx[count < self.time_cut]
                    count = count[count < self.time_cut]
                else:
                    pass

            if len(hidx) == 0:  # Once everyone has lost, the fun stops
                break

        results /= nbonds

        return results

    def solve(self, p_guess=None):
        """Fit results to an multi exponential decay and integrate to find
        characteristic time

        Parameters
        ----------
        p_guess : tuple of floats, optional
            Initial guess for the leastsq fit, must match the shape of the
            expected coefficients


        Continuous defition results are fitted to a double exponential with
        :func:`scipy.optimize.leastsq`, intermittent definition are fit to a
        triple exponential.

        The results of this fitting procedure are saved into the *fit*,
        *tau* and *estimate* keywords in the solution dict.

         - *fit* contains the coefficients, (A1, tau1, tau2) or
           (A1, A2, tau1, tau2, tau3)
         - *tau* contains the calculated lifetime in ps for the hydrogen
           bonding
         - *estimate* contains the estimate provided by the fit of the time
           autocorrelation function

        In addition, the output of the :func:`~scipy.optimize.leastsq` function
        is saved into the solution dict

         - *infodict*
         - *mesg*
         - *ier*

        """

        if self.solution["results"] is None:
            raise ValueError(
                "Results have not been generated use, the run method first"
            )

        # Prevents an odd bug with leastsq where it expects
        # double precision data sometimes...
        time = self.solution["time"].astype(np.float64)
        results = self.solution["results"].astype(np.float64)

        def within_bounds(p):
            """Returns True/False if boundary conditions are met or not.
            Uses length of p to detect whether it's handling continuous /
            intermittent

            Boundary conditions are:
             0 < A_x < 1
             sum(A_x) < 1
             0 < tau_x
            """
            if len(p) == 3:
                A1, tau1, tau2 = p
                return (A1 > 0.0) & (A1 < 1.0) & (tau1 > 0.0) & (tau2 > 0.0)
            elif len(p) == 5:
                A1, A2, tau1, tau2, tau3 = p
                return (
                    (A1 > 0.0)
                    & (A1 < 1.0)
                    & (A2 > 0.0)
                    & (A2 < 1.0)
                    & ((A1 + A2) < 1.0)
                    & (tau1 > 0.0)
                    & (tau2 > 0.0)
                    & (tau3 > 0.0)
                )

        def err(p, x, y):
            """Custom residual function, returns real residual if all
            boundaries are met, else returns a large number to trick the
            leastsq algorithm
            """
            if within_bounds(p):
                return y - self._my_solve(x, *p)
            else:
                return np.full_like(y, 100000)

        def double(x, A1, tau1, tau2):
            """Sum of two exponential functions"""
            A2 = 1 - A1
            return A1 * np.exp(-x / tau1) + A2 * np.exp(-x / tau2)

        def triple(x, A1, A2, tau1, tau2, tau3):
            """Sum of three exponential functions"""
            A3 = 1 - (A1 + A2)
            return (
                A1 * np.exp(-x / tau1)
                + A2 * np.exp(-x / tau2)
                + A3 * np.exp(-x / tau3)
            )

        if self.bond_type == "continuous":
            self._my_solve = double

            if p_guess is None:
                p_guess = (0.5, 10 * self.sample_time, self.sample_time)

                p, cov, infodict, mesg, ier = scipy.optimize.leastsq(
                    err, p_guess, args=(time, results), full_output=True
                )
            self.solution["fit"] = p
            A1, tau1, tau2 = p
            A2 = 1 - A1
            self.solution["tau"] = A1 * tau1 + A2 * tau2
        else:
            self._my_solve = triple

            if p_guess is None:
                p_guess = (
                    0.33,
                    0.33,
                    10 * self.sample_time,
                    self.sample_time,
                    0.1 * self.sample_time,
                )

            p, cov, infodict, mesg, ier = scipy.optimize.leastsq(
                err, p_guess, args=(time, results), full_output=True
            )
            self.solution["fit"] = p
            A1, A2, tau1, tau2, tau3 = p
            A3 = 1 - A1 - A2
            self.solution["tau"] = A1 * tau1 + A2 * tau2 + A3 * tau3

        self.solution["infodict"] = infodict
        self.solution["mesg"] = mesg
        self.solution["ier"] = ier

        if ier in [1, 2, 3, 4]:  # solution found if ier is one of these values
            self.solution["estimate"] = self._my_solve(
                self.solution["time"], *p
            )
        else:
            warnings.warn("Solution to results not found", RuntimeWarning)

    def __repr__(self):
        return (
            "<MDAnalysis HydrogenBondAutoCorrel analysis measuring the "
            "{btype} lifetime of {n} different hydrogens>"
            "".format(btype=self.bond_type, n=len(self.h))
        )
