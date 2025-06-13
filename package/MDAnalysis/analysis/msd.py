######## Documentation
#
# Issue: The default '_conclude_simple()' function previously calculated time-averaged mean squared displacement considering linear timespacing between the frames even when a dump file with non-linear timestep gap is provided
#
#   Modified by Sirsha Ganguly
#
# Modification:
# '_conclude_modified()' is the new function added to calculate time-averaged mean squared displacement of non-linear dumps
    # Note: This new function is generalized for non-linear dumps.
    #       Moreover, if you use this function to calculate for linear timedump this gives same result as '_conclude_simple()' which is implemented in the default version.
    #       So, I believe the algorithm in '_conclude_modified()' can potentially be used to calculate time-averaged msd in future.
    #       For now it only executes when the dump is non-linear. [ This is determined by the added function is_linear() ]
    #       It is tested with a .dump file of a Coarse Grained system of single type 
#
# General Algorithm of '_conclude_modified()' to calculate time-averaged msd:
    # It creates a dictionary/hash-map with the key being the time difference between the frames (Delta t) and value being a list containing the frame to frame MSD values [msd1, msd2, .....] for that particular Delta t
        # Dictionary to collect MSDs: {Δt: [msd1, msd2, ...]}
    # The reference frame is changed with each iteration and each time a new Key (i.e, Delta t) is found it is appended in the dictionary/hash-map
    # If a duplicate Delta_t is found it is added to the value (msd list) of that specefic Key in the dictionary/hash-map
    # Lastly in the dictionary/hash-map for each Delta_t the msd values inside the list are averaged
#

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

r"""
Mean Squared Displacement --- :mod:`MDAnalysis.analysis.msd`
==============================================================

:Authors: Hugo MacDermott-Opeskin
:Year: 2020
:Copyright: Lesser GNU Public License v2.1+

This module implements the calculation of Mean Squared Displacements (MSDs)
by the Einstein relation. MSDs can be used to characterize the speed at
which particles move and has its roots in the study of Brownian motion.
For a full explanation of the theory behind MSDs and the subsequent calculation
of self-diffusivities the reader is directed to :footcite:p:`Maginn2019`.
MSDs can be computed from the following expression, known as the
**Einstein formula**:

.. math::

   MSD(r_{d}) = \bigg{\langle} \frac{1}{N} \sum_{i=1}^{N} |r_{d}
    - r_{d}(t_0)|^2 \bigg{\rangle}_{t_{0}}

where :math:`N` is the number of equivalent particles the MSD is calculated
over, :math:`r` are their coordinates and :math:`d` the desired dimensionality
of the MSD. Note that while the definition of the MSD is universal, there are
many practical considerations to computing the MSD that vary between
implementations. In this module, we compute a "windowed" MSD, where the MSD
is averaged over all possible lag-times :math:`\tau \le \tau_{max}`,
where :math:`\tau_{max}` is the length of the trajectory, thereby maximizing
the number of samples.

The computation of the MSD in this way can be computationally intensive due to
its :math:`N^2` scaling with respect to :math:`\tau_{max}`. An algorithm to
compute the MSD with :math:`N log(N)` scaling based on a Fast Fourier
Transform is known and can be accessed by setting ``fft=True`` [Calandri2011]_
[Buyl2018]_. The FFT-based approach requires that the
`tidynamics <https://github.com/pdebuyl-lab/tidynamics>`_ package is
installed; otherwise the code will raise an :exc:`ImportError`.

Please cite [Calandri2011]_ [Buyl2018]_ if you use this module in addition to
the normal MDAnalysis citations.

.. warning::
    To correctly compute the MSD using this analysis module, you must supply
    coordinates in the **unwrapped** convention. That is, when atoms pass
    the periodic boundary, they must not be **wrapped** back into the primary
    simulation cell. MDAnalysis does not currently offer this functionality in
    the ``MDAnalysis.transformations`` API despite having functions with
    similar names. We plan to implement the appropriate transformations in the
    future. In the meantime, various simulation packages provide utilities to
    convert coordinates to the unwrapped convention. In GROMACS for example,
    this can be done using ``gmx trjconv`` with the ``-pbc nojump`` flag.

Computing an MSD
----------------
This example computes a 3D MSD for the movement of 100 particles undergoing a
random walk. Files provided as part of the MDAnalysis test suite are used
(in the variables :data:`~MDAnalysis.tests.datafiles.RANDOM_WALK` and
:data:`~MDAnalysis.tests.datafiles.RANDOM_WALK_TOPO`)

First load all modules and test data

.. code-block:: python

    import MDAnalysis as mda
    import MDAnalysis.analysis.msd as msd
    from MDAnalysis.tests.datafiles import RANDOM_WALK_TOPO, RANDOM_WALK

Given a universe containing trajectory data we can extract the MSD
analysis by using the class :class:`EinsteinMSD`

.. code-block:: python

    u = mda.Universe(RANDOM_WALK_TOPO, RANDOM_WALK)
    MSD = msd.EinsteinMSD(u, select='all', msd_type='xyz', fft=True)
    MSD.run()

The MSD can then be accessed as

.. code-block:: python

    msd =  MSD.results.timeseries

Visual inspection of the MSD is important, so let's take a look at it with a
 simple plot.

.. code-block:: python

    import matplotlib.pyplot as plt
    nframes = MSD.n_frames
    timestep = 1 # this needs to be the actual time between frames
    lagtimes = np.arange(nframes)*timestep # make the lag-time axis
    fig = plt.figure()
    ax = plt.axes()
    # plot the actual MSD
    ax.plot(lagtimes, msd, lc="black", ls="-", label=r'3D random walk')
    exact = lagtimes*6
    # plot the exact result
    ax.plot(lagtimes, exact, lc="black", ls="--", label=r'$y=2 D\tau$')
    plt.show()

This gives us the plot of the MSD with respect to lag-time (:math:`\tau`).
We can see that the MSD is approximately linear with respect to :math:`\tau`.
This is a numerical example of a known theoretical result that the MSD of a
random walk is linear with respect to lag-time, with a slope of :math:`2d`.
In this expression :math:`d` is the dimensionality of the MSD. For our 3D MSD,
this is 3. For comparison we have plotted the line :math:`y=6\tau` to which an
ensemble of 3D random walks should converge.

.. _figure-msd:

.. figure:: /images/msd_demo_plot.png
    :scale: 100 %
    :alt: MSD plot

Note that a segment of the MSD is required to be linear to accurately
determine self-diffusivity. This linear segment represents the so called
"middle" of the MSD plot, where ballistic trajectories at short time-lags are
excluded along with poorly averaged data at long time-lags. We can select the
"middle" of the MSD by indexing the MSD and the time-lags. Appropriately
linear segments of the MSD can be confirmed with a log-log plot as is often
reccomended :footcite:p:`Maginn2019` where the "middle" segment can be identified
as having a slope of 1.

.. code-block:: python

    plt.loglog(lagtimes, msd)
    plt.show()

Now that we have identified what segment of our MSD to analyse, let's compute
a self-diffusivity.

Computing Self-Diffusivity
--------------------------------
Self-diffusivity is closely related to the MSD.

.. math::

   D_d = \frac{1}{2d} \lim_{t \to \infty} \frac{d}{dt} MSD(r_{d})

From the MSD, self-diffusivities :math:`D` with the desired dimensionality
:math:`d` can be computed by fitting the MSD with respect to the lag-time to
a linear model. An example of this is shown below, using the MSD computed in
the example above. The segment between :math:`\tau = 20` and :math:`\tau = 60`
is used to demonstrate selection of a MSD segment.

.. code-block:: python

    from scipy.stats import linregress
    start_time = 20
    start_index = int(start_time/timestep)
    end_time = 60
    linear_model = linregress(lagtimes[start_index:end_index],
    				  		  msd[start_index:end_index])
    slope = linear_model.slope
    error = linear_model.stderr
    # dim_fac is 3 as we computed a 3D msd with 'xyz'
    D = slope * 1/(2*MSD.dim_fac)

We have now computed a self-diffusivity!

Combining Multiple Replicates
--------------------------------
It is common practice to combine replicates when calculating MSDs. An example
of this is shown below using MSD1 and MSD2.

.. code-block:: python

    u1 = mda.Universe(RANDOM_WALK_TOPO, RANDOM_WALK)
    MSD1 = msd.EinsteinMSD(u1, select='all', msd_type='xyz', fft=True)
    MSD1.run()

    u2 = mda.Universe(RANDOM_WALK_TOPO, RANDOM_WALK)
    MSD2 = msd.EinsteinMSD(u2, select='all', msd_type='xyz', fft=True)
    MSD2.run()

    combined_msds = np.concatenate((MSD1.results.msds_by_particle,
                                    MSD2.results.msds_by_particle), axis=1)
    average_msd = np.mean(combined_msds, axis=1)

The same cannot be achieved by concatenating the replicas in a single run as
the jump between the last frame of the first trajectory and frame 0 of the
next trajectory will lead to an artificial inflation of the MSD and hence
any subsequent diffusion coefficient calculated.

Notes
_____

There are several factors that must be taken into account when setting up and
processing trajectories for computation of self-diffusivities.
These include specific instructions around simulation settings, using
unwrapped trajectories and maintaining a relatively small elapsed time between
saved frames. Additionally, corrections for finite size effects are sometimes
employed along with various means of estimating errors
:footcite:p:`Yeh2004,Bulow2020` The reader is directed to the following review,
which describes many of the common pitfalls :footcite:p:`Maginn2019`. There are
other ways to compute self-diffusivity, such as from a Green-Kubo integral. At
this point in time, these methods are beyond the scope of this module.


Note also that computation of MSDs is highly memory intensive. If this is
proving a problem, judicious use of the ``start``, ``stop``, ``step`` keywords
to control which frames are incorporated may be required.

References
----------

.. footbibliography::


Classes
-------

.. autoclass:: EinsteinMSD
    :members:
    :inherited-members:

"""

import numpy as np
import logging
from ..due import due, Doi
from .base import AnalysisBase
from ..core import groups
from tqdm import tqdm

logger = logging.getLogger("MDAnalysis.analysis.msd")

due.cite(
    Doi("10.21105/joss.00877"),
    description="Mean Squared Displacements with tidynamics",
    path="MDAnalysis.analysis.msd",
    cite_module=True,
)
due.cite(
    Doi("10.1051/sfn/201112010"),
    description="FCA fast correlation algorithm",
    path="MDAnalysis.analysis.msd",
    cite_module=True,
)
del Doi



class EinsteinMSD(AnalysisBase):

    @staticmethod
    def is_linear(sequence): # This function is to check if the timesteps in the array is linear or not! If it's linear default MDA code will carry-out other wise the modified code carries!
        if len(sequence) < 2:
            return True
        step = sequence[1] - sequence[0]
        return all(sequence[i+1] - sequence[i] == step for i in range(len(sequence) - 1))


    r"""Class to calculate Mean Squared Displacement by the Einstein relation.

    Parameters
    ----------
    u : Universe or AtomGroup
        An MDAnalysis :class:`Universe` or :class:`AtomGroup`.
        Note that :class:`UpdatingAtomGroup` instances are not accepted.
    select : str
        A selection string. Defaults to "all" in which case
        all atoms are selected.
    msd_type : {'xyz', 'xy', 'yz', 'xz', 'x', 'y', 'z'}
        Desired dimensions to be included in the MSD. Defaults to 'xyz'.
    fft : bool
        If ``True``, uses a fast FFT based algorithm for computation of
        the MSD. Otherwise, use the simple "windowed" algorithm.
        The tidynamics package is required for `fft=True`.
        Defaults to ``True``.

    Attributes
    ----------
    dim_fac : int
        Dimensionality :math:`d` of the MSD.
    results.timeseries : :class:`numpy.ndarray`
        The averaged MSD over all the particles with respect to lag-time.
    results.msds_by_particle : :class:`numpy.ndarray`
        The MSD of each individual particle with respect to lag-time.
    ag : :class:`AtomGroup`
        The :class:`AtomGroup` resulting from your selection
    n_frames : int
        Number of frames included in the analysis.
    n_particles : int
        Number of particles MSD was calculated over.


    .. versionadded:: 2.0.0
    """

    def __init__(self, u, select="all", msd_type="xyz", fft=True, **kwargs):
        r"""
        Parameters
        ----------
        u : Universe or AtomGroup
            An MDAnalysis :class:`Universe` or :class:`AtomGroup`.
        select : str
            A selection string. Defaults to "all" in which case
            all atoms are selected.
        msd_type : {'xyz', 'xy', 'yz', 'xz', 'x', 'y', 'z'}
            Desired dimensions to be included in the MSD.
        fft : bool
            If ``True``, uses a fast FFT based algorithm for computation of
            the MSD. Otherwise, use the simple "windowed" algorithm.
            The tidynamics package is required for `fft=True`.
        """
        if isinstance(u, groups.UpdatingAtomGroup):
            raise TypeError(
                "UpdatingAtomGroups are not valid for MSD " "computation"
            )

        super(EinsteinMSD, self).__init__(u.universe.trajectory, **kwargs)

        # args
        self.select = select
        self.msd_type = msd_type
        self._parse_msd_type()
        self.fft = fft

        # local
        self.ag = u.select_atoms(self.select)
        self.n_particles = len(self.ag)
        self._position_array = None

        # result
        self.results.msds_by_particle = None
        self.results.timeseries = None

    def _prepare(self):
        # self.n_frames only available here
        # these need to be zeroed prior to each run() call
        self.results.msds_by_particle = np.zeros(
            (self.n_frames, self.n_particles)
        )
        self._position_array = np.zeros(
            (self.n_frames, self.n_particles, self.dim_fac)
        )
        # self.results.timeseries not set here

    def _parse_msd_type(self):
        r"""Sets up the desired dimensionality of the MSD."""
        keys = {
            "x": [0],
            "y": [1],
            "z": [2],
            "xy": [0, 1],
            "xz": [0, 2],
            "yz": [1, 2],
            "xyz": [0, 1, 2],
        }

        self.msd_type = self.msd_type.lower()

        try:
            self._dim = keys[self.msd_type]
        except KeyError:
            raise ValueError(
                "invalid msd_type: {} specified, please specify one of xyz, "
                "xy, xz, yz, x, y, z".format(self.msd_type)
            )

        self.dim_fac = len(self._dim)

    def _single_frame(self):
        r"""Constructs array of positions for MSD calculation."""
        # shape of position array set here, use span in last dimension
        # from this point on
        self._position_array[self._frame_index] = self.ag.positions[
            :, self._dim
        ]

    def _conclude(self):
        if self.fft:
            self._conclude_fft()
        else:
            if self.is_linear([ts.time for ts in self._trajectory]):
                self._conclude_simple()
            else:
                self._conclude_modified() # New modified code

    def _conclude_simple(self):
        r"""Calculates the MSD via the simple "windowed" algorithm."""
        lagtimes = np.arange(1, self.n_frames)
        positions = self._position_array.astype(np.float64)
        for lag in tqdm(lagtimes):
            disp = positions[:-lag, :, :] - positions[lag:, :, :]
            sqdist = np.square(disp).sum(axis=-1)
            self.results.msds_by_particle[lag, :] = np.mean(sqdist, axis=0)
        self.results.timeseries = self.results.msds_by_particle.mean(axis=1)

    def _conclude_fft(self):  # with FFT, np.float64 bit prescision required.
        r"""Calculates the MSD via the FCA fast correlation algorithm."""
        try:
            import tidynamics
        except ImportError:
            raise ImportError(
                """ERROR --- tidynamics was not found!

                tidynamics is required to compute an FFT based MSD (default)

                try installing it using pip eg:

                    pip install tidynamics

                or set fft=False"""
            )

        positions = self._position_array.astype(np.float64)
        for n in tqdm(range(self.n_particles)):
            self.results.msds_by_particle[:, n] = tidynamics.msd(
                positions[:, n, :]
            )
        self.results.timeseries = self.results.msds_by_particle.mean(axis=1)
    

    def _conclude_modified(self):
        from collections import defaultdict
        print("The Dump file has non-linear time spacing")
        
        dump_times = [ts.time for ts in self._trajectory]
        # print(dump_times)
        n_frames = len(dump_times)
        n_atoms = self.n_particles
        positions = np.zeros((n_frames, n_atoms, 3))
        positions = self._position_array.astype(np.float64)

        msd_dict = defaultdict(list) # Dictionary to collect MSDs: {Δt: [msd1, msd2, ...]}

        # Looping over all the frames as if the referenced gets shifted frame to frame
        for i in range(n_frames):
            for j in range(i + 1, n_frames):
                delta_t = dump_times[j] - dump_times[i]

                # Compute displacement and squared displacement
                disp = positions[j] - positions[i]
                squared_disp = np.sum(disp ** 2, axis=1)
                msd = np.mean(squared_disp)

                # Store MSD under corresponding Δt
                msd_dict[delta_t].append(msd)

        msd_dict[0] = [0]

        # Prepare averaged results
        delta_t_values = []
        avg_msds = []
        for delta_t in sorted(msd_dict.keys()):
            msd_list = msd_dict[delta_t]
            avg_msd = np.mean(msd_list)
            delta_t_values.append(delta_t)
            avg_msds.append(avg_msd)
        
        self.results.timeseries = delta_t_values
        self.results.msds_by_particle = avg_msds

        # self.results.timeseries = self.results.msds_by_particle.mean(axis=1)