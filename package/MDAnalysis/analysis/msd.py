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

r"""
Mean Squared Displacement --- :mod:`MDAnalysis.analysis.msd`
==============================================================

:Authors: Hugo MacDermott-Opeskin
:Year: 2020
:Copyright: GNU Public License v2

This module implements the calculation of Mean Squared Displacements (MSDs)
by the Einstein relation. MSDs can be used to characterize the speed at
which particles move and has its roots in the study of Brownian motion.
For a full explanation of the theory behind MSDs and the subsequent
calculation of self-diffusivities the reader is directed to [Maginn2019]_.
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
    from MDAnalysis.tests.datafiles import RANDOM_WALK, RANDOM_WALK_TOPO

Given a universe containing trajectory data we can extract the MSD
analysis by using the class :class:`EinsteinMSD`

.. code-block:: python

    u = mda.Universe(RANDOM_WALK, RANDOM_WALK_TOPO)
    MSD = msd.EinsteinMSD(u, select='all', msd_type='xyz', fft=True)
    MSD.run()

The MSD can then be accessed as

.. code-block:: python

    msd =  MSD.timeseries

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
reccomended [Maginn2019]_ where the "middle" segment can be identified as
having a slope of 1.

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
    error = linear_model.rvalue
    # dim_fac is 3 as we computed a 3D msd with 'xyz'
    D = slope * 1/(2*MSD.dim_fac)

We have now computed a self-diffusivity!


Notes
_____

There are several factors that must be taken into account when setting up and
processing trajectories for computation of self-diffusivities.
These include specific instructions around simulation settings, using
unwrapped trajectories and maintaining a relatively small elapsed time between
saved frames. Additionally, corrections for finite size effects are sometimes
employed along with various means of estimating errors [Yeh2004]_ [Bulow2020]_.
The reader is directed to the following review, which describes many of the
common pitfalls [Maginn2019]_. There are other ways to compute
self-diffusivity, such as from a Green-Kubo integral. At this point in time,
these methods are beyond the scope of this module.


Note also that computation of MSDs is highly memory intensive. If this is
proving a problem, judicious use of the ``start``, ``stop``, ``step`` keywords to control which frames are incorporated may be required.

References
----------

.. [Maginn2019] Maginn, E. J., Messerly, R. A., Carlson, D. J.; Roe, D. R.,
                Elliott, J. R. Best Practices for Computing Transport
                Properties 1. Self-Diffusivity and Viscosity from Equilibrium
                Molecular Dynamics [Article v1.0]. Living J. Comput. Mol. Sci.
                2019, 1 (1).

.. [Yeh2004] Yeh, I. C.; Hummer, G. System-Size Dependence of Diffusion
                Coefficients and Viscosities from Molecular Dynamics
                Simulations with Periodic Boundary Conditions.
                J. Phys. Chem. B 2004, 108 (40), 15873–15879.

.. [Bulow2020] von Bülow, S.; Bullerjahn, J. T.; Hummer, G. Systematic
                Errors in Diffusion Coefficients from Long-Time Molecular
                Dynamics Simulations at Constant Pressure. 2020.
                arXiv:2003.09205 [Cond-Mat, Physics:Physics].



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

logger = logging.getLogger('MDAnalysis.analysis.msd')

due.cite(Doi("10.21105/joss.00877"),
         description="Mean Squared Displacements with tidynamics",
         path="MDAnalysis.analysis.msd",
         cite_module=True)
due.cite(Doi("10.1051/sfn/201112010"),
         description="FCA fast correlation algorithm",
         path="MDAnalysis.analysis.msd",
         cite_module=True)
del Doi


class EinsteinMSD(AnalysisBase):
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
    timeseries : :class:`numpy.ndarray`
        The averaged MSD over all the particles with respect to lag-time.
    msds_by_particle : :class:`numpy.ndarray`
        The MSD of each individual particle with respect to lag-time.
    ag : :class:`AtomGroup`
        The :class:`AtomGroup` resulting from your selection
    n_frames : int
        Number of frames included in the analysis.
    n_particles : int
        Number of particles MSD was calculated over.
    """

    def __init__(self, u, select='all', msd_type='xyz', fft=True, **kwargs):
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
            raise TypeError("UpdatingAtomGroups are not valid for MSD "
                            "computation")

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
        self.msds_by_particle = None
        self.timeseries = None

    def _prepare(self):
        # self.n_frames only available here
        # these need to be zeroed prior to each run() call
        self.msds_by_particle = np.zeros((self.n_frames, self.n_particles))
        self._position_array = np.zeros(
            (self.n_frames, self.n_particles, self.dim_fac))
        # self.timeseries not set here

    def _parse_msd_type(self):
        r""" Sets up the desired dimensionality of the MSD.

        """
        keys = {'x': [0], 'y': [1], 'z': [2], 'xy': [0, 1],
                'xz': [0, 2], 'yz': [1, 2], 'xyz': [0, 1, 2]}

        self.msd_type = self.msd_type.lower()

        try:
            self._dim = keys[self.msd_type]
        except KeyError:
            raise ValueError(
                'invalid msd_type: {} specified, please specify one of xyz, '
                'xy, xz, yz, x, y, z'.format(self.msd_type))

        self.dim_fac = len(self._dim)

    def _single_frame(self):
        r""" Constructs array of positions for MSD calculation.

        """
        # shape of position array set here, use span in last dimension
        # from this point on
        self._position_array[self._frame_index] = (
            self.ag.positions[:, self._dim])

    def _conclude(self):
        if self.fft:
            self._conclude_fft()
        else:
            self._conclude_simple()

    def _conclude_simple(self):
        r""" Calculates the MSD via the simple "windowed" algorithm.

        """
        lagtimes = np.arange(1, self.n_frames)
        positions = self._position_array.astype(np.float64)
        for lag in lagtimes:
            disp = positions[:-lag, :, :] - positions[lag:, :, :]
            sqdist = np.square(disp).sum(axis=-1)
            self.msds_by_particle[lag, :] = np.mean(sqdist, axis=0)
        self.timeseries = self.msds_by_particle.mean(axis=1)

    def _conclude_fft(self):  # with FFT, np.float64 bit prescision required.
        r""" Calculates the MSD via the FCA fast correlation algorithm.

        """
        try:
            import tidynamics
        except ImportError:
            raise ImportError("""ERROR --- tidynamics was not found!

                tidynamics is required to compute an FFT based MSD (default)

                try installing it using pip eg:

                    pip install tidynamics

                or set fft=False""")

        positions = self._position_array.astype(np.float64)
        for n in range(self.n_particles):
            self.msds_by_particle[:, n] = tidynamics.msd(
                positions[:, n, :])
        self.timeseries = self.msds_by_particle.mean(axis=1)
