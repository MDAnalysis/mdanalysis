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
# This module written by Hugo MacDermott-Opeskin 2020

r"""
Mean Squared Displacement --- :mod:`MDAnalysis.analysis.msd`
==============================================================

This module implements the calculation of Mean Squared Displacmements (MSDs) by the Einstein relation.
MSDs can be used to characterise the speed at which particles move and has its roots
in the study of Brownian motion. For a full explanation of the theory behind MSDs and the subsequent calculation of self diffusivities the reader is directed to [Maginn2019]_.
MSDs can be computed from the following expression, known as "Einstein" formula:

.. math::

   MSD(r_{d}) = \bigg{\langle} \frac{1}{N} \sum_{i=1}^{N} |r_{d}  - r_{d}(t_0)|^2 \bigg{\rangle}_{t_{0}}

Where :math:`N` is the number of equivalent particles the MSD is calculated over, :math:`r` are their coordinates and :math:`d` the desired
dimensionality of the MSD. Note that while the definition of the MSD is universal, there are many practical considerations to computing the MSD
that vary between implementations. In this module, we compute a "windowed" MSD, where the MSD is averaged over all possible lag times :math:`\tau \le \tau_{max}`,
where :math:`\tau_{max}` is the length of the trajectory, thereby maximising the number of samples.

The computation of the MSD in this way can be computationally intensive due to it's :math:`N^2` scaling with respect to :math:`\tau_{max}`. 
An algorithm to compute the MSD with :math:`N log(N)` scaling based on a Fast Fourier Transform is known and can be accessed by setting fft=True [Calandri2011]_.

Computing an MSD
----------------
This example computes a 3D MSD for the movement of 100 particles undergoing a random walk.
Files provided as part of the MDAnalysis test suite are used
(in the variables :data:`~MDAnalysis.tests.datafiles.RANDOM_WALK` and
:data:`~MDAnalysis.tests.datafiles.RANDOM_WALK_TOPO`)

First load all modules and test data

    >>> import MDAnalysis as mda
    >>> import MDAnalysis.analysis.msd as msd
    >>> from MDAnalysis.tests.datafiles import RANDOM_WALK, RANDOM_WALK_TOPO

Given a universe containing trajectory data we can extract the MSD
Analyis by using the class :class:`EinsteinMSD` 

    >>> u = mda.Universe(RANDOM_WALK, RANDOM_WALK_TOPO)
    >>> MSD = msd.EinsteinMSD(u, 'all', msd_type='xyz', fft=True)
    >>> MSD.run()

The MSD can then be accessed as

    >>> msd =  MSD.timeseries

Visual inspection of the MSD is important, so lets take a look at it with a simple plot.

    >>> import matplotlib.pyplot as plt
    >>> nframes = MSD.n_frames
    >>> timestep = 1 # this needs to be the actual time between frames in your trajectory
    >>> lagtimes = np.arange(nframes)*timestep # make the lag time axis
    >>> plt.plot(msd, lagtimes)
    >>> plt.show()

We can see that the MSD is approximately linear, this is a numerical proof of a known theoretical result that the MSD of a random walk is approximately linear with respect to lagtime, where the slope is approximately :math:`2*D`.
A segment of the MSD is required to be linear to accurately determine self-diffusivity.
This linear segment represents the so called "middle" of the MSD plot, where ballistic trajectories at short time-lags are excluded along with poorly averaged data at long time-lags.
We can select the "middle" of the MSD by indexing the MSD and the time-lags. Appropriately linear segments of the MSD can be confirmed with a log-log plot as is often reccomended [Maginn2019]_ where the "middle" segment can be identified as having a slope of 1.

    >>> plt.loglog(msd, lagtimes)
    >>> plt.show()

Now that we have identified what segment of our MSD to analyse, lets compute a self diffusivity.

Computing Self Diffusivity
--------------------------------
Self diffusivity is closely related to the MSD.

.. math::

   D_d = \frac{1}{2d} \lim_{t \to \infty} \frac{d}{dt} MSD(r_{d}) 

From the MSD, self-diffusivities :math:`D` with the desired dimensionality :math:`d` can be computed by fitting the MSD with respect to the lag time to a linear model. 
An example of this is shown below, using the MSD computed in the example above. The segment between :math:`\tau = 20` and :math:`\tau = 80` is used to demonstrate selection of an MSD segment.

    >>> from scipy.stats import linregress as lr
    >>> start_time = 20
    >>> start_index = start_time/timestep
    >>> end_time = 80
    >>> end_index = end_time/timestep
    >>> linear_model = lr(lagtimes[start_index:end_index], msd[start_index:end_index])
    >>> slope = linear_model.slope
    >>> error = linear_model.rvalue
    >>> D = slope * 1/(2*MSD.dim_fac) #dim_fac is 3 as we computed a 3D msd ('xyz')

We have now computed a self diffusivity!


Notes
_____

There are several factors that must be taken into account when setting up and processing trajectories for computation of self-diffusivities.
These include specific instructions around simulation settings, using unwrapped trajectories and maintaining relatively small elapsed time between saved frames. 
Additionally corrections for finite size effects are sometimes employed along with varied means of estimating errors.
The reader is directed to the following review, which describes many of the common pitfalls [Maginn2019]_. There are other ways to compute self diffusivity including from a Green-Kubo integral. At this point in time these methods are beyond the scope of this module


References
----------

.. [Maginn2019] Maginn, E. J.; Messerly, R. A.; Carlson, D. J.; Roe, D. R.; Elliott, J. R. Best Practices for Computing Transport Properties 1. Self-Diffusivity and Viscosity from Equilibrium Molecular Dynamics [Article v1.0]. Living J. Comput. Mol. Sci. 2019, 1 (1).
.. [Calandri2011] Calandrini, V.; Pellegrini, E.; Calligari, P.; Hinsen, K.; Kneller, G. R. NMoldyn-Interfacing Spectroscopic Experiments, Molecular Dynamics Simulations and Models for Time Correlation Functions. Collect. SFN 2011, 12, 201–232.

Classes and Functions
---------------------

.. autoclass:: EinsteinMSD
    :members:

"""

from __future__ import division, absolute_import

import numpy as np
import logging
import tidynamics
from ..due import due, Doi
from .base import AnalysisBase

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
    
    Attributes
    ----------
    dim_fac : int
        Dimensionality :math:`d` of the MSD.
    timeseries : :class:`np.ndarray`
        The averaged MSD with respect to lag-time.
    msd_per_particle : :class:`np.ndarray`
        The MSD of each individual particle with respect to lag-time.
    n_frames : int
        Number of frames in trajectory.
    n_particles : int
        Number of particles MSD was calculated over.

    """

    def __init__(self, u, selection=None, msd_type='xyz', fft=True, **kwargs):
        r"""
        Parameters
        ----------
        u : Universe
            An MDAnalysis :class:`Universe`.
        selection : str
            An MDAnalysis selection string. Defaults to `None`.
        msd_type : {'xyz', 'xy', 'yz', 'xz', 'x', 'y', 'z'}
            Desired dimensions to be included in the MSD. Defaults to 'xyz'.
        fft : bool
            Use a fast FFT based algorithm for computation of the MSD.
            Otherwise, use the naive or "simple" algorithm. Defaults to ``True``.

        """
        self.u = u
        
        super(EinsteinMSD, self).__init__(self.u.trajectory, **kwargs)

        #args
        self.selection = selection
        self.msd_type = msd_type
        self.fft = fft
        
        #local
        self._dim = None
        self._position_array = None
        self._atoms = None
        
        #indexing
        self._n_frames = 0
        self._n_particles = 0

        #result
        self.dim_fac = 0
        self.timeseries = None
        self.msds_by_particle = None

    def _prepare(self):
        self._parse_msd_type()
        self._atoms = self.u.select_atoms(self.selection)
        self._n_frames = len(self.u.trajectory)
        self._n_particles = len(self._atoms)
        self._position_array = np.zeros((self._n_frames, self._n_particles, 3))
        self.msds_by_particle = np.zeros((self._n_frames, self._n_particles))
        # self.timeseries not set here
        
    def _parse_msd_type(self):
        r""" Sets up the desired dimensionality of the MSD.
        
        Parameters
        ----------
        self.msd_type : {'xyz', 'xy', 'yz', 'xz', 'x', 'y', 'z'}
            Dimensions to be included in the MSD.

        Returns
        -------
        self._dim : list
            Array-like used to slice the positions to obtain desired dimensionality
        self.dim_fac : int
            Dimension factor :math:`d` of the MSD.

        """
        keys = {'x':[0], 'y':[1], 'z':[2], 'xy':[0,1], 'xz':[0,2], 'yz':[1,2], 'xyz':[0,1,2]}

        self.msd_type = self.msd_type.lower()

        try:
            self._dim = keys[self.msd_type]
        except KeyError:
            raise ValueError('invalid msd_type: {} specified, please specify one of xyz, xy, xz, yz, x, y, z'.format(self.msd_type)

        self.dim_fac = len(self._dim)

    def _single_frame(self):
        r""" Constructs array of positions for MSD calculation.
        
        Parameters
        ----------
        self.u : :class:`Universe`
            MDAnalysis Universe
        
        Returns
        -------
        self._position_array : :class:`np.ndarray`
            Array of particle positions with respect to time shape = (n_frames, n_particles, 3)
        """

        self._position_array[self._frame_index,:,:] = self._atoms.positions
            
    def _conclude(self):
        if self.fft == True:
            self.timeseries = self._conclude_fft()
        else:
            self.timeseries = self._conclude_simple()

    def _conclude_simple(self): 
        r""" Calculates the MSD via the simple "windowed" algorithm.
        
        Parameters
        ----------
        self._position_array : :class:`np.ndarray`
            Array of particle positions with respect to time shape (n_frames, n_particles, 3).

        Returns
        -------
        self.timeseries : :class:`np.ndarray`
            The MSD as a function of lag-time.

        """
        lagtimes = np.arange(1,self.n_frames)
        self.msds_by_particle[0,:] = np.zeros(self._n_particles) # preset the zero lagtime so we don't have to iterate through
        for n in range(self._n_particles):
            for lag in lagtimes:
                disp = self._position_array[:-lag,n,self._dim if lag else None] - self._position_array[lag:,n,self._dim]
                sqdist = np.square(disp, dtype=np.float64).sum(axis=1, dtype=np.float64) #accumulation in np.float64 required
                self.msds_by_particle[lag,n] = np.mean(sqdist, dtype=np.float64)
        msd = self.msds_by_particle.mean(axis=1, dtype=np.float64)
        return msd

    def _conclude_fft(self): #with FFT, np.float64 bit prescision required.
        r""" Calculates the MSD via the FCA fast correlation algorithm.
        
        Parameters
        ----------
        self._position_array : :class:`np.ndarray`
            Array of particle positions with respect to time shape (n_frames, n_particles, 3).

        Returns
        -------
        self.timeseries : :class:`np.ndarray`
            The MSD as a function of lagtime.

        """
        reshape_positions = self._position_array[:,:,self._dim].astype(np.float64)
        for n in range(self._n_particles):
            self.msds_by_particle[:,n] = tidynamics.msd(reshape_positions[:,n,:])
        msd = self.msds_by_particle.mean(axis=1, dtype=np.float64)
        return msd
