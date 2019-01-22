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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


r"""
Dielectric --- :mod:`MDAnalysis.analysis.dielectric`
===========================================================

:Authors: Mattia Felice Palermo, Philip Loche
:Year: 2018
:Copyright: GNU Public License v2

A tool to compute the average dipole moment and the static dielectric constant

.. math::

   \varepsilon = 1 + \frac{\langle M^2 \rangle - \langle M \rangle^2}
                            {3 \varepsilon_ 0 V k_B T}

for a system simulated in tin foil boundary conditions. For a distance

References
----------
Neumann, Martin.
Dipole Moment Fluctuation Formulas in Computer Simulations of Polar Systems. 
Molecular Physics 50, no. 4 (November 1983), 841–858, doi:10.1080/00268978300102721
"""

from __future__ import division

import numpy as np
from numpy.testing import assert_almost_equal

from MDAnalysis.units import constants, convert
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.lib.mdamath import make_whole

class DielectricConstant(AnalysisBase):
    r"""Dielectric constant

    Parameters
    ----------
    selection : AtomGroup
          any atomgroup
    start : int
          The frame to start at [0]
    stop : int
          The frame to end at [-1]
    step : int
          The step size through the trajectory in frames [0]
    temperature : float
        Temperature (Kelvin) at which the system has been simulated [298.15]
    verbose : bool (optional)
          Show detailed progress of the calculation if set to ``True``; the
          default is ``False``.
        
    Attributes
    ----------
    results : dict
          * M: Directional dependant dipole moment 
            :math:`\langle \boldsymbol M \rangle`.
          * M2: Directional dependant squared dipole moment 
            :math:`\langle \boldsymbol M^2 \rangle`.
          * fluct: Directional dependant dipole moment fluctuation 
            :math:`\langle \boldsymbol M^2 \rangle - \langle \boldsymbol M \rangle^2`.
          * eps: Directional dependant static dielectric constant
          * eps_mean: Static dielectric constant

    Example
    -------
    Create a DielectricConstant object by supplying a selection,
    then use the :meth:`run` method::
    
      diel = DielectricConstant(selection)
      diel.run()

    The static dielectric constant of the system will be returned
    within the `results` variable of the object.
    
    
    .. versionadded:: 0.19.1

    """
    def __init__(self, selection, temperature=300, **kwargs):
        super(DielectricConstant, self).__init__(selection.universe.trajectory,
                                                 **kwargs)
        self.selection = selection
        self.temperature = temperature
        self.charges = selection.charges
        self.volume = 0
        
        # Dictionary containing results
        self.results = {"M":  np.zeros(3),
                        "M2": np.zeros(3),
                        "fluct": np.zeros(3),
                        "eps": np.zeros(3),
                        "eps_mean": 0}

    def _prepare(self):
        error = "ERROR: Total charge of the selection is not zero."
        assert_almost_equal(0, np.sum(self.charges), err_msg=error)
        
    def _single_frame(self):
        # Make molecules whole
        for frag in self.selection.fragments:
            make_whole(frag)
        
        self.volume += self.selection.universe.trajectory.ts.volume
        
        M = np.dot(self.charges, self.selection.positions)
        self.results["M"] += M
        self.results["M2"] += M * M

    def _conclude(self):
        self.results["M"] /= self.n_frames
        self.results["M2"] /= self.n_frames
        self.volume /= self.n_frames
        
        self.results["fluct"] = self.results["M2"] - self.results["M"] * self.results["M"]
        
        self.results["eps"] = 1 + self.results["fluct"] / (
                convert(constants["Boltzman_constant"], "kJ/mol", "eV") * 
                self.temperature * self.volume * constants["electric_constant"])
                
        self.results["eps_mean"] = 1 + self.results["fluct"].mean() / (
                convert(constants["Boltzman_constant"], "kJ/mol", "eV")  * 
                self.temperature * self.volume * constants["electric_constant"])