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
:Year: 2019
:Copyright: GNU Public License v2

A tool to compute the average dipole moment :math:`M` and
the static dielectric constant

.. math::

   \varepsilon = 1 + \frac{\langle M^2 \rangle - \langle M \rangle^2}
                            {3 \varepsilon_ 0 V k_B T}

for a system simulated in tin foil boundary conditions.
"""

from __future__ import absolute_import, division

import numpy as np
from numpy.testing import assert_almost_equal

from MDAnalysis.units import constants, convert
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.lib.mdamath import make_whole
from MDAnalysis.due import due, Doi
from MDAnalysis.exceptions import NoDataError

due.cite(Doi("10.1080/00268978300102721"),
         description="Dielectric analysis",
         path="MDAnalysis.analysis.dielectric",
         cite_module=True)
del Doi

class DielectricConstant(AnalysisBase):
    r"""Dielectric constant

    Parameters
    ----------
    selection : AtomGroup
          any atomgroup
    temperature : float
        Temperature (Kelvin) at which the system has been simulated [298.15]
    make_whole : bool
        Make molecules whole; If the input already contains whole molecules
        this can be disabled to gain speedup; the default is ``True``
    verbose : bool
          Show detailed progress of the calculation if set to ``True``; the
          default is ``False``

    Attributes
    ----------
    results : dict
          * M: Directional dependant dipole moment
            :math:`\langle \boldsymbol M \rangle` in :math:`eÅ`.
          * M2: Directional dependant squared dipole moment
            :math:`\langle \boldsymbol M^2 \rangle` in :math:`(eÅ)^2`
          * fluct: Directional dependant dipole moment fluctuation
            :math:`\langle \boldsymbol M^2 \rangle - \langle \boldsymbol M \rangle^2`
            in :math:`(eÅ)^2`
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


    .. versionadded:: 0.20.0

    """
    def __init__(self, selection, temperature=300, make_whole=True, **kwargs):
        super(DielectricConstant, self).__init__(selection.universe.trajectory,
                                                 **kwargs)
        self.selection = selection
        self.temperature = temperature
        self.make_whole = make_whole
        self.volume = 0
        
        try:
            self.charges = selection.charges
        except:
            raise NoDataError("No charges defined for selection.")
            
        if np.sum(self.charges) != 0:
            raise NotImplementedError("Analysis for non neutral systems not available!")

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
        if self.make_whole:
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

        self.results["eps"] = self.results["fluct"] / (
                convert(constants["Boltzman_constant"], "kJ/mol", "eV") *
                self.temperature * self.volume * constants["electric_constant"])

        self.results["eps_mean"] = self.results["eps"].mean()

        self.results["eps"] += 1
        self.results["eps_mean"] += 1
