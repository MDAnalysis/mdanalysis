# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
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
Dielectric --- :mod:`MDAnalysis.analysis.dielectric`
===========================================================

:Authors: Mattia Felice Palermo, Philip Loche
:Year: 2022
:Copyright: Lesser GNU Public License v2.1+
"""

import numpy as np

from MDAnalysis.units import constants, convert
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.due import due, Doi
from MDAnalysis.exceptions import NoDataError

due.cite(
    Doi("10.1080/00268978300102721"),
    description="Dielectric analysis",
    path="MDAnalysis.analysis.dielectric",
    cite_module=True,
)
del Doi


class DielectricConstant(AnalysisBase):
    r"""
    Computes the average dipole moment

    .. math::
        \boldsymbol M = \sum_i q_i \boldsymbol r_i


    where :math:`q_i` is the charge and :math:`\boldsymbol r_i` the position of
    atom :math:`i` in the given :class:`MDAnalysis.core.groups.AtomGroup`.
    Also, the static dielectric constant

    .. math::
        \varepsilon = 1 + \frac{\langle M^2 \rangle - \langle M \rangle^2}
                              {3 \varepsilon_ 0 V k_B T}


    is calculated for a system in tin foil boundary conditions, which is
    the usual case if electrostatics are handled with a Ewald summation
    technique. See [Neumann1983]_ for details on the derivation.

    .. warning::
      Applying this class requires that no free charges, such as ions or
      charged fragments, are present in the simulation.

    Parameters
    ----------
    atomgroup : MDAnalysis.core.groups.AtomGroup
      Atomgroup on which the analysis is executed
    temperature : float
      Temperature (Kelvin) at which the system has been simulated
    make_whole : bool
      Make molecules whole; If the input already contains whole molecules
      this can be disabled to gain speedup
    verbose : bool
      Show detailed progress of the calculation

    Attributes
    ----------
    results.M : numpy.ndarray
      Directional dependant dipole moment
      :math:`\langle \boldsymbol M \rangle` in :math:`eÅ`.
    results.M2 : numpy.ndarray
      Directional dependant squared dipole moment
      :math:`\langle \boldsymbol M^2 \rangle` in :math:`(eÅ)^2`
    results.fluct : float
      Directional dependant dipole moment fluctuation
      :math:`\langle \boldsymbol M^2 \rangle - \langle \boldsymbol M \rangle^2`
      in :math:`(eÅ)^2`
    results.eps : numpy.ndarray
      Directional dependant static dielectric constant
    results.eps_mean : float
      Static dielectric constant

    Example
    -------
    Create a :class:`DielectricConstant` instance by supplying an
    :class:`~MDAnalysis.core.groups.AtomGroup`,
    then use the :meth:`run` method::

      import MDAnalysis as mda
      from MDAnalysis.analysis.dielectric import DielectricConstant
      from MDAnalysisTests.datafiles import PSF_TRICLINIC, DCD_TRICLINIC

      # Load a pure water system
      universe = mda.Universe(PSF_TRICLINIC, DCD_TRICLINIC)

      diel = DielectricConstant(universe.atoms)
      diel.run()
      print(diel.results)

    The static dielectric constant of the provided atomgroup is saved
    within the :class:`~MDAnalysis.analysis.base.Results` attribute.


    .. versionadded:: 2.1.0
    """

    def __init__(self, atomgroup, temperature=300, make_whole=True, **kwargs):
        super(DielectricConstant, self).__init__(
            atomgroup.universe.trajectory, **kwargs
        )
        self.atomgroup = atomgroup
        self.temperature = temperature
        self.make_whole = make_whole

    def _prepare(self):
        if not hasattr(self.atomgroup, "charges"):
            raise NoDataError("No charges defined given atomgroup.")

        if not np.allclose(
            self.atomgroup.total_charge(compound="fragments"), 0.0, atol=1e-5
        ):
            raise NotImplementedError(
                "Analysis for non-neutral systems or"
                " systems with free charges are not"
                " available."
            )

        self.volume = 0
        self.results.M = np.zeros(3)
        self.results.M2 = np.zeros(3)
        self.results.fluct = np.zeros(3)
        self.results.eps = np.zeros(3)
        self.results.eps_mean = 0

    def _single_frame(self):
        if self.make_whole:
            self.atomgroup.unwrap()

        self.volume += self.atomgroup.universe.trajectory.ts.volume

        M = np.dot(self.atomgroup.charges, self.atomgroup.positions)
        self.results.M += M
        self.results.M2 += M * M

    def _conclude(self):
        self.results.M /= self.n_frames
        self.results.M2 /= self.n_frames
        self.volume /= self.n_frames

        self.results.fluct = self.results.M2 - self.results.M * self.results.M

        self.results.eps = self.results.fluct / (
            convert(constants["Boltzmann_constant"], "kJ/mol", "eV")
            * self.temperature
            * self.volume
            * constants["electric_constant"]
        )

        self.results.eps_mean = self.results.eps.mean()

        self.results.eps += 1
        self.results.eps_mean += 1
