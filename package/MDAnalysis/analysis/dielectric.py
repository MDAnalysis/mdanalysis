# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
Dielectric --- :mod:`MDAnalysis.analysis.dielectric`
===========================================================

A tool to compute the trajectory average static dielectric
constant for a system simulated in tin foil boundary conditions.

References
----------
  Mol. Phys. 1983, vol. 50, No. 4, p. 848
  DOI: 10.1080/00268978300102721
"""
from __future__ import division

import numpy as np
from numpy.testing import assert_almost_equal

from MDAnalysis.analysis.base import AnalysisBase

class DielectricConstant(AnalysisBase):
    """Dielectric constant computation
    DielectricConstant(selection, temperature=298.15)

    Parameters
    ----------
    selection : `AtomGroup` object
      Any atomgroup

    Keywords
    --------
    temperature : float
          Temperature (Kelvin) at which the system has been simulated [298.15]

    start : int
          The frame to start at [0]
    stop : int
          The frame to end at [-1]
    step : int
          The step size through the trajectory in frames [0]

    Notes
    -----
    The dielectric constant computation works only for selections
    whose total charge equals to zero.

    Example
    -------
    Create a DielectricConstant object by supplying a selection,
    then use the `run` method:
      diel = DielectricConstant(selection)
      diel.run()

    The static dielectric constant of the system will be returned
    within the `results` variable of the object.

    diel_val = diel.results['dielectric']

    """
    def __init__(self, selection, temperature=298.15,
                 start=None, stop=None, step=None):
        self._ags = [selection]
        self._universe = selection.universe
        self._setup_frames(self._universe.trajectory, start, stop, step)

        self.charges = selection.charges
        self.temperature = temperature

        self.inst_diel = []
        self.results = {}

        error = "ERROR: Total charge of the selection is not zero."
        assert_almost_equal(0, np.sum(self.charges), err_msg=error)

    def _single_frame(self):
        positions = self._ags[0].positions
        volume = self._universe.trajectory.ts.volume

        self.inst_diel.append(compute_dielectric(positions, self.charges,
                                                volume, self.temperature))

    def _conclude(self):
        self.results['dielectric'] = np.array(self.inst_diel).mean()

def compute_dielectric(positions, charges, volume, temperature):
    """ Single frame dielectric constant computation

    This function computes the dielectric constant given atomic position
    and charges, configuration volume and system temperature.

    Parameters
    ----------
    positions : array
        List or numpy array containing the positions of the atoms in the
        form [[x1, y1, z1], [x2, y2, z2], ..., [xn, yn, zn]].

    charges : array
        List or numpy array containing atomic charges in the form
        [q1, q2, ..., 1n].

    volume : float
        Volume of the system in Angstrom**3.

    temperature : float
        Temperature of the simulate system in kelvin

    Returns
    -------
    float
        Value of the isotropic static dielectric constant

    """
    # Dipole moment of the selection
    # 4.8032 converts to debye
    dipole = np.add.reduce(np.multiply(positions, charges[:, None])) * 4.8032
    dipole2 = np.outer(dipole, dipole)

    # alpha = Debye_to_CA**2/(eps0*boltz*temp*ang_to_m)
    alpha = 91017.8315508915/temperature
    diel = 1 + alpha*np.diagonal(dipole2/volume)

    return diel
