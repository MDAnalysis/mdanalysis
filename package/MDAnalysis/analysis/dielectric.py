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

A tool to compute the static dielectric constant for a system simulated in 
tin foil boundary conditions.

See Mol. Phys. 1983, vol. 50, No. 4, p. 848 
DOI: 10.1080/00268978300102721
"""
from __future__ import division

import numpy as np

from base import AnalysisBase
#from MDAnalysis.analysis.base import AnalysisBase

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
          Temperature (in Kelvin) at which the system has been simulated [298.15]

    start : int
          The frame to start at [0]
    stop : int
          The frame to end at [-1]
    step : int
          The step size through the trajectory in frames [0]

    Example
    -------
    Create a DielectricConstant object by supplying a selection,
    then use the `run` method:
      diel = DielectricConstant(selection)
      diel.run()

    The static dielectric constant of the system will be returned within the `results`
    variable of the object.

    diel_val = diel.results['dielectric']
    """
    def __init__(self, selection, temperature=298.15, start=None, stop=None, step=None):
        self._ags = [selection]
        self._universe = selection.universe
        self._setup_frames(self._universe.trajectory, start, stop, step)

        self.charges = selection.charges
        self.natoms = len(self._ags[0])
        self.temperature = temperature
        
        self.inst_volume = []
        self.inst_M2V2 = []

    def _single_frame(self):
        positions = self._ags[0].positions

        # Dipole moment of the selection
        dipole = np.add.reduce(np.multiply(positions, self.charges[:, None]))
        volume = self._universe.trajectory.ts.volume

        # 4.8032 converts dipole to debye
        MV = dipole / volume * 4.8032

        # Compute self tensor product of the dipole-volume ratio
        self.inst_M2V2.append(np.outer(MV, MV))
        self.inst_volume.append(volume)

    def _conclude(self):
        avevol = np.array(self.inst_volume).mean()
        aveM2V2 = np.mean(self.inst_M2V2, axis=0)

        # alpha = Debye_to_CA**2/(3.*eps0*boltz*temp*ang_to_m)
        alpha = 30339.2771836305/self.temperature
        diel = 1 + 3 * alpha*np.diagonal(aveM2V2*avevol)
        self.results['dielectric'] = np.sum(diel)/3
