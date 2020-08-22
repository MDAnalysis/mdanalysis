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

import numpy as np

from . import base


class OpenMMSimulationReader(base.ReaderBase):
    """Reader for the OpenMM Simulation objects

    """
    format = 'OPENMMSIMULATION'
    units = {'time': 'ps', 'length': 'nm', 'velocity': 'nm/ps',
            'force': 'kJ/(mol*nm)'}

    @staticmethod
    def _format_hint(thing):
        """Can this reader read *thing*?
        .. versionadded:: 1.0.0
        """
        try:
            from simtk.openmm.app import Simulation
        except ImportError:
            return False
        else:
            return isinstance(thing, Simulation) 

    def __init__(self, filename, **kwargs):
        super(OpenMMSimulationReader, self).__init__(filename, **kwargs)
        self.filename = filename
        self.n_atoms = self.filename.system.getNumParticles()
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        self.ts.dt = (self.filename.context.getState(1).getTime()._value - 
            self.filename.context.getState(0).getTime()._value
        )
        print('first {}'.format(self.ts.dt))

        self._frame = 1
        self._read_frame(self._frame)

    @property
    def n_frames(self):
        return self.filename.currentStep

    def _read_frame(self, frame):
        print('called read frame on {}'.format(frame))
        self._frame = frame
        self.ts.frame = self._frame

        state = self.filename.context.getState(self._frame, getVelocities=True, getForces=True, getEnergy=True)
        self.ts.dt = (self.filename.context.getState(self._frame).getTime()._value - 
            self.filename.context.getState(self._frame - 1).getTime()._value
        )

        self.ts.data['time'] = state.getTime()
        self.ts.data['potential_energy'] = state.getPotentialEnergy()
        self.ts.data['kinetic_energy'] = state.getKineticEnergy()
        self.ts.triclinic_dimensions = state.getPeriodicBoxVectors(asNumpy=True)._value
        self.ts.positions = state.getPositions(asNumpy=True)._value
        self.ts.velocities = state.getVelocities(asNumpy=True)._value
        self.ts.forces = state.getForces(asNumpy=True)._value

        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)
            self.ts.triclinic_dimensions = self.convert_pos_from_native(
                self.ts.triclinic_dimensions, inplace=False
            )
            self.convert_velocities_from_native(self.ts._velocities)
            self.convert_forces_from_native(self.ts._forces)
            self.convert_time_from_native(self.ts.dt)

        return self.ts
    
    def _read_next_timestep(self):
        return self._read_frame(self._frame + 1)

    #def _reopen(self):
        #self.ts.frame = 0 
        #self._frame = 0
    

