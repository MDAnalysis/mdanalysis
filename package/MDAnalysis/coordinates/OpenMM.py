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

"""OpenMM structure I/O --- :mod:`MDAnalysis.coordinates.OpenMM`
================================================================


Read coordinates data from a 
`OpenMM <http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.simulation.Simulation.html#simtk.openmm.app.simulation.Simulation>`_ 
:class:`simtk.openmm.app.simulation.Simulation` with :class:`OpenMMReader` 
into a MDAnalysis Universe. 

Also converts other objects within the 
`OpenMM Application Layer <http://docs.openmm.org/latest/api-python/app.html>`_:

    - `simtk.openmm.app.pdbfile.PDBFile <http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.pdbfile.PDBFile.html#simtk.openmm.app.pdbfile.PDBFile>`_ 
    - `simtk.openmm.app.modeller.Modeller <http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.modeller.Modeller.html#simtk.openmm.app.modeller.Modeller>`_
    - `simtk.openmm.app.pdbxfile.PDBxFile <http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.pdbxfile.PDBxFile.html#simtk.openmm.app.pdbxfile.PDBxFile>`_


Classes
-------

.. autoclass:: OpenMMSimulationReader
   :members:

.. autoclass:: OpenMMAppReader
   :members:


"""


import numpy as np

from . import base


class OpenMMSimulationReader(base.SingleFrameReaderBase):
    """Reader for OpenMM Simulation objects


    .. versionadded:: 2.0.0
    """

    format = "OPENMMSIMULATION"
    units = {"time": "ps", "length": "nm", "velocity": "nm/ps", 
            "force": "kJ/(mol*nm)"}

    @staticmethod
    def _format_hint(thing):
        """Can this reader read *thing*?
        """
        try:
            from simtk.openmm.app import Simulation
        except ImportError:
            return False
        else:
            return isinstance(thing, Simulation)

    def _read_first_frame(self):
        self.n_atoms = self.filename.topology.getNumAtoms()

        self.ts = self._mda_timestep_from_omm_context()

        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)
            self.ts.triclinic_dimensions = self.convert_pos_from_native(
                self.ts.triclinic_dimensions, inplace=False
            )
            self.convert_velocities_from_native(self.ts._velocities)
            self.convert_forces_from_native(self.ts._forces)
            self.convert_time_from_native(self.ts.dt)

    def _mda_timestep_from_omm_context(self):
        """ Construct Timestep object from Openmm context """

        state = self.filename.context.getState(-1, getVelocities=True, 
                getForces=True, getEnergy=True)

        n_atoms = self.filename.context.getSystem().getNumParticles()

        ts = self._Timestep(n_atoms, **self._ts_kwargs)
        ts.frame = 0
        ts.data["time"] = state.getTime()._value
        ts.data["potential_energy"] = state.getPotentialEnergy()
        ts.data["kinetic_energy"] = state.getKineticEnergy()
        ts.triclinic_dimensions = state.getPeriodicBoxVectors(
                asNumpy=True)._value
        ts.positions = state.getPositions(asNumpy=True)._value
        ts.velocities = state.getVelocities(asNumpy=True)._value
        ts.forces = state.getForces(asNumpy=True)._value

        return ts


class OpenMMAppReader(base.SingleFrameReaderBase):
    """Reader for OpenMM Application layer objects 

    See also `OpenMM Application layer <http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.simulation.Simulation.html#simtk.openmm.app.simulation.Simulation>`_ 

    .. versionadded:: 2.0.0
    """

    format = "OPENMMAPP"
    units = {"time": "ps", "length": "nm"}

    @staticmethod
    def _format_hint(thing):
        """Can this reader read *thing*?
        """
        try:
            from simtk.openmm import app
        except ImportError:
            return False
        else:
            return isinstance(thing, (app.PDBFile, app.Modeller,
                app.PDBxFile))

    def _read_first_frame(self):
        self.n_atoms = self.filename.topology.getNumAtoms()

        self.ts = self._mda_timestep_from_omm_app()

        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)
            self.ts.triclinic_dimensions = self.convert_pos_from_native(
                self.ts.triclinic_dimensions, inplace=False
            )

    def _mda_timestep_from_omm_app(self):
        """ Construct Timestep object from Openmm Application object """

        omm_object = self.filename
        n_atoms = omm_object.topology.getNumAtoms()

        ts = self._Timestep(n_atoms, **self._ts_kwargs)
        ts.frame = 0
        if omm_object.topology.getPeriodicBoxVectors() is not None:
            ts.triclinic_dimensions = np.array(
                omm_object.topology.getPeriodicBoxVectors()._value
            )
        ts.positions = np.array(omm_object.getPositions()._value)

        return ts


