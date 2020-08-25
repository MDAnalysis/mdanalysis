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

Read coordinates data from a `OpenMM <http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.simulation.Simulation.html#simtk.openmm.app.simulation.Simulation>`_ :class:`simtk.openmm.app.simulation.Simulation` with :class:`OpenMMReader` 
into a MDAnalysis Universe. 

Also converts:

    - `simtk.openmm.app.pdbfile.PDBFile <http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.pdbfile.PDBFile.html#simtk.openmm.app.pdbfile.PDBFile>`_ 
    - `simtk.openmm.app.modeller.Modeller <http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.modeller.Modeller.html#simtk.openmm.app.modeller.Modeller>`_
    - `simtk.openmm.app.pdbxfile.PDBxFile <http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.pdbxfile.PDBxFile.html#simtk.openmm.app.pdbxfile.PDBxFile>`_


Classes
-------

.. autoclass:: OpenMMSimulationReader
   :members:

.. autoclass:: OpenMMPDBFileReader
   :members:

.. autoclass:: OpenMMPDBxFileReader
   :members:

.. autoclass:: OpenMMModellerReader
   :members:


"""


import numpy as np

from . import base


class OpenMMSimulationReader(base.SingleFrameReaderBase):
    """Reader for the OpenMM Simulation objects

    """

    format = "OPENMMSIMULATION"
    units = {"time": "ps", "length": "nm", "velocity": "nm/ps", "force": "kJ/(mol*nm)"}

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

    def _read_first_frame(self):
        self.n_atoms = self.filename.topology.getNumAtoms()

        self.ts = _mda_timestep_from_omm_context(
            self.filename.context, self._Timestep, **self._ts_kwargs
        )

        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)
            self.ts.triclinic_dimensions = self.convert_pos_from_native(
                self.ts.triclinic_dimensions, inplace=False
            )
            self.convert_velocities_from_native(self.ts._velocities)
            self.convert_forces_from_native(self.ts._forces)
            self.convert_time_from_native(self.ts.dt)


class OpenMMPDBFileReader(base.SingleFrameReaderBase):
    """Reader for the OpenMM PDBFile objects

    """

    format = "OPENMMPDBFILE"
    units = {"time": "ps", "length": "nm"}

    @staticmethod
    def _format_hint(thing):
        """Can this reader read *thing*?
        .. versionadded:: 1.0.0
        """
        try:
            from simtk.openmm.app import PDBFile
        except ImportError:
            return False
        else:
            return isinstance(thing, PDBFile)

    def _read_first_frame(self):
        self.n_atoms = self.filename.topology.getNumAtoms()

        self.ts = _mda_timestep_from_omm_modeller_pdbfile(
            self.filename, self._Timestep, **self._ts_kwargs
        )

        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)
            self.ts.triclinic_dimensions = self.convert_pos_from_native(
                self.ts.triclinic_dimensions, inplace=False
            )


class OpenMMPDBxFileReader(base.SingleFrameReaderBase):
    """Reader for the OpenMM PDBxFile objects

    """

    format = "OPENMMPDBXFILE"
    units = {"time": "ps", "length": "nm"}

    @staticmethod
    def _format_hint(thing):
        """Can this reader read *thing*?
        .. versionadded:: 1.0.0
        """
        try:
            from simtk.openmm.app import PDBxFile
        except ImportError:
            return False
        else:
            return isinstance(thing, PDBxFile)

    def _read_first_frame(self):
        self.n_atoms = self.filename.topology.getNumAtoms()

        self.ts = _mda_timestep_from_omm_modeller_pdbfile(
            self.filename, self._Timestep, **self._ts_kwargs
        )

        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)
            self.ts.triclinic_dimensions = self.convert_pos_from_native(
                self.ts.triclinic_dimensions, inplace=False
            )


class OpenMMModellerReader(base.SingleFrameReaderBase):
    """Reader for the OpenMM Modeller objects

    """

    format = "OPENMMMODELLER"
    units = {"time": "ps", "length": "nm"}

    @staticmethod
    def _format_hint(thing):
        """Can this reader read *thing*?
        .. versionadded:: 1.0.0
        """
        try:
            from simtk.openmm.app import Modeller
        except ImportError:
            return False
        else:
            return isinstance(thing, Modeller)

    def _read_first_frame(self):
        self.n_atoms = self.filename.topology.getNumAtoms()

        self.ts = _mda_timestep_from_omm_modeller_pdbfile(
            self.filename, self._Timestep, **self._ts_kwargs
        )

        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)
            self.ts.triclinic_dimensions = self.convert_pos_from_native(
                self.ts.triclinic_dimensions, inplace=False
            )


def _mda_timestep_from_omm_context(omm_context, timestep_module, **ts_kwargs):
    """ Construct Timestep object from Openmm context 

    Parameters
    ----------
    omm_context: simtk.openmm.context
    timestep_module: MDAnalysis.coordinates.base.timestep 
        This is the module, but the object gets created within this function """

    state = omm_context.getState(-1, getVelocities=True, getForces=True, getEnergy=True)

    n_atoms = omm_context.getSystem().getNumParticles()

    ts = timestep_module(n_atoms, **ts_kwargs)
    ts.frame = 0
    ts.data["time"] = state.getTime()._value
    ts.data["potential_energy"] = state.getPotentialEnergy()
    ts.data["kinetic_energy"] = state.getKineticEnergy()
    ts.triclinic_dimensions = state.getPeriodicBoxVectors(asNumpy=True)._value
    ts.positions = state.getPositions(asNumpy=True)._value
    ts.velocities = state.getVelocities(asNumpy=True)._value
    ts.forces = state.getForces(asNumpy=True)._value

    return ts


def _mda_timestep_from_omm_modeller_pdbfile(omm_object, timestep_module, **ts_kwargs):
    """ Construct Timestep object from Openmm context 

    Parameters
    ----------
    omm_object: simtk.openmm.app.PDBFile, simtk.openmm.app.Modeller, or simtk.openmm.app.PDBxFile
    timestep_module: MDAnalysis.coordinates.base.timestep 
        This is the module, but the object gets created within this function """

    n_atoms = omm_object.topology.getNumAtoms()

    ts = timestep_module(n_atoms, **ts_kwargs)
    ts.frame = 0
    if omm_object.topology.getPeriodicBoxVectors() is not None:
        ts.triclinic_dimensions = np.array(
            omm_object.topology.getPeriodicBoxVectors()._value
        )
    ts.positions = np.array(omm_object.getPositions()._value)

    return ts
