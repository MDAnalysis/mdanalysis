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

"""OpenMM structure I/O --- :mod:`MDAnalysis.converters.OpenMM`
================================================================


Read coordinates data from a
`OpenMM <http://docs.openmm.org/latest/api-python/generated/openmm.app.simulation.Simulation.html#openmm.app.simulation.Simulation>`_
:class:`openmm.app.simulation.Simulation` with :class:`OpenMMReader`
into a MDAnalysis Universe.

Also converts other objects within the
`OpenMM Application Layer <http://docs.openmm.org/latest/api-python/app.html>`_:

- `openmm.app.pdbfile.PDBFile <http://docs.openmm.org/latest/api-python/generated/openmm.app.pdbfile.PDBFile.html#openmm.app.pdbfile.PDBFile>`_
- `openmm.app.modeller.Modeller <http://docs.openmm.org/latest/api-python/generated/openmm.app.modeller.Modeller.html#openmm.app.modeller.Modeller>`_
- `openmm.app.pdbxfile.PDBxFile <http://docs.openmm.org/latest/api-python/generated/openmm.app.pdbxfile.PDBxFile.html#openmm.app.pdbxfile.PDBxFile>`_

Example
-------
OpenMM can read various file formats into OpenMM objects.
MDAnalysis can then convert some of these OpenMM objects into MDAnalysis Universe objects.

    >>> import openmm.app as app
    >>> import MDAnalysis as mda
    >>> from MDAnalysis.tests.datafiles import PDBX
    >>> pdbxfile = app.PDBxFile(PDBX)
    >>> mda.Universe(pdbxfile)
    <Universe with 60 atoms>



Classes
-------

.. autoclass:: OpenMMSimulationReader
   :members:

.. autoclass:: OpenMMAppReader
   :members:


"""


import numpy as np

from ..coordinates import base


class OpenMMSimulationReader(base.SingleFrameReaderBase):
    """Reader for OpenMM Simulation objects


    .. versionadded:: 2.0.0
    """

    format = "OPENMMSIMULATION"
    units = {"time": "ps", "length": "nm", "velocity": "nm/ps",
             "force": "kJ/(mol*nm)", "energy": "kJ/mol"}

    @staticmethod
    def _format_hint(thing):
        """Can this reader read *thing*?
        """
        try:
            from openmm.app import Simulation
        except ImportError:
            try:  # pragma: no cover
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
            self.ts.dimensions[3:] = _sanitize_box_angles(self.ts.dimensions[3:])
            self.convert_velocities_from_native(self.ts._velocities)
            self.convert_forces_from_native(self.ts._forces)
            self.convert_time_from_native(self.ts.dt)

    def _mda_timestep_from_omm_context(self):
        """ Construct Timestep object from OpenMM context """
        try:
            import openmm.unit as u
        except ImportError:  # pragma: no cover
            import simtk.unit as u

        state = self.filename.context.getState(-1, getVelocities=True,
                getForces=True, getEnergy=True)

        n_atoms = self.filename.context.getSystem().getNumParticles()

        ts = self._Timestep(n_atoms, **self._ts_kwargs)
        ts.frame = 0
        ts.data["time"] = state.getTime()._value
        ts.data["potential_energy"] = (
            state.getPotentialEnergy().in_units_of(u.kilojoule/u.mole)._value
        )
        ts.data["kinetic_energy"] = (
            state.getKineticEnergy().in_units_of(u.kilojoule/u.mole)._value
        )
        ts.triclinic_dimensions = state.getPeriodicBoxVectors(
                asNumpy=True)._value
        ts.dimensions[3:] = _sanitize_box_angles(ts.dimensions[3:])
        ts.positions = state.getPositions(asNumpy=True)._value
        ts.velocities = state.getVelocities(asNumpy=True)._value
        ts.forces = state.getForces(asNumpy=True)._value

        return ts


class OpenMMAppReader(base.SingleFrameReaderBase):
    """Reader for OpenMM Application layer objects

    See also `the object definition in the OpenMM Application layer <http://docs.openmm.org/latest/api-python/generated/openmm.app.simulation.Simulation.html#openmm.app.simulation.Simulation>`_

    .. versionadded:: 2.0.0
    """

    format = "OPENMMAPP"
    units = {"time": "ps", "length": "nm"}

    @staticmethod
    def _format_hint(thing):
        """Can this reader read *thing*?
        """
        try:
            from openmm import app
        except ImportError:
            try:  # pragma: no cover
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
            if self.ts.dimensions is not None:
                self.ts.triclinic_dimensions = self.convert_pos_from_native(
                    self.ts.triclinic_dimensions, inplace=False
                )
                self.ts.dimensions[3:] = _sanitize_box_angles(self.ts.dimensions[3:])

    def _mda_timestep_from_omm_app(self):
        """ Construct Timestep object from OpenMM Application object """

        omm_object = self.filename
        n_atoms = omm_object.topology.getNumAtoms()

        ts = self._Timestep(n_atoms, **self._ts_kwargs)
        ts.frame = 0
        if omm_object.topology.getPeriodicBoxVectors() is not None:
            ts.triclinic_dimensions = np.array(
                omm_object.topology.getPeriodicBoxVectors()._value
            )
            ts.dimensions[3:] = _sanitize_box_angles(ts.dimensions[3:])
        ts.positions = np.array(omm_object.getPositions()._value)

        return ts


def _sanitize_box_angles(angles):
    """ Ensure box angles correspond to first quadrant

    See `discussion on unitcell angles <https://github.com/MDAnalysis/mdanalysis/pull/2917/files#r620558575>`_
    """
    inverted = 180 - angles

    return np.min(np.array([angles, inverted]), axis=0)
