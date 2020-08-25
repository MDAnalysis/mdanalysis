# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
import pytest
from numpy.testing import assert_almost_equal

from MDAnalysisTests.datafiles import CONECT
from MDAnalysisTests.coordinates.reference import (RefAdKSmall,
                                                   RefAdK)

from MDAnalysisTests.coordinates.base import _SingleFrameReader

import MDAnalysis as mda

mm = pytest.importorskip('simtk.openmm')
unit = pytest.importorskip('simtk.unit')
app = pytest.importorskip('simtk.openmm.app')
pmd = pytest.importorskip('parmed')

class TestOpenMMPDBFileReader(_SingleFrameReader):
    __test__ = True
    def setUp(self):
        self.universe = mda.Universe(app.PDBFile(RefAdKSmall.filename))
        self.ref = mda.Universe(RefAdKSmall.filename)
        self.prec = 3

    def test_dimensions(self):
        assert_almost_equal(
            # Angles seem to parse differently when openmm reads the pdb file
            self.universe.trajectory.ts.dimensions[0:3],
            self.ref.trajectory.ts.dimensions[0:3],
            self.prec,
            "OpenMMPDBFileReader failed to get unitcell dimensions from OpenMMPDBFile")
    
    def test_coordinates(self):
        up = self.universe.atoms.positions
        rp = self.ref.atoms.positions
        assert_almost_equal(up, rp, decimal=3)

class TestOpenMMModellerReader(_SingleFrameReader):
    __test__ = True
    def setUp(self):
        pdb_obj = app.PDBFile(RefAdKSmall.filename)
        modeller = app.Modeller(pdb_obj.topology, pdb_obj.positions)
        self.universe = mda.Universe(modeller)
        self.ref = mda.Universe(RefAdKSmall.filename)
        self.prec = 3

    def test_dimensions(self):
        assert_almost_equal(
            # Angles seem to parse differently when openmm reads the pdb file
            self.universe.trajectory.ts.dimensions[0:3],
            self.ref.trajectory.ts.dimensions[0:3],
            self.prec,
            "OpenMMModellerReader failed to get unitcell dimensions from OpenMMModeller")
    
    def test_coordinates(self):
        up = self.universe.atoms.positions
        rp = self.ref.atoms.positions
        assert_almost_equal(up, rp, decimal=3)


class TestOpenMMSimulationReader(_SingleFrameReader):
    __test__ = True

    def setUp(self):
        pdb = app.PDBFile(RefAdKSmall.filename)
        forcefield = app.ForceField('amber99sbildn.xml')
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds)
        integrator = mm.LangevinIntegrator(
                300*unit.kelvin, 1/unit.picoseconds, 2.0*unit.femtoseconds
        )
        sim = app.Simulation(pdb.topology, system, integrator)
        sim.context.setPositions(pdb.positions)
        self.universe = mda.Universe(sim)
        self.ref = mda.Universe(RefAdKSmall.filename)
        self.prec = 3

    def test_dimensions(self):
        assert_almost_equal(
            # Angles seem to parse differently when openmm reads the pdb file
            self.universe.trajectory.ts.dimensions[0:3],
            self.ref.trajectory.ts.dimensions[0:3],
            self.prec,
            "OpenMMSimulationReader failed to get unitcell dimensions from OpenMMSimulation")
    
    def test_coordinates(self):
        up = self.universe.atoms.positions
        rp = self.ref.atoms.positions
        assert_almost_equal(up, rp, decimal=3)

    def test_pickle_singleframe_reader(self): pass

