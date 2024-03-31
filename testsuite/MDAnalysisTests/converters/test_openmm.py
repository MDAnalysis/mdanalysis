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
import numpy as np
from numpy.testing import assert_allclose

from MDAnalysisTests.datafiles import PDBX
from MDAnalysisTests.coordinates.reference import RefAdKSmall

from MDAnalysisTests.coordinates.base import _SingleFrameReader

import MDAnalysis as mda

try:
    import openmm as mm
    from openmm import unit, app
except ImportError:
    try:
        from simtk import openmm as mm
        from simtk import unit
        from simtk.openmm import app
    except ImportError:
        pytest.skip(allow_module_level=True)


class TestOpenMMBasicSimulationReader():
    @pytest.fixture
    def omm_sim_uni(self):
        system = mm.System()
        topology = app.Topology()
        chain = topology.addChain("CHAIN")
        hydrogen = app.element.Element.getByAtomicNumber(1)
        residue = topology.addResidue("RES", chain)
        for i in range(5):
            system.addParticle(1.0)
            topology.addAtom(hydrogen.symbol, hydrogen, residue)
        positions = np.ones((5, 3)) * unit.angstrom
        integrator = mm.LangevinIntegrator(
            273 * unit.kelvin,
            1.0 / unit.picoseconds, 2.0 * unit.femtoseconds,
        )
        simulation = app.Simulation(topology, system, integrator)
        simulation.context.setPositions(positions)

        return mda.Universe(simulation)

    def test_dimensions(self, omm_sim_uni):
        assert_allclose(
            omm_sim_uni.trajectory.ts.dimensions,
            np.array([20., 20., 20., 90., 90., 90.]),
            rtol=0,
            atol=1e-3,
            err_msg=("OpenMMBasicSimulationReader failed to get unitcell "
                     "dimensions from OpenMM Simulation Object"),
        )

    def test_coordinates(self, omm_sim_uni):
        up = omm_sim_uni.atoms.positions
        reference = np.ones((5, 3))
        assert_allclose(up, reference, rtol=0, atol=1e-3)

    def test_basic_topology(self, omm_sim_uni):
        assert omm_sim_uni.atoms.n_atoms == 5
        assert omm_sim_uni.residues.n_residues == 1
        assert omm_sim_uni.residues.resnames[0] == "RES"
        assert omm_sim_uni.segments.n_segments == 1
        assert omm_sim_uni.segments.segids[0] == '0'
        assert len(omm_sim_uni.bonds.indices) == 0

    def test_data(self, omm_sim_uni):
        data = omm_sim_uni.trajectory.ts.data
        assert isinstance(data["kinetic_energy"], float)
        assert isinstance(data["potential_energy"], float)


class TestOpenMMPDBFileReader(_SingleFrameReader):
    __test__ = True

    def setUp(self):
        self.universe = mda.Universe(app.PDBFile(RefAdKSmall.filename))
        self.ref = mda.Universe(RefAdKSmall.filename)
        self.prec = 3

    def test_dimensions(self):
        assert_allclose(
            self.universe.trajectory.ts.dimensions,
            self.ref.trajectory.ts.dimensions,
            rtol=0,
            atol=1e-3,
            err_msg=("OpenMMPDBFileReader failed to get unitcell dimensions "
                     "from OpenMMPDBFile"),
        )

    def test_coordinates(self):
        up = self.universe.atoms.positions
        rp = self.ref.atoms.positions
        assert_allclose(up, rp, rtol=0, atol=1e-3)


class TestOpenMMModellerReader(_SingleFrameReader):
    __test__ = True

    def setUp(self):
        pdb_obj = app.PDBFile(RefAdKSmall.filename)
        modeller = app.Modeller(pdb_obj.topology, pdb_obj.positions)
        self.universe = mda.Universe(modeller)
        self.ref = mda.Universe(RefAdKSmall.filename)
        self.prec = 3

    def test_dimensions(self):
        assert_allclose(
            self.universe.trajectory.ts.dimensions,
            self.ref.trajectory.ts.dimensions,
            rtol=0,
            atol=1e-3,
            err_msg=("OpenMMModellerReader failed to get unitcell dimensions "
                     "from OpenMMModeller"),
        )

    def test_coordinates(self):
        up = self.universe.atoms.positions
        rp = self.ref.atoms.positions
        assert_allclose(up, rp, rtol=0, atol=1e-3)


class TestOpenMMSimulationReader(_SingleFrameReader):
    __test__ = True

    def setUp(self):
        pdb = app.PDBFile(RefAdKSmall.filename)
        forcefield = app.ForceField("amber99sbildn.xml")
        system = forcefield.createSystem(
            pdb.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds
        )
        integrator = mm.LangevinIntegrator(
            300 * unit.kelvin, 1 / unit.picoseconds, 2.0 * unit.femtoseconds
        )
        sim = app.Simulation(pdb.topology, system, integrator)
        sim.context.setPositions(pdb.positions)
        self.universe = mda.Universe(sim)
        self.ref = mda.Universe(RefAdKSmall.filename)
        self.prec = 3

    def test_dimensions(self):
        assert_allclose(
            self.universe.trajectory.ts.dimensions,
            self.ref.trajectory.ts.dimensions,
            rtol=0,
            atol=1e-3,
            err_msg=("OpenMMSimulationReader failed to get unitcell "
                     "dimensions from OpenMMSimulation"),
        )

    def test_coordinates(self):
        up = self.universe.atoms.positions
        rp = self.ref.atoms.positions
        assert_allclose(up, rp, rtol=0, atol=1e-3)

    @pytest.mark.xfail(reason='OpenMM pickling not supported yet')
    def test_pickle_singleframe_reader(self):
        """
        See `OpenMM SwigPyObject serialisation discussion <https://github.com/MDAnalysis/mdanalysis/issues/2887>`_
        """
        super().test_pickle_singleframe_reader()


@pytest.fixture
def PDBX_U():
    return mda.Universe(app.PDBxFile(PDBX))


def test_pdbx_coordinates(PDBX_U):
    ref_pos = 10 * np.array(
        [
            [0.6537, 3.7168, 1.4751],
            [0.6141, 3.7813, 1.3496],
            [0.4851, 3.8633, 1.3566],
            [0.4428, 3.914, 1.2532],
            [0.6074, 3.6749, 1.2371],
            [0.7366, 3.5948, 1.2334],
            [0.8485, 3.6434, 1.1667],
            [0.7494, 3.477, 1.3064],
            [0.9699, 3.5741, 1.1709],
            [0.8706, 3.4079, 1.3107],
            [0.9803, 3.4575, 1.2437],
            [0.4265, 3.8799, 1.476],
            [0.3019, 3.9561, 1.4964],
            [0.3076, 4.1011, 1.4488],
            [0.2061, 4.1541, 1.4068],
            [0.2633, 3.9555, 1.6453],
            [0.1999, 3.8275, 1.6944],
            [0.1911, 3.7306, 1.6147],
            [0.1619, 3.8224, 1.8133],
            [0.424, 4.1656, 1.4575],
            [0.4388, 4.3071, 1.4222],
            [0.5294, 4.3324, 1.3023],
            [0.5726, 4.4457, 1.2828],
            [0.4853, 4.387, 1.5449],
            [0.3731, 4.4096, 1.6449],
            [0.4257, 4.4404, 1.7833],
            [0.3124, 4.4673, 1.8792],
            [0.5552, 4.2303, 1.2196],
            [0.6382, 4.2504, 1.1005],
            [0.555, 4.3159, 0.9885],
            [0.4442, 4.2696, 0.9581],
            [0.7097, 4.119, 1.0495],
            [0.7974, 4.0487, 1.1586],
            [0.7916, 4.1436, 0.9212],
            [0.9078, 4.133, 1.2278],
            [0.6112, 4.4208, 0.9256],
            [0.5533, 4.4859, 0.8069],
            [0.6548, 4.4763, 0.6937],
            [0.6177, 4.4543, 0.5787],
            [0.5139, 4.6313, 0.8343],
            [0.3789, 4.6433, 0.904],
            [0.7844, 4.4865, 0.7289],
            [0.8992, 4.4777, 0.6384],
            [0.9414, 4.3301, 0.6208],
            [1.0515, 4.2927, 0.6608],
            [1.0151, 4.5643, 0.6923],
            [0.8528, 4.2473, 0.5599],
            [0.8763, 4.1044, 0.5361],
            [1.0028, 4.0746, 0.4537],
            [1.0652, 3.9719, 0.478],
            [0.7534, 4.0355, 0.4762],
            [0.6391, 4.0245, 0.5727],
            [0.5337, 4.1106, 0.5851],
            [0.6204, 3.9232, 0.6731],
            [0.4486, 4.0675, 0.6847],
            [0.5002, 3.9534, 0.7413],
            [0.695, 3.8109, 0.714],
            [0.4534, 3.8758, 0.8486],
            [0.6468, 3.7325, 0.8175],
            [0.5277, 3.7648, 0.8837],
        ]
    )
    rp = PDBX_U.atoms.positions
    assert_allclose(ref_pos, rp, rtol=0, atol=1e-3)
