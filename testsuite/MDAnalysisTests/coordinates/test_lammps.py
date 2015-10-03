def test_datareader_VE():
    from MDAnalysis.coordinates.LAMMPS import DATAReader
    assert_raises(ValueError, DATAReader, 'filename')


class _TestLammpsData_Coords(TestCase):
    """Tests using a .data file for loading single frame.

    All topology loading from MDAnalysisTests.data is done in test_topology
    """

    def setUp(self):
        self.u = MDAnalysis.Universe(self.filename)

    def tearDown(self):
        del self.u

    def test_n_atoms(self):
        assert_equal(self.u.atoms.n_atoms, self.n_atoms)

    def test_coords(self):
        assert_equal(self.u.atoms[0].pos, self.pos_atom1)

    def test_velos(self):
        assert_equal(self.u.atoms[0].velocity, self.vel_atom1)

    def test_dimensions(self):
        assert_equal(self.u.dimensions, self.dimensions)

    def test_singleframe(self):
        assert_raises(IOError, self.u.trajectory.next)

    def test_seek(self):
        assert_raises(IndexError, self.u.trajectory.__getitem__, 1)

    def test_seek_2(self):
        ts = self.u.trajectory[0]
        assert_equal(type(ts), MDAnalysis.coordinates.base.Timestep)

    def test_iter(self):
        # Check that iterating works, but only gives a single frame
        assert_equal(len(list(iter(self.u.trajectory))), 1)


class TestLammpsData_Coords(_TestLammpsData_Coords, RefLAMMPSData):
    pass


class TestLammpsDataMini_Coords(_TestLammpsData_Coords, RefLAMMPSDataMini):
    pass
