import pytest

import numpy as np
from numpy.testing import (assert_allclose,
                           assert_equal)

from gridData import Grid, mrc

from . import datafiles

@pytest.fixture(scope="module")
def g1():
    return Grid(datafiles.CCP4, file_format="MRC")

@pytest.fixture(scope="module")
def g2():
    data = mrc.MRC()
    data.read(datafiles.CCP4)
    grid, edges = data.histogramdd()
    return Grid(grid=grid, edges=edges)

def test_ccp4_Grid(g1):
    _test_ccp4(g1)

def test_ccp4_mrc(g2):
    _test_ccp4(g2)

def _test_ccp4(g):
    POINTS = 192
    assert_equal(g.grid.flat, np.arange(1, POINTS+1))
    assert_equal(g.grid.size, POINTS)
    assert_allclose(g.delta, [3./4, .5, 2./3])
    assert_equal(g.origin, np.zeros(3))



@pytest.fixture(scope="module")
def ccp4data():
    return mrc.MRC(datafiles.CCP4_1JZV)

@pytest.mark.parametrize('name,value', [
    # nx, ny, nz are named nc, nr, ns in the CCP4 module
    ('nx', 96),
    ('ny', 76),
    ('nz', 70),
    ('mode', 2),
    ('nxstart', -4),
    ('nystart', -23),
    ('nzstart', 102),
    ('mx', 84),
    ('my', 84),
    ('mz', 160),
    ('cella', np.rec.array((45.8, 45.8, 89.65),
                        dtype=[('x', '<f4'), ('y', '<f4'), ('z', '<f4')])),
    ('cellb', np.rec.array((90., 90., 90.),
                           dtype=[('alpha', '<f4'), ('beta', '<f4'), ('gamma', '<f4')])),
    ('mapc', 2),
    ('mapr', 1),
    ('maps', 3),
    ('dmin', -0.9930942058563232),
    ('dmax', 9.050403594970703),
    ('dmean', -0.0005801090155728161),
    ('ispg', 92),
    ('nsymbt', 640),
    ('nversion', 0),
    ('origin', np.rec.array((0., 0., 0.),
                            dtype=[('x', '<f4'), ('y', '<f4'), ('z', '<f4')])),
    ('map', b'MAP '),
    ('machst', np.array([68, 65,  0,  0], dtype=np.uint8)),
    ('rms', 0.40349153),
    ('nlabl', 1),
    ('label', np.array([
        b' Map from fft                                                                   ',
        b'                                                                                ',
        b'                                                                                ',
        b'                                                                                ',
        b'                                                                                ',
        b'                                                                                ',
        b'                                                                                ',
        b'                                                                                ',
        b'                                                                                ',
        b'                                                                                '],
                    dtype='|S80'))
])
def test_ccp4_read_header(ccp4data, name, value):
    if type(value) is float:
        assert_allclose(ccp4data.header[name], value, rtol=1e-06)
    else:
        assert_equal(ccp4data.header[name], value)

def test_axes_orientation(ccp4data):
    # correctly interpret mapc, mapr, maps = 2, 1, 3
    # for nx, ny, nz = 96, 76, 70.
    # see also #76
    assert_equal(ccp4data.shape, (76, 96, 70))

def test_delta(ccp4data):
    assert_allclose(ccp4data.delta, np.array(
        [[0.5452381, 0.       , 0.       ],
         [0.       , 0.5452381, 0.       ],
         [0.       , 0.       , 0.5603125]], dtype=np.float32))

def test_origin(ccp4data):
    # shift with nxstart, nystart, nzstart and delta
    #
    # (visual comparison of CCP4 and DX file in ChimeraX at same
    # level shows full agreement)
    assert_allclose(ccp4data.origin, [-12.5404758,  -2.1809523,  57.151876 ])

def test_triclinic_ValueError():
    with pytest.raises(ValueError,
                       match="Only orthorhombic unitcells are currently "
                       "supported, not"):
        Grid(datafiles.MRC_EMD3001, file_format="MRC")

def test_mrcfile_volume_check():
    with pytest.raises(ValueError, match="is not a volumetric density"):
        Grid(datafiles.ISPG_0)

def test_mrcfile_volume_force():
    grid = Grid(datafiles.ISPG_0, assume_volumetric=True)
    assert_allclose(np.sum(grid.grid), 829.925)


class TestGridMRC():
    @pytest.fixture(scope="class")
    def grid(self):
        return Grid(datafiles.CCP4_1JZV)

    def test_shape(self, grid, ccp4data):
        assert_equal(grid.grid.shape, ccp4data.shape)

    def test_mrc_header(self, grid, ccp4data):
        # undocumented MRC header in Grid
        assert grid._mrc_header == ccp4data.header

    def test_delta(self, grid, ccp4data):
        assert_allclose(grid.delta, np.diag(ccp4data.delta))

    def test_origin(self, grid, ccp4data):
        assert_allclose(grid.origin, ccp4data.origin)

    def test_data(self, grid, ccp4data):
        assert_allclose(grid.grid, ccp4data.array)

class TestMRCWrite:
    """Tests for MRC write functionality"""
    
    def test_mrc_write_roundtrip(self, ccp4data, tmpdir):
        """Test writing and reading back preserves data"""
        outfile = str(tmpdir / "roundtrip.mrc")
        
        # Write the file
        ccp4data.write(outfile)
        
        # Read it back
        mrc_read = mrc.MRC(outfile)
        
        # Check data matches
        assert_allclose(mrc_read.array, ccp4data.array)
        assert_allclose(mrc_read.origin, ccp4data.origin, rtol=1e-5, atol=1e-3)
        assert_allclose(mrc_read.delta, ccp4data.delta, rtol=1e-5)
    
    def test_mrc_write_header_preserved(self, ccp4data, tmpdir):
        """Test that header fields are preserved"""
        outfile = str(tmpdir / "header.mrc")
        
        ccp4data.write(outfile)
        mrc_read = mrc.MRC(outfile)
        
        # Check axis ordering preserved
        assert mrc_read.header.mapc == ccp4data.header.mapc
        assert mrc_read.header.mapr == ccp4data.header.mapr
        assert mrc_read.header.maps == ccp4data.header.maps
        
        # Check offsets preserved
        assert mrc_read.header.nxstart == ccp4data.header.nxstart
        assert mrc_read.header.nystart == ccp4data.header.nystart
        assert mrc_read.header.nzstart == ccp4data.header.nzstart
    
    def test_mrc_write_new_file(self, tmpdir):
        """Test creating new MRC file from scratch"""
        outfile = str(tmpdir / "new.mrc")
        
        # Create new MRC object
        mrc_new = mrc.MRC()
        mrc_new.array = np.arange(24).reshape(2, 3, 4).astype(np.float32)
        mrc_new.delta = np.diag([1.0, 2.0, 3.0])
        mrc_new.origin = np.array([5.0, 10.0, 15.0])
        mrc_new.rank = 3
        
        # Write and read back
        mrc_new.write(outfile)
        mrc_read = mrc.MRC(outfile)
        
        # Verify
        assert_allclose(mrc_read.array, mrc_new.array, rtol=1e-5)
        assert_allclose(mrc_read.origin, mrc_new.origin, rtol=1e-4)
        assert_allclose(np.diag(mrc_read.delta), np.diag(mrc_new.delta), rtol=1e-5)
    
    def test_mrc_write_zero_voxel_raises(self, tmpdir):
        """Test that zero voxel size raises ValueError"""
        outfile = str(tmpdir / "invalid.mrc")
        
        mrc_obj = mrc.MRC()
        mrc_obj.array = np.ones((2, 2, 2), dtype=np.float32)
        mrc_obj.delta = np.diag([0.0, 1.0, 1.0])
        mrc_obj.origin = np.array([0.0, 0.0, 0.0])
        mrc_obj.rank = 3
        
        with pytest.raises(ValueError, match="Voxel size must be positive"):
            mrc_obj.write(outfile)


class TestGridMRCWrite:
    """Tests for Grid.export() with MRC format"""
    
    def test_grid_export_mrc(self, tmpdir):
        """Test Grid.export() with file_format='mrc'"""
        outfile = str(tmpdir / "grid.mrc")
        
        # Create simple grid
        data = np.arange(60).reshape(3, 4, 5).astype(np.float32)
        g = Grid(data, origin=[0, 0, 0], delta=[1.0, 1.0, 1.0])
        
        # Export and read back
        g.export(outfile, file_format='mrc')
        g_read = Grid(outfile)
        
        # Verify
        assert_allclose(g_read.grid, g.grid, rtol=1e-5)
        assert_allclose(g_read.origin, g.origin, rtol=1e-4)
        assert_allclose(g_read.delta, g.delta, rtol=1e-5)
    
    def test_grid_export_mrc_roundtrip(self, tmpdir):
        """Test MRC → Grid → export → Grid preserves data"""
        outfile = str(tmpdir / "roundtrip_grid.mrc")
        
        # Load original
        g_orig = Grid(datafiles.CCP4_1JZV)
        
        # Export and reload
        g_orig.export(outfile, file_format='mrc')
        g_read = Grid(outfile)
        
        # Verify
        assert_allclose(g_read.grid, g_orig.grid, rtol=1e-5)
        assert_allclose(g_read.origin, g_orig.origin, rtol=1e-4)
        assert_allclose(g_read.delta, g_orig.delta, rtol=1e-5)
        assert_equal(g_read.grid.shape, g_orig.grid.shape)
    
    def test_grid_export_mrc_preserves_header(self, tmpdir):
        """Test that Grid preserves MRC header through export"""
        outfile = str(tmpdir / "header_grid.mrc")
        
        g_orig = Grid(datafiles.CCP4_1JZV)
        orig_mapc = g_orig._mrc_header.mapc
        orig_mapr = g_orig._mrc_header.mapr
        orig_maps = g_orig._mrc_header.maps
        
        # Export and check
        g_orig.export(outfile, file_format='mrc')
        g_read = Grid(outfile)
        
        assert g_read._mrc_header.mapc == orig_mapc
        assert g_read._mrc_header.mapr == orig_mapr
        assert g_read._mrc_header.maps == orig_maps

    def test_mrc_write_4x4x4_with_header(self, tmpdir):
        """Test writing 4x4x4 MRC file with custom header values."""
        
        # Create 4x4x4 data
        data = np.arange(64, dtype=np.float32).reshape((4, 4, 4))
        outfile = str(tmpdir / "test_with_header.mrc")
        
        # Create and write MRC
        m = mrc.MRC()
        m.array = data
        m.delta = np.diag([1.5, 2.0, 2.5])
        m.origin = np.array([10.0, 20.0, 30.0])
        m.rank = 3
        m.write(outfile)
        
        # Read back and verify
        m_read = mrc.MRC(outfile)
        assert_allclose(m_read.array, data)
        assert_allclose(np.diag(m_read.delta), [1.5, 2.0, 2.5])
        assert_allclose(m_read.origin, [10.0, 20.0, 30.0], rtol=1e-4, atol=1.0)


    def test_mrc_write_4x4x4_without_header(self, tmpdir):
        """Test writing 4x4x4 MRC file with default header."""
        
        # Create 4x4x4 random data
        np.random.seed(42)
        data = np.random.rand(4, 4, 4).astype(np.float32)
        outfile = str(tmpdir / "test_without_header.mrc")
        
        # Create and write MRC
        m = mrc.MRC()
        m.array = data
        m.delta = np.diag([1.0, 1.0, 1.0])
        m.origin = np.array([0.0, 0.0, 0.0])
        m.rank = 3
        m.write(outfile)
        
        # Read back and verify
        m_read = mrc.MRC(outfile)
        assert_allclose(m_read.array, data)
        assert_allclose(np.diag(m_read.delta), [1.0, 1.0, 1.0])
        assert_allclose(m_read.origin, [0.0, 0.0, 0.0], rtol=1e-4, atol=1.0)
