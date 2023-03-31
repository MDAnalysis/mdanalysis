import numpy as np
from numpy.testing import (assert_array_equal, assert_array_almost_equal,
                           assert_almost_equal)

import pytest

from gridData import Grid

def f_arithmetic(g):
    return g + g - 2.5 * g / (g + 5.3)

@pytest.fixture(scope="class")
def data():
    d = dict(
        griddata=np.arange(1, 28).reshape(3, 3, 3),
        origin=np.zeros(3),
        delta=np.ones(3))
    d['grid'] = Grid(d['griddata'], origin=d['origin'],
                     delta=d['delta'])
    return d

class TestGrid(object):
    @pytest.fixture
    def pklfile(self, data, tmpdir):
        g = data['grid']
        fn = tmpdir.mkdir('grid').join('grid.dat')
        g.save(fn)  # always saves as pkl
        return fn

    def test_init(self, data):
        g = Grid(data['griddata'], origin=data['origin'],
                 delta=1)
        assert_array_equal(g.delta, data['delta'])

    def test_init_wrong_origin(self, data):
        with pytest.raises(TypeError):
            Grid(data['griddata'], origin=np.ones(4), delta=data['delta'])

    def test_init_wrong_delta(self, data):
        with pytest.raises(TypeError):
            Grid(data['griddata'], origin=data['origin'], delta=np.ones(4))

    def test_empty_Grid(self):
        g = Grid()
        assert isinstance(g, Grid)

    def test_init_missing_delta_ValueError(self, data):
        with pytest.raises(ValueError):
            Grid(data['griddata'], origin=data['origin'])

    def test_init_missing_origin_ValueError(self, data):
        with pytest.raises(ValueError):
            Grid(data['griddata'], delta=data['delta'])

    def test_init_wrong_data_exception(self):
        with pytest.raises(IOError):
            Grid("__does_not_exist__")

    def test_load_wrong_fileformat_ValueError(self):
        with pytest.raises(ValueError):
            Grid(grid=True, file_format="xxx")

    def test_equality(self, data):
        assert data['grid'] == data['grid']
        assert data['grid'] != 'foo'
        g = Grid(data['griddata'], origin=data['origin'] + 1, delta=data['delta'])
        assert data['grid'] != g

    def test_addition(self, data):
        g = data['grid'] + data['grid']
        assert_array_equal(g.grid.flat, (2 * data['griddata']).flat)
        g = 2 + data['grid']
        assert_array_equal(g.grid.flat, (2 + data['griddata']).flat)
        g = g + data['grid']
        assert_array_equal(g.grid.flat, (2 + (2 * data['griddata'])).flat)

    def test_subtraction(self, data):
        g = data['grid'] - data['grid']
        assert_array_equal(g.grid.flat, np.zeros(27))
        g = 2 - data['grid']
        assert_array_equal(g.grid.flat, (2 - data['griddata']).flat)

    def test_multiplication(self, data):
        g = data['grid'] * data['grid']
        assert_array_equal(g.grid.flat, (data['griddata'] ** 2).flat)
        g = 2 * data['grid']
        assert_array_equal(g.grid.flat, (2 * data['griddata']).flat)

    def test_division(self, data):
        g = data['grid'] / data['grid']
        assert_array_equal(g.grid.flat, np.ones(27))
        g = 2 / data['grid']
        assert_array_equal(g.grid.flat, (2 / data['griddata']).flat)

    def test_floordivision(self, data):
        g = data['grid'].__floordiv__(data['grid'])
        assert_array_equal(g.grid.flat, np.ones(27, dtype=np.int64))
        g = 2 // data['grid']
        assert_array_equal(g.grid.flat, (2 // data['griddata']).flat)

    def test_power(self, data):
        g = data['grid'] ** 2
        assert_array_equal(g.grid.flat, (data['griddata'] ** 2).flat)
        g = 2 ** data['grid']
        assert_array_equal(g.grid.flat, (2 ** data['griddata']).flat)

    def test_compatibility_type(self, data):
        assert data['grid'].check_compatible(data['grid'])
        assert data['grid'].check_compatible(3)
        g = Grid(data['griddata'], origin=data['origin'] - 1, delta=data['delta'])
        assert data['grid'].check_compatible(g)

    def test_wrong_compatibile_type(self, data):
        with pytest.raises(TypeError):
            data['grid'].check_compatible("foo")

    def test_non_orthonormal_boxes(self, data):
        delta = np.eye(3)
        with pytest.raises(NotImplementedError):
            Grid(data['griddata'], origin=data['origin'], delta=delta)

    def test_centers(self, data):
        # this only checks the edges. If you know an alternative
        # algorithm that isn't an exact duplicate of the one in
        # g.centers to test this please implement it.
        g = Grid(data['griddata'], origin=np.ones(3), delta=data['delta'])
        centers = np.array(list(g.centers()))
        assert_array_equal(centers[0], g.origin)
        assert_array_equal(centers[-1] - g.origin,
                           (np.array(g.grid.shape) - 1) * data['delta'])

    def test_resample_factor_failure(self, data):
        pytest.importorskip('scipy')

        with pytest.raises(ValueError):
            g = data['grid'].resample_factor(0)

    def test_resample_factor(self, data):
        pytest.importorskip('scipy')

        g = data['grid'].resample_factor(2)
        assert_array_equal(g.delta, np.ones(3) * .5)
        # zooming in by a factor of 2. Each subinterval is
        # split in half, so 3 gridpoints (2 subintervals)
        # becomes 5 gridpoints (4 subintervals)
        assert_array_equal(g.grid.shape, np.ones(3) * 5)
        # check that the values are identical with the
        # correct stride.
        assert_array_almost_equal(g.grid[::2, ::2, ::2],
                                  data['grid'].grid)

    def test_load_pickle(self, data, tmpdir):
        g = data['grid']
        fn = str(tmpdir.mkdir('grid').join('grid.pkl'))
        g.save(fn)

        h = Grid()
        h.load(fn)

        assert h == g

    def test_init_pickle_pathobjects(self, data, tmpdir):
        g = data['grid']
        fn = tmpdir.mkdir('grid').join('grid.pickle')
        g.save(fn)

        h = Grid(fn)

        assert h == g

    @pytest.mark.parametrize("fileformat", ("pkl", "PKL", "pickle", "python"))
    def test_load_fileformat(self, data, pklfile, fileformat):
        h = Grid(pklfile, file_format="pkl")
        assert h == data['grid']

    # At the moment, reading the file with the wrong parser does not give
    # good error messages.
    @pytest.mark.xfail
    @pytest.mark.parametrize("fileformat", ("ccp4", "plt", "dx"))
    def test_load_wrong_fileformat(self, data, pklfile, fileformat):
        with pytest.raises('ValueError'):
            Grid(pklfile, file_format=fileformat)

    # just check that we can export without stupid failures; detailed
    # format checks in separate tests
    @pytest.mark.parametrize("fileformat", ("dx", "pkl"))
    def test_export(self, data, fileformat, tmpdir):
        g = data['grid']
        fn = tmpdir.mkdir('grid_export').join("grid.{}".format(fileformat))
        g.export(fn)   # check that path objects work
        h = Grid(fn)   # use format autodetection
        assert g == h

    @pytest.mark.parametrize("fileformat", ("ccp4", "plt"))
    def test_export_not_supported(self, data, fileformat, tmpdir):
        g = data['grid']
        fn = tmpdir.mkdir('grid_export').join("grid.{}".format(fileformat))
        with pytest.raises(ValueError):
            g.export(fn)


def test_inheritance(data):
    class DerivedGrid(Grid):
        pass

    dg = DerivedGrid(data['griddata'], origin=data['origin'],
                     delta=data['delta'])
    result = f_arithmetic(dg)

    assert isinstance(result, DerivedGrid)

    ref = f_arithmetic(data['grid'])
    assert_almost_equal(result.grid, ref.grid)

def test_anyarray(data):
    ma = np.ma.MaskedArray(data['griddata'])
    mg = Grid(ma, origin=data['origin'], delta=data['delta'])

    assert isinstance(mg.grid, ma.__class__)

    result = f_arithmetic(mg)
    ref = f_arithmetic(data['grid'])

    assert_almost_equal(result.grid, ref.grid)
