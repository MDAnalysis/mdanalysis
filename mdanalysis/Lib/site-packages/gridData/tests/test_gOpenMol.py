import numpy as np
from numpy.testing import (assert_almost_equal,
                           assert_equal)

from gridData import Grid

from . import datafiles

def test_gOpenMol():
    g = Grid(datafiles.gOpenMol)
    POINTS = 192
    assert g.grid.size == 165048
    assert_almost_equal(g.delta, [1.0, 1.0, 1.0])
    assert_equal(g.grid.shape, (46, 46, 78))
    assert_almost_equal(g.origin, [0.5995016, 0.5995016, 0.5919984])
    assert_almost_equal(g.grid[::20, ::20, ::30],
                        np.array([[[1.02196848, 0.        , 0.88893718],
                                   [0.99051529, 0.        , 0.95906246],
                                   [0.96112466, 0.        , 0.88996845]],

                                  [[0.97247058, 0.        , 0.91574967],
                                   [1.00237465, 1.34423399, 0.87810922],
                                   [0.97917157, 0.        , 0.84717268]],

                                  [[0.99103099, 0.        , 0.86521846],
                                   [0.96421844, 0.        , 0.        ],
                                   [0.98432779, 0.        , 0.8817184 ]]])
    )
    assert_almost_equal(g.grid.mean(), 0.5403224581733577)
