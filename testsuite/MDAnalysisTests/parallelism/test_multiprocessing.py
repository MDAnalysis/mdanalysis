"""Test that MDAnalysis plays nicely with multiprocessing

"""
import multiprocessing
import numpy as np
import pytest

import MDAnalysis as mda
from MDAnalysisTests.datafiles import (
    PSF, DCD
)

from numpy.testing import assert_equal


@pytest.fixture
def u():
    return mda.Universe(PSF, DCD)


def cog(u, ag, frame_id):
    u.trajectory[frame_id]

    return ag.center_of_geometry()


def test_multiprocess_COM(u):
    ag = u.atoms[10:20]

    ref = np.array([cog(u, ag, i)
                    for i in range(4)])

    p = multiprocessing.Pool(2)

    res = np.array([p.apply(cog, args=(u, ag, i))
                    for i in range(4)])

    assert_equal(ref, res)
