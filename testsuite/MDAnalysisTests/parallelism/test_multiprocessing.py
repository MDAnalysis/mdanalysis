"""Test that MDAnalysis plays nicely with multiprocessing

"""
import multiprocessing
import numpy as np
import pytest

import MDAnalysis as mda
from MDAnalysisTests.datafiles import (
    PSF, DCD,
    GRO, XTC,
    PDB,
    XYZ,
)

from numpy.testing import assert_equal


@pytest.fixture(params=[
    (PSF, DCD),
    (GRO, XTC),
    (PDB,),
    (XYZ,),
])
def u(request):
    if len(request.param) == 1:
        f = request.param
        return mda.Universe(f)
    else:
        top, trj = request.param
        return mda.Universe(top, trj)

# Define target functions here
# inside test functions doesn't work
def cog(u, ag, frame_id):
    u.trajectory[frame_id]

    return ag.center_of_geometry()


def getnames(u, ix):
    # Check topology stuff works
    return u.atoms[ix].name


def test_multiprocess_COG(u):
    ag = u.atoms[10:20]

    ref = np.array([cog(u, ag, i)
                    for i in range(4)])

    p = multiprocessing.Pool(2)
    res = np.array([p.apply(cog, args=(u, ag, i))
                    for i in range(4)])
    p.close()
    assert_equal(ref, res)


def test_multiprocess_names(u):
    ref = [getnames(u, i)
           for i in range(10)]

    p = multiprocessing.Pool(2)
    res = [p.apply(getnames, args=(u, i))
                   for i in range(10)]
    p.close()

    assert_equal(ref, res)
