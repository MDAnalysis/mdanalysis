import numpy as np
import pytest

from MDAnalysis.topology.MinimalParser import MinimalParser

from MDAnalysisTests.datafiles import (
    DCD,
    INPCRD,
    LAMMPSdcd2,
    NCDF,
    TRJ,
    TRJncdf,
    TRR,
    TRZ,
    XTC,
)


@pytest.mark.parametrize(
    'filename,expected_n_atoms', [
        (DCD, 3341),
        (INPCRD, 5),
        (LAMMPSdcd2, 12421),
        (NCDF, 2661),
        (TRR, 47681),
        (XTC, 47681),
        (np.zeros((10, 3)), 10),  # memory reader
    ])
def test_minimal_parser(filename, expected_n_atoms):
    with MinimalParser(filename) as p:
        top = p.parse()
    assert top.n_atoms == expected_n_atoms


@pytest.mark.parametrize('filename', [
    TRJ,
    TRJncdf,
    TRZ,
])
def test_minimal_parser_fail(filename):
    with MinimalParser(filename) as p:
        with pytest.raises(NotImplementedError):
            p.parse()
