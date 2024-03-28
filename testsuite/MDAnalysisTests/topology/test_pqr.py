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
from io import StringIO
from numpy.testing import assert_equal, assert_almost_equal
import MDAnalysis as mda

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import (
    PQR,
    PQR_icodes,
)


class TestPQRParser(ParserBase):
    parser = mda.topology.PQRParser.PQRParser
    ref_filename = PQR
    expected_attrs = ['ids', 'names', 'charges', 'radii', 'record_types',
                      'resids', 'resnames', 'icodes',
                      'segids']
    guessed_attrs = ['masses', 'types']
    expected_n_atoms = 3341
    expected_n_residues = 214
    expected_n_segments = 1

    def test_attr_size(self, top):
        assert len(top.ids) == top.n_atoms
        assert len(top.names) == top.n_atoms
        assert len(top.charges) == top.n_atoms
        assert len(top.radii) == top.n_atoms
        assert len(top.resids) == top.n_residues
        assert len(top.resnames) == top.n_residues
        assert len(top.segids) == top.n_segments


class TestPQRParser2(TestPQRParser):
    ref_filename = PQR_icodes
    expected_n_atoms = 5313
    expected_n_residues = 474


def test_record_types():
    u = mda.Universe(PQR_icodes)

    assert u.atoms[4052].record_type == 'ATOM'
    assert u.atoms[4053].record_type == 'HETATM'

    assert_equal(u.atoms[:10].record_types, 'ATOM')
    assert_equal(u.atoms[4060:4070].record_types, 'HETATM')


GROMACS_PQR = '''
REMARK    The B-factors in this file hold atomic radii
REMARK    The occupancy in this file hold atomic charges
TITLE     system
REMARK    THIS IS A SIMULATION BOX
CRYST1   35.060   34.040   38.990  90.00  90.00  90.00 P 1           1
MODEL        1
ATOM      1  O    ZR     1      15.710  17.670  23.340 -0.67  1.48           O
TER
ENDMDL
'''

def test_gromacs_flavour():
    u = mda.Universe(StringIO(GROMACS_PQR), format='PQR')

    assert len(u.atoms) == 1
    # topology things
    assert u.atoms[0].type == 'O'
    assert u.atoms[0].segid == 'SYSTEM'
    assert not u._topology.types.is_guessed
    assert_almost_equal(u.atoms[0].radius, 1.48, decimal=5)
    assert_almost_equal(u.atoms[0].charge, -0.67, decimal=5)
    # coordinatey things
    assert_almost_equal(u.atoms[0].position, [15.710, 17.670, 23.340], decimal=4)
