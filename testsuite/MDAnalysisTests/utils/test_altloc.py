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
from MDAnalysis import Universe
from numpy.testing import assert_equal
from MDAnalysisTests.datafiles import PDB_full


@pytest.fixture()
def u():
    return Universe(PDB_full, guess_bonds=True)


def test_atomgroups(u):
    segidB0 = len(u.select_atoms("segid B and (not altloc B)"))
    segidB1 = len(u.select_atoms("segid B and (not altloc A)"))
    assert_equal(segidB0, segidB1)
    altlocB0 = len(u.select_atoms("segid B and (altloc A)"))
    altlocB1 = len(u.select_atoms("segid B and (altloc B)"))
    assert_equal(altlocB0, altlocB1)
    sum = len(u.select_atoms("segid B"))
    assert_equal(sum, segidB0 + altlocB0)


def test_bonds(u):
    # need to force topology to load before querying individual atom bonds
    bonds0 = u.select_atoms("segid B and (altloc A)")[0].bonds
    bonds1 = u.select_atoms("segid B and (altloc B)")[0].bonds
    assert_equal(len(bonds0), len(bonds1))


def test_write_read(u, tmpdir):
    outfile = str(tmpdir.join('test.pdb'))
    u.select_atoms("all").write(outfile)
    u2 = Universe(outfile)
    assert len(u.atoms) == len(u2.atoms)
