# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import absolute_import
from numpy.testing import (
    assert_,
    assert_array_equal,
    assert_equal,
    assert_raises,
    assert_warns,
)
import MDAnalysis as mda

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import (
    PDB,
    PDB_small,
    PDB_conect,
    PDB_conect2TER,
    PDB_singleconect,
    PDB_chainidnewres,
    PDB_sameresid_diffresname
)
from MDAnalysis.topology.PDBParser import PDBParser


_PDBPARSER = mda.topology.PDBParser.PDBParser


class TestPDBParser(ParserBase):
    """This one has neither chainids or segids"""
    parser = mda.topology.PDBParser.PDBParser
    filename = PDB
    expected_attrs = ['ids', 'names',
                      'resids', 'resnames']
    guessed_attrs = ['types', 'masses']
    expected_n_atoms = 47681
    expected_n_residues = 11302
    expected_n_segments = 1


class TestPDBParserSegids(ParserBase):
    """Has segids"""
    parser = mda.topology.PDBParser.PDBParser
    filename = PDB_small
    expected_attrs = ['ids', 'names',
                      'resids', 'resnames', 'segids']
    guessed_attrs = ['types', 'masses']
    expected_n_atoms = 3341
    expected_n_residues = 214
    expected_n_segments = 1

class TestPDBConect(object):
    """Testing PDB topology parsing (PDB)"""
    @staticmethod
    def test_conect_parser():
        lines = ("CONECT1233212331",
                 "CONECT123331233112334",
                 "CONECT123341233312335",
                 "CONECT123351233412336",
                 "CONECT12336123271233012335",
                 "CONECT12337 7718 84081234012344",
                 "CONECT1233812339123401234112345")
        results = ((12332, [12331]),
                   (12333, [12331, 12334]),
                   (12334, [12333, 12335]),
                   (12335, [12334, 12336]),
                   (12336, [12327, 12330, 12335]),
                   (12337, [7718, 8408, 12340, 12344]),
                   (12338, [12339, 12340, 12341, 12345]))
        for line, res in zip(lines, results):
            bonds = mda.topology.PDBParser._parse_conect(line)
            assert_equal(bonds[0], res[0])
            for bond, res_bond in zip(bonds[1], res[1]):
                assert_equal(bond, res_bond)

    @staticmethod
    def test_conect_parser_runtime():
        assert_raises(RuntimeError,
                      mda.topology.PDBParser._parse_conect,
                      'CONECT12337 7718 84081234012344123')

    @staticmethod
    def test_conect_topo_parser():
        """Check that the parser works as intended,
        and that the returned value is a dictionary
        """
        with _PDBPARSER(PDB_conect) as p:
            top = p.parse()
            assert_(isinstance(top, mda.core.topology.Topology))


def test_conect2ter():
    def parse():
        with PDBParser(PDB_conect2TER) as p:
            struc = p.parse()
        return struc
    assert_warns(UserWarning, parse)
    struc = parse()

    assert_(hasattr(struc, 'bonds'))
    assert_(len(struc.bonds.values) == 4)


def test_single_conect():
    def parse():
        with PDBParser(PDB_singleconect) as p:
            struc = p.parse()
        return struc
    assert_warns(UserWarning, parse)
    struc = parse()
    assert_(hasattr(struc, 'bonds'))
    assert_(len(struc.bonds.values) == 2)


def test_new_chainid_new_res():
    # parser must start new residue when chainid starts

    u = mda.Universe(PDB_chainidnewres)

    assert_(len(u.residues) == 4)
    assert_array_equal(u.residues.resids, [1, 2, 3, 3])
    assert_(len(u.segments) == 4)
    assert_array_equal(u.segments.segids, ['A', 'B', 'C', 'D'])
    assert_(len(u.segments[0].atoms) == 5)
    assert_(len(u.segments[1].atoms) == 5)
    assert_(len(u.segments[2].atoms) == 5)
    assert_(len(u.segments[3].atoms) == 7)

def test_sameresid_diffresname():
    with _PDBPARSER(PDB_sameresid_diffresname) as p:
        top = p.parse()
    resids = [9, 9]
    resnames = ['GLN', 'POPC']
    for i, (resid, resname) in enumerate(zip(resids, resnames)):
        assert_(top.resids.values[i] == resid)
        assert_(top.resnames.values[i] == resname)
