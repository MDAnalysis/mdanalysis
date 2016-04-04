# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from numpy.testing import (
    assert_,
)
import MDAnalysis as mda

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import (
    PDB,
    PDB_small,
    PDB_conect
)


class TestPDBParser(ParserBase):
    """This one has neither chainids or segids"""
    parser = mda.topology.PrimitivePDBParser.PrimitivePDBParser
    filename = PDB
    expected_attrs = ['ids', 'names',
                      'resids', 'resnames']
    expected_n_atoms = 47681
    expected_n_residues = 11302
    expected_n_segments = 1


class TestPDBParserSegids(ParserBase):
    """Has segids"""
    parser = mda.topology.PrimitivePDBParser.PrimitivePDBParser
    filename = PDB_small
    expected_attrs = ['ids', 'names',
                      'resids', 'resnames', 'segids']
    expected_n_atoms = 3341
    expected_n_residues = 214
    expected_n_segments = 1

class TestPDBConect(object):
    """Testing PDB topology parsing (PrimitivePDB)"""
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
            bonds = mda.topology.PrimitivePDBParser._parse_conect(line)
            assert_equal(bonds[0], res[0])
            for bond, res_bond in zip(bonds[1], res[1]):
                assert_equal(bond, res_bond)

        assert_raises(RuntimeError,
                      mda.topology.PrimitivePDBParser._parse_conect,
                      'CONECT12337 7718 84081234012344123')

    @staticmethod
    def test_conect_topo_parser():
        """Check that the parser works as intended,
        and that the returned value is a dictionary
        """
        with mda.topology.PrimitivePDBParser(PDB_conect) as p:
            top = p.parse()
            assert_(isinstance(top, mda.core.topology.Topology))
