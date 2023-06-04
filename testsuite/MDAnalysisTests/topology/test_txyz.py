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
import MDAnalysis as mda
import pytest
import numpy as np

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import TXYZ, ARC, ARC_PBC
from numpy.testing import assert_equal

class TestTXYZParser(ParserBase):
    parser = mda.topology.TXYZParser.TXYZParser
    guessed_attrs = ['masses']
    expected_attrs = ['ids', 'names', 'bonds', 'types', 'elements']
    expected_n_residues = 1
    expected_n_atoms = 9
    expected_n_segments = 1

    def test_number_of_bonds(self, top):
        assert len(top.bonds.values) == 8

    @pytest.fixture(params=(TXYZ, ARC))
    def filename(self, request):
        return request.param

    def test_atom_type_type(self, top):
        """
        ``u.atoms.types`` contains strings.

        See `issue #2435 <https://github.com/MDAnalysis/mdanalysis/issues/2435>
        """
        types = top.types.values
        type_is_str = [isinstance(atom_type, str) for atom_type in types]
        assert all(type_is_str)

def test_TXYZ_elements():
    """The test checks whether elements attribute are assigned
    properly given a TXYZ file with valid elements record.
    """
    u = mda.Universe(TXYZ, format='TXYZ')
    element_list = np.array(['C', 'H', 'H', 'O', 'H', 'C', 'H', 'H', 'H'], dtype=object)
    assert_equal(u.atoms.elements, element_list)
    
    
def test_missing_elements_noattribute():
    """Check that:

    1) a warning is raised if elements are missing
    2) the elements attribute is not set
    """
    wmsg = ("Element information is missing, elements attribute will not be "
            "populated")
    with pytest.warns(UserWarning, match=wmsg):
        u = mda.Universe(ARC_PBC)
    with pytest.raises(AttributeError):
        _ = u.atoms.elements

