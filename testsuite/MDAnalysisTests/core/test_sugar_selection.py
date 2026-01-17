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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis of BioMolecular Simulations. J. Comput. Chem. 32 (2011), 2319-2327,
# doi:10.1002/jcc.21787
#

import pytest
import numpy as np
import MDAnalysis as mda


def test_sugar_selection():
    # 0GA = Glucose, ALA = Protein, 4LB = Galactose
    resnames = np.array(["0GA", "ALA", "4LB"])
    # 3 atoms per residue
    resindices = np.repeat([0, 1, 2], 3)
    # 1 segment
    segindices = np.repeat([0], 3)

    u = mda.Universe.empty(
        n_atoms=9,
        n_residues=3,
        n_segments=1,
        atom_resindex=resindices,
        residue_segindex=segindices,
        trajectory=True,
    )
    u.add_TopologyAttr("resname", resnames)

    sugar = u.select_atoms("sugar")

    assert len(sugar) == 6
    assert np.all(np.isin(sugar.resnames, ["0GA", "4LB"]))
    assert "ALA" not in sugar.resnames


def test_sugar_selection_empty():
    resnames = np.array(["ALA", "VAL"])
    resindices = np.repeat([0, 1], 3)
    segindices = np.repeat([0], 2)

    u = mda.Universe.empty(
        n_atoms=6,
        n_residues=2,
        n_segments=1,
        atom_resindex=resindices,
        residue_segindex=segindices,
        trajectory=True,
    )
    u.add_TopologyAttr("resname", resnames)

    sugar = u.select_atoms("sugar")

    assert len(sugar) == 0
