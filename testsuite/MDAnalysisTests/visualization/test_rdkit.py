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

import os
from tempfile import NamedTemporaryFile
from io import BytesIO
import numpy as np
from numpy.testing import assert_equal
import pytest
from PIL import Image
import MDAnalysis as mda
from MDAnalysisTests.util import import_not_available
try:
    from MDAnalysis.visualization.RDKit import RDKitDrawer
except ImportError:
    class RDKitDrawer:
        pass


requires_rdkit = pytest.mark.skipif(import_not_available("rdkit"),
                                    reason="requires RDKit")


@pytest.mark.skipif(not import_not_available("rdkit"),
                    reason="only for min dependencies build")
class TestRequiresRDKit:
    def test_requires_rdkit(self):
        with pytest.raises(ImportError, match="RDKit is needed"):
            from MDAnalysis.visualization.RDKit import RDKitDrawer


@requires_rdkit
class TestRDKitDrawer:
    @pytest.fixture
    def u(self):
        return mda.Universe.from_smiles("CCO", numConfs=3)

    @pytest.fixture
    def drawer(self):
        return RDKitDrawer()

    def test_init(self, drawer):
        assert mda._FORMATTERS["RDKIT"][(
            mda.AtomGroup, "image/png")] == drawer._repr_atomgroup

    def test_repr_n_atoms_limit(self, u):
        d = RDKitDrawer(max_atoms=1)
        assert d._repr_atomgroup(u.atoms) is None

    def test_prepare_no3D(self, u, drawer):
        mol = drawer._prepare_atomgroup_for_drawing(u.atoms, keep_3D=False)
        assert_equal(mol.GetConformer().GetPositions()[:, 2], 
                     [0]*mol.GetNumAtoms())

    def test_ag_to_img_png(self, u, drawer):
        png = drawer.atomgroup_to_image(u.atoms)
        img = Image.open(BytesIO(png))
        assert img.size == drawer.size
        assert len(png) == 6349

    def test_ag_to_img_svg(self, u, drawer):
        svg = drawer.atomgroup_to_image(u.atoms, useSVG=True)
        assert len(svg) == 1275

    def test_ag_to_gif(self, u, drawer):
        with NamedTemporaryFile() as f:
            drawer.atomgroup_to_gif(u.atoms, output=f.name)
            assert os.path.isfile(f.name)
