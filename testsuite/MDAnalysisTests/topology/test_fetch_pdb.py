# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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


# Note need to make test not run when there is no internet


from pathlib import Path

#import pytest
import MDAnalysis as mda

class TestDocstringExamples:
    """This class test the example found in fetch_pdb's docstring"""

    def test_one_file_download(self, tmp_path):
        """Tests docstring's mda.fetch_pdb("1AKE", file_format="cif") """
        assert isinstance(mda.fetch_pdb("1AKE", cache_path=tmp_path, file_format="cif"), str)

    def test_multiple_files_download(self, tmp_path):
        """Tests docstring's mda.fetch_pdb(["1AKE", "4BWZ"], progressbar=True) """
        list_of_path_string = mda.fetch_pdb(["1AKE", "4BWZ"], cache_path=tmp_path, progressbar=True)
        assert all(isinstance(PDB_ID, str) for PDB_ID in list_of_path_string)

    def test_one_file_to_universe(self, tmp_path):
        """Tests docstring's mda.Universe(mda.fetch_pdb("1AKE"), file_format="pdb.gz") """
        assert isinstance(mda.Universe(mda.fetch_pdb("1AKE"), file_format="pdb.gz", cache_path=tmp_path, progressbar=True), mda.Universe)

    def test_multiple_files_to_universe(self, tmp_path):
        """Tests docstring's [mda.Universe(mda.fetch_pdb(PDB_ID), file_format="pdb.gz") for PDB_ID in ("1AKE", "4BWZ")] """
        list_of_path_string = [mda.Universe(mda.fetch_pdb(PDB_ID), file_format="pdb.gz") for PDB_ID in ("1AKE", "4BWZ")]
        assert all(isinstance(PDB_ID, mda.Universe) for PDB_ID in list_of_path_string)




