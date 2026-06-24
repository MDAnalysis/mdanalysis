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

import pytest
import MDAnalysis as mda

import re

from MDAnalysis.fetch.fetchers import HAS_POOCH
from MDAnalysis.fetch.pdb import SUPPORTED_FILE_FORMATS_DOWNLOADER
from urllib import request
from pathlib import Path

try:
    request.urlopen("https://files.wwpdb.org/", timeout=2)
    HAS_ACCESS_TO_WWPDB = True
except request.URLError:
    HAS_ACCESS_TO_WWPDB = False


@pytest.mark.skipif(not HAS_POOCH, reason="Pooch is not installed.")
@pytest.mark.skipif(
    not HAS_ACCESS_TO_WWPDB,
    reason="Can not connect to https://files.wwpdb.org/",
)
def test_download_one_file(tmp_path):

    path = mda.fetch.from_PDB(["1AKE"], cache_path=tmp_path)
    assert path.exists()
    assert path.name == "1AKE.cif.gz"

@pytest.mark.skipif(not HAS_POOCH, reason="Pooch is not installed.")
@pytest.mark.skipif(
    not HAS_ACCESS_TO_WWPDB,
    reason="Can not connect to https://files.wwpdb.org/",
)
def test_download_multiple_files(tmp_path):

    paths = mda.fetch.from_PDB(["1AKE", "4AKE"], cache_path=tmp_path)
    assert all(isinstance(path, Path) for path in paths)
    assert [path.name for path in paths] == list(
        ["1AKE.cif.gz", "4AKE.cif.gz"]
    )

@pytest.mark.skipif(not HAS_POOCH, reason="Pooch is not installed.")
@pytest.mark.skipif(
    not HAS_ACCESS_TO_WWPDB,
    reason="Can not connect to https://files.wwpdb.org/",
)
def test_download_file_format(tmp_path):

    path = mda.fetch.from_PDB(["1AKE"], cache_path=tmp_path, file_format="pdb")
    assert path.exists()
    assert path.name == "1AKE.pdb"


@pytest.mark.skipif(not HAS_POOCH, reason="Pooch is not installed.")
@pytest.mark.skipif(
    not HAS_ACCESS_TO_WWPDB,
    reason="Can not connect to https://files.wwpdb.org/",
)
def test_invalid_file_format(tmp_path):
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Invalid file format. Supported file formats "
            f"are {SUPPORTED_FILE_FORMATS_DOWNLOADER}"
        ),
    ):
        mda.fetch.from_PDB(
            pdb_ids="1AKE", cache_path=tmp_path, file_format="barfoo"
        )
