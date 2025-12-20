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
from MDAnalysis.fetch.pdb import (
    DEFAULT_CACHE_NAME_DOWNLOADER,
    HAS_POOCH,
    SUPPORTED_FILE_FORMATS_DOWNLOADER,
)

import re
from urllib import request
from shutil import rmtree
from pathlib import Path

if HAS_POOCH:
    from requests.exceptions import HTTPError
    import pooch

try:
    request.urlopen("https://files.wwpdb.org/", timeout=2)
    HAS_ACCESS_TO_WWPDB = True
except request.URLError:
    HAS_ACCESS_TO_WWPDB = False


@pytest.mark.skipif(
    HAS_POOCH,
    reason="Pooch is installed.",
)
def test_pooch_installation(tmp_path):
    with pytest.raises(
        ModuleNotFoundError,
        match="pooch is needed as a dependency for from_PDB()",
    ):
        mda.fetch.from_PDB("1AKE", cache_path=tmp_path, file_format="cif")


@pytest.mark.skipif(not HAS_POOCH, reason="Pooch is not installed.")
@pytest.mark.skipif(
    not HAS_ACCESS_TO_WWPDB,
    reason="Can not connect to https://files.wwpdb.org/",
)
class TestDocstringExamples:
    """This class tests all the examples found in from_PDB's docstring"""

    @pytest.mark.parametrize("pdb_id", ["1AKE", "4BWZ"])
    def test_one_file_download(self, tmp_path, pdb_id):
        path = mda.fetch.from_PDB(
            pdb_id, cache_path=tmp_path, file_format="cif"
        )
        assert isinstance(path, Path)
        assert Path(path).name == f"{pdb_id}.cif"

    def test_multiple_files_download(self, tmp_path):
        list_of_path_strings = mda.fetch.from_PDB(
            ["1AKE", "4BWZ"], cache_path=tmp_path, progressbar=True
        )
        assert all(isinstance(pdb_id, Path) for pdb_id in list_of_path_strings)
        assert all(
            [
                Path(path).name == f"{name}.cif.gz"
                for path, name in zip(
                    list_of_path_strings, ["1AKE", "4BWZ"], strict=True
                )
            ]
        )

    @pytest.mark.parametrize(
        "pdb_id, n_atoms", [("1AKE", 3816), ("4BWZ", 2824)]
    )
    def test_files_to_universe(self, tmp_path, pdb_id, n_atoms):
        u = mda.Universe(
            mda.fetch.from_PDB(
                pdb_id,
                file_format="pdb.gz",
                cache_path=tmp_path,
                progressbar=True,
            )
        )
        assert isinstance(u, mda.Universe) and (len(u.atoms) == n_atoms)


@pytest.fixture()
def clean_up_default_cache():
    rmtree(pooch.os_cache(DEFAULT_CACHE_NAME_DOWNLOADER), ignore_errors=True)
    yield
    rmtree(pooch.os_cache(DEFAULT_CACHE_NAME_DOWNLOADER))


@pytest.mark.skipif(not HAS_POOCH, reason="Pooch is not installed.")
@pytest.mark.skipif(
    not HAS_ACCESS_TO_WWPDB,
    reason="Can not connect to https://files.wwpdb.org/",
)
class TestExpectedBehaviors:

    def test_no_cache_path(self, clean_up_default_cache):
        assert isinstance(mda.fetch.from_PDB("1AKE", cache_path=None), Path)

    def test_str_input_gives_path_output(self, tmp_path):
        assert isinstance(
            mda.fetch.from_PDB(
                pdb_ids="1AKE", cache_path=tmp_path, file_format="cif"
            ),
            Path,
        )

    def test_list_input_gives_list_output(self, tmp_path):
        assert isinstance(
            mda.fetch.from_PDB(pdb_ids=["1AKE"], cache_path=tmp_path), list
        )


@pytest.mark.skipif(not HAS_POOCH, reason="Pooch is not installed.")
@pytest.mark.skipif(
    not HAS_ACCESS_TO_WWPDB,
    reason="Can not connect to https://files.wwpdb.org/",
)
class TestExpectedErrors:

    def test_invalid_pdb(self, tmp_path):
        with pytest.raises(HTTPError):
            mda.fetch.from_PDB(pdb_ids="foobar", cache_path=tmp_path)

    def test_invalid_file_format(self, tmp_path):
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
