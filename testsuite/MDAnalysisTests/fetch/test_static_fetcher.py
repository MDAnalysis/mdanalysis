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
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, N. M. Melo, S. L. Seyler,
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

from pathlib import Path
import re
from shutil import rmtree
from urllib import request

import pytest

from MDAnalysis.fetch.fetchers import (
    DEFAULT_CACHE_NAME_DOWNLOADER,
    HAS_POOCH,
    StaticFetcher,
)

if HAS_POOCH:
    import pooch

try:
    request.urlopen("https://files.wwpdb.org/", timeout=2)
    HAS_ACCESS_TO_WWPDB = True
except request.URLError:
    HAS_ACCESS_TO_WWPDB = False


BASE_URL = "https://files.wwpdb.org/download/"
SINGLE_PDB = "1AKE.pdb"
MULTIPLE_PDBS = ("1AKE.pdb", "4AKE.pdb")
REGISTRY_NAME = "test_db.txt"


@pytest.fixture()
def clean_up_default_cache():
    rmtree(pooch.os_cache(DEFAULT_CACHE_NAME_DOWNLOADER), ignore_errors=True)
    yield
    rmtree(pooch.os_cache(DEFAULT_CACHE_NAME_DOWNLOADER), ignore_errors=True)


@pytest.mark.skipif(HAS_POOCH, reason="Pooch is installed.")
def test_pooch_installation():
    with pytest.raises(
        ModuleNotFoundError,
        match="pooch is needed as a dependency for Fetchers",
    ):
        StaticFetcher()


@pytest.mark.skipif(not HAS_POOCH, reason="Pooch is not installed.")
class TestExpectedErrors:

    def test_missing_base_url(self, tmp_path):
        downloader = StaticFetcher(cache_path=tmp_path)

        with pytest.raises(
            ValueError, match=re.escape("base_url is not defined in fetch()")
        ):
            downloader.fetch(file_name=SINGLE_PDB)

    def test_invalid_downloader(self, tmp_path):
        downloader = StaticFetcher(cache_path=tmp_path)

        with pytest.raises(
            ValueError,
            match=re.escape(
                "Invalid downloader 'barfoo'. Valid options are "
                "'HTTP', 'FTP', 'SFTP', 'DOI'."
            ),
        ):
            downloader.fetch(
                base_url=BASE_URL,
                file_name=SINGLE_PDB,
                downloader="barfoo",
            )

    def test_invalid_hash(self, tmp_path):
        with pytest.raises(ValueError, match='Invalid hash "barfoo"'):
            StaticFetcher(cache_path=tmp_path, hash="barfoo")


@pytest.mark.skipif(not HAS_POOCH, reason="Pooch is not installed.")
@pytest.mark.skipif(
    not HAS_ACCESS_TO_WWPDB,
    reason="Can not connect to https://files.wwpdb.org/",
)
class TestExpectedBehaviors:

    def test_default_cache_path(self, clean_up_default_cache):
        downloader = StaticFetcher()
        path = downloader.fetch(base_url=BASE_URL, file_name=SINGLE_PDB)

        assert isinstance(path, Path)
        assert path.name == SINGLE_PDB
        assert path.exists()

    def test_create_database(self, tmp_path):
        downloader = StaticFetcher(cache_path=tmp_path)
        path = downloader.fetch(
            base_url=BASE_URL, file_name=SINGLE_PDB, db_name=REGISTRY_NAME
        )

        assert isinstance(path, Path)
        assert path.name == SINGLE_PDB
        assert path.exists()
        assert (tmp_path / REGISTRY_NAME).exists()

    def test_existing_database(self, tmp_path):
        downloader = StaticFetcher(cache_path=tmp_path)
        downloader.fetch(
            base_url=BASE_URL, file_name=SINGLE_PDB, db_name=REGISTRY_NAME
        )

        downloader = StaticFetcher(cache_path=tmp_path)
        path = downloader.fetch(
            base_url=BASE_URL, file_name=SINGLE_PDB, db_name=REGISTRY_NAME
        )

        assert isinstance(path, Path)
        assert path.name == SINGLE_PDB
        assert path.exists()

    def test_no_database(self, tmp_path):
        downloader = StaticFetcher(cache_path=tmp_path)
        path = downloader.fetch(
            base_url=BASE_URL, file_name=SINGLE_PDB, db_name=None
        )

        assert isinstance(path, Path)
        assert path.name == SINGLE_PDB
        assert downloader.db_path is None
        assert path.exists()
        assert not (tmp_path / REGISTRY_NAME).exists()

    def test_multiple_downloads_create_database(self, tmp_path):
        downloader = StaticFetcher(cache_path=tmp_path)
        paths = downloader.fetch(
            base_url=BASE_URL, file_name=MULTIPLE_PDBS, db_name=REGISTRY_NAME
        )

        assert all(isinstance(path, Path) for path in paths)
        assert [path.name for path in paths] == list(MULTIPLE_PDBS)
        assert all(path.exists() for path in paths)
        assert (tmp_path / REGISTRY_NAME).exists()

    def test_multiple_downloads_existing_database(self, tmp_path):
        downloader = StaticFetcher(cache_path=tmp_path)
        downloader.fetch(
            base_url=BASE_URL, file_name=MULTIPLE_PDBS, db_name=REGISTRY_NAME
        )

        assert (tmp_path / REGISTRY_NAME).exists()

        downloader = StaticFetcher(cache_path=tmp_path)
        paths = downloader.fetch(
            base_url=BASE_URL, file_name=MULTIPLE_PDBS, db_name=REGISTRY_NAME
        )

        assert all(isinstance(path, Path) for path in paths)
        assert [path.name for path in paths] == list(MULTIPLE_PDBS)
        assert all(path.exists() for path in paths)
