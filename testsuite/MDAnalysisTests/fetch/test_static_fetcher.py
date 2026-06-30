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

from servers import temporary_http_server

import hashlib
import pytest


from MDAnalysis.fetch.fetchers import (
    DEFAULT_CACHE_NAME_DOWNLOADER,
    HAS_POOCH,
    StaticFetcher,
)

if HAS_POOCH:
    import pooch

REGISTRY_NAME = "hashes.txt"


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
        with temporary_http_server() as (host, port, temp_folder):
            downloader = StaticFetcher(cache_path=tmp_path)

            with pytest.raises(
                ValueError,
                match=re.escape("base_url is not defined in fetch()"),
            ):
                downloader.fetch(file_name="TEST_FILE1.txt")

    def test_invalid_downloader(self, tmp_path):
        with temporary_http_server() as (host, port, temp_folder):
            base_url = f"http://{host}:{port}/"
            downloader = StaticFetcher(cache_path=tmp_path)

            with pytest.raises(
                ValueError,
                match=re.escape(
                    "Invalid downloader 'barfoo'. Valid options are "
                    "'HTTP', 'FTP', 'SFTP', 'DOI'."
                ),
            ):
                downloader.fetch(
                    base_url=base_url,
                    file_name="TEST_FILE1.txt",
                    downloader="barfoo",
                )

    def test_invalid_hash(self, tmp_path):
        hash = "foo"

        with temporary_http_server() as (host, port, temp_folder):
            base_url = f"http://{host}:{port}/"

            with pytest.raises(
                ValueError,
                match=re.escape(
                    f'Invalid hash "{hash}". Valid hashes algorithms are {hashlib.algorithms_available}.'
                ),
            ):
                downloader = StaticFetcher(cache_path=tmp_path, hash=hash)


@pytest.mark.skipif(not HAS_POOCH, reason="Pooch is not installed.")
class TestExpectedBehaviors:

    def test_default_cache_path(self, clean_up_default_cache):

        with temporary_http_server() as (host, port, temp_folder):
            base_url = f"http://{host}:{port}/"
            downloader = StaticFetcher()
            path = downloader.fetch(
                base_url=base_url, file_name="TEST_FILE1.txt"
            )

            assert isinstance(path, Path)
            assert path.name == "TEST_FILE1.txt"
            assert (
                path.read_text()
                == "The USA is going to win the 2026 World Cup!\nU-S-A! U-S-A! U-S-A!"
            )
            assert path.exists()

    def test_create_database(self, tmp_path):
        with temporary_http_server() as (host, port, temp_folder):
            base_url = f"http://{host}:{port}/"
            downloader = StaticFetcher(cache_path=tmp_path)
            path = downloader.fetch(
                base_url=base_url,
                file_name="TEST_FILE1.txt",
                db_name=REGISTRY_NAME,
            )

            assert isinstance(path, Path)
            assert path.name == "TEST_FILE1.txt"
            assert (
                path.read_text()
                == "The USA is going to win the 2026 World Cup!\nU-S-A! U-S-A! U-S-A!"
            )
            assert path.exists()

            assert (downloader.db_path).exists()
            assert (
                downloader.db_path
            ).read_text() == "TEST_FILE1.txt sha256:c4bdb6ba200a917b8384ffeffa4999bf05bd4e479f6580d795aca509c9122dc4\n"

    def test_different_hashes(self, tmp_path):
        with temporary_http_server() as (host, port, temp_folder):
            base_url = f"http://{host}:{port}/"
            downloader = StaticFetcher(cache_path=tmp_path, hash="md5")
            path = downloader.fetch(
                base_url=base_url,
                file_name="TEST_FILE1.txt",
                db_name=REGISTRY_NAME,
            )

            assert (
                downloader.db_path
            ).read_text() == "TEST_FILE1.txt md5:b2f138521297db74b6b280feeb14f9f6\n"

    def test_existing_database(self, tmp_path):
        with temporary_http_server() as (host, port, temp_folder):
            base_url = f"http://{host}:{port}/"
            downloader = StaticFetcher(cache_path=tmp_path)

            p1 = downloader.fetch(
                base_url=base_url,
                file_name="TEST_FILE1.txt",
                db_name=REGISTRY_NAME,
            )

            downloader = StaticFetcher(cache_path=tmp_path)
            p2 = downloader.fetch(
                base_url=base_url,
                file_name="TEST_FILE1.txt",
                db_name=REGISTRY_NAME,
            )

            assert p1.stat().st_mtime == p2.stat().st_mtime

    def test_no_database(self, tmp_path):
        with temporary_http_server() as (host, port, temp_folder):
            base_url = f"http://{host}:{port}/"
            downloader = StaticFetcher(cache_path=tmp_path)
            path = downloader.fetch(
                base_url=base_url, file_name="TEST_FILE1.txt", db_name=None
            )

            assert isinstance(path, Path)
            assert path.exists()
            assert path.name == "TEST_FILE1.txt"

            assert downloader.db_path is None
            assert not (tmp_path / REGISTRY_NAME).exists()

    def test_multiple_downloads_no_database(self, tmp_path):

        with temporary_http_server() as (host, port, temp_folder):
            base_url = f"http://{host}:{port}/"
            downloader = StaticFetcher(cache_path=tmp_path)
            paths = downloader.fetch(
                base_url=base_url,
                file_name=("TEST_FILE1.txt", "TEST_FILE2.txt"),
                db_name=None,
            )

            assert all(isinstance(path, Path) for path in paths)
            assert all(path.exists() for path in paths)
            assert [path.name for path in paths] == list(
                ("TEST_FILE1.txt", "TEST_FILE2.txt")
            )

            assert downloader.db_path is None
            assert not (tmp_path / REGISTRY_NAME).exists()

    def test_multiple_downloads_create_database(self, tmp_path):

        with temporary_http_server() as (host, port, temp_folder):
            base_url = f"http://{host}:{port}/"
            downloader = StaticFetcher(cache_path=tmp_path)
            paths = downloader.fetch(
                base_url=base_url,
                file_name=("TEST_FILE1.txt", "TEST_FILE2.txt"),
                db_name=REGISTRY_NAME,
            )

            assert downloader.db_path.read_text() == (
                "TEST_FILE1.txt sha256:c4bdb6ba200a917b8384ffeffa4999bf05bd4e479f6580d795aca509c9122dc4\n"
                "TEST_FILE2.txt sha256:0ec192c0f90d1332f2abca4398596d3978434ecbae6abea8ffd989412b592458\n"
            )

    def test_multiple_downloads_existing_database(self, tmp_path):
        with temporary_http_server() as (host, port, temp_folder):
            base_url = f"http://{host}:{port}/"
            downloader = StaticFetcher(cache_path=tmp_path)
            paths1 = downloader.fetch(
                base_url=base_url,
                file_name=("TEST_FILE1.txt", "TEST_FILE2.txt"),
                db_name=REGISTRY_NAME,
            )
            mtime1 = [path.stat().st_mtime for path in paths1]

            downloader = StaticFetcher(cache_path=tmp_path)
            paths2 = downloader.fetch(
                base_url=base_url,
                file_name=("TEST_FILE1.txt", "TEST_FILE2.txt"),
                db_name=REGISTRY_NAME,
            )

            mtime2 = [path.stat().st_mtime for path in paths2]

            assert mtime1 == mtime2
