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

import re
from pathlib import Path
from shutil import rmtree
from unittest.mock import Mock

from servers import temporary_http_server

import hashlib
import pytest


from MDAnalysis.fetch.fetchers import (
    DEFAULT_CACHE_NAME_DOWNLOADER,
    HAS_POOCH,
    StaticFetcher,
    pooch,
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

    def test_invalid_downloader(self, tmp_path):
        with temporary_http_server() as (host, port, temp_folder):
            base_url = f"http://{host}:{port}/"
            downloader = StaticFetcher(cache_path=tmp_path)

            with pytest.raises(
                ValueError,
                match=re.escape(
                    "Invalid downloader 'barfoo'. Valid options "
                    + "are ('auto', 'http', 'https', 'ftp', 'sftp', 'doi')"
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

    def test_append_error(self, tmp_path):

        with temporary_http_server() as (host, port, temp_folder):
            base_url = f"http://{host}:{port}/"
            fetcher = StaticFetcher(cache_path=tmp_path)

            file1 = fetcher.fetch(
                base_url=base_url,
                file_name="TEST_FILE1.txt",
                db_name=REGISTRY_NAME,
            )

            # match later
            with pytest.raises(ValueError):
                file2 = fetcher.fetch(
                    base_url=base_url,
                    file_name="TEST_FILE2.txt",
                    db_name=REGISTRY_NAME,
                )

    def test_invalid_auto_downloader(self, tmp_path):

        base_url = "foo"
        with pytest.raises(ValueError):

            fetcher = StaticFetcher(
                cache_path=tmp_path,
            )
            fetcher.fetch(
                base_url=base_url, file_name="bar", downloader="auto"
            )

    def test_invalid_manual_downloader(self, tmp_path):

        base_url = "foo"
        with pytest.raises(ValueError):

            fetcher = StaticFetcher(
                cache_path=tmp_path,
            )
            fetcher.fetch(base_url=base_url, file_name="bar", downloader="FOO")


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
            assert path.read_text() == "Sally sells seashells by the seashore"
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
            assert path.read_text() == "Sally sells seashells by the seashore"
            assert path.exists()

            assert Path(downloader.cache_path / REGISTRY_NAME).exists()
            assert (
                downloader.cache_path / REGISTRY_NAME
            ).read_text() == "TEST_FILE1.txt sha256:a625aaf4ca5e2d358b216165cee3247a93a40e699bb864193499d230ab7aad7e\n"

    def test_append_database(self, tmp_path):
        with temporary_http_server() as (host, port, temp_folder):
            base_url = f"http://{host}:{port}/"
            downloader = StaticFetcher(cache_path=tmp_path)
            path = downloader.fetch(
                base_url=base_url,
                file_name="TEST_FILE1.txt",
                db_name=REGISTRY_NAME,
            )

            path = downloader.fetch(
                base_url=base_url,
                file_name="TEST_FILE2.txt",
                db_name=REGISTRY_NAME,
                append_db=True,
            )

        assert Path(downloader.cache_path / REGISTRY_NAME).read_text() == (
            "TEST_FILE1.txt sha256:a625aaf4ca5e2d358b216165cee3247a93a40e699bb864193499d230ab7aad7e\n"
            "TEST_FILE2.txt sha256:eff7c015c379263afdf464bd1baf266909d0e4d4af7cccb722dd4994ff4e998c\n"
        )

    # Add paramertization
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
                Path(downloader.cache_path / REGISTRY_NAME).read_text()
                == "TEST_FILE1.txt md5:adf1020ce2ffe073600990e1c9c72ce8\n"
            )

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

            assert not Path(downloader.cache_path / REGISTRY_NAME).exists()
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

            assert not Path(downloader.cache_path / REGISTRY_NAME).exists()
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

            assert Path(downloader.cache_path / REGISTRY_NAME).read_text() == (
                "TEST_FILE1.txt sha256:a625aaf4ca5e2d358b216165cee3247a93a40e699bb864193499d230ab7aad7e\n"
                "TEST_FILE2.txt sha256:eff7c015c379263afdf464bd1baf266909d0e4d4af7cccb722dd4994ff4e998c\n"
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

    def test_environment_variable_override(self, tmp_path, monkeypatch):
        monkeypatch.setenv("MDANALYSIS_FETCHER_DATA", str(tmp_path))

        with temporary_http_server() as (host, port, temp_folder):
            base_url = f"http://{host}:{port}/"
            downloader = StaticFetcher()

            path = downloader.fetch(
                base_url=base_url,
                file_name="TEST_FILE1.txt",
            )

        assert Path(tmp_path / "TEST_FILE1.txt").exists()


@pytest.mark.skipif(not HAS_POOCH, reason="Pooch is not installed.")
class TestRegistry:
    def test_append_registry(self, tmp_path):

        with temporary_http_server() as (host, port, temp_folder):
            base_url = f"http://{host}:{port}/"
            fetcher = StaticFetcher(cache_path=tmp_path)

            file1 = fetcher.fetch(
                base_url=base_url,
                file_name="TEST_FILE1.txt",
                db_name="db_hash1.txt",
            )

            file2 = fetcher.fetch(
                base_url=base_url,
                file_name="TEST_FILE2.txt",
                db_name="db_hash2.txt",
            )

            registry = file1.parent / "db_hash1.txt"
            fetcher.append_registry(registry, ["TEST_FILE2.txt"])

            assert Path(registry).read_text() == (
                "TEST_FILE1.txt sha256:a625aaf4ca5e2d358b216165cee3247a93a40e699bb864193499d230ab7aad7e\n"
                "TEST_FILE2.txt sha256:eff7c015c379263afdf464bd1baf266909d0e4d4af7cccb722dd4994ff4e998c\n"
            )

    def test_write_registry(self, tmp_path):
        # This can't call StaticFetcher directly for an effective test
        # Maybe refactor the file creation into a function handle or fixture
        files = {
            "file1.txt": "Molecular \n",
            "file2.txt": "Dynamics. \n",
        }
        for filename, content in files.items():
            with open(tmp_path / filename, "w") as f:
                f.write(content)

        fetcher = StaticFetcher(cache_path=tmp_path)

        fetcher.write_registry(
            tmp_path / "file_1_hash.txt", [tmp_path / "file1.txt"]
        )
        assert (
            tmp_path / "file_1_hash.txt"
        ).read_text() == "file1.txt sha256:2da169c5aae36a823c202da49fb11935b76277efcb5cd42a4cf238ddda2a9b20\n"

        fetcher.write_registry(
            tmp_path / "file_2_hash.txt", [tmp_path / "file2.txt"]
        )
        assert (
            tmp_path / "file_2_hash.txt"
        ).read_text() == "file2.txt sha256:3a0dbd9e2abc4a7bbae6adfe92e2858218135926dacd4a7d3fb4ca2dbdbe457a\n"

        fetcher.write_registry(
            tmp_path / "file_1_and_2_hash.txt",
            [tmp_path / "file1.txt", tmp_path / "file2.txt"],
        )
        assert (tmp_path / "file_1_and_2_hash.txt").read_text() == (
            "file1.txt sha256:2da169c5aae36a823c202da49fb11935b76277efcb5cd42a4cf238ddda2a9b20\n"
            "file2.txt sha256:3a0dbd9e2abc4a7bbae6adfe92e2858218135926dacd4a7d3fb4ca2dbdbe457a\n"
        )

    ## Test works, but doesn't work on github action. Idk why and need to find out at some point.
    def test_check_registry(self, tmp_path):

        files = {
            "file1.txt": "Molecular \\n",
            "file2.txt": "Dynamics. \\n",
            "file3.txt": "Analysis. \\n",
        }

        for filename, content in files.items():
            with open(tmp_path / filename, "w") as f:
                f.write(content)

        fetcher = StaticFetcher(cache_path=tmp_path)

        # Write only file1
        fetcher.write_registry(
            tmp_path / "file_1_2_and_3_hash.txt",
            [tmp_path / "file1.txt"],
        )

        # Show that file2 and file3 are missing
        assert sorted(
            fetcher.check_registry(tmp_path / "file_1_2_and_3_hash.txt")
        ) == sorted(
            [
                tmp_path / "file3.txt",
                tmp_path / "file2.txt",
            ]
        )

    def test_check_registry_ignore(self, tmp_path):

        files = {
            "file1.txt": "Molecular \\n",
            "file2.txt": "Dynamics. \\n",
            "file3.txt": "Analysis. \\n",
        }

        for filename, content in files.items():
            with open(tmp_path / filename, "w") as f:
                f.write(content)

        fetcher = StaticFetcher(cache_path=tmp_path)

        # Write only file1
        fetcher.write_registry(
            tmp_path / "file_1_2_and_3_hash.txt",
            [tmp_path / "file1.txt"],
        )

        assert fetcher.check_registry(
            tmp_path / "file_1_2_and_3_hash.txt",
            ignore=[tmp_path / "file2.txt"],
        ) == [
            tmp_path / "file3.txt",
        ]


@pytest.mark.skipif(not HAS_POOCH, reason="Pooch is not installed.")
class TestDownloaders:

    @pytest.mark.parametrize(
        ("downloader", "constructor_name"),
        [
            ("ftp", "FTPDownloader"),
            ("sftp", "SFTPDownloader"),
            ("doi", "DOIDownloader"),
        ],
    )
    def test_monkeypatch_other_downloaders(
        self,
        monkeypatch,
        tmp_path,
        downloader,
        constructor_name,
    ):
        downloaded_path = tmp_path / "example.dat"

        # Monkeypatches the pooch.create to return a fake path
        pooch_fetch = Mock(return_value=str(downloaded_path))
        pooch_instance = Mock(fetch=pooch_fetch)

        monkeypatch.setattr(
            pooch,
            "create",
            Mock(return_value=pooch_instance),
        )

        # Monkekpatch a generic pooch Downloader
        # This will be used to mock the behavior of a generic pooch Downloader
        downloader_instance = object()
        downloader_constructor = Mock(return_value=downloader_instance)

        monkeypatch.setattr(
            pooch,
            constructor_name,
            downloader_constructor,
        )

        # Since pooch.create() is monkeypatched to return a fake path. This won't raise an exception
        # And the fake downloader will be passed into the fetch method to mock its behavior.
        result = StaticFetcher(cache_path=tmp_path).fetch(
            base_url="https://example.com/",
            file_name="example.dat",
            db_name=None,
            downloader=downloader,
        )

        assert result == downloaded_path

    def test_downloader_error(self, tmp_path):
        with pytest.raises(
            ValueError,
        ):

            with temporary_http_server() as (host, port, temp_folder):
                base_url = f"http://{host}:{port}/"
                downloader = StaticFetcher()

                downloader = StaticFetcher(
                    cache_path=tmp_path,
                )
                path = downloader.fetch(
                    base_url=base_url,
                    file_name="TEST_FILE1.txt",
                    downloader="FOO",
                )
