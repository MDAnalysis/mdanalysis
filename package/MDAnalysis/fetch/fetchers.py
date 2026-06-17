# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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

import hashlib

from pathlib import Path
from abc import ABC, abstractmethod

try:
    import pooch
except ImportError:
    HAS_POOCH = False
else:
    HAS_POOCH = True

DEFAULT_CACHE_NAME_DOWNLOADER = "MDAnalysis_pdbs"
DEFAULT_TIMEOUT = 10
DEFAULT_RETRIES = 2


class BaseFetcher(ABC):
    """Blueprint Class for all Fetchers"""

    def __init__(
        self,
    ):
        pass

    @abstractmethod
    def fetch(self, base_url, verbose, timeout, retries):
        # Starts file retrieval workflow
        #
        # All fetchers should call _check_pooch()
        self._check_pooch()

        # Global Variable attributes
        self.base_url = base_url
        self.timeout = timeout  # timeout
        self.retries = retries  # number of retries

    def _check_pooch(
        self,
    ):
        # Note that requests is a major dependency of pooch and is guaranteed to be installed
        if not HAS_POOCH:
            raise ModuleNotFoundError(
                "pooch is needed as a dependency for Fetchers"
            )

    def _validate_fetch_args(self, args):
        """Checks to see if @abstractmethod fetch parameters are initalized correctly"""

        if "base_url" not in args:
            raise ValueError("base_url is not defined in fetch()")

        args.setdefault("progressbar", False)
        args.setdefault("timeout", DEFAULT_TIMEOUT)
        args.setdefault("retries", DEFAULT_RETRIES)

        return args


class StaticFetcher(BaseFetcher):
    """Fetcher automatically downloads file in entirety and cache it to disk"""

    def __init__(self, cache_path=None, hash="sha256"):

        self._check_pooch()

        self.cache_path = self._check_cache_path_input(cache_path)
        self.hash = self._check_hash_input(hash)
        self.db_path = None

    def fetch(
        self,
        file_name=None,
        verbose=False,
        db_name="hashes.txt",
        downloader="HTTP",
        **kwargs,
    ):
        kwargs = self._validate_fetch_args(kwargs)

        registry_dictionary = {}
        LOAD_FROM_CACHE = False
        CREATE_DATABASE = False

        if db_name is not None:
            self.db_path = self.cache_path / Path(db_name)

            if self.db_path.exists():
                LOAD_FROM_CACHE = True
            else:
                CREATE_DATABASE = True

        if LOAD_FROM_CACHE:
            # Reads pooch registry file format
            # https://www.fatiando.org/pooch/latest/registry-files.html#registry-file-format
            with open(self.db_path, mode="r") as f:
                for line in f:
                    key, value = line.strip().split()
                    registry_dictionary[key] = value

            # Adds files not in cache
            for file in file_name:
                if file not in registry_dictionary:
                    registry_dictionary[file] = None

        else:  # No Database (just download)
            # This block of code allows file_name to be a tuple instead of a string
            if isinstance(file_name, str):
                _file_name = (file_name,)
            else:
                _file_name = file_name

            registry_dictionary = {name: None for name in _file_name}

        main_downloader = pooch.create(
            path=self.cache_path,
            base_url=kwargs["base_url"],
            registry=registry_dictionary,
            retry_if_failed=kwargs["retries"],
        )

        download_kwargs = kwargs.copy()
        download_kwargs.pop("base_url")
        download_kwargs.pop("retries")

        match downloader:
            case "HTTP":
                fetch_downloader = pooch.HTTPDownloader(**download_kwargs)
            case "FTP":
                fetch_downloader = pooch.FTPDownloader(**download_kwargs)
            case "SFTP":
                fetch_downloader = pooch.SFTPDownloader(**download_kwargs)
            case "DOI":
                fetch_downloader = pooch.DOIDownloader(**download_kwargs)
            case _:
                raise ValueError(
                    f"Invalid downloader '{downloader}'. Valid options are 'HTTP', 'FTP', 'SFTP', 'DOI'."
                )

        paths = [
            Path(
                main_downloader.fetch(
                    fname=file_name,
                    progressbar=verbose,
                    downloader=fetch_downloader,
                )
            )
            for file_name in registry_dictionary.keys()
        ]

        if CREATE_DATABASE:
            pooch.make_registry(directory=self.cache_path, output=self.db_path)

        if len(paths) == 1:
            return paths[0]
        else:
            return paths

    def _check_cache_path_input(self, cache_path):
        if cache_path is None:
            return Path(pooch.os_cache(DEFAULT_CACHE_NAME_DOWNLOADER))
        else:
            return Path(cache_path)

    def _check_hash_input(self, hash):
        if hash in hashlib.algorithms_available:
            return hash
        else:
            raise ValueError(
                f'Invalid hash "{hash}". Valid hashes algorithms are {hashlib.algorithms_available}.'
            )


class DynamicFetcher(BaseFetcher):
    """Fetcher yields a Python Generator for dynamic downloading and analysis"""

    pass
