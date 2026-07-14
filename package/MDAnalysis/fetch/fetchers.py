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

"""
Fetchers --- :mod:`MDAnalysis.fetch.fetchers`
=============================================

This module contains the Fetchers classes that can be used to retrieve or fetch files
from remote servers. These classes uses the third party library :mod:`pooch` as
a dependency.

Classes
-------

.. autoclass:: StaticFetcher
    :members:
    :inherited-members:

Variables
---------

These are global submodule level variables that affect the runtime behavior across
all Fetcher Classes. Changing these values will affect all Fetchers!


.. autodata:: DEFAULT_CACHE_NAME_DOWNLOADER
.. autodata:: DEFAULT_TIMEOUT
.. autodata:: DEFAULT_RETRIES

"""

import hashlib

from pathlib import Path
from abc import ABC, abstractmethod

try:
    import pooch
except ImportError:
    HAS_POOCH = False
else:
    HAS_POOCH = True


#: Name of the :mod:`pooch` cache directory
#: ``pooch.os_cache(DEFAULT_CACHE_NAME_DOWNLOADER)``;
#:
#: see :func:`pooch.os_cache` for further details.'
#:
#: .. versionadded:: 2.11.0
DEFAULT_CACHE_NAME_DOWNLOADER = "MDAnalysis_pdbs"

#: Default time in seconds to wait for a response from the server before timing out.
#:
#: .. versionadded:: 2.11.0
DEFAULT_TIMEOUT = 10

#: Default number of attempts to retry a download if it fails.
#:
#: .. versionadded:: 2.11.0
DEFAULT_RETRIES = 2


class _BaseFetcher(ABC):
    """Blueprint Class for all Fetchers

    This shouldn't be initalized directly but should be inherited by other
    Fetchers classes.

    """

    def __init__(
        self,
    ):
        pass

    @abstractmethod
    def fetch(self, base_url, verbose, timeout, retries):
        # Starts file retrieval workflow
        # All fetchers should call _check_pooch()
        #
        # These arguments should be implemented by all child Fetchers.
        pass

    def _check_pooch(
        self,
    ):
        if not HAS_POOCH:
            raise ModuleNotFoundError(
                "pooch is needed as a dependency for Fetchers"
            )

    def _validate_fetch_args(self, args):
        """This initalized the fetcher library variables"""

        args.setdefault("timeout", DEFAULT_TIMEOUT)
        args.setdefault("retries", DEFAULT_RETRIES)

        return args


class StaticFetcher(_BaseFetcher):
    """
    Downloads files from a URL to disk and caches them to a local directory.

    Parameters
    ----------
    cache_path : str or pathlib.Path, optional
        Path to the cache directory. If set to ``None``, the default cache
        directory will be used as specified by
        :data:`DEFAULT_CACHE_NAME_DOWNLOADER`.

        If the directory does not exist, it will attempted to be created.

    hash : str
        Hash algorithm to use for verifying the integrity of downloaded files.
        The default is ``sha256``. Valid options are any hash algorithm
        available in the :mod:`hashlib` module.


    Attributes
    ----------
    cache_path : pathlib.Path or ``None``
        Path to the cache directory.

    db_path : pathlib.Path or ``None``
        Path to the database file used for caching. Created after calling
        fetch().

    hash : str or ``None``
        Hash algorithm used for verifying the integrity of downloaded files.

    Notes
    -----
    The download directory can be overridden by setting the environment
    variable ``MDANALYSIS_FETCHER_DATA`` to a valid path. This class uses
    :mod:`pooch` as a backend for downloading and caching files.

    """

    def __init__(self, cache_path=None, hash="sha256"):

        self._check_pooch()

        self.cache_path = self._check_cache_path_input(cache_path)
        self.hash = self._check_hash_input(hash)
        self.db_path = None

    def fetch(
        self,
        base_url,
        file_name,
        verbose=False,
        db_name="hashes.txt",
        append_db=False,
        downloader="HTTP",
        **kwargs,
    ):
        """
        Download one or more files from a static base URL and cache them
        locally.

        Parameters
        ----------
        base_url : str
            Base URL from which to download the file(s). This should be a valid
            URL pointing to the directory containing the files to be downloaded.
        file_name : str or sequence of str
            Name of the file or files to download.
            Note that the request is phrased as {base_url}/{file_name}.
        verbose : bool, optional
            If True, shows fetcher progress.
            Default is False.
        db_name : str, optional
            Name of the local hash database file used to verify cached downloads.
            Default is "hashes.txt". If None, no registry database is read or written.
        timeout : float, optional
            Time in seconds to wait for a response from the server before timing out.
            Default is :data:`DEFAULT_TIMEOUT`.
        retries : int, optional
            Number of attempts to retry a download if it fails. Default is :data:`DEFAULT_RETRIES`.
        downloader : str, optional
            Downloader backend to use. If a string is provided, it must identify
            a supported downloader such as "HTTP".
            Default is "HTTP".

        Returns
        -------
        pathlib.Path or list of pathlib.Path
            The downloaded file path for a single file, or a list of paths for
            multiple files.

        Example
        -------

        A script to download a protein from the RCSB Protein Data Bank.

        .. code-block:: python

            from MDAnalysis.fetch import StaticFetcher

            fetcher = StaticFetcher(cache_path=cache_path)

            # Download a single file from the RCSB Protein Data Bank
            path = fetcher.fetch(
                file_name="1AKE.cif",
                base_url="https://files.wwpdb.org/download/",
            )

            # Download multiple files from the RCSB Protein Data Bank
            path = fetcher.fetch(
                file_name=["1AKE.cif", "4AKE.cif"],
                base_url="https://files.wwpdb.org/download/",
            )

        Notes
        -----
        The download directory can be overridden by setting the environment
        variable ``MDANALYSIS_FETCHER_DATA`` to a valid path. This class uses
        :mod:`pooch` as a backend for downloading and caching files. The
        cache database is created on demand when ``db_name`` does not
        exist.

        """
        # Keywords arguments are reserved for common _BaseFetcher.fetch() arguements.
        kwargs = self._validate_fetch_args(kwargs)

        LOAD_FROM_CACHE = False
        CREATE_DATABASE = False
        MISSING_FILES = False
        APPEND_DATABASE = append_db

        registry_dictionary = {}

        if db_name is not None:
            self.db_path = self.cache_path / Path(db_name)

            if self.db_path.exists():
                LOAD_FROM_CACHE = True
            else:
                CREATE_DATABASE = True

        if LOAD_FROM_CACHE:
            registry_dictionary = self.read_registry(self.db_path)
            missing_files_list = self.check_registry(self.db_path)

            if len(missing_files_list) != 0:
                MISSING_FILES = True

        if MISSING_FILES and not APPEND_DATABASE:
            raise ValueError(
                f"There are unknown files in the registry! The missing files are {missing_files_list}. To fix this, please set append_db=True to append the database"
            )

        # Code to process non-registry files
        # One-liner that forces strings into tuple
        no_db_files = (file_name,) if isinstance(file_name, str) else file_name
        for file in no_db_files:
            if file not in registry_dictionary:
                registry_dictionary[file] = None

        ## Pooch setup
        main_downloader = pooch.create(
            path=self.cache_path,
            base_url=base_url,
            registry=registry_dictionary,
            retry_if_failed=kwargs["retries"],
            env="MDANALYSIS_FETCHER_DATA",
        )

        download_kwargs = kwargs.copy()
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
        ##

        if CREATE_DATABASE:
            self.write_registry(self.db_path, paths)

        if APPEND_DATABASE and LOAD_FROM_CACHE:
            self.fix_registry(self.db_path, registry_dictionary)

        return paths[0] if len(paths) == 1 else paths

    # Reads pooch registry file format
    # https://www.fatiando.org/pooch/latest/registry-files.html#registry-file-format

    def fix_registry(self, db_path, file_dict):
        """Append newly downloaded files to an existing registry."""

        new_files = [
            self.cache_path / file_name
            for file_name, file_hash in file_dict.items()
            if file_hash is None
        ]
        self.write_registry(db_path, new_files, mode="a")




    def check_registry(self, db_path):
        """
        Return cache files that are missing from the registry database.

        Reads the registry at ``db_path`` and compares its recorded filenames
        against the files currently present in ``self.cache_path``. The registry
        database file itself is ignored.

        Args:
            db_path: Path to the registry database to read.

        Returns:
            list[pathlib.Path]: A list of cache file paths whose filenames are not
            present in the registry.
        """

        registry_dictionary = self.read_registry(db_path)
        database_files = set(registry_dictionary.keys())

        cache_files = (
            path
            for path in self.cache_path.rglob("*")
            if path != self.db_path and path.is_file()
        )

        return [
            path for path in cache_files if path.name not in database_files
        ]

    def read_registry(self, db_path):
        """
        Read a registry file into a dictionary of filenames and hashes.

        Each line in the registry file is expected to contain a filename and its
        corresponding hash value, separated by whitespace.

        Args:
            db_path: Path to the registry file to read.

        Returns:
            dict[str, str]: A dictionary where each key is a filename and each value
            is the file's stored hash.
        """

        hash_dict = {}

        with open(db_path, mode="r") as f:
            for line in f:
                key, value = line.strip().split()
                hash_dict[key] = value

        return hash_dict

    def write_registry(self, db_path, files, mode="w"):
        """Method to exclusively write pooch registry"""

        with open(db_path, mode=mode) as f:
            for file in files:
                digest = pooch.file_hash(file, alg=self.hash)
                f.write(f"{file.name} {self.hash}:{digest}\n")

    ### Arugment Validation Methods
    def _check_cache_path_input(self, cache_path):
        if cache_path is None:
            path = Path(pooch.os_cache(DEFAULT_CACHE_NAME_DOWNLOADER))
        else:
            path = Path(cache_path)

        Path(path).mkdir(parents=True, exist_ok=True)
        return path

    def _check_hash_input(self, hash):
        if hash in hashlib.algorithms_available:
            return hash
        else:
            raise ValueError(
                f'Invalid hash "{hash}". Valid hashes algorithms are {hashlib.algorithms_available}.'
            )


# class DynamicFetcher(_BaseFetcher):
#     """Fetcher yields a Python Generator for dynamic downloading and analysis"""

#     raise NotImplementedError
