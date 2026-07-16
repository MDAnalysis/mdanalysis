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

This module contains fetcher classes that retrieve files from remote servers.
These classes use the third-party library :mod:`pooch` as a dependency.

Classes
-------

.. autoclass:: StaticFetcher
    :members:
    :inherited-members:

Variables
---------

These module-level variables affect runtime behavior across all fetcher classes.
Changing these values affects all fetchers.


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
#: See :func:`pooch.os_cache` for further details.
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
    """Base class for all fetchers.

    This class should not be initialized directly; fetcher implementations
    should inherit from it.

    """

    def __init__(
        self,
    ):
        pass

    @abstractmethod
    def fetch(self, base_url, verbose, timeout, retries):
        """Retrieve files from a remote server."""
        # Starts file retrieval workflow
        # All fetchers should call _check_pooch()
        #
        # These arguments should be implemented by all child Fetchers.
        pass

    def _check_pooch(
        self,
    ):
        """Raise an error if :mod:`pooch` is not installed."""
        if not HAS_POOCH:
            raise ModuleNotFoundError(
                "pooch is needed as a dependency for Fetchers"
            )

    def _validate_fetch_args(self, args):
        """Add default timeout and retry values to fetch arguments."""

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
        If the directory does not exist, it will be created.

    hash : str, optional
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

    hash : str
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
            The requested URL has the form ``{base_url}/{file_name}``.
        verbose : bool, optional
            If ``True``, show download progress. The default is ``False``.
        db_name : str or None, optional
            Name of the local hash database file used to verify cached downloads.
            The default is ``"hashes.txt"``. If ``None``, no registry database
            is read or written.
        append_db : bool, optional
            If ``True``, add downloaded files that are missing from an existing
            registry to that registry. If ``False``, missing registry entries
            raise a :class:`ValueError`. The default is ``False``.
        timeout : float, optional
            Time in seconds to wait for a response from the server before timing
            out. The default is :data:`DEFAULT_TIMEOUT`.
        retries : int, optional
            Number of times to retry a failed download. The default is
            :data:`DEFAULT_RETRIES`.
        downloader : str, optional
            Downloader backend to use. Supported values are ``"HTTP"``,
            ``"FTP"``, ``"SFTP"``, and ``"DOI"``. The default is ``"HTTP"``.

        Returns
        -------
        pathlib.Path or list of pathlib.Path
            The downloaded file path for a single file, or a list of paths for
            multiple files.

        Examples
        --------

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
        # Keywords arguments that are reserved for common _BaseFetcher.fetch() arguments.
        kwargs = self._validate_fetch_args(kwargs)

        LOAD_FROM_CACHE = False
        CREATE_DATABASE = False
        MISSING_FILES = False
        APPEND_DATABASE = append_db

        ## Reading from Registry
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
                "There are unknown files in the registry! The missing files are"
                + f" {missing_files_list}. To fix this, please set append_db=True"
                + " to append the database"
            )

        # Ensure fetch() only get requested files
        requested_files = (file_name,) if isinstance(file_name, str) else tuple(file_name)
        for name in requested_files:
            registry_dictionary.setdefault(name, None)
        
        ##

        ## Download code using pooch
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
                    f"Invalid downloader '{downloader}'. Valid options "
                    + "are 'HTTP', 'FTP', 'SFTP', 'DOI'."
                )
        
        paths = [
            Path(
                main_downloader.fetch(
                    fname=file_name,
                    progressbar=verbose,
                    downloader=fetch_downloader,
                )
            )
            for file_name in requested_files
        ]
        
        ##

        ## Registry write code
        if CREATE_DATABASE:
            self.write_registry(self.db_path, paths)

        if APPEND_DATABASE and LOAD_FROM_CACHE:
            self.fix_registry(self.db_path, registry_dictionary)
        
        ##

        return paths[0] if len(paths) == 1 else paths

    def fix_registry(self, db_path, file_dict):
        """
        Append newly downloaded files to an existing Pooch registry.

        Parameters
        ----------
        db_path : str or path-like
            Path to the registry file to update.
        file_dict : dict
            Dictionary mapping filenames to hash values. Files with a hash value of
            ``None`` are treated as newly downloaded files and appended to the
            registry.

        Returns
        -------
        None
            This method updates the registry file in place and does not return a
            value.

        Notes
        -----
        For each entry in ``file_dict`` with a value of ``None``, this method builds
        the corresponding file path relative to :attr:`cache_path` and appends its
        hash to the registry using :meth:`write_registry`.
        """
        new_files = [
            self.cache_path / file_name
            for file_name, file_hash in file_dict.items()
            if file_hash is None
        ]
        self.write_registry(db_path, new_files, mode="a")

    def check_registry(self, db_path):
        """
        Return cache files that are missing from the registry.

        Parameters
        ----------
        db_path : str or path-like
            Path to the registry file to read.

        Returns
        -------
        missing_files : list of pathlib.Path
            Cache file paths whose filenames are not present in the registry.
            The registry database file itself is excluded from the result.

        Notes
        -----
        This method compares filenames from the registry against files found
        recursively under :attr:`cache_path`. A cache file is considered missing
        from the registry when ``path.name`` is not a key in the registry
        dictionary.
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
        Read a Pooch registry file into a dictionary.

        Parameters
        ----------
        db_path : str or path-like
            Path to the registry file to read.

        Returns
        -------
        hash_dict : dict
            Dictionary mapping each filename in the registry to its stored hash
            value. Hash values are expected to include the hash algorithm prefix,
            for example ``"sha256:<digest>"``.

        Notes
        -----
        Each line in the registry file is expected to have the format::

            <filename> <hash_algorithm>:<digest>
        """
        hash_dict = {}

        with open(db_path, mode="r") as f:
            for line in f:
                key, value = line.strip().split()
                hash_dict[key] = value

        return hash_dict

    def write_registry(self, db_path, files, mode="w"):
        """
        Write a Pooch registry file with hashes for the given files.

        Parameters
        ----------
        db_path : str or path-like
            Path to the registry file to write.
        files : iterable of path-like
            Files to include in the registry. Each file must provide a ``name``
            attribute and be readable by ``pooch.file_hash``.
        mode : str, optional
            File opening mode used when writing the registry. Default is ``"w"``.

        Returns
        -------
        None
            This method writes the registry to disk and does not return a value.

        Notes
        -----
        Each registry line is written in the format::

            <filename> <hash_algorithm>:<digest>
        """
        with open(db_path, mode=mode) as f:
            for file in files:
                digest = pooch.file_hash(file, alg=self.hash)
                f.write(f"{file.name} {self.hash}:{digest}\n")

    # Argument validation methods
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
                f'Invalid hash "{hash}". Valid hashes algorithms'
                + f" are {hashlib.algorithms_available}."
            )


# class DynamicFetcher(_BaseFetcher):
#     """Fetcher yields a Python Generator for dynamic downloading and analysis"""

#     raise NotImplementedError
