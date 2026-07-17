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
import re

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
        """Add the default global timeout and retry variables to fetch arguments."""

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

    def fetch(
        self,
        base_url,
        file_name,
        verbose=False,
        db_name="hashes.txt",
        append_db=False,
        downloader="auto",
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
            Downloader backend to use. Supported values are ``"auto"``, ``"http"``,
            ``"ftp"``, ``"sftp"``, and ``"doi"``. The default is ``"auto"``.

        Returns
        -------
        pathlib.Path or list of pathlib.Path
            The downloaded file path for a single file, or a list of paths for
            multiple files.

        Examples
        --------
        Download a single CIF file from the RCSB Protein Data Bank.

        >>> StaticFetcher.fetch(file_name="1AKE.cif",
            base_url="https://files.wwpdb.org/download/")
        './MDAnalysis_pdbs/1AKE.cif'

        Download multiple CIF files from the RCSB Protein Data Bank.
        >>> StaticFetcher.fetch(file_name=["1AKE.cif", "4AKE.cif"],
            base_url="https://files.wwpdb.org/download/")
        ['./MDAnalysis_pdbs/1AKE.cif', './MDAnalysis_pdbs/4AKE.cif']

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

        registry_dictionary = {}

        # Process file names
        requested_files = (
            (file_name,) if isinstance(file_name, str) else tuple(file_name)
        )
        requested_files_abs = [
            self.cache_path / name for name in requested_files
        ]

        ## Reading from Registry
        if db_name is not None:
            db_path = self.cache_path / Path(db_name)

            if db_path.exists():
                LOAD_FROM_CACHE = True
            else:
                CREATE_DATABASE = True

        if LOAD_FROM_CACHE:
            registry_dictionary = self.read_registry(db_path)
            missing_files_list = self.check_registry(
                db_path, files=list(requested_files_abs)
            )

            if len(missing_files_list) != 0:
                MISSING_FILES = True

        if MISSING_FILES and not APPEND_DATABASE:
            raise ValueError(
                "fetch() is requesting files not found in the registry. "
                + f"The missing files are {missing_files_list}. "
                + "To fix this, please set append_db=True to append the registry."
            )

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
        fetch_downloader = self._set_downloader(
            base_url, downloader, **download_kwargs
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
            self.write_registry(db_path, paths)

        if APPEND_DATABASE and LOAD_FROM_CACHE:
            self.append_registry(db_path, requested_files)

        ##

        return paths[0] if len(paths) == 1 else paths

    def append_registry(self, db_path, files):
        """
        Append cached files to an existing Pooch registry.

        Each entry in ``files`` is resolved relative to :attr:`cache_path`. The
        file hash is computed using the fetcher's configured hash algorithm :attr:`hash` and a
        new registry line is appended to ``db_path``.

        Parameters
        ----------
        db_path : str or path-like
            Path to the registry file to update.
        files : iterable of str or path-like
            File names or paths for cached files to append to the registry. Relative
            paths are interpreted relative to :attr:`cache_path`.

        Returns
        -------
        None
            This method updates the registry file in place and does not return a
            value.

        Example
        -------
        >>> from MDAnalysis.fetch.fetchers import StaticFetcher
        >>> fetcher = StaticFetcher()
        >>> file1 = fetcher.fetch(
        ...     file_name="1AKE.cif",
        ...     base_url="https://files.wwpdb.org/download/",
        ...     db_name="db_hash1.txt",
        ... )
        >>> file2 = fetcher.fetch(
        ...     file_name="4AKE.cif",
        ...     base_url="https://files.wwpdb.org/download/",
        ...     db_name="db_hash2.txt",
        ... )
        >>> registry = file1.parent / "db_hash1.txt"
        >>> fetcher.append_registry(registry, ["4AKE.cif"])
        >>> registry.read_text()
        1AKE.cif sha256:01f41b1b42318a1a5df7f650dbab881677aa0e8d825f7c42dd26ae16a94c0948
        4AKE.cif sha256:fcb2ff49a3e255797fee277ce28e0acace67f6e6ddf432841f8451f00cbde9e9

        Notes
        -----
        Existing registry entries are preserved. This method does not check for or
        remove duplicate file entries.

        Each appended registry line has the format::

            <filename> <hash_algorithm>:<digest>

        """
        new_files = [self.cache_path / file_name for file_name in files]
        self.write_registry(Path(db_path), new_files, mode="a")

    def check_registry(self, db_path, files=[], ignore=[]):
        """
        Return cache files that are missing from the registry.

        This method compares filenames within the registry against files found
        recursively under :attr:`cache_path`. A cache file is considered missing
        when it is recorded in the registry, but not physically present on disk.

        Parameters
        ----------
        db_path : str or path-like
            Path to the registry file to read.
        ignore : list of str or path-like
            Files to be ignored

        Returns
        -------
        missing_files : list of pathlib.Path
            Cache file paths whose filenames are not present in the registry.
            The registry database file itself is excluded from the result.

        Example
        -------
        .. code-block:: python
        >>> files = {
        ...     "file1.txt": "Molecular \\n",
        ...     "file2.txt": "Dynamics. \\n",
        ...     "file3.txt": "Analysis. \\n"
        ... }
        >>> for filename, content in files.items():
        ...     with open(tmp_path / filename, "w") as f:
        ...         f.write(content)
        ...
        >>> fetcher = StaticFetcher(cache_path=tmp_path)
        >>> fetcher.write_registry(
        ...     "file_1_2_and_3_hash.txt",
        ...     files=["file1.txt"],
        ... )
        >>> fetcher.check_registry("file_1_2_and_3_hash.txt")
        [Path('./MDAnalysis_pdbs/file3.txt'), Path('./MDAnalysis_pdbs/file2.txt')]
        >>> fetcher.check_registry("file_1_2_and_3_hash.txt", ignore=["file2.txt"])
        [Path('./MDAnalysis_pdbs/file3.txt')]


        Notes
        -----
        Each line in the registry file is expected to have the format::

            <filename> <hash_algorithm>:<digest>


        """
        registry_dictionary = self.read_registry(db_path)
        database_files = set(registry_dictionary.keys())

        cache_files = [
            path
            for path in self.cache_path.rglob("*")
            if path != db_path and path.is_file()
        ] + files

        return [
            path
            for path in cache_files
            if (path.name not in database_files) and (path not in ignore)
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

        Example
        -------
        .. code-block:: python

        >>> files = {
        ...     "file1.txt": "Molecular \\n",
        ...     "file2.txt": "Dynamics. \\n",
        ... }
        >>> for filename, content in files.items():
        ...     with open(tmp_path / filename, "w") as f:
        ...         f.write(content)
        ...
        >>> fetcher = StaticFetcher(cache_path=tmp_path)
        >>> fetcher.write_registry(
        ...     "file_1_and_2_hash.txt",
        ...     ["file1.txt", "file2.txt"],
        ... )
        >>> fetcher.read_registry("file_1_and_2_hash.txt")
        {'file1.txt': 'sha256:2da169c5aae36a823c202da49fb11935b76277efcb5cd42a4cf238ddda2a9b20',
        'file2.txt': 'sha256:3a0dbd9e2abc4a7bbae6adfe92e2858218135926dacd4a7d3fb4ca2dbdbe457a'}
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
        files : iterable of str or path-like
            Files to include in the registry. Each file must provide a ``name``
            attribute and be readable by ``pooch.file_hash``.
        mode : str, optional
            File opening mode used when writing the registry. Default is ``"w"``.

        Returns
        -------
        None
            This method writes the registry to disk and does not return a value.

        Example
        -------
        .. code-block:: python
        >>> files = {
        ...     "file1.txt": "Molecular \\n",
        ...     "file2.txt": "Dynamics. \\n",
        ... }
        >>> for filename, content in files.items():
        ...     with open(tmp_path / filename, "w") as f:
        ...         f.write(content)
        ...
        >>> fetcher = StaticFetcher(cache_path=tmp_path)
        >>> fetcher.write_registry(
        ...     "file_1_and_2_hash.txt",
        ...     ["file1.txt", "file2.txt"],
        ... )
        >>> Path("file_1_and_2_hash.txt").read_text()
        file1.txt sha256:2da169c5aae36a823c202da49fb11935b76277efcb5cd42a4cf238ddda2a9b20
        file2.txt sha256:3a0dbd9e2abc4a7bbae6adfe92e2858218135926dacd4a7d3fb4ca2dbdbe457a

        Notes
        -----
        Each registry line is written in the format::

            <filename> <hash_algorithm>:<digest>
        """
        with open(db_path, mode=mode) as f:
            for file in files:
                file = Path(file)
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

    def _set_downloader(self, base_url, downloader, **kwargs):
        """Sets Downloader in fetch() by matching a regex against the download link"""

        SUPPORTED_DOWNLOADERS = ("auto", "http", "https", "ftp", "sftp", "doi")

        if downloader not in SUPPORTED_DOWNLOADERS:
            raise ValueError(
                f"Invalid downloader '{downloader}'. Valid options "
                + f"are {SUPPORTED_DOWNLOADERS}"
            )

        if downloader == "auto":
            # The regex below is AI generated, but the overall idea of using regular expressions
            # to check the first bit of the url was thought up by me.
            #
            # I thought of these four examples and prompt AI to come with a regex that capture the
            # susbtring before "://""
            #
            # doi://10.6084/m9.figshare.14763051.v1/tiny-data.txt
            # https://www.example.com/page
            # ftp://ftp.example.com/files/document.txt
            # sftp://username@example.com/path/to/folder

            regex = r"^([^:]+)://"
            match = re.match(regex, base_url)

            if match:
                _downloader = match.group(1)
        else:
            _downloader = downloader

        match _downloader:
            case "http" | "https":
                return pooch.HTTPDownloader(**kwargs)
            case "ftp":
                return pooch.FTPDownloader(**kwargs)
            case "sftp":
                return pooch.SFTPDownloader(**kwargs)
            case "doi":
                return pooch.DOIDownloader(**kwargs)


# class DynamicFetcher(_BaseFetcher):
#     """Fetcher yields a Python Generator for dynamic downloading and analysis"""

#     raise NotImplementedError
