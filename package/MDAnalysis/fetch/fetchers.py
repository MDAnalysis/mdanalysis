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

import csv
import sqlite3
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
ALLOWED_EXTENSIONS_DATABASE = {".csv", ".db"}

class BaseFetcher(ABC):
    """Blueprint Class for all Fetchers"""

    def __init__(
        self,
        reuse_connection=True,
    ):

        self.pooled = reuse_connection # Connection Pooling

  
    @abstractmethod
    def fetch(self, base_url, progressbar, timeout, retries):
        # Starts file retrieval workflow
        #
        # All fetchers should call _check_pooch()
        self._check_pooch()

        # Global Variable attributes  
        self.base_url = base_url
        self.verbose = progressbar  # Progressbar
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


class StaticFetcher(BaseFetcher):
    """Fetcher automatically downloads file in entirety and cache it to disk"""

    def __init__(self, cache_path=None, hash="sha256", **kwargs):

        ## TODO put guard parameter from ABC Fetcher
        #super().__init__(kwargs["reuse_connection"])
        super().__init__()

        self.cache_path = self._set_cache_path(cache_path)
        self.hash = self._check_hash_input(hash)
        self.hash = self._check_hash_input(hash)


    def fetch(
        self,
        filename=None,
        force=False,
        ignore_hash=False,
        db_name="hashes.db",
        **kwargs,
    ):


        ###
        # ## All variable to used for methods
        # # ABC Fetcher variables (guaranteed to exist)
        # self.base_url = kwargs["base_url"]
        # self.verbose = kwargs["progressbar"]  # Progressbar
        # self.timeout = kwargs["timeout"]  # timeout
        # self.retries = kwargs["retries"]  # number of retries

        # # Static Fetcher Specific variables
        # self.filename = (
        #     filename  # If not none, then user can change otherwise use default
        # )

        # self.override = force  # Boolean to override files (download despite being present)
        # self._ignore_hash = (
        #     ignore_hash  # If true, ignore hash and keep downloading
        # )
        # ##
        # ###


        ## Pseudocode
        self._check_pooch() # Check dependencies

        
        if db_name is not None:

            db_name = Path(db_name)
            
            if db_name.suffix.lower() not in ALLOWED_EXTENSIONS_DATABASE:
                raise ValueError(
                    f"Database name should have one of these extensions: {ALLOWED_EXTENSIONS_DATABASE}"
            )

        

        ## TODO  Workflow
        ## Prequel: Start Connection Pooling with Server 
        #
        #
        # 1. Check filename or get file name (content-deposition) via HTTP GET
        # 2. Check against database:
        # 2a. Write database if doesn't exist (cancel with db_name=None)
        # 2b. Check against database -> (_read_cache):
            # Check header for file_name and hash (ONLY SUPPORT ONE TYPE OF HASH per DB FILE for maintainability sake)
            # If mismatch with hash, toss exception (override with ignore_hash -- PUT BIG WARNING IN THIS)
            # If matchs, skip download and just return pathlib.Path() (override with force)
            # If empty, contuine with download and write hash to database (_write_cache())

        # 



            


    def _set_cache_path(self, cache_path):

        if cache_path is None:
            return pooch.os_cache(DEFAULT_CACHE_NAME_DOWNLOADER)
        else:
            return cache_path

    def _check_hash_input(self, hash):
        if hash in hashlib.algorithms_available:
            return hash
        else:
            raise ValueError(
                f"Invalid hash \"{hash}\". Valid hashes algorithms are {hashlib.algorithms_available}. See 'hashlib.algorithms_available'"
            )

    def _write_cache(self, db_path):
        # Create/query hash file (either a csv or database file)
        #
        # Not using Pooch.make_registry as that implements MD5 checksum which is not secure!

        pass

    def _read_cache(self, db_path):
        # Check and loads hash (either a csv or database file)

        breakpoint()
        db_extension = db_path.suffix.lower()

        if db_extension == ".csv":
            with open(db_extension, newline='') as csvfile:
                file = csv.reader(csvfile, dialect='unix')



        elif db_extension == ".db":
            pass



class DynamicFetcher(BaseFetcher):
    """Fetcher yields a Python Generator for dynamic downloading and analysis"""

    pass
