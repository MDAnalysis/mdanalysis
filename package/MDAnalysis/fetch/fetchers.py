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
    import requests
except ImportError:
    HAS_POOCH = False
else:
    HAS_POOCH = True

DEFAULT_CACHE_NAME_DOWNLOADER = "MDAnalysis_pdbs"

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

        self.cache_path = self._check_cache_path_input(cache_path)
        self.hash = self._check_hash_input(hash)

    def fetch(
        self,
        base_url,
        filename=None,
        force=False,
        ignore_hash=False,
        db_name="hashes.txt",
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
            HAS_DATABASE = True
            self.db_path = self.cache_path / Path(db_name)
        else:
            HAS_DATABASE = False
            self.db_path = None

        if HAS_DATABASE:
            #self.db_path.parent.mkdir(parents=True, exist_ok=True)
            #self._file_extension = self.db_path.suffix

            if not self.db_path.exists():
                CREATE_DATABASE = True
            else:
                CREATE_DATABASE = False


            CREATE_DATABASE = True
            if CREATE_DATABASE: # Load a None registry dictionary
                registry_dictionary = {
                    '1AKE.pdb': None
                }

            else: # Loads from file
                #registry_dictionary = pooch.Pooch.load_registry(fname = (self.cache_path / 'test.txt'))
                pass

            import ipdb; ipdb.set_trace()


        
            ## Should stilll be ok
            downloader = pooch.create(
                path=self.cache_path,
                base_url=base_url,
                registry=registry_dictionary,
            )

            downloader.load_registry(fname = (self.cache_path / 'test.txt'))
            paths = [
                Path(downloader.fetch(fname=file_name, progressbar=True))
                for file_name in registry_dictionary.keys()
            ]
            
            print(self.cache_path)

            ## Add guard block here ro make it work
            pooch.make_registry(directory=self.cache_path, output=(self.cache_path / 'test.txt'))
            if len(paths) == 1:
                return paths[0]
            else:
                return paths

            ## SAVE to registry

            

        





            
        ## TODO  Workflow
        ## Prequel: Start Connection Pooling with Server 
        #
        #
        # 1. Check filename or get file name (content-deposition) via HTTP GET
        # 2. Check against database:
        # 2a. Write database if doesn't exist (override with db_name=None)
        # 2b. Check against database -> (_read_cache):
            # Check header for file_name and hash (ONLY SUPPORT ONE TYPE OF HASH per DB FILE for maintainability sake)
            # If mismatch with hash, toss exception (override with ignore_hash -- PUT BIG WARNING IN THIS)
            # If matchs, skip download and just return pathlib.Path() (override with force)
            # If empty, contuine with download and write hash to database (_write_cache())

            # Note replace with pooch instead

        # 


    def _check_cache_path_input(self, cache_path):

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


    def _create_database(self):
   
        self.db_path.parent.mkdir(parents=True,exist_ok=True)
        ## CSV
        with self.db_path.open(mode='x') as f:
            writer = csv.writer(f)
            writer.writerow(['File', f'Hash:{self.hash}'])        


    def _read_database(self, db_path):
        # Check and loads hash (either a csv or database file)
        pass


    def _write_database(self, db_path):
        # Create/query hash file (either a csv or database file)
        #
        # Not using Pooch.make_registry as that implements MD5 checksum which is not secure!

        pass

    



class DynamicFetcher(BaseFetcher):
    """Fetcher yields a Python Generator for dynamic downloading and analysis"""

    pass
