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
        self.db_path = None

    def fetch(
        self,
        base_url,
        file_name=None,
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
        # self.file_name = (
        #     file_name  # If not none, then user can change otherwise use default
        # )

        # self.override = force  # Boolean to override files (download despite being present)
        # self._ignore_hash = (
        #     ignore_hash  # If true, ignore hash and keep downloading
        # )
        # ##
        # ###

        self._check_pooch() 

        CREATE_DATABASE = False
        if db_name is not None: # HAS DATABASE
            self.db_path = self.cache_path / Path(db_name)
            if not self.db_path.exists():
                CREATE_DATABASE = True

        if CREATE_DATABASE: # Load a None registry dictionary (bc of no cache)
            print('creating')
            registry_dictionary = {
                file_name : None
            }

            registry=registry_dictionary

            downloader = pooch.create(
                path=self.cache_path,
                base_url=base_url,
                registry=registry_dictionary
            )

        else: # Loads from file
            print('loading')
            registry=open(self.db_path, mode='r')

            registry_dictionary = {}

            with open(self.db_path, mode='r') as f:
                for line in f:
                    key, value = line.strip().split()
                    registry_dictionary[key] = value

            downloader = pooch.create(
                path=self.cache_path,
                base_url=base_url,
                registry=registry_dictionary
            )
            
            import ipdb; ipdb.set_trace()


        
        
        paths = [
            Path(downloader.fetch(fname=file_name, progressbar=True))
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
            return pooch.os_cache(DEFAULT_CACHE_NAME_DOWNLOADER)
        else:
            return cache_path

    def _check_hash_input(self, hash):
        if hash in hashlib.algorithms_available:
            return hash
        else:
            raise ValueError(
                f"Invalid hash \"{hash}\". Valid hashes algorithms are {hashlib.algorithms_available}."
            )
    
class DynamicFetcher(BaseFetcher):
    """Fetcher yields a Python Generator for dynamic downloading and analysis"""

    pass
