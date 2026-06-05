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
        base_url,
        progressbar,
        keep_session,
    ):

        self.base_url = base_url
        self.progressbar = progressbar  # Progressbar
        self.keep_session  # Connection Spooling

    @abstractmethod
    def fetch(self, timeout, retries):
        # Starts file retrieval workflow
        #
        # All fetchers should call _check_pooch

        self.timeout = timeout  # timeout
        self.retries = retries  # number of retries
        self._check_pooch()

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

    def __init__(self, cache_path):
        super().__init__(base_url, progressbar)
        self.cache_path = _set_cache_path(cache_path)

    def fetch(self):
        # Starts file retrieval workflow
        pass

    def _set_cache_path(self, cache_path):

        if cache_path is None:
            return pooch.os_cache(DEFAULT_CACHE_NAME_DOWNLOADER)
        else:
            return cache_path

    def _write_cache(self):
        # Create/query hash file (either a csv or database file)
        pass

    def _read_cache(self):
        # Check and loads hash (either a csv or database file)
        pass


class DynamicFetcher(BaseFetcher):
    """Fetcher yields a Python Generator for dynamic downloading and analysis"""

    pass
