# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
#

"""
AlphaFold Fetchers --- :mod:`MDAnalysis.fetch.alphafold`
========================================================

This suite of functions download predicted structure files from the AlphaFold
Protein Structure Database.

.. _AlphaFold DB: https://alphafold.ebi.ac.uk/

Functions
---------

.. autofunction:: from_AlphaFold

"""
from pathlib import Path

try:
    import pooch
except ImportError:
    HAS_POOCH = False
else:
    HAS_POOCH = True

DEFAULT_CACHE_NAME_DOWNLOADER = "MDAnalysis_alphafold"

SUPPORTED_FILE_FORMATS_DOWNLOADER = ("pdb", "cif")


def from_AlphaFold(
    uniprot_ids,
    cache_path=None,
    progressbar=False,
    file_format="pdb",
):
    """
    Download one or more AlphaFold predicted structures and cache them locally.

    Given one or multiple UniProt IDs, downloads the corresponding structure files
    and stores them in a local cache directory.

    Parameters
    ----------
    uniprot_ids : str or sequence of str
        A single UniProt ID or a sequence of IDs to fetch.
    cache_path : str or pathlib.Path
        Directory where downloaded file(s) will be cached.
    file_format : str
        The file extension/format to download ("pdb" or "cif").
    progressbar : bool
        If True, display a progress bar during file downloads.

    Returns
    -------
    :class:`~pathlib.Path` or list of :class:`~pathlib.Path`
        The path(s) to the downloaded file(s).
    """

    if not HAS_POOCH:
        raise ModuleNotFoundError(
            "pooch is needed as a dependency for from_AlphaFold()"
        )
    elif file_format not in SUPPORTED_FILE_FORMATS_DOWNLOADER:
        raise ValueError(
            f"Invalid file format. Supported: {SUPPORTED_FILE_FORMATS_DOWNLOADER}"
        )

    if isinstance(uniprot_ids, str):
        _ids = (uniprot_ids,)
    else:
        _ids = uniprot_ids

    if cache_path is None:
        cache_path = pooch.os_cache(DEFAULT_CACHE_NAME_DOWNLOADER)

    registry_dictionary = {
        f"AF-{uid}-F1-model_v4.{file_format}": None for uid in _ids
    }

    downloader = pooch.create(
        path=cache_path,
        base_url="https://alphafold.ebi.ac.uk/files/",
        registry=registry_dictionary,
    )

    paths = [
        Path(downloader.fetch(fname=file_name, progressbar=progressbar))
        for file_name in registry_dictionary.keys()
    ]

    return paths if not isinstance(uniprot_ids, str) else paths[0]
