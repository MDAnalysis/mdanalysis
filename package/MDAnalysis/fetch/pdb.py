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
PDB Fetchers --- :mod:`MDAnalysis.fetch.pdb`
============================================

This suite of functions download structure files from the Research Collaboratory for
Structural Bioinformatics (RCSB) `Protein Data Batabank`_ (PDB).

.. _Protein Data Batabank: https://www.rcsb.org/

Variables
---------

.. autodata:: DEFAULT_CACHE_NAME_DOWNLOADER


Functions
---------

.. autofunction:: from_PDB

"""
from pathlib import Path

try:
    import pooch
except ImportError:
    HAS_POOCH = False
else:
    HAS_POOCH = True

#: Name of the :mod:`pooch` cache directory ``pooch.os_cache(DEFAULT_CACHE_NAME_DOWNLOADER)``;
#: see :func:`pooch.os_cache` for further details.
#:
#: .. versionadded:: 2.11.0
DEFAULT_CACHE_NAME_DOWNLOADER = "MDAnalysis_pdbs"

# These file formats are here https://www.rcsb.org/docs/programmatic-access/file-download-services#pdb-entry-files"
SUPPORTED_FILE_FORMATS_DOWNLOADER = (
    "cif",
    "cif.gz",
    "bcif",
    "bcif.gz",
    "xml",
    "xml.gz",
    "pdb",
    "pdb.gz",
    "pdb1",
    "pdb1.gz",
)


def from_PDB(
    pdb_ids,
    cache_path=None,
    progressbar=False,
    file_format="cif.gz",
):
    """
    Download one or more PDB files from the RCSB Protein Data Bank and cache
    them locally.

    Given one or multiple PDB IDs, downloads the corresponding structure files
    format and stores them in a local cache directory. If files are cached on
    disk, *from_PDB* will skip the download and use the cached version instead.

    Returns the path(s) as a :class:`~pathlib.Path` to the downloaded file(s).

    Parameters
    ----------
    pdb_ids : str or sequence of str
        A single PDB ID as a string, or a sequence of PDB IDs to fetch.
    cache_path : str or pathlib.Path
        Directory where downloaded file(s) will be cached.
        The default ``None`` argument uses the :mod:`pooch` default cache with
        project name :data:`DEFAULT_CACHE_NAME_DOWNLOADER`.
    file_format : str
        The file extension/format to download (e.g., "cif", "pdb").
        See the Notes section below for a list of all supported file formats.
    progressbar : bool
        If True, display a progress bar during file downloads. Default is False.

    Returns
    -------
    :class:`~pathlib.Path` or list of :class:`~pathlib.Path`
        The path(s) to the downloaded file(s). Returns a single
        :class:`~pathlib.Path` if a single pdb id is given, or a list of
        :class:`~pathlib.Path` if multiple pdb ids are provided.

    Raises
    ------
    ValueError
        For an invalid file format. Supported file formats are under Notes.

    :class:`requests.exceptions.HTTPError`
        If an invalid PDB code is specified.

    Notes
    -----
    This function uses the `RCSB File Download Services`_ for directly downloading
    structure files via https.

    .. _`RCSB File Download Services`:
       https://www.rcsb.org/docs/programmatic-access/file-download-services

    The RCSB currently provides data in ``'cif'`` , ``'cif.gz'`` , ``'bcif'`` ,
    ``'bcif.gz'`` , ``'xml'`` , ``'xml.gz'`` , ``'pdb'`` , ``'pdb.gz'``,
    ``'pdb1'``, ``'pdb1.gz'`` file formats and can therefore be downloaded.
    Not all of these formats can be currently read with MDAnalysis.

    Caching, controlled by the `cache_path` parameter, is handled internally by
    :mod:`pooch`. The default cache name is taken from
    :data:`DEFAULT_CACHE_NAME_DOWNLOADER`. To clear cache (and subsequently force
    re-fetching), it is required to delete the cache folder as specified by
    `cache_path`.

    Examples
    --------
    Download a single PDB file:

    >>> import MDAnalysis as mda
    >>> from MDAnalysis.fetch import from_PDB
    >>> from_PDB("1AKE", file_format="cif")
    './MDAnalysis_pdbs/1AKE.cif'

    Download multiple PDB files with a progress bar:

    >>> mda.fetch.from_PDB(["1AKE", "4BWZ"], progressbar=True)
    ['./MDAnalysis_pdbs/1AKE.pdb.gz', './MDAnalysis_pdbs/4BWZ.pdb.gz']

    Download a single PDB file and convert it to a universe:

    >>> mda.Universe(mda.fetch.from_PDB("1AKE"), file_format="pdb.gz")
    <Universe with 3816 atoms>

    Download multiple PDB files and convert each of them into a universe:

    >>> [mda.Universe(pdb) for pdb in mda.fetch.from_PDB(["1AKE", "4BWZ"], progressbar=True)]
    [<Universe with 3816 atoms>, <Universe with 2824 atoms>]


    .. versionadded:: 2.11.0
    """

    if not HAS_POOCH:
        raise ModuleNotFoundError(
            "pooch is needed as a dependency for from_PDB()"
        )
    elif file_format not in SUPPORTED_FILE_FORMATS_DOWNLOADER:
        raise ValueError(
            "Invalid file format. Supported file formats "
            f"are {SUPPORTED_FILE_FORMATS_DOWNLOADER}"
        )

    if isinstance(pdb_ids, str):
        _pdb_ids = (pdb_ids,)
    else:
        _pdb_ids = pdb_ids

    if cache_path is None:
        cache_path = pooch.os_cache(DEFAULT_CACHE_NAME_DOWNLOADER)

    # Have to do this dictionary approach instead of using pooch.retrieve in order
    # to prevent the hardcoded known_hash warning from showing up.
    registry_dictionary = {
        f"{pdb_id}.{file_format}": None for pdb_id in _pdb_ids
    }

    downloader = pooch.create(
        path=cache_path,
        base_url="https://files.wwpdb.org/download/",
        registry=registry_dictionary,
    )

    paths = [
        Path(downloader.fetch(fname=file_name, progressbar=progressbar))
        for file_name in registry_dictionary.keys()
    ]

    return paths if not isinstance(pdb_ids, str) else paths[0]
