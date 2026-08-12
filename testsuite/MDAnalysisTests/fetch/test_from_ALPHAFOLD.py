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


import pytest
import hashlib
import re

from MDAnalysis.fetch.pdb import from_ALPHAFOLD
from MDAnalysis.fetch.fetchers import HAS_POOCH
from MDAnalysis.fetch.pdb import _SUPPORTED_FILE_FORMATS_ALPHAFOLD
from urllib import request
from pathlib import Path

try:
    request.urlopen("https://alphafold.ebi.ac.uk/", timeout=2)
    HAS_ACCESS_TO_ALPHAFOLD = True
except request.URLError:
    HAS_ACCESS_TO_ALPHAFOLD = False


@pytest.mark.skipif(
    not HAS_ACCESS_TO_ALPHAFOLD,
    reason="Can not connect to Alphafold",
)
@pytest.mark.skipif(not HAS_POOCH, reason="Pooch is not installed.")
@pytest.mark.parametrize(
    "id, expected_sha256",
    [
        (
            "Q9I1F6",
            "378a7aaa35d7a63fbb3cecef6f46483e2037863e83baf9a8b07a0ca223c947d1",
        ),
        (
            "P04637",
            "4332bff49301d04434b6e3b43c439aa5bd8238760e88c7142f7d0c0851a68e5d",
        ),
    ],
)
def test_download_one_file_default_settings(tmp_path, id, expected_sha256):
    p1 = from_ALPHAFOLD(
        id=id,
        cache_path=tmp_path,
    )

    assert p1.exists()
    assert p1.suffix == ".cif"
    assert isinstance(p1, Path)
    assert hashlib.sha256(p1.read_bytes()).hexdigest() == expected_sha256


@pytest.mark.skipif(
    not HAS_ACCESS_TO_ALPHAFOLD,
    reason="Can not connect to Alphafold",
)
@pytest.mark.skipif(not HAS_POOCH, reason="Pooch is not installed.")
@pytest.mark.parametrize("file_format", ["bcif", "cif", "pdb"])
def test_different_format(tmp_path, file_format):

    p1 = from_ALPHAFOLD(
        id="Q9I1F6", cache_path=tmp_path, file_format=file_format
    )

    assert p1.suffix == f".{file_format}"


@pytest.mark.skipif(
    not HAS_ACCESS_TO_ALPHAFOLD,
    reason="Can not connect to Alphafold",
)
@pytest.mark.skipif(not HAS_POOCH, reason="Pooch is not installed.")
def test_invalid_format(tmp_path):
    with pytest.raises(
        ValueError,
        match=re.escape(
            f"Invalid file format: boo. "
            + f"Supported formats are: {list(_SUPPORTED_FILE_FORMATS_ALPHAFOLD.keys())}"
        ),
    ):
        from_ALPHAFOLD(id="Q9I1F6", cache_path=tmp_path, file_format="boo")
