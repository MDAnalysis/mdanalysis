# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
from MDAnalysis.fetch.pdb import from_DOI
import hashlib

from MDAnalysis.fetch.fetchers import HAS_POOCH
from urllib import request

from pathlib import Path

try:
    request.urlopen("https://figshare.com", timeout=2)
    HAS_ACCESS_TO_FIGSHARE = True
except request.URLError:
    HAS_ACCESS_TO_FIGSHARE = False
try:
    request.urlopen("https://zenodo.org", timeout=2)
    HAS_ACCESS_TO_ZENODO = True
except request.URLError:
    HAS_ACCESS_TO_ZENODO = False


@pytest.mark.skipif(not HAS_POOCH, reason="Pooch is not installed.")
@pytest.mark.skipif(
    not HAS_ACCESS_TO_FIGSHARE,
    reason="Can not connect to https://figshare.com",
)
@pytest.mark.parametrize(
    "doi, file_name, expected_md5",
    [
        (
            "https://doi.org/10.6084/m9.figshare.5108170",
            "adk4AKE.psf",
            "c6f6b0711709cdfebe4620191b6d8e23",
        ),
        (
            "doi.org/10.6084/m9.figshare.5108170",
            "adk4AKE.psf",
            "c6f6b0711709cdfebe4620191b6d8e23",
        ),
    ],
)
def test_url_prefix(tmp_path, doi, file_name, expected_md5):
    p1 = from_DOI(
        doi,
        file_name=file_name,
        cache_path=tmp_path,
    )

    assert p1.exists()
    assert isinstance(p1, Path)
    assert hashlib.md5(p1.read_bytes()).hexdigest() == expected_md5


@pytest.mark.skipif(not HAS_POOCH, reason="Pooch is not installed.")
@pytest.mark.skipif(
    not HAS_ACCESS_TO_ZENODO,
    reason="Can not connect to https://zenodo.org",
)
@pytest.mark.parametrize(
    "doi, file_name, expected_md5",
    [
        (
            "https://doi.org/10.6084/m9.figshare.5108170",
            "adk4AKE.psf",
            "c6f6b0711709cdfebe4620191b6d8e23",
        ),
        (
            "https://doi.org/10.5281/zenodo.10028894",
            "MDAKits UGM.pdf",
            "8b2816dc971916dda3c77c0df65e833d",
        ),
    ],
)
def test_different_sources(tmp_path, doi, file_name, expected_md5):

    p1 = from_DOI(
        doi,
        file_name=file_name,
        cache_path=tmp_path,
    )

    assert p1.exists()
    assert isinstance(p1, Path)
    assert hashlib.md5(p1.read_bytes()).hexdigest() == expected_md5
