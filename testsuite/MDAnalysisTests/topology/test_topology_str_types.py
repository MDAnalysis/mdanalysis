# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
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
from __future__ import division, absolute_import


from six import string_types
import MDAnalysis
import pytest

from MDAnalysis.tests.datafiles import (
    CRD,
    LAMMPSdata,
    DLP_CONFIG_minimal,
    GRO,
    TPR,
    DMS,
    GMS_SYMOPT,
    MMTF,
    mol2_molecule,
    PRM7,
    PDB_small,
    PDBQT_input,
    PQR,
    PRM,
    PSF,
    PRM12,
    HoomdXMLdata,
    XPDB_small,
    XYZ_mini,
    DLP_HISTORY_minimal, )


@pytest.mark.parametrize('prop', [
    'name',
    'resname',
    'type',
    'segid',
    'moltype',
])
# topology formats curated from values available in
# MDAnalysis._PARSERS
@pytest.mark.parametrize( 'top_format, top', [
    ('CONFIG', DLP_CONFIG_minimal),
    ('CRD', CRD),
    ('DATA', LAMMPSdata),
    ('DMS', DMS),
    ('GMS', GMS_SYMOPT),
    ('GRO', GRO),
    ('HISTORY', DLP_HISTORY_minimal),
    ('MMTF', MMTF),
    ('MOL2', mol2_molecule),
    ('PARM7', PRM7),
    ('PDB', PDB_small),
    ('PDBQT', PDBQT_input),
    ('PQR', PQR),
    ('PRMTOP', PRM),
    ('PSF', PSF),
    ('TOP', PRM12),
    ('TPR', TPR),
    ('XML', HoomdXMLdata),
    ('XPDB', XPDB_small),
    ('XYZ', XYZ_mini)
])
def test_str_types(top_format, top, prop):
    # Python 2/3 topology string type checking
    # Related to Issue #1336
    u = MDAnalysis.Universe(top, format=top_format)
    if hasattr(u.atoms[0], prop):
        assert isinstance(getattr(u.atoms[0], prop), string_types)
