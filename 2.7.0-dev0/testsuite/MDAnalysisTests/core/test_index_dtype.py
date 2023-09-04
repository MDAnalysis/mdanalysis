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

"""32 bit compat tests

Tests for making sure that integer arrays used for indexing use `np.intp`.
This dtype is important for platform independent indexing of other arrays.

"""
import numpy as np
import pytest
from MDAnalysisTests import make_Universe


@pytest.fixture()
def u():
    return make_Universe()


def test_ag_ix(u):
    assert u.atoms.ix.dtype == np.intp


def test_rg_ix(u):
    assert u.residues.ix.dtype == np.intp


def test_sg_ix(u):
    assert u.segments.ix.dtype == np.intp


def test_atom_ix_array(u):
    assert u.atoms[0].ix_array.dtype == np.intp


def test_residue_ix_array(u):
    assert u.residues[0].ix_array.dtype == np.intp


def test_segment_ix_array(u):
    assert u.segments[0].ix_array.dtype == np.intp


def test_atomgroup_indices(u):
    assert u.atoms.indices.dtype == np.intp


def test_atomgroup_residue_upshift(u):
    assert u.atoms.resindices.dtype == np.intp


def test_atomgroup_segment_upshift(u):
    assert u.atoms.segindices.dtype == np.intp


def test_residuegroup_atom_downshift(u):
    # downshift arrays are a list (one for each residue)
    assert all((arr.dtype == np.intp)
               for arr in u.residues.indices)


def test_residuegroup_resindices(u):
    assert u.residues.resindices.dtype == np.intp


def test_residuegroup_segment_upshift(u):
    assert u.residues.segindices.dtype == np.intp


def test_segmentgroup_atom_downshift(u):
    assert all((arr.dtype == np.intp)
               for arr in u.segments.indices)


def test_segmentgroup_residue_downshift(u):
    assert all((arr.dtype == np.intp)
               for arr in u.segments.resindices)


def test_segmentgroup_segindices(u):
    assert u.segments.segindices.dtype == np.intp
