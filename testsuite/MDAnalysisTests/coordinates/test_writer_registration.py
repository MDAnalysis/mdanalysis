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
from __future__ import absolute_import

import pytest

import MDAnalysis as mda
from MDAnalysis.coordinates.base import WriterBase


class TestWriterCreation(object):
    class MagicWriter(WriterBase):
        # this writer does the 'magic' format
        format = 'MAGIC'

        def __init__(self, filename, n_atoms=None):
            self.filename = filename
            self.n_atoms = n_atoms

    class MultiMagicWriter(MagicWriter):
        # this writer does the 'magic' and 'magic2' formats
        # but *only* supports multiframe writing.
        format = ['MAGIC', 'MAGIC2']
        multiframe = True
        singleframe = False

    def test_default_multiframe(self):
        assert isinstance(mda.Writer('this.magic'), self.MultiMagicWriter)

    def test_singleframe(self):
        # check that singleframe=False has been respected
        assert isinstance(mda.Writer('this.magic', multiframe=False), self.MagicWriter)

    def test_multiframe_magic2(self):
        # this will work as we go for multiframe
        assert isinstance(mda.Writer('that.magic2'), self.MultiMagicWriter)

    def test_singleframe_magic2(self):
        # this should fail, there isn't a singleframe writer for magic2
        with pytest.raises(TypeError):
            mda.Writer('that.magic2', multiframe=False)
