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

import MDAnalysis as mda
from MDAnalysis.auxiliary.base import AuxReader

import pytest


def test_get_auxreader_for_none():
    with pytest.raises(ValueError, match="Must provide either auxdata or format"):
        mda.auxiliary.core.get_auxreader_for()


def test_get_auxreader_for_wrong_auxdata():
    with pytest.raises(ValueError, match="Unknown auxiliary data format for auxdata:"):
        mda.auxiliary.core.get_auxreader_for(auxdata="test.none")


def test_get_auxreader_for_wrong_format():
    with pytest.raises(ValueError, match="Unknown auxiliary data format"):
        mda.auxiliary.core.get_auxreader_for(format="none")


def test_notimplemented_read_next_timestep():
    with pytest.raises(NotImplementedError, match="BUG: Override "
                       "_read_next_step()"):
        reader = mda.auxiliary.base.AuxReader()


class No_Memory_Usage(AuxReader):
    def _read_next_step(self):
        pass


def test_notimplemented_memory_usage():
    with pytest.raises(NotImplementedError, match="BUG: Override "
                       "_memory_usage()"):
        reader = No_Memory_Usage()
        reader._memory_usage()
