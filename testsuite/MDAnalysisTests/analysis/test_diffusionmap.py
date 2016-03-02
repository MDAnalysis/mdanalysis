# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import print_function

import MDAnalysis

import MDAnalysis.analysis.align as align
import MDAnalysis.analysis.rms as rms
import MDAnalysis.lib.qcprot as qcp
import MDAnalysis.analysis.diffusionmap as diffusionmap
from MDAnalysis import SelectionError

from numpy.testing import (TestCase, dec,
                           assert_almost_equal, assert_raises, assert_equal)
import numpy as np
from nose.plugins.attrib import attr

import tempdir
from os import path

from MDAnalysisTests.datafiles import PSF, DCD, FASTA
from MDAnalysisTests import executable_not_found, parser_not_found


class TestDiffusionmap(object):
    def __init__(self):
        u = MDAnalysis.Universe(PSF,DCD)
        eg,ev=diffusionmap.diffusionmap(u)
        self.eg=eg
        self.ev=ev
       

    def test_eg(self):
        assert_equal(self.eg.shape, (98,))
        assert_almost_equal(self.eg[0], 1.0)
        assert_almost_equal(self.eg[-1],0.03174182)


    def test_ev(self):
        assert_equal(self.ev.shape, (98,98))
        assert_almost_equal(self.ev[0,0], .095836037343022831)



