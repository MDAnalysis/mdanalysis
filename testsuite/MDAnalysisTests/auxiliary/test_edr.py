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

import pytest
import numpy as np

import os

import MDAnalysis as mda

from MDAnalysisTests.datafiles import (AUX_EDR,
                                       AUX_EDR_TPR,
                                       AUX_EDR_XTC,
                                       AUX_EDR_RAW)
from MDAnalysisTests.auxiliary.base import (BaseAuxReaderTest,
                                            BaseAuxReference)


def read_auxstep_data(step):
    """parses raw edr data (plain text) and returns dictionary with auxdata for
    step"""
    with open(AUX_EDR_RAW) as f:
        rawdata = f.readlines()
    n_entries = 52  # number of aux terms per step
    stepdata = rawdata[step * n_entries: (step + 1) * n_entries]
    aux_dict = {}
    for line in stepdata:
        # .join() necessary for entries with multi-word keys
        aux_dict[" ".join(line.split()[:-1])] = float(line.split()[-1])
    return(aux_dict)


class EDRReference(BaseAuxReference):
    def __init__(self):
        super(EDRReference, self).__init__()
        self.testdata = AUX_EDR
        self.reader = mda.auxiliary.EDR.EDRReader
        self.n_steps = 4
        self.dt = 0.02

        def reference_auxstep(i):
            # create a reference AuxStep for step i
            auxstep = mda.auxiliary.base.AuxStep(dt=self.dt,
                                                initial_time=self.initial_time)
            auxstep.step = i
            auxstep._data = read_auxstep_data(i)
            return auxstep

        self.auxsteps = [reference_auxstep(i) for i in range(self.n_steps)]

        # add the auxdata and format for .xvg to the reference description
        self.description['auxdata'] = os.path.abspath(self.testdata)
        self.description['format'] = self.reader.format

        # for testing the selection of data/time
        self.time_selector = "Time"
        self.select_time_ref = [step._data[self.time_selector]
                                for step in self.auxsteps]
        self.data_selector = None  # selects all data
        self.select_data_ref = [step._data for step in self.auxsteps]

        # testing __getitem__ with slice and list. Should allow us to iterate
        # through the specified auxiliary steps...
        self.iter_list = [0, -2]
        self.iter_list_auxsteps = [self.auxsteps[0], self.auxsteps[2]]
        self.iter_slice = slice(None, None, 2)  # every second step
        self.iter_slice_auxsteps = [self.auxsteps[0], self.auxsteps[2]]


class TestEDRReader(BaseAuxReaderTest):
    @staticmethod
    @pytest.fixture()
    def ref():
        return EDRReference()

    @staticmethod
    @pytest.fixture
    def ref_universe(ref):
        u = mda.Universe(AUX_EDR_TPR, AUX_EDR_XTC)
        u.trajectory.add_auxiliary('test', ref.testdata)
        return u

    @staticmethod
    @pytest.fixture()
    def reader(ref):
        return ref.reader(
            ref.testdata,
            initial_time=ref.initial_time,
            dt=ref.dt, auxname=ref.name,
            time_selector=None,
            data_selector=None
        )
