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
from numpy.testing import assert_almost_equal

from pathlib import Path

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

# The EDRReader behaves differently from the XVGReader, creating dummy test
# data similar to what is done for XVG is not possible. The EDRReference and
# some tests in TestEDRReader had to be changed.


class EDRReference(BaseAuxReference):
    def __init__(self):
        super(EDRReference, self).__init__()
        self.testdata = AUX_EDR
        self.reader = mda.auxiliary.EDR.EDRReader
        self.n_steps = 4
        self.dt = 0.02
        self.description = {'dt': self.dt, 'represent_ts_as': 'closest',
                            'initial_time': self.initial_time,
                            'time_selector': "Time", 'data_selector': None,
                            'constant_dt': True, 'cutoff': -1,
                            'auxname': self.name}

        def reference_auxstep(i):
            # create a reference AuxStep for step i
            t_init = self.initial_time
            auxstep = mda.auxiliary.base.AuxStep(dt=self.dt,
                                                 initial_time=t_init)
            auxstep.step = i
            auxstep._data = read_auxstep_data(i)
            return auxstep

        self.auxsteps = [reference_auxstep(i) for i in range(self.n_steps)]

        # add the auxdata and format for .xvg to the reference description
        self.description['auxdata'] = Path(self.testdata).resolve()
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

        def reference_timestep(dt=0.02, offset=0):
            # return a trajectory timestep with specified dt, offset + move to
            # frame 1; for use in auxiliary reading of different timesteps
            ts = mda.coordinates.base.Timestep(0, dt=dt,
                                               time_offset=offset)
            ts.frame = 1
            return ts

        self.bond_step1_val = 1426.2252197265625
        self.bond_step2_val = 1482.0098876953125
        self.bond_step1_2_avg = (self.bond_step1_val + self.bond_step2_val) / 2
        # test reading a timestep with lower frequency. Auxiliary steps with
        # times between [0.02 ps, 0.06 ps) will be assigned to this timestep,
        # i.e. step 1 (0.02 ps) and step 2 (0.04 ps).
        self.lower_freq_ts = reference_timestep(dt=0.04, offset=0)
        # 'closest' representative value will match step 2
        # Pick value for "Bond" to check against.
        self.lowf_closest_rep = self.bond_step2_val
        # 'average' representative value
        self.lowf_avg_time = 0.03
        self.lowf_average_rep = [self.lowf_avg_time, self.bond_step1_2_avg]

        # test reading a timestep with higher frequency. Auxiliart steps with
        # times between [0.25ps, 0.75ps) will be assigned to this timestep,
        # i.e. no auxiliary steps
        self.higher_freq_ts = reference_timestep(dt=0.01, offset=0)
        self.highf_rep = self.format_data(np.nan)

        # test reading a timestep that is offset. Auxiliary steps with
        # times between [0.015ps, 0.035ps) will be assigned to this timestep,
        # i.e step 1 (0.02 ps)
        self.offset_ts = reference_timestep(dt=0.02, offset=0.005)
        # 'closest' representative value will match step 1 data
        self.offset_closest_rep = self.bond_step1_val

        # testing cutoff for representative values
        self.cutoff = 0
        # for 'average': use low frequenct timestep, only step 2 within 0ps
        # cutoff
        self.lowf_cutoff_average_rep = self.bond_step2_val
        # for 'closest': use offset timestep; no timestep within 0ps cutoff
        # changed here to match requirements for EDR data
        self.offset_cutoff_closest_rep = np.array(np.nan)


class TestEDRReader(BaseAuxReaderTest):
    @staticmethod
    @pytest.fixture
    def ref():
        return EDRReference()

    @staticmethod
    @pytest.fixture
    def ref_universe(ref):
        u = mda.Universe(AUX_EDR_TPR, AUX_EDR_XTC)
        u.trajectory.add_auxiliary('test', ref.testdata, "Bond")
        return u

    @staticmethod
    @pytest.fixture
    def reader(ref):
        return ref.reader(
            ref.testdata,
            initial_time=ref.initial_time,
            dt=ref.dt, auxname=ref.name,
            time_selector="Time",
            data_selector=None
        )

    def test_time_non_constant_dt(self, reader):
        reader.constant_dt = False
        reader.time_selector = None
        with pytest.raises(ValueError, match="If dt is not constant, "
                                             "must have a valid time "
                                             "selector"):
            reader.time

    def test_iterate_through_trajectory(self, ref, ref_universe):
        # check the representative values of aux for each frame are as expected
        # trajectory here has same dt, offset; so there's a direct
        # correspondence between frames and steps
        for i, ts in enumerate(ref_universe.trajectory):
            assert_almost_equal(ts.aux.test, ref.auxsteps[i].data["Bond"])

    def test_iterate_as_auxiliary_from_trajectory(self, ref, ref_universe):
        # check representative values of aux for each frame are as expected
        # trajectory here has same dt, offset, so there's a direct
        # correspondence between frames and steps, and iter_as_aux will run
        # through all frames
        for i, ts in enumerate(ref_universe.trajectory.iter_as_aux('test')):
            assert_almost_equal(ts.aux.test, ref.auxsteps[i].data["Bond"])

    def test_step_to_frame_time_diff(self, reader, ref):
        # Timestep is 0.002 (10 %) longer than auxiliary data
        ts = mda.coordinates.base.Timestep(0, dt=ref.dt + 0.002)

        # Test all 4 frames
        for idx in range(4):

            frame, time_diff = reader.step_to_frame(idx, ts,
                                                    return_time_diff=True)

            assert frame == idx
            np.testing.assert_almost_equal(time_diff, idx * 0.002)

    def test_read_lower_freq_timestep(self, ref, reader):
        # test reading a timestep with lower frequency
        ts = ref.lower_freq_ts
        reader.update_ts(ts)
        # check the value set in ts is as we expect
        assert_almost_equal(ts.aux.test["Bond"], ref.lowf_closest_rep,
                            err_msg="Representative value in ts.aux "
                                    "does not match")

    def test_read_higher_freq_timestep(self, ref, reader):
        # try reading a timestep with higher frequency
        ts = ref.higher_freq_ts
        reader.update_ts(ts)
        assert_almost_equal(ts.aux.test, ref.highf_rep,
                            err_msg="Representative value in ts.aux "
                                    "does not match")

    def test_read_offset_timestep(self, ref, reader):
        # try reading a timestep offset from auxiliary
        ts = ref.offset_ts
        reader.update_ts(ts)
        assert_almost_equal(ts.aux.test["Bond"], ref.offset_closest_rep,
                            err_msg="Representative value in ts.aux "
                                    "does not match")

    def test_represent_as_average(self, ref, reader):
        # test the 'average' option for 'represent_ts_as'
        # reset the represent method to 'average'...
        reader.represent_ts_as = 'average'
        # read timestep; use the low freq timestep
        ts = ref.lower_freq_ts
        reader.update_ts(ts)
        # check the representative value set in ts is as expected
        test_value = [ts.aux.test["Time"], ts.aux.test["Bond"]]
        assert_almost_equal(test_value, ref.lowf_average_rep,
                            err_msg="Representative value does not match when "
                                    "using with option 'average'")

    def test_represent_as_average_with_cutoff(self, ref, reader):
        # test the 'represent_ts_as' 'average' option when we have a cutoff set
        # set the cutoff...
        reader.cutoff = ref.cutoff
        # read timestep; use the low frequency timestep
        ts = ref.lower_freq_ts
        reader.update_ts(ts)
        # check representative value set in ts is as expected
        assert_almost_equal(ts.aux.test["Bond"], ref.lowf_cutoff_average_rep,
                            err_msg="Representative value does not match when "
                                    "applying cutoff")

    def test_add_all_terms_from_file(self, ref, ref_universe):
        ref_universe.trajectory.add_auxiliary("*", ref.testdata)
        # adding "test" manually to match above addition of test term
        ref_terms = ["test"] + [key for key in read_auxstep_data(0).keys()]
        terms = [key for key in ref_universe.trajectory._auxs]
        assert ref_terms == terms

    def test_add_all_terms_from_reader(self, ref_universe, reader):
        ref_universe.trajectory.add_auxiliary("*", reader)
        ref_terms = ["test"] + [key for key in read_auxstep_data(0).keys()]
        terms = [key for key in ref_universe.trajectory._auxs]
        assert ref_terms == terms

    def test_add_term_list_custom_names_from_file(self, ref, ref_universe):
        ref_universe.trajectory.add_auxiliary(["bond", "temp"], ref.testdata,
                                              ["Bond", "Temperature"])
        ref_dict = read_auxstep_data(0)
        assert ref_universe.trajectory.ts.aux.bond == ref_dict["Bond"]
        assert ref_universe.trajectory.ts.aux.temp == ref_dict["Temperature"]

    def test_add_term_list_custom_names_from_reader(self, ref_universe,
                                                    reader):
        ref_universe.trajectory.add_auxiliary(["bond", "temp"], reader,
                                              ["Bond", "Temperature"])
        ref_dict = read_auxstep_data(0)
        assert ref_universe.trajectory.ts.aux.bond == ref_dict["Bond"]
        assert ref_universe.trajectory.ts.aux.temp == ref_dict["Temperature"]

    def test_add_single_term_custom_name_from_file(self, ref, ref_universe):
        ref_universe.trajectory.add_auxiliary("temp", ref.testdata,
                                              "Temperature")
        ref_dict = read_auxstep_data(0)
        assert ref_universe.trajectory.ts.aux.temp == ref_dict["Temperature"]

    def test_add_single_term_custom_name_from_reader(self, ref_universe,
                                                     reader):
        ref_universe.trajectory.add_auxiliary("temp", reader, "Temperature")
        ref_dict = read_auxstep_data(0)
        assert ref_universe.trajectory.ts.aux.temp == ref_dict["Temperature"]

    def test_terms_update_on_iter(self, ref_universe, reader):
        ref_universe.trajectory.add_auxiliary(["bond", "temp"], reader,
                                              ["Bond", "Temperature"])
        ref_dict = read_auxstep_data(0)
        assert ref_universe.trajectory.ts.aux.bond == ref_dict["Bond"]
        assert ref_universe.trajectory.ts.aux.temp == ref_dict["Temperature"]
        ref_dict = read_auxstep_data(1)
        ref_universe.trajectory.next()
        assert ref_universe.trajectory.ts.aux.bond == ref_dict["Bond"]
        assert ref_universe.trajectory.ts.aux.temp == ref_dict["Temperature"]

    def test_invalid_data_selector(self, ref, ref_universe):
        with pytest.raises(KeyError, match="'Nonsense' is not a key"):
            ref_universe.trajectory.add_auxiliary("Nonsense", AUX_EDR)

    def test_read_all_times(self, reader):
        all_times_expected = np.array([0., 0.02, 0.04, 0.06])
        assert np.all(all_times_expected == reader.read_all_times())
