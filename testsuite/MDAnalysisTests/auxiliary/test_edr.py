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
from numpy.testing import assert_allclose
import pickle

from pathlib import Path

import MDAnalysis as mda
from MDAnalysis import units
from MDAnalysis.auxiliary.EDR import HAS_PYEDR

from MDAnalysisTests.datafiles import (AUX_EDR,
                                       AUX_EDR_TPR,
                                       AUX_EDR_XTC,
                                       AUX_EDR_RAW,
                                       AUX_EDR_SINGLE_FRAME)
from MDAnalysisTests.auxiliary.base import (BaseAuxReaderTest,
                                            BaseAuxReference,
                                            assert_auxstep_equal)


def read_raw_data_file(step):
    """parses raw edr data (plain text) and returns dictionary with auxdata for
    step"""
    with open(AUX_EDR_RAW) as f:
        rawdata = f.readlines()
    n_entries = 52  # number of aux terms per step
    stepdata = rawdata[step * n_entries: (step + 1) * n_entries]
    aux_dict = {}
    edr_units = {}
    for line in stepdata:
        aux_dict[line.split(",")[0]] = float(line.split(",")[1])
        edr_units[line.split(",")[0]] = line.split(",")[2].strip()
    return aux_dict, edr_units


def get_auxstep_data(step):
    return read_raw_data_file(step)[0]


def get_edr_unit_dict(step):
    return read_raw_data_file(step)[1]


@pytest.mark.skipif(HAS_PYEDR, reason="pyedr present")
def test_pyedr_not_present_raises():
    with pytest.raises(ImportError, match="please install pyedr"):
        aux = mda.auxiliary.EDR.EDRReader(AUX_EDR)


@pytest.mark.skipif(not HAS_PYEDR, reason="pyedr not installed")
class EDRReference(BaseAuxReference):
    """
    Class to hold reference values for :class:`TestEDRReader` to test against.

    Because of the binary file format, it is not possible to create a test file
    of the format that is assumed for :class:`BaseAuxReference`. As such, some
    adaptations were necessary for the EDRReader.
    """

    def __init__(self):
        super(EDRReference, self).__init__()
        self.testdata = AUX_EDR
        self.reader = mda.auxiliary.EDR.EDRReader
        self.n_steps = 4
        self.dt = 0.02
        self.description = {'dt': self.dt, 'represent_ts_as': 'closest',
                            'initial_time': self.initial_time,
                            'time_selector': "Time", 'data_selector': None,
                            'constant_dt': True, 'cutoff': None,
                            'auxname': self.name}

        def reference_auxstep(i):
            # create a reference AuxStep for step i
            t_init = self.initial_time
            auxstep = mda.auxiliary.EDR.EDRStep(dt=self.dt,
                                                initial_time=t_init)
            auxstep.step = i
            auxstep._data = get_auxstep_data(i)
            return auxstep

        self.auxsteps = [reference_auxstep(i) for i in range(self.n_steps)]

        # add the auxdata and format for .edr to the reference description
        self.description['auxdata'] = Path(self.testdata).resolve()
        self.description['format'] = self.reader.format

        # for testing the selection of data/time
        self.time_selector = "Time"
        self.select_time_ref = [step._data[self.time_selector]
                                for step in self.auxsteps]
        self.data_selector = "Bond"  # selects all data
        self.select_data_ref = [step._data["Bond"] for step in self.auxsteps]

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

        # for testing EDRReader.get_data()
        self.times = np.array([0., 0.02, 0.04, 0.06])
        self.bonds = np.array([1374.82324219, 1426.22521973,
                               1482.0098877, 1470.33752441])
        self.angles = np.array([3764.52734375, 3752.83032227,
                                3731.59179688, 3683.40942383])


@pytest.mark.skipif(not HAS_PYEDR, reason="pyedr not installed")
class TestEDRReader(BaseAuxReaderTest):
    """ Class to conduct tests for the auxiliary EDRReader

    Normally, it would be desirable to use the tests from
    :class:`BaseAuxReaderTest`, but this is not possible for some of the tests
    there because of the differences between :class:`BaseAuxReference` and
    :class:`EDRReference` which are ultimately due to the different file
    format. Therefore, some tests had to be overridden here. New tests for EDR-
    specific functionality were also added.
    """

    @staticmethod
    @pytest.fixture
    def ref():
        return EDRReference()

    @staticmethod
    @pytest.fixture
    def ref_universe(ref):
        u = mda.Universe(AUX_EDR_TPR, AUX_EDR_XTC)
        u.trajectory.add_auxiliary(ref.testdata, {"test": "Bond"})
        return u

    @staticmethod
    @pytest.fixture
    def reader(ref):
        reader = ref.reader(
            ref.testdata,
            initial_time=ref.initial_time,
            dt=ref.dt, auxname=ref.name,
            time_selector="Time",
            data_selector=None
        )
        ref_units = get_edr_unit_dict(0)
        if reader.unit_dict != ref_units:
            for term, ref_unit in ref_units.items():
                data = reader.data_dict[term]
                reader_unit = reader.unit_dict[term]
                try:
                    reader.data_dict[term] = units.convert(data,
                                                           reader_unit,
                                                           ref_unit)
                except ValueError:
                    continue  # some units not supported yet
        reader.rewind()
        return reader

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
            assert_allclose(ts.aux.test, ref.auxsteps[i].data["Bond"])

    def test_iterate_as_auxiliary_from_trajectory(self, ref, ref_universe):
        # check representative values of aux for each frame are as expected
        # trajectory here has same dt, offset, so there's a direct
        # correspondence between frames and steps, and iter_as_aux will run
        # through all frames
        for i, ts in enumerate(ref_universe.trajectory.iter_as_aux('test')):
            assert_allclose(ts.aux.test, ref.auxsteps[i].data["Bond"])

    def test_step_to_frame_time_diff(self, reader, ref):
        # Timestep is 0.002 (10 %) longer than auxiliary data
        ts = mda.coordinates.base.Timestep(0, dt=ref.dt + 0.002)

        # Test all 4 frames
        for idx in range(4):

            frame, time_diff = reader.step_to_frame(idx, ts,
                                                    return_time_diff=True)

            assert frame == idx
            assert_allclose(time_diff, idx * 0.002)

    def test_read_lower_freq_timestep(self, ref, reader):
        # test reading a timestep with lower frequency
        ts = ref.lower_freq_ts
        reader.update_ts(ts)
        # check the value set in ts is as we expect
        assert_allclose(ts.aux.test["Bond"], ref.lowf_closest_rep,
                        err_msg="Representative value in ts.aux "
                                "does not match")

    def test_read_higher_freq_timestep(self, ref, reader):
        # try reading a timestep with higher frequency
        ts = ref.higher_freq_ts
        reader.update_ts(ts)
        assert_allclose(ts.aux.test, ref.highf_rep,
                        err_msg="Representative value in ts.aux "
                                "does not match")

    def test_read_offset_timestep(self, ref, reader):
        # try reading a timestep offset from auxiliary
        ts = ref.offset_ts
        reader.update_ts(ts)
        assert_allclose(ts.aux.test["Bond"], ref.offset_closest_rep,
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
        assert_allclose(test_value, ref.lowf_average_rep,
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
        assert_allclose(ts.aux.test["Bond"], ref.lowf_cutoff_average_rep,
                        err_msg="Representative value does not match when "
                                "applying cutoff")

    def test_add_all_terms_from_file(self, ref, ref_universe):
        ref_universe.trajectory.add_auxiliary(auxdata=ref.testdata)
        # adding "test" manually to match above addition of test term
        ref_terms = ["test"] + [key for key in get_auxstep_data(0).keys()]
        terms = [key for key in ref_universe.trajectory._auxs]
        assert ref_terms == terms

    def test_add_all_terms_from_reader(self, ref_universe, reader):
        ref_universe.trajectory.add_auxiliary(auxdata=reader)
        ref_terms = ["test"] + [key for key in get_auxstep_data(0).keys()]
        terms = [key for key in ref_universe.trajectory._auxs]
        assert ref_terms == terms

    def test_add_term_list_custom_names_from_file(self, ref, ref_universe):
        ref_universe.trajectory.add_auxiliary(ref.testdata,
                                              {"bond": "Bond",
                                               "temp": "Temperature"},
                                              )
        ref_dict = get_auxstep_data(0)
        assert ref_universe.trajectory.ts.aux.bond == ref_dict["Bond"]
        assert ref_universe.trajectory.ts.aux.temp == ref_dict["Temperature"]

    def test_add_term_list_custom_names_from_reader(self, ref_universe,
                                                    reader):
        ref_universe.trajectory.add_auxiliary(reader,
                                              {"bond": "Bond",
                                               "temp": "Temperature"},
                                              )
        ref_dict = get_auxstep_data(0)
        assert ref_universe.trajectory.ts.aux.bond == ref_dict["Bond"]
        assert ref_universe.trajectory.ts.aux.temp == ref_dict["Temperature"]

    def test_raise_error_if_auxname_already_assigned(self, ref_universe,
                                                     reader):
        with pytest.raises(ValueError, match="Auxiliary data with name"):
            ref_universe.trajectory.add_auxiliary(reader, "test", "Bond")

    def test_add_single_term_custom_name_from_file(self, ref, ref_universe):
        ref_universe.trajectory.add_auxiliary(ref.testdata,
                                              {"temp": "Temperature"},
                                              )
        ref_dict = get_auxstep_data(0)
        assert ref_universe.trajectory.ts.aux.temp == ref_dict["Temperature"]

    def test_add_single_term_custom_name_from_reader(self, ref_universe,
                                                     reader):
        ref_universe.trajectory.add_auxiliary(reader, {"temp": "Temperature"})
        ref_dict = get_auxstep_data(0)
        assert ref_universe.trajectory.ts.aux.temp == ref_dict["Temperature"]

    def test_terms_update_on_iter(self, ref_universe, reader):
        ref_universe.trajectory.add_auxiliary(reader,
                                              {"bond": "Bond",
                                               "temp": "Temperature"},
                                              )
        ref_dict = get_auxstep_data(0)
        assert ref_universe.trajectory.ts.aux.bond == ref_dict["Bond"]
        assert ref_universe.trajectory.ts.aux.temp == ref_dict["Temperature"]
        ref_dict = get_auxstep_data(1)
        ref_universe.trajectory.next()
        assert ref_universe.trajectory.ts.aux.bond == ref_dict["Bond"]
        assert ref_universe.trajectory.ts.aux.temp == ref_dict["Temperature"]

    def test_invalid_data_selector(self, ref, ref_universe):
        with pytest.raises(KeyError, match="'Nonsense' is not a key"):
            ref_universe.trajectory.add_auxiliary(AUX_EDR,
                                                  {"something": "Nonsense"},
                                                  )

    def test_read_all_times(self, reader):
        all_times_expected = np.array([0., 0.02, 0.04, 0.06])
        assert_allclose(all_times_expected, reader.read_all_times())

    def test_get_data_from_string(self, ref, reader):
        returned = reader.get_data("Bond")
        assert isinstance(returned, dict)
        assert_allclose(ref.times, returned["Time"])
        assert_allclose(ref.bonds, returned["Bond"])

    def test_get_data_from_list(self, ref, reader):
        returned = reader.get_data(["Bond", "Angle"])
        assert isinstance(returned, dict)
        assert_allclose(ref.times, returned["Time"])
        assert_allclose(ref.bonds, returned["Bond"])
        assert_allclose(ref.angles, returned["Angle"])

    def test_get_data_everything(self, ref, reader):
        returned = reader.get_data()
        returned_asterisk = reader.get_data()
        assert returned.keys() == returned_asterisk.keys()
        ref_terms = [key for key in get_auxstep_data(0).keys()]
        assert ref_terms == reader.terms
        assert_allclose(ref.bonds, returned["Bond"])

    @pytest.mark.parametrize("get_data_input", (42,
                                                "Not a valid term",
                                                ["Bond", "Not a valid term"]))
    def test_get_data_invalid_selections(self, reader, get_data_input):
        with pytest.raises(KeyError, match="data selector"):
            reader.get_data(get_data_input)

    def test_warning_when_space_in_aux_spec(self, ref_universe, reader):
        with pytest.warns(UserWarning, match="Auxiliary name"):
            ref_universe.trajectory.add_auxiliary(reader,
                                                  {"Pres. DC": "Pres. DC"},
                                                  )

    def test_warn_too_much_memory_usage(self, ref_universe, reader):
        with pytest.warns(UserWarning, match="AuxReader: memory usage "
                          "warning! Auxiliary data takes up 3[0-9.]*e-06 GB of"
                          r" memory \(Warning limit: 1e-08 GB\)"):
            ref_universe.trajectory.add_auxiliary(reader,
                                                  {"temp": "Temperature"},
                                                  memory_limit=10,
                                                  )

    def test_auxreader_picklable(self, reader):
        new_reader = pickle.loads(pickle.dumps(reader))
        for step_index, auxstep in enumerate(reader):
            assert_auxstep_equal(new_reader[step_index], auxstep)

    def test_units_are_converted_by_EDRReader(self, reader):
        original_units = get_edr_unit_dict(0)
        reader_units = reader.unit_dict
        # so far, lengths and speeds are converted
        for term in ["Box-X", "Box-Vel-XX"]:
            assert original_units[term] != reader_units[term]

    def test_warning_when_unknown_unit(self, ref_universe, reader):
        with pytest.warns(UserWarning, match="Could not find"):
            ref_universe.trajectory.add_auxiliary(reader,
                                                  {"temp": "Temperature"},
                                                  )

    def test_unit_conversion_is_optional(self, ref):
        reader = ref.reader(
            ref.testdata,
            initial_time=ref.initial_time,
            dt=ref.dt, auxname=ref.name,
            time_selector="Time",
            data_selector=None,
            convert_units=False
        )
        ref_units = get_edr_unit_dict(0)
        # The units from AUX_EDR match the ones from the reference
        # data file AUX_EDR_RAW. If the EDRReader does not convert the units on
        # reading the file, then the two unit dictionaries should be identical.
        assert reader.unit_dict == ref_units


@pytest.mark.skipif(not HAS_PYEDR, reason="pyedr not installed")
def test_single_frame_input_file():
    """Previously, EDRReader could not handle EDR input files with only one
       frame. See Issue #3999."""
    reader = mda.auxiliary.EDR.EDRReader(AUX_EDR_SINGLE_FRAME,
                                         convert_units=False)
    ref_dict = get_auxstep_data(0)
    reader_data_dict = reader.auxstep.data
    assert ref_dict == reader_data_dict
