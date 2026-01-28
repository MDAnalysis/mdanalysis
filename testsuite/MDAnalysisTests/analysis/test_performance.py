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
"""Tests for optional performance timing in AnalysisBase"""

import MDAnalysis as mda
import pytest
from MDAnalysis.analysis import base

from MDAnalysisTests.datafiles import DCD, PSF


class DummyAnalysis(base.AnalysisBase):
    """A simple analysis class for testing timing functionality"""

    def __init__(self, trajectory, **kwargs):
        super(DummyAnalysis, self).__init__(trajectory, **kwargs)

    def _prepare(self):
        self.results.data = []

    def _single_frame(self):
        # Simple operation to record
        self.results.data.append(self._ts.frame)

    def _conclude(self):
        # Convert to list to ensure results are set
        pass


class TestPerformanceTiming:
    """Test suite for performance timing feature"""

    @pytest.fixture
    def universe(self):
        """Fixture to load test trajectory"""
        return mda.Universe(PSF, DCD)

    def test_timing_disabled_by_default(self, universe):
        """Test that timing is disabled by default"""
        analysis = DummyAnalysis(universe.trajectory)
        analysis.run(stop=2)
        
        # When timing is disabled, results.timing should not exist
        assert not hasattr(analysis.results, 'timing')

    def test_timing_disabled_explicit(self, universe):
        """Test that timing can be explicitly disabled"""
        analysis = DummyAnalysis(universe.trajectory, enable_timing=False)
        analysis.run(stop=2)
        
        # When timing is disabled, results.timing should not exist
        assert not hasattr(analysis.results, 'timing')

    def test_timing_enabled(self, universe):
        """Test that timing is recorded when enabled"""
        analysis = DummyAnalysis(universe.trajectory, enable_timing=True)
        analysis.run(stop=2)
        
        # When timing is enabled, results.timing should exist
        assert hasattr(analysis.results, 'timing')
        assert isinstance(analysis.results.timing, dict)

    def test_timing_has_prepare_key(self, universe):
        """Test that 'prepare' timing is recorded"""
        analysis = DummyAnalysis(universe.trajectory, enable_timing=True)
        analysis.run(stop=2)
        
        assert 'prepare' in analysis.results.timing
        assert isinstance(analysis.results.timing['prepare'], float)

    def test_timing_has_frame_iteration_key(self, universe):
        """Test that 'frame_iteration' timing is recorded"""
        analysis = DummyAnalysis(universe.trajectory, enable_timing=True)
        analysis.run(stop=2)
        
        assert 'frame_iteration' in analysis.results.timing
        assert isinstance(analysis.results.timing['frame_iteration'], float)

    def test_timing_has_total_key(self, universe):
        """Test that 'total' timing is recorded"""
        analysis = DummyAnalysis(universe.trajectory, enable_timing=True)
        analysis.run(stop=2)
        
        assert 'total' in analysis.results.timing
        assert isinstance(analysis.results.timing['total'], float)

    def test_timing_values_are_floats(self, universe):
        """Test that all timing values are floats"""
        analysis = DummyAnalysis(universe.trajectory, enable_timing=True)
        analysis.run(stop=5)
        
        for key, value in analysis.results.timing.items():
            assert isinstance(value, float), f"Timing value for '{key}' is not a float"

    def test_timing_values_are_positive(self, universe):
        """Test that all timing values are positive"""
        analysis = DummyAnalysis(universe.trajectory, enable_timing=True)
        analysis.run(stop=5)
        
        for key, value in analysis.results.timing.items():
            assert value >= 0, f"Timing value for '{key}' is negative"

    def test_timing_with_frame_selection(self, universe):
        """Test that timing works with frame selection"""
        analysis = DummyAnalysis(universe.trajectory, enable_timing=True)
        analysis.run(start=1, stop=5, step=2)
        
        assert hasattr(analysis.results, 'timing')
        assert 'prepare' in analysis.results.timing
        assert 'frame_iteration' in analysis.results.timing
        assert 'total' in analysis.results.timing

    def test_timing_with_frames_list(self, universe):
        """Test that timing works with frames list"""
        analysis = DummyAnalysis(universe.trajectory, enable_timing=True)
        analysis.run(frames=[0, 2, 4])
        
        assert hasattr(analysis.results, 'timing')
        assert 'prepare' in analysis.results.timing
        assert 'frame_iteration' in analysis.results.timing
        assert 'total' in analysis.results.timing

    def test_timing_prepare_less_than_total(self, universe):
        """Test that prepare time is less than total time"""
        analysis = DummyAnalysis(universe.trajectory, enable_timing=True)
        analysis.run(stop=5)
        
        assert analysis.results.timing['prepare'] < analysis.results.timing['total']

    def test_timing_frame_iteration_less_than_total(self, universe):
        """Test that frame_iteration time is less than total time"""
        analysis = DummyAnalysis(universe.trajectory, enable_timing=True)
        analysis.run(stop=5)
        
        assert analysis.results.timing['frame_iteration'] < analysis.results.timing['total']

    def test_timing_consistency_across_runs(self, universe):
        """Test that timing measurements are consistent in repeated runs"""
        analysis1 = DummyAnalysis(universe.trajectory, enable_timing=True)
        analysis1.run(stop=3)
        
        # Reload universe for clean state
        universe2 = mda.Universe(PSF, DCD)
        analysis2 = DummyAnalysis(universe2.trajectory, enable_timing=True)
        analysis2.run(stop=3)
        
        # Both should have timing info
        assert hasattr(analysis1.results, 'timing')
        assert hasattr(analysis2.results, 'timing')
        
        # Same keys should be present
        assert set(analysis1.results.timing.keys()) == set(analysis2.results.timing.keys())

    def test_timing_empty_frames(self, universe):
        """Test timing behavior with empty frame selection"""
        analysis = DummyAnalysis(universe.trajectory, enable_timing=True)
        analysis.run(stop=0)
        
        # Should still record timing even with no frames
        assert hasattr(analysis.results, 'timing')
        assert 'prepare' in analysis.results.timing

    def test_timing_does_not_affect_analysis_results(self, universe):
        """Verify that enabling timing doesn't change analysis results"""
        analysis_no_timing = DummyAnalysis(universe.trajectory, enable_timing=False)
        analysis_no_timing.run(stop=5)
        
        # Reload universe for clean state
        universe2 = mda.Universe(PSF, DCD)
        analysis_with_timing = DummyAnalysis(universe2.trajectory, enable_timing=True)
        analysis_with_timing.run(stop=5)
        
        # Results should be identical regardless of timing enabled
        assert analysis_no_timing.results.data == analysis_with_timing.results.data

    def test_timing_logical_relationships_with_tolerance(self, universe):
        """Test that timing values have correct logical relationships"""
        analysis = DummyAnalysis(universe.trajectory, enable_timing=True)
        analysis.run(stop=5)
        
        timing = analysis.results.timing
        # Allow 1ms floating point tolerance due to measurement precision
        tolerance = 0.001
        
        # prepare should not exceed total
        assert timing['prepare'] <= timing['total'] + tolerance, \
            "Prepare time should not exceed total time"
        # frame_iteration should not exceed total
        assert timing['frame_iteration'] <= timing['total'] + tolerance, \
            "Frame iteration time should not exceed total time"

    def test_timing_values_no_negative(self, universe):
        """Ensure timing values are never negative (catches clock issues)"""
        analysis = DummyAnalysis(universe.trajectory, enable_timing=True)
        analysis.run(stop=5)
        
        for key, value in analysis.results.timing.items():
            assert value >= 0, f"Timing value for '{key}' is negative: {value}"
