#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
"""
Test optimized histogram implementation --- :mod:`MDAnalysisTests.lib.test_histogram_opt`
=========================================================================================

Tests the optimized histogram implementation against numpy.histogram for correctness
and performance improvements.

"""
import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

try:
    from MDAnalysis.lib.histogram_opt import optimized_histogram, HAS_NUMBA
except ImportError:
    HAS_NUMBA = False
    optimized_histogram = None


@pytest.mark.skipif(not HAS_NUMBA, reason="Numba not available")
class TestOptimizedHistogram:
    """Test the optimized histogram implementation."""
    
    def test_correctness_uniform(self):
        """Test correctness with uniform distribution."""
        np.random.seed(42)
        data = np.random.uniform(0, 15, 10000).astype(np.float64)
        
        np_hist, np_edges = np.histogram(data, bins=75, range=(0.0, 15.0))
        opt_hist, opt_edges = optimized_histogram(data, bins=75, range=(0.0, 15.0))
        
        assert_allclose(np_hist, opt_hist, rtol=1e-14, atol=1)
        assert_allclose(np_edges, opt_edges, rtol=1e-14)
    
    def test_correctness_normal(self):
        """Test correctness with normal distribution."""
        np.random.seed(42)
        data = np.random.normal(7.5, 2, 10000).clip(0, 15).astype(np.float64)
        
        np_hist, np_edges = np.histogram(data, bins=75, range=(0.0, 15.0))
        opt_hist, opt_edges = optimized_histogram(data, bins=75, range=(0.0, 15.0))
        
        assert_allclose(np_hist, opt_hist, rtol=1e-14, atol=1)
        assert_allclose(np_edges, opt_edges, rtol=1e-14)
    
    def test_edge_cases_zeros(self):
        """Test with all zeros."""
        data = np.zeros(1000, dtype=np.float64)
        
        np_hist, np_edges = np.histogram(data, bins=75, range=(0.0, 15.0))
        opt_hist, opt_edges = optimized_histogram(data, bins=75, range=(0.0, 15.0))
        
        assert_allclose(np_hist, opt_hist, rtol=1e-14, atol=1)
        assert_allclose(np_edges, opt_edges, rtol=1e-14)
    
    def test_edge_cases_max_values(self):
        """Test with values at maximum range."""
        data = np.ones(1000, dtype=np.float64) * 14.999
        
        np_hist, np_edges = np.histogram(data, bins=75, range=(0.0, 15.0))
        opt_hist, opt_edges = optimized_histogram(data, bins=75, range=(0.0, 15.0))
        
        assert_allclose(np_hist, opt_hist, rtol=1e-14, atol=1)
        assert_allclose(np_edges, opt_edges, rtol=1e-14)
    
    def test_boundary_values(self):
        """Test with boundary values."""
        data = np.array([0.0, 14.999, 15.0, 7.5] * 250, dtype=np.float64)
        
        np_hist, np_edges = np.histogram(data, bins=75, range=(0.0, 15.0))
        opt_hist, opt_edges = optimized_histogram(data, bins=75, range=(0.0, 15.0))
        
        # Allow for small differences at boundaries due to floating point precision
        assert_allclose(np_hist, opt_hist, rtol=1e-14, atol=1)
        assert_allclose(np_edges, opt_edges, rtol=1e-14)
    
    def test_bin_edges_values(self):
        """Test with values exactly at bin edges."""
        data = np.linspace(0, 15, 1001, dtype=np.float64)
        
        np_hist, np_edges = np.histogram(data, bins=75, range=(0.0, 15.0))
        opt_hist, opt_edges = optimized_histogram(data, bins=75, range=(0.0, 15.0))
        
        # Allow for small differences at boundaries due to floating point precision
        assert_allclose(np_hist, opt_hist, rtol=1e-14, atol=1)
        assert_allclose(np_edges, opt_edges, rtol=1e-14)
    
    def test_serial_parallel_consistency(self):
        """Test that serial and parallel versions produce identical results."""
        np.random.seed(42)
        data = np.random.random(100000).astype(np.float64) * 15.0
        
        hist_serial, edges_serial = optimized_histogram(
            data, bins=75, range=(0.0, 15.0), use_parallel=False
        )
        hist_parallel, edges_parallel = optimized_histogram(
            data, bins=75, range=(0.0, 15.0), use_parallel=True
        )
        
        assert_allclose(hist_serial, hist_parallel, rtol=1e-14, atol=1)
        assert_array_equal(edges_serial, edges_parallel)
    
    def test_different_bin_counts(self):
        """Test with different bin counts."""
        np.random.seed(42)
        data = np.random.random(10000).astype(np.float64) * 15.0
        
        for bins in [10, 50, 100, 200]:
            np_hist, np_edges = np.histogram(data, bins=bins, range=(0.0, 15.0))
            opt_hist, opt_edges = optimized_histogram(data, bins=bins, range=(0.0, 15.0))
            
            assert_allclose(np_hist, opt_hist, rtol=1e-14, atol=1,
                          err_msg=f"Failed for bins={bins}")
            assert_allclose(np_edges, opt_edges, rtol=1e-14,
                          err_msg=f"Failed for bins={bins}")
    
    def test_different_ranges(self):
        """Test with different histogram ranges."""
        np.random.seed(42)
        data = np.random.random(10000).astype(np.float64) * 20.0
        
        for range_val in [(0.0, 10.0), (0.0, 20.0), (5.0, 15.0)]:
            np_hist, np_edges = np.histogram(data, bins=75, range=range_val)
            opt_hist, opt_edges = optimized_histogram(data, bins=75, range=range_val)
            
            assert_allclose(np_hist, opt_hist, rtol=1e-14, atol=1,
                          err_msg=f"Failed for range={range_val}")
            assert_allclose(np_edges, opt_edges, rtol=1e-14,
                          err_msg=f"Failed for range={range_val}")
    
    def test_non_contiguous_array(self):
        """Test with non-contiguous array."""
        np.random.seed(42)
        # Create a non-contiguous array by slicing
        data = np.random.random(20000).astype(np.float64) * 15.0
        data_non_contig = data[::2]  # Every other element
        
        assert not data_non_contig.flags['C_CONTIGUOUS']
        
        np_hist, np_edges = np.histogram(data_non_contig, bins=75, range=(0.0, 15.0))
        opt_hist, opt_edges = optimized_histogram(data_non_contig, bins=75, range=(0.0, 15.0))
        
        assert_allclose(np_hist, opt_hist, rtol=1e-14, atol=1)
        assert_allclose(np_edges, opt_edges, rtol=1e-14)
    
    @pytest.mark.parametrize("size", [100, 1000, 10000, 100000])
    def test_scaling(self, size):
        """Test correctness at different scales."""
        np.random.seed(42)
        data = np.random.random(size).astype(np.float64) * 15.0
        
        np_hist, np_edges = np.histogram(data, bins=75, range=(0.0, 15.0))
        opt_hist, opt_edges = optimized_histogram(data, bins=75, range=(0.0, 15.0))
        
        assert_allclose(np_hist, opt_hist, rtol=1e-14, atol=1,
                      err_msg=f"Failed for size={size}")
        assert_allclose(np_edges, opt_edges, rtol=1e-14,
                      err_msg=f"Failed for size={size}")


@pytest.mark.skipif(HAS_NUMBA, reason="Testing fallback when Numba not available")
class TestHistogramFallback:
    """Test that the module falls back gracefully when Numba is not available."""
    
    def test_import_without_numba(self):
        """Test that import works without Numba."""
        # This test will only run if Numba is not installed
        # The import should still work but HAS_NUMBA should be False
        from MDAnalysis.lib import histogram_opt
        assert not histogram_opt.HAS_NUMBA