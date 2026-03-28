# ASV Benchmarking Guide for MDAnalysis

**For**: Developers and contributors
**Purpose**: Using ASV to benchmark and compare performance improvements
**Last Updated**: March 28, 2026

---

## Table of Contents
1. [Quick Start](#quick-start)
2. [Writing Benchmarks](#writing-benchmarks)
3. [Running Benchmarks](#running-benchmarks)
4. [Comparing Results](#comparing-results)
5. [Best Practices](#best-practices)
6. [Troubleshooting](#troubleshooting)

---

## Quick Start

### Installation

```bash
# ASV should already be installed
pip install asv line_profiler snakeviz
```

### First Benchmark Run

```bash
cd benchmarks_repo
python -m asv run --conf asv_local.conf.json --quick
```

This runs all benchmarks with reduced iterations. Full runs take much longer.

---

## Writing Benchmarks

### Basic Benchmark Class

```python
class TimeSortingAlgorithm:
    """Benchmark for sorting performance"""
    
    # Parameters to vary across runs
    param_names = ['array_size', 'data_type']
    params = [
        [100, 1000, 10000],                    # array_size values
        ['random', 'sorted', 'reverse_sorted'] # data_type values
    ]
    
    def setup(self, array_size, data_type):
        """Called before each benchmark iteration"""
        import numpy as np
        
        if data_type == 'random':
            self.array = np.random.random(array_size)
        elif data_type == 'sorted':
            self.array = np.arange(array_size, dtype=float)
        else:  # reverse sorted
            self.array = np.arange(array_size, dtype=float)[::-1]
    
    def time_sort(self, array_size, data_type):
        """Timed function - measures wall-clock time"""
        import numpy as np
        return np.sort(self.array.copy())
    
    def peakmem_sort(self, array_size, data_type):
        """Memory tracking function"""
        # If function is memory intensive, track peak usage
        import numpy as np
        return np.sort(self.array.copy())
    
    def teardown(self, array_size, data_type):
        """Called after each benchmark iteration"""
        self.array = None
```

### Available Timing Functions

| Function | Measures | Use Case |
|----------|----------|----------|
| `time_*` | Wall-clock time | Most common - total execution time |
| `timeraw_*` | Raw timing (no warmup) | Low-level operations |
| `mem_*` | Memory allocation | Track memory overhead |
| `peakmem_*` | Peak memory usage | Identify memory spikes |
| `track_*` | Arbitrary metrics | Custom tracking |

---

## Running Benchmarks

### Configuration Files

**For local development** (`asv_local.conf.json`):
```bash
asv run --conf asv_local.conf.json
```

**For CI/CD** (`asv.conf.json`):
```bash
asv run --conf config/asv.conf.json
```

### Common Commands

```bash
# Quick benchmark (reduced iterations)
asv run --quick

# Full benchmark suite
asv run

# Benchmark specific function
asv run --bench ^TrajectoryLoading

# Generate HTML results
asv run && asv publish

# Preview in browser
asv preview
```

### Understanding Output

```
Running benchmarks for commit abc1234...

TrajectoryLoading.time_load_dcd          178.42ms
TrajectoryLoading.time_load_xtc          192.15ms  
AtomSelection.time_simple_selection      2.34ms
AtomSelection.time_complex_selection     8.91ms
DistanceCalculation.time_pairwise_500    11.30ms
DistanceCalculation.time_pairwise_1000   35.61ms

Traceback (most recent call last):
  # Some benchmarks might fail if test data unavailable
  # This is normal - they'll show <notrun>
```

---

## Comparing Results

### Before/After Comparison

```bash
# Stash current results
asv run -b <benchmark_name> -o baseline.json

# Apply optimization
git apply optimization.patch

# Run after optimization
asv run -b <benchmark_name> -o optimized.json

# Compare
asv compare baseline.json optimized.json
```

### Interpreting Comparison Output

```
benchmarks.DistanceCalculation.time_pairwise

              before      after    ratio
---  ---  ---  -----  ----------  -------
     500  1.0  363ms   11.3ms    0.03x  ⭐ (33x faster!)
    1000  1.0 1456ms   35.6ms    0.02x  ⭐ (41x faster!)
```

**Ratio < 1.0** = After is faster (improvement!)
**Ratio > 1.0** = After is slower (regression!)

---

## Performance Profiling

### Generate Profiles

```bash
# Create cProfile data
python -m cProfile -o profile.prof your_script.py

# Visualize with snakeviz
snakeviz profile.prof
```

### Identify Hotspots

```bash
# Line-level profiling
kernprof -l script.py
python -m line_profiler script.py.lprof
```

---

## Best Practices

### ✅ Do This

1. **Setup and cleanup**: Use `setup()` and `teardown()` methods
2. **Parameter ranges**: Test across realistic data sizes
3. **Multiple runs**: ASV averages multiple iterations automatically
4. **Clear names**: Use descriptive function names (e.g., `time_rmsd_10k_atoms`)
5. **Document parameters**: Explain what each parameter tests

### ❌ Don't Do This

1. **Global state**: Avoid modifying module-level variables in benchmarks
2. **External files**: Don't rely on files that might not exist
3. **randomness**: Set `np.random.seed()` for reproducibility
4. **Background tasks**: Don't benchmark code using other processes
5. **Unbounded loops**: Always specify exact iteration count

### Example: Good vs Bad Benchmark

```python
# ❌ BAD: Unpredictable, unreproducible
class UntestableBenchmark:
    def time_analysis(self):
        import random
        # Nondeterministic!
        u = mda.Universe("random_universe.pdb")
        for i in range(random.randint(10, 1000)):
            u.trajectory.next()

# ✅ GOOD: Predictable, reproducible, parameterized
class TestableBenchmark:
    param_names = ['n_frames', 'n_atoms']
    params = [[100, 1000], [1000, 10000]]
    
    def setup(self, n_frames, n_atoms):
        import numpy as np
        np.random.seed(42)  # Reproducible
        self.coords = np.random.random((n_atoms, 3))
        self.n_frames = n_frames
    
    def time_distance_calc(self, n_frames, n_atoms):
        # Predictable, fully controlled
        diff = self.coords[:, np.newaxis, :] - self.coords[np.newaxis, :, :]
        return np.sqrt(np.sum(diff**2, axis=2))
```

---

## Common Patterns

### Pattern 1: Scaling Analysis

```python
class ScalingBenchmark:
    """Test how performance scales with system size"""
    
    param_names = ['n_atoms']
    params = [[100, 1000, 10000, 100000]]
    
    def setup(self, n_atoms):
        import numpy as np
        self.coords = np.random.random((n_atoms, 3))
    
    def time_distance_matrix(self, n_atoms):
        import numpy as np
        diff = self.coords[:, np.newaxis, :] - self.coords[np.newaxis, :, :]
        return np.sqrt(np.sum(diff**2, axis=2))
```

**Use case**: Identify O(n), O(n²), O(n³) behavior

### Pattern 2: Multiple Implementations

```python
class AlgorithmComparison:
    """Compare different implementations"""
    
    param_names = ['implementation']
    params = [['python', 'numpy', 'cython']]
    
    def setup(self, implementation):
        self.impl = implementation
    
    def time_compute(self, implementation):
        if implementation == 'python':
            return compute_python()
        elif implementation == 'numpy':
            return compute_numpy()
        else:
            return compute_cython()
```

### Pattern 3: Caching Impact

```python
class CachingBenchmark:
    """Measure caching effects"""
    
    param_names = ['cache_enabled']
    params = [[True, False]]
    
    def setup(self, cache_enabled):
        self.cache_enabled = cache_enabled
        if cache_enabled:
            enable_cache()
        else:
            disable_cache()
    
    def time_repeated_selection(self, cache_enabled):
        u.select_atoms("protein")  # First call
        u.select_atoms("protein")  # Second call (cached?)
        u.select_atoms("protein")  # Third call (cached?)
```

---

## System Information

ASV automatically captures:
- Python version
- Installed packages and versions
- CPU model
- RAM amount
- Operating system

This allows comparing results across machines and time.

---

## Troubleshooting

### Issue: "Benchmark did not complete"

**Cause**: Infinite loops, missing dependencies, exceptions
**Solution**: 
```python
def setup(self, param):
    try:
        # Your setup code
    except Exception as e:
        # Handle gracefully
        pass

def time_operation(self, param):
    if not hasattr(self, 'data'):
        return  # Skip if setup failed
```

### Issue: "Results inconsistent across runs"

**Cause**: Non-deterministic code, system load, interference
**Solution**:
```python
def setup(self, param):
    import numpy as np
    np.random.seed(42)  # ← Fix random seed!
```

### Issue: ASV takes forever to run

**Cause**: Too many iterations, too-large parameters
**Solution**:
```bash
asv run --quick  # ← Use quick mode for development
asv run --sample=5  # ← Reduce sample count
```

### Issue: Memory errors on large benchmarks

**Solution**:
```python
def teardown(self, param):
    # Clean up large objects
    self.data = None
    gc.collect()
```

---

## Continuous Integration

### GitHub Actions Example

```yaml
- name: Run ASV benchmarks
  run: |
    cd benchmarks_repo
    python -m asv run --conf asv_local.conf.json --quick
    python -m asv preview
```

---

## Additional Resources

- [ASV Documentation](https://asv.readthedocs.io/)
- [MDAnalysis Performance Guide](https://docs.mdanalysis.org/)
- [NumPy Performance Tips](https://numpy.org/doc/stable/reference/internals.html)

---

## Example: Full Benchmark Module

See `benchmarks/bench_core_operations.py` for real examples

```python
"""
Complete example with multiple benchmark classes
"""

class CompleteExample:
    """Template for complete benchmarks"""
    
    param_names = ['param1', 'param2']
    params = [
        [1, 10, 100],
        ['option_a', 'option_b']
    ]
    timeout = 60  # Max 60 seconds per test
    
    def setup_cache(self):
        """Cache that persists across parameter combinations"""
        # Download data once
        return {'data': load_large_file()}
    
    def setup(self, param1, param2):
        """Setup for this parameter combination"""
        pass
    
    def time_operation(self, param1, param2):
        """The timed operation"""
        pass
    
    def teardown(self, param1, param2):
        """Cleanup"""
        pass
```

---

**Next**: See `performance.md` for optimization results and `PHASE_5_OPTIMIZATION_RESULTS.md` for real benchmarking examples.
