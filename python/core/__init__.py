# $Id$
"""Core functions of MDAnalysis: The basic class is an AtomGroup; the
whole simulation is called the Universe. Selections are computed on an
AtomGroup and return another AtomGroup. Timeseries are a convenient way to analyse trajectories.

To get started, load the Universe:

  u = Universe(psffilename,dcdfilename)  
"""
__all__ = ['AtomGroup', 'Selection', 'Timeseries', 
           'distances', 'rms_fitting']

# NOTE: KDTree routines are significantly faster for some distance
#       selections. However, they cannot deal with periodic boxes and thus
#       ignore periodicity; if periodicity is crucial, disable KDTree
#       routines with
#use_KDTree_routines = False
use_KDTree_routines = 'fast' 
# True, 'fast'   - only use KDTree routines that are typically faster than others
# 'always'       - always use KDTree routines where available (eg for benchmarking)
# False, 'never' - always use alternatives

 
