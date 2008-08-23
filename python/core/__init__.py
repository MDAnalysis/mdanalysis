# $Id$
"""Core functions of MDAnalysis: The basic class is an AtomGroup; the
whole simulation is called the Universe. Selections are computed on an
AtomGroup and return another AtomGroup. Timeseries are a convenient way to analyse trajectories.

To get started, load the Universe:

  u = Universe(psffilename,dcdfilename)  
"""
__all__ = ['AtomGroup', 'Selection', 'Timeseries', 
           'distances', 'rms_fitting']
