"""
"""
from MDAnalysis import *
from MDAnalysis import _numtools
from pylab import *
import Numeric

system = AtomGroup.Universe("data/kalp16.psf", "data/kalp16.dcd")
skip = 1000
asel = system.selectAtoms("segid KALP")
data = system._dcd.timeseries(asel, skip)
