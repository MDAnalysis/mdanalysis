import numpy
from MDAnalysis import *
from pylab import *

universe = Universe(...)

group = universe.selectAtoms("resname DMPC and ( name N or name P )")

skip = 5  # probably don't have enough memory to load the entire trajectory
data = universe.dcd.timeseries(group, skip=skip, format="afc")

# The data is structured as (N*2,numframes,3), where N is the number of lipids
# since you are pulling out the coordinates of the phosphate and nitrogen for each lipid
# v will contain the vector connecting each phosphate and nitrogen of the headgroup
v = numpy.subtract(data[0::2], data[1::2])
len_v = numpy.sqrt(numpy.add.reduce(numpy.power(v, 2), axis=-1))
# Normalize v
v /= len_v[..., numpy.newaxis]
# this is only necessary if you have a lot of lipids and/or a lot of frames - it never hurts to free up memory
del len_v

lagtime = numpy.arange(1, len(data[0]) / 2 + 1)
C_t = numpy.array(
    [numpy.average(1.5 * numpy.power(numpy.add.reduce(v[:, :-lag] * v[:, lag:], axis=-1), 2) - 0.5) for lag in lagtime])

plot(C_t)

