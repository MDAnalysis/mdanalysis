import numpy
from MDAnalysis import *

universe = Universe("system.gro")

tail_carbons = numpy.arange(2, 15)

skip = 5
order_param = numpy.zeros(len(tail_carbons))

for i, carbon in enumerate(tail_carbons):
    selection = "resname DMPC and ( name C2%d or name H%dR or name H%dS or name C3%d or name H%dX or name H%dY )" % \
                ((carbon,) * 6)
    group = universe.selectAtoms(selection)

    data = universe.dcd.timeseries(group, format="afc", skip=skip)

    # There are two deuteriums/carbon atom position in each acyl chain
    cd = numpy.concatenate((data[1::3] - data[0::3], data[2::3] - data[0::3]), axis=0)
    del data
    cd_r = numpy.sqrt(numpy.sum(numpy.power(cd, 2), axis=-1))

    # Dot product with the z axis
    cos_theta = cd[..., 2] / cd_r
    S_cd = -0.5 * (3. * numpy.square(cos_theta) - 1)

    # Depending on how you want to treat each deuterium, you can skip the next step
    S_cd.shape = (S_cd.shape[0], S_cd.shape[1] / 4, -1)  # 4 deuterium order parameters/lipid carbon position
    order_param[i] = numpy.average(S_cd)
