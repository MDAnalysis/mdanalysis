import numpy as np
from MDAnalysis import *

universe = Universe("system.gro")

tail_carbons = np.arange(2, 15)

skip = 5
order_param = np.zeros(len(tail_carbons))

for i, carbon in enumerate(tail_carbons):
    selection = "resname DMPC and ( name C2%d or name H%dR or name H%dS or name C3%d or name H%dX or name H%dY )" % \
                ((carbon,) * 6)
    group = universe.select_atoms(selection)

    data = universe.dcd.timeseries(group, format="afc", skip=skip)

    # There are two deuteriums/carbon atom position in each acyl chain
    cd = np.concatenate((data[1::3] - data[0::3], data[2::3] - data[0::3]), axis=0)
    del data
    cd_r = np.sqrt(np.sum(np.power(cd, 2), axis=-1))

    # Dot product with the z axis
    cos_theta = cd[..., 2] / cd_r
    S_cd = -0.5 * (3. * np.square(cos_theta) - 1)

    # Depending on how you want to treat each deuterium, you can skip the next step
    S_cd.shape = (S_cd.shape[0], S_cd.shape[1] / 4, -1)  # 4 deuterium order parameters/lipid carbon position
    order_param[i] = np.average(S_cd)
