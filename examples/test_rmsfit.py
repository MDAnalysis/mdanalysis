from MDAnalysis import *
import Numeric
from MDAnalysis import rms_fitting
path = "kalp16"
universe = Universe(path+".psf", path+".dcd")

asel = universe.selectAtoms("segid KALP and name CA")
ref = asel.coordinates()
ref -= asel.centerOfMass().astype(Numeric.Float32)

dout = DCD.DCDWriter("temp.dcd", universe.numberOfAtoms())
ts_transpose = Numeric.transpose(universe.dcd.ts._pos)
universe.dcd.skip = 500

for ts in universe.dcd:
    print ts.frame
    coor = asel.coordinates()
    com = asel.centerOfMass().astype(Numeric.Float32)
    coor -= com
    r = rms_fitting.rms_rotation_matrix(coor, ref, asel.masses())
    # apply to entire system
    newts = Numeric.array(Numeric.transpose(Numeric.matrixmultiply(ts_transpose-com, r))).astype(Numeric.Float32)
    dout.write_next_timestep(DCD.Timestep(newts))

del(dout)
