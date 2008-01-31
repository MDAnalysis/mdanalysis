from Scientific.Geometry import *
from Scientific.Geometry import Transformation
from MDAnalysis import *
import Numeric
from MDAnalysis import _rms_matrix
kalp = 23
system = AtomGroup.System("data/kalp%d.psf"%kalp, "data/kalp%d.dcd"%kalp)
asel = system.selectAtoms("segid KALP and resid 3 and not backbone")
ref = asel.coordinates()
ref -= asel.centerOfMass().astype(Numeric.Float32)

dout = DCD.DCDWriter("temp.dcd", len(system._atoms))
ts_transpose = Numeric.transpose(system._dcd.ts._pos)

for ts in system._dcd:
    print ts.frame
    coor = asel.coordinates()
    com = asel.centerOfMass().astype(Numeric.Float32)
    coor -= com
    r = _rms_matrix._rms_rotation_matrix(coor, ref, asel.masses())
    # apply to entire system
    newts = Numeric.array(Numeric.transpose(Numeric.matrixmultiply(ts_transpose-com, r))).astype(Numeric.Float32)
    dout.write_next_timestep(DCD.Timestep(newts))

del(dout)
