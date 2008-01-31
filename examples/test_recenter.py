from Scientific.Geometry import *
from Scientific.Geometry import Transformation
from MDAnalysis import *
import Numeric
from MDAnalysis import _rms_matrix
system = AtomGroup.Universe("data/water1728.psf", "data/water_dyn.dcd.large_pconst_drag")

#ts_transpose = Numeric.transpose(system.dcd.ts._pos)

for a in system._atoms:
    a.mass = 0.72

def recenter1(ref, coord_set, cellsize):
    new_coor = coord_set - ref.astype(Numeric.Float32)
    signs = Numeric.sign(new_coor)
    p = Numeric.greater(Numeric.absolute(new_coor), cellsize/2)
    v = Numeric.not_equal(Numeric.sum(p, axis=1), 0.)
    new_coor -= (p*signs*v[...,Numeric.NewAxis]*cellsize).astype(Numeric.Float32)
    return new_coor

def recenter2(ref, coord_set, cellsize):
    ref = ref[...,Numeric.NewAxis].astype(Numeric.Float32)
    cellsize = cellsize[...,Numeric.NewAxis]
    new_coor = coord_set - ref
    signs = Numeric.sign(new_coor)
    p = Numeric.greater(Numeric.absolute(new_coor), cellsize/2)
    v = Numeric.not_equal(Numeric.sum(p, axis=1), 0.)
    new_coor -= (p*signs*v[...,Numeric.NewAxis]*cellsize).astype(Numeric.Float32)
    return new_coor


def timing_test():
    import timeit
    s1 = """\
from MDAnalysis import AtomGroup
import Numeric
from __main__ import recenter1
system = AtomGroup.Universe("data/water1728.psf", "data/water_dyn.dcd.large_pconst_drag")
ts = system.dcd.ts
cellsize = Numeric.array(ts.dimensions[:3])
ts_transpose = Numeric.transpose(ts._pos)
ref = ts[0]
        """
    t1 = timeit.Timer('recenter1(ref, ts_transpose, cellsize)', s1)
    num = 10000
    print "Recenter1: %.2f usec/pass" % (num*t1.timeit(number=num)/num)
    s2 = """\
from MDAnalysis import AtomGroup
import Numeric
from __main__ import recenter2
system = AtomGroup.Universe("data/water1728.psf", "data/water_dyn.dcd.large_pconst_drag")
ts = system.dcd.ts
cellsize = Numeric.array(ts.dimensions[:3])
ref = ts[0]
        """
    t2 = timeit.Timer('recenter2(ref, ts._pos, cellsize)', s2)
    print "Recenter2: %.2f usec/pass" % (num*t2.timeit(number=num)/num)

def temp():
    radius = min(system.dcd.ts.dimensions[0:3])/2

    from Scientific.Statistics import Histogram
    a = Numeric.zeros(1, Numeric.Float32)

    rdf_sq = Histogram.Histogram(a, 100, (0, radius*radius))
    for i, a in enumerate(system._atoms[:20]):
        print i
        ref = a.pos
        for ts in system.dcd:
            cellsize = Numeric.array(ts.dimensions[0:3])
            new_coor = recenter2(ref, ts._pos, cellsize)
            dist_sq = Numeric.add.reduce(Numeric.power(new_coor,2),axis=0)
            rdf_sq._addData(dist_sq)

    #for ts in system.dcd:
    #    print ts.frame
    #    cellsize = Numeric.array(ts.dimensions[0:3])
    #    for i, ref in enumerate(ts):
    #        new_coor = recenter2(ref, ts._pos, cellsize)
    #        dist = Numeric.sqrt(Numeric.add.reduce(Numeric.power(new_coor,2), axis=0))
            #p = Numeric.less_equal(dist, radius)
    #        rdf.addData(dist)
        #new_coor = recenter2(ref, ts._pos, cellsize)
        #new_coor = Numeric.transpose(recenter(ref, ts_transpose, cellsize))

    #del(dout)

if __name__ == "__main__":
    ts = system.dcd.ts
    cellsize = Numeric.array(ts.dimensions[:3])
    ts_transpose = Numeric.transpose(ts._pos)
    ref = ts[0]
    new_coor1 = recenter1(ref, ts_transpose, cellsize)
    new_coor2 = recenter2(ref, ts._pos, cellsize)
