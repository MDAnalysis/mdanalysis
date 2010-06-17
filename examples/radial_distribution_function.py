# Calculating G(r) (radial distribution function) of a water box, 
# taking into account periodic boundaries

from MDAnalysis import *
import numpy

universe = Universe(...)
group = universe.selectAtoms("resname TIP3 and name OH2")

dmin, dmax = 1.0, 5.0
rdf, edges = numpy.histogramdd([0], bins=100, range=[(dmin, dmax)])
rdf *= 0

for ts in universe.dcd:
    box = ts.dimensions[:3]
    coor = group.coordinates()
    dist = distances.distance_array(coor, coor, box)
    new_rdf, edges = numpy.histogramdd(numpy.ravel(dist), bins=100, range=[(dmin, dmax)])
    rdf += new_rdf

numframes = universe.dcd.numframes / universe.dcd.skip

#Normalize RDF
density = 1
radii = edges[0]
vol = (4./3.)*numpy.pi*density*(numpy.power(radii[1:],3)-numpy.power(radii[:-1], 3))
rdf /= vol*numframes


