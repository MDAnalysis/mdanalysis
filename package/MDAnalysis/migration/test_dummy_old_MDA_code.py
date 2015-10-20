import MDAnalysis
from MDAnalysis.tests.datafiles import GRO, XTC
universe = MDAnalysis.Universe(GRO, XTC)

#old selection
all_selection = universe.selectAtoms('all')

#additional old selectAtoms selection (this comment shouldn't be modified despite containing the method name)
all_selection.selectAtoms('bynum 1:10')

#testing atomgroup methods to properties (and exclusion of comments from conversion):

#all_selection.residues()
all_selection.residues()
#all_selection.charges()
all_selection.charges()
#all_selection.indices()
all_selection.indices()
#all_selection.masses()
all_selection.masses()
#all_selection.names()
all_selection.names()
#all_selection.types()
all_selection.types()
#all_selection.radii()
all_selection.radii()
#all_selection.resids()
all_selection.resids()
#all_selection.resnames()
all_selection.resnames()
#all_selection.resnums()
all_selection.resnums()
#all_selection.segids()
all_selection.segids()

#similarly for atomgroup count method renaming:

#all_selection.numberOfAtoms()
all_selection.numberOfAtoms()

#all_selection.numberOfResidues()
all_selection.numberOfResidues()

#all_selection.numberOfSegments()
all_selection.numberOfSegments()

#for old import statements:

#import MDAnalysis.KDTree
import MDAnalysis.KDTree
#from MDAnalysis import KDTree
from MDAnalysis import KDTree

#import MDAnalysis.core.transformations
import MDAnalysis.core.transformations
#from MDAnalysis.core import transformations
from MDAnalysis.core import transformations

#import MDAnalysis.core.util
import MDAnalysis.core.util
#from MDAnalysis.core import util
from MDAnalysis.core import util

#import MDAnalysis.core.log
import MDAnalysis.core.log
#from MDAnalysis.core import log
from MDAnalysis.core import log

#import MDAnalysis.core.units
import MDAnalysis.core.units
#from MDAnalysis.core import units
from MDAnalysis.core import units

#import MDAnalysis.core.distances
import MDAnalysis.core.distances
#from MDAnalysis.core import distances
from MDAnalysis.core import distances

#import MDAnalysis.core.parallel
import MDAnalysis.core.parallel
#from MDAnalysis.core import parallel
from MDAnalysis.core import parallel

# These methods are now properties returning an object
#AtomGroup.bond() -> AtomGroup.bond.value()
AtomGroup.bond()
#AtomGroup.angle() -> AtomGroup.angle.value()
AtomGroup.angle()
#AtomGroup.torsion() -> AtomGroup.dihedral.value()
AtomGroup.torsion()
#AtomGroup.improper() -> AtomGroup.improper.value()
AtomGroup.improper()

#atomgroup, atom and universe torsion to dihedral conversions
#AtomGroup.torsions -> AtomGroup.dihedrals
AtomGroup.torsions
#Atom.torsions -> Atom.dihedrals
Atom.torsions
#Universe.torsions -> Universe.dihedrals
Universe.torsions

#camelcase fixes

# from core.AtomGroup
#totalMass -> total_mass
ag.totalMass
#totalCharge -> total_charge
ag.totalCharge
#centerOfGeometry -> center_of_geometry
ag.centerOfGeometry
#centerOfMass -> center_of_mass
ag.centerOfMass
#radiusOfGyration -> radius_of_gyration
ag.radiusOfGyration
#shapeParameter -> shape_parameter
ag.shapeParameter
#momentOfInertia -> moment_of_inertia
ag.momentOfInertia
#principalAxes -> principal_axes
ag.principalAxes
#packIntoBox -> pack_into_box
ag.packIntoBox
#asUniverse -> as_universe
ag.asUniverse
#align_principalAxis -> align_principal_axis
ag.align_principalAxis

# from lib.distances
#applyPBC -> apply_PBC
lib.distances.applyPBC

#frame_count = universe.trajectory.numframes
frame_count = universe.trajectory.numframes
traj = universe.trajectory
#frame_count = traj.numframes
frame_count = traj.numframes

# From MDAnalysis.lib.distances
#calc_torsions() -> calc_dihedrals()
#from MDAnalysis.lib.distances import calc_torsions
from MDAnalysis.lib.distances import calc_torsions
#MDAnalysis.lib.distances.calc_torsions()
MDAnalysis.lib.distances.calc_torsions()
result = MDAnalysis.lib.distances.calc_torsions()
#dist.calc_torsions()
dist.calc_torsions()

#atomgroup method pluralizations
#set_mass(new) --> set_masses(new)
ag.set_mass(new)

#set_charge(new) --> set_charges(new)
ag.set_charge(new)

#set_name(new) --> set_names(new)
ag.set_name(new)

#set_type(new) --> set_types(new)
ag.set_type(new)

#set_radius(new) --> set_radii(new)
ag.set_radius(new)

#set_bfactor(new) --> set_bfactors(new)
ag.set_bfactor(new)

#set_altloc(new) --> set_altlocs(new)
ag.set_altloc(new)

#set_serial(new) --> set_serials(new)
ag.set_serial(new)

#set_resid(new) --> set_resids(new)
ag.set_resid(new)

#set_resname(new) --> set_resnames(new)
ag.set_resname(new)

#set_resnum(new) --> set_resnums(new)
ag.set_resnum(new)

#set_segid(new) --> set_segids(new)
ag.set_segid(new)

#this test case has caused issues:
g.set_resid(resid * np.ones(len(g)))

#frame numbering is now 0-based:
#ts.frame - 1 -> ts.frame - 0
ts.frame - 1

#ts.frame + 2 -> ts.frame + 3
ts.frame + 2

#ts.frame == 3 -> ts.frame == 2
ts.frame == 3

#ts.frame != 5 -> ts.frame != 4
ts.frame != 5

#another
ts.frame = 9

#+1
[ts.frame for ts in self.trajectory[2:9:3]]

#+1
[ts.frame for ts in self.trajectory]

assert_equal(self.ts.frame, 1, "rewinding to frame 1")

#decoy comment
assert_almost_equal(ts.frame, 544)

assert_almost_equal(ts.dummy, 544)

#frame warning with indentation complexity:
class Dummy(object):
    assert_almost_equal(ts.frame, 544)

    ts.frame = 77

#numatoms to n_atoms keyword argument conversion while preserving the conversion from numberOfAtoms() to n_atoms as well:
with MDAnalysis.Writer(pdbtrj, multiframe=True, bonds=False, numatoms=u.atoms.numberOfAtoms()) as PDB:
    pass

#alternative call syntax:
with MDAnalysis.coordinates.core.writer(pdbtrj, multiframe=True, bonds=False, numatoms=u.atoms.numberOfAtoms()) as PDB:
    pass

#the above fix should be specific to .writer or .Writer, so the following should not be recognized (as a probe for specificity) from the keyword argument standpoint [method replacement is ok]:
with MDAnalysis.coordinates.core.writerr(pdbtrj, multiframe=True, bonds=False, numatoms=u.atoms.numberOfAtoms()) as PDB:
    pass

#however, the fixer should be sufficiently flexible to recognize a different input filename, the omission of default arguments, spacing between 'numatoms' and '=', and an explicit integer value for numatoms, along with some additional kwargs:
with MDAnalysis.Writer(other_filename, numatoms = 55, start = 0, step = 2) as GRO:
    pass
