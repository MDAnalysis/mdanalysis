# Build a all-atom alpha helix for sequence AAAAAAPAAAAAA
# Oliver Beckstein 2012
# Placed into the Public Domain
#
# Based on
# http://salilab.org/modeller/wiki/Make%20alpha%20helix
# Example for model.build_sequence(), secondary_structure.alpha()

from modeller import *
from modeller.optimizers import conjugate_gradients

log.verbose()

# Set up environment
e = environ()
# use all-hydrogen topology:
e.libs.topology.read('${LIB}/top_allh.lib')
e.libs.parameters.read('${LIB}/par.lib')
e.io.hydrogen = True

# Build an extended chain model from primary sequence
m = model(e)
m.build_sequence('AAAAAAPAAAAAA')

# Make stereochemical restraints on all atoms
allatoms = selection(m)
m.restraints.make(allatoms, restraint_type='STEREO', spline_on_site=False)

# Constrain all residues to be alpha-helical
# (Could also use m.residue_range() rather than m.residues here.)
m.restraints.add(secondary_structure.alpha(m.residues))

# Get an optimized structure with CG, and write it out
cg = conjugate_gradients()
cg.optimize(allatoms, max_iterations=1000)
m.write(file='A6PA6_alpha.pdb')
