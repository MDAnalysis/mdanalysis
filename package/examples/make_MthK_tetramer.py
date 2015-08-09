"""Example: building a MthK K-channel tetramer

Rotate the monomer in pdb:3LDD around the 4-fold axis formed by the
K+-ions in the filter.

.. Note::

   This example is not exact; use the BIOMT record in the PDB to
   create the real tetramer.
"""
import MDAnalysis

# need permissive to read HETATM (apparently...)
P = MDAnalysis.Universe('./data/3ldd.pdb', permissive=True)
filterK = P.select_atoms('resname K and resid 1:4')
monomer = P.select_atoms('protein')
axis = (filterK[0], filterK[-1])  # first to last filter ion
monomer.write('A.pdb')
monomer.rotateby(90, axis, filterK)
monomer.write('B.pdb')
monomer.rotateby(90, axis, filterK)
monomer.write('C.pdb')
monomer.rotateby(90, axis, filterK)
monomer.write('D.pdb')

