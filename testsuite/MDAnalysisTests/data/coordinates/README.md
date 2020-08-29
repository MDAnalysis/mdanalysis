Coordinate tests can have special small test-files that cover everything that
needs to be tested. That can mean mail formated files or files that don't
fullfill their standard completely but are still read by most tools (eq. PDB).

Here follows a list of each file with a short description how it was generated
and how the content can be validated.

text.xyz
--------
## Creation
Written with MDAnalysis using the 'create_data.py' script

## Validation
Manually examining the file in a text editor of your choice. The atom names
should be the same as in 'test_topology.pdb'.

test.xyz.bz2
------------
## Creation

    bzip2 test.xyz

## Validation
unpack and check that the content is the same as test.xyz

test.xtc
--------
## Creation
Written with MDAnalysis using the 'create_data.py script

## Validation
With `gmx dump -f test.xtc` you can look at the content of the file in
plain text using Gromacs utilities.

test.trr
--------
## Creation
Written with MDAnalysis using the 'create_data.py script

## Validation
With `gmx dump -f test.trr` you can look at the content of the file in
plain text using Gromacs utilities.
