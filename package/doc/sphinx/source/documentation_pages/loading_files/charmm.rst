.. _load_charmm:

####################################
Loading CHARMM files with MDAnalysis
####################################

.. _load_psf:

Loading PSF files
-----------------

PSF_ files from CHARMM, NAMD (whitespace separated)
and XPLOR are supported as topology files.
Both standard and extended (EXT) PSF formats are supported.
Universes created from a PSF file will have the following attributes:
 - ids
 - atom names
 - atom types
 - masses
 - charges
 - resids
 - resnames
 - segids
 - bonds
 - angles
 - dihedrals
 - impropers

For implementation details, see :mod:`MDAnalysis.topology.PSFParser`.

.. _PSF: http://www.charmm.org/documentation/c35b1/struct.html

.. _load_dcd:

Loading DCD files
-----------------

A DCD file can be used as a minimal topology file, where only the number of atoms
is read.

Generally, DCD trajectories produced by any code can be read to provide coordinates
although there can be issues with the unitcell (simulation box) representation
(see :attr:`MDAnalysis.coordinates.DCD.DCDReader.dimensions`).

DCDs can also be written but the :class:`MDAnalysis.coordinates.DCD.DCDWriter`
follows recent NAMD/VMD convention for the unitcell but still writes AKMA time.
Reading and writing these trajectories within MDAnalysis will work seamlessly but
if you process those trajectories
with other tools you might need to watch out that time and unitcell dimensions
are correctly interpreted.

For full implementation details see :mod:`MDAnalysis.coordinates.DCD`.

.. _load_crd:

Loading CRD files
-----------------
