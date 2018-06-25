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

.. _load_crd:

Loading CRD files
-----------------
