.. Contains the formatted docstrings from the coordinates modules located in 'mdanalysis/MDAnalysis/coordinates'

**************************
Topology modules
**************************

The topology module contains the functions to read topology
files. MDAnalysis uses topology files to identify atoms and bonds
between the atoms. It can use topology files from MD packages such as
CHARMM's and NAMD's PSF format or Amber's PRMTOP files. In addition,
it can also glean atom information from single frame coordinate files
such the PDB, CRD, or PQR formats.

Typically, MDAnalysis recognizes formats by the file extension and
hence most users probably do not need to concern themselves with
classes and functions described here. However, if MDAnalysis does not
properly recognize a file format then a user can explicitly set the
topology file format in the *topology_format* keyword argument to
:class:`~MDAnalysis.core.AtomGroup.Universe`.

.. rubric:: Contents

.. toctree::
   :maxdepth: 1

   topology/init
   topology/core
   topology/tables
   topology/PSFParser
   topology/TOPParser
   topology/CRDParser
   topology/GROParser
   topology/PDBParser
   topology/PrimitivePDBParser
   topology/PDBQTParser
   topology/PQRParser
   topology/DMSParser

