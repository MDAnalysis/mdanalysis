.. Contains the formatted docstrings from the topology modules located in 'mdanalysis/MDAnalysis/topology'

**************************
Topology modules
**************************

The topology module contains the functions to read topology
files. MDAnalysis uses topology files to identify atoms and bonds
between the atoms. It can use topology files from MD packages such as
CHARMM's and NAMD's PSF format or Amber's PRMTOP files. In addition,
it can also glean atom information from single frame coordinate files
such the PDB, CRD, or PQR formats (see the :ref:`Supported topology
formats`).

Typically, MDAnalysis recognizes formats by the file extension and
hence most users probably do not need to concern themselves with
classes and functions described here. However, if MDAnalysis does not
properly recognize a file format then a user can explicitly set the
topology file format in the *topology_format* keyword argument to
:class:`~MDAnalysis.core.universe.Universe`.

.. rubric:: Topology formats

.. toctree::
   :maxdepth: 1

   topology/init
   topology/CRDParser
   topology/DLPolyParser
   topology/DMSParser
   topology/FHIAIMSParser
   topology/GMSParser
   topology/GROParser
   topology/GSDParser
   topology/HoomdXMLParser
   topology/ITPParser
   topology/LAMMPSParser
   topology/MinimalParser
   topology/MMTFParser
   topology/MOL2Parser
   topology/PDBParser
   topology/ExtendedPDBParser
   topology/PDBQTParser
   topology/PQRParser
   topology/PSFParser
   topology/TOPParser
   topology/TPRParser
   topology/TXYZParser
   topology/XYZParser

.. rubric:: Topology core modules

The remaining pages are primarily of interest to developers as they
contain functions and classes that are used in the implementation of
the topology readers.

.. toctree::
   :maxdepth: 1

   topology/base
   topology/core
   topology/tpr_util
