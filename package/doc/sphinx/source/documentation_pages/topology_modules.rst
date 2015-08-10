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
:class:`~MDAnalysis.core.AtomGroup.Universe`.

.. rubric:: Topology formats

.. toctree::
   :maxdepth: 1

   topology/init
   topology/PSFParser
   topology/TOPParser
   topology/CRDParser
   topology/GROParser
   topology/PDBParser
   topology/PrimitivePDBParser
   topology/ExtendedPDBParser
   topology/PDBQTParser
   topology/PQRParser
   topology/DLPolyParser
   topology/DMSParser
   topology/TPRParser
   topology/MOL2Parser
   topology/XYZParser
   topology/LAMMPSParser
   topology/GMSParser
   topology/HoomdXMLParser

.. rubric:: Topology core modules

The remaining pages are primarily of interest to developers as they
contain functions and classes that are used in the implementation of
the topology readers.

.. toctree::
   :maxdepth: 1

   topology/base
   topology/core
   topology/tables
   topology/tpr_util
