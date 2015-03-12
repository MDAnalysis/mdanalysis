NumPy access to Gromacs XDR trajectories --- :mod:`MDAnalysis.coordinates.xdrfile`
==================================================================================

Gromacs_ trajectory formats are based on the XDR_ (eXternal Data
Representation) format. Within this module, classes are provided to
access frames of Gromacs :mod:`~MDAnalysis.coordinates.XTC` or
:mod:`~MDAnalysis.coordinates.TRR` trajectories. Relevant data
(positions, velocities, forces, box dimensions) are made available as
NumPy_ arrays. The module also contains low-level code in
:mod:`~MDAnalysis.coordinates.xdrfile.libxdrfile2` that performs the
actual file access and includes some additions such as trajectory
indices that provide random access for the Gromacs trajectories (which
is not available with current standard Gromacs tools).

.. _Gromacs: http://www.gromacs.org
.. _XDR: http://en.wikipedia.org/wiki/External_Data_Representation
.. _numpy: http://www.numpy.org/

.. rubric:: Contents

.. toctree::
   :maxdepth: 1

   xdrfile/XTC
   xdrfile/TRR
   xdrfile/core
   xdrfile/libxdrfile2
   xdrfile/statno


