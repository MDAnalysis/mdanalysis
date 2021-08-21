.. module:: MDAnalysis.converters

.. _`Converter modules`:

**************************
Converter modules
**************************

.. versionadded:: 2.0.0

The :mod:`MDAnalysis.converters` module contains the Converter classes that
MDAnalysis uses to convert MDAnalysis  structures to and from other Python
packages. 

If you are converting *to* MDAnalysis, you can use the normal syntax for 
creating a Universe from files. Typically MDAnalysis will recognise which 
library it is dealing with, so you will not need to pass in a ``format`` keyword.

For example::

    import MDAnalysis as mda
    from MDAnalysis import GRO
    import parmed as pmd

    pgro = pmd.load_file(GRO)  # creates parmed structure
    ugro = mda.Universe(pgro)  # you can just pass it in

If you are converting *from* MDAnalysis, use the 
:func:`~MDAnalysis.core.groups.AtomGroup.convert_to` method. With this, 
you will have to specify a package name (case-insensitive). ::

    pgro2 = ugro.atoms.convert_to('PARMED')  # converts back to parmed structure

Another syntax is also available for tab-completion support::

    pgro2 = ugro.atoms.convert_to.parmed()


.. rubric:: Available converters

.. toctree::
   :maxdepth: 1

   converters/ParmEd
   converters/RDKit
   converters/OpenMM

.. rubric:: Converter functionalities

.. toctree::
    :maxdepth: 1

    core/accessors

