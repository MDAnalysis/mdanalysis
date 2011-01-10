.. Contains the formatted docstrings from the KDTree modules located in 'mdanalysis/MDAnalysis/KDTree'

**************************
KDTree module
**************************

:Author: Thomas Hamelryck
:Year:   2002
:Licence: Biopython

The KD tree data structure can be used for all kinds of searches that
involve N-dimensional vectors. For example, neighbor searches (find all points
within a radius of a given point) or finding all point pairs in a set
that are within a certain radius of each other. See "Computational Geometry: 
Algorithms and Applications" (Mark de Berg, Marc van Kreveld, Mark Overmars, 
Otfried Schwarzkopf).

MDAnalysis uses Biopython's KDTree module for distance selections
between atoms but the user is of course welcome to use the module for
any purpose.

.. toctree::
   :maxdepth: 1

   KDTree/KDTree
   KDTree/NeighborSearch

