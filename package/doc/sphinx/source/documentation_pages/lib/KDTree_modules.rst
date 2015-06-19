.. Contains the formatted docstrings from the KDTree modules located in 'mdanalysis/MDAnalysis/lib/KDTree'

**************************
KDTree module
**************************

:Author: Thomas Hamelryck
:Year:   2002
:Licence: Biopython

The KD tree data structure can be used for all kinds of searches that
involve N-dimensional vectors. For example, neighbor searches (find all points
within a radius of a given point) or finding all point pairs in a set
that are within a certain radius of each other.

MDAnalysis uses `Biopython`_'s KDTree module for distance selections
between atoms but its use goes far beyond this narrow application.

.. SeeAlso:: The algorithms are based on [deBerg2000]_ with improvements
             suggested by [Bentley1990]_. See the documentation of the
             individual classes for details.

.. rubric:: References

.. [deBerg2000]  Mark de Berg, Marc van Kreveld, Mark Overmars,  
                 Otfried Schwarzkopf. 
                 `Computational Geometry: Algorithms and Applications`_.
                 Springer, 2nd edition. 2000.
.. [Bentley1990] J.L. Bentley, "Kd trees for semidynamic point sets," in Sixth
                 Annual ACM Symposium on Computational Geometry, vol. 91. San
                 Francisco, 1990 


.. _`Computational Geometry: Algorithms and Applications`:
   http://books.google.co.uk/books?id=C8zaAWuOIOcC
.. _Biopython: http://biopython.org

.. rubric:: Contents

.. toctree::
   :maxdepth: 1

   ./KDTree/KDTree
   ./KDTree/NeighborSearch

