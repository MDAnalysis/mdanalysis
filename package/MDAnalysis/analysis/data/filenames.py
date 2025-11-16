# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
r"""Analysis data files
===================

:mod:`MDAnalysis.analysis.data` contains data files that are used as part of
analysis. These can be experimental or theoretical data. Files are stored
inside the package and made accessible via variables in 
:mod:`MDAnalysis.analysis.data.filenames`. These variables are documented
below, including references to the literature and where they are used inside
:mod:`MDAnalysis.analysis`.

Data files
----------

.. data:: Rama_ref

   Reference Ramachandran histogram for
   :class:`MDAnalysis.analysis.dihedrals.Ramachandran`.  The data were
   calculated on a data set of 500 PDB structures taken from
   :footcite:p:`Lovell2003`. This is a numpy array in the :math:`\phi` and
   :math:`\psi` backbone dihedral angles.

   Load and plot it with ::

      import numpy as np
      import matplotlib.pyplot as plt
      from MDAnalysis.analysis.data.filenames import Rama_ref
      X, Y = np.meshgrid(np.arange(-180, 180, 4), np.arange(-180, 180, 4))
      Z = np.load(Rama_ref)
      ax.contourf(X, Y, Z, levels=[1, 17, 15000])

   The given levels will draw contours that contain 90% and 99% of the data
   points.The reference data are shown in 
   :ref:`Ramachandran reference plot figure <figure-rama-ref-plot>`. An example
   of analyzed data together with the reference data are shown in
   :ref:`Ramachandran plot figure <figure-ramachandran>` as an example.
   
   
.. _figure-rama-ref-plot:

.. figure:: /images/rama_ref_plot.png
   :scale: 80%
   :alt: Ramachandran Ref Plot
   
   Reference Ramachandran plot, with contours that contain 90% 
   ("allowed region") and 99% ("generously allowed region") of the data points 
   from the reference data set.

.. data:: Janin_ref

   Reference Janin histogram for :class:`MDAnalysis.analysis.dihedrals.Janin`.
   The data were calculated on a data set of 500 PDB structures taken from
   :footcite:p:`Lovell2003`. This is a numpy array in the :math:`\chi_1` and
   :math:`\chi_2` sidechain dihedral angles.

   Load and plot it with ::

      import numpy as np
      import matplotlib.pyplot as plt
      from MDAnalysis.analysis.data.filenames import Janin_ref
      X, Y = np.meshgrid(np.arange(0, 360, 6), np.arange(0, 360, 6))
      Z = np.load(Janin_ref)
      ax.contourf(X, Y, Z, levels=[1, 6, 600])

   The given levels will draw contours that contain 90% and 98% of the
   data. The reference data are shown in 
   :ref:`Janin reference plot figure <figure-janin-ref-plot>`. An example of 
   analyzed data together with the reference data are shown in
   :ref:`Janin plot figure<figure-janin-ref-plot>` as an example.

.. _figure-janin-ref-plot:

.. figure:: /images/janin_ref_plot.png
   :scale: 85 %
   :alt: Janin Ref Plot

   Janin reference plot with contours that contain 90% ("allowed region") and 
   98% ("generously allowed region") of the data points from the reference data 
   set.

"""


__all__ = [
    "Rama_ref",
    "Janin_ref",
    # reference plots for Ramachandran and Janin classes
]


from importlib import resources

_base_ref = resources.files("MDAnalysis.analysis.data")
Rama_ref = (_base_ref / "rama_ref_data.npy").as_posix()
Janin_ref = (_base_ref / "janin_ref_data.npy").as_posix()

# This should be the last line: clean up namespace
del resources
