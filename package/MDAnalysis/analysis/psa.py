# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
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

r"""
Calculating path similarity --- :mod:`MDAnalysis.analysis.psa`
==========================================================================

:Author: Sean Seyler
:Year: 2015
:Copyright: GNU Public License v3

.. versionadded:: 0.10.0

.. deprecated:: 2.8.0

  This module is deprecated in favour of the mdakit
  `pathsimanalysis <https://github.com/MDAnalysis/PathSimAnalysis>`_ and
  will be removed in MDAnalysis 3.0.0.


See Also
--------
:mod:`pathsimanalysis.psa`

"""

import warnings

from pathsimanalysis import (
    get_path_metric_func,
    sqnorm,
    get_msd_matrix,
    reshaper,
    get_coord_axes,
    hausdorff,
    hausdorff_wavg,
    hausdorff_avg,
    hausdorff_neighbors,
    discrete_frechet,
    dist_mat_to_vec,
    Path,
    PSAPair,
    PSAnalysis,
)


wmsg = ('Deprecation in version 2.8.0:\n'
        'MDAnalysis.analysis.psa is deprecated in favour of the MDAKit '
        'PathSimAnalysis (https://github.com/MDAnalysis/PathSimAnalysis) '
        'and will be removed in MDAnalysis version 3.0.0')
warnings.warn(wmsg, category=DeprecationWarning)
