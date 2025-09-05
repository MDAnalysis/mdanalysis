# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2023 The MDAnalysis Development Team and contributors
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

"""Water dynamics analysis --- :mod:`MDAnalysis.analysis.waterdynamics`
=======================================================================

.. warning::
    This module is deprecated and will be removed in MDAnalysis 3.0.0.
    Please use the dedicated `waterdynamics MDAKit <https://www.mdanalysis.org/waterdynamics/>`_
    instead.

:Author: Alejandro Bernardin
:Year: 2014-2015
:Copyright: GNU Lesser General Public License v2.1 or later (LGPLv2.1+)

.. versionadded:: 0.11.0
.. deprecated:: 2.8.0
    This module is deprecated in favor of the `waterdynamics MDAKit
    <https://www.mdanalysis.org/waterdynamics/>`_ and will be removed in
    MDAnalysis 3.0.0.

This module provides analysis tools for studying water dynamics in molecular
simulations. It has been moved to a separate package for better maintenance
and development.

Migration Guide
--------------
To migrate to the new package:

1. Install the waterdynamics MDAKit::

    pip install waterdynamics

2. Update your imports::

    # Old
    from MDAnalysis.analysis.waterdynamics import WaterOrientationalRelaxation
    
    # New
    from waterdynamics.waterdynamics import WaterOrientationalRelaxation

See Also
--------
* `waterdynamics MDAKit Documentation <https://www.mdanalysis.org/waterdynamics/>`_
* `MDAnalysis Documentation <https://www.mdanalysis.org/>`_
* `MDAnalysis GitHub Repository <https://github.com/MDAnalysis/mdanalysis>`_
"""

import warnings
from typing import Any, Optional, Type, TypeVar, Union

# Type variable for the analysis classes
T = TypeVar('T', bound='AnalysisBase')

# Deprecation warning for the entire module
warnings.warn(
    "The 'MDAnalysis.analysis.waterdynamics' module is deprecated and will be "
    "removed in MDAnalysis 3.0.0. Please install and use the 'waterdynamics' "
    "MDAKit instead (https://www.mdanalysis.org/waterdynamics/).",
    category=DeprecationWarning,
    stacklevel=2
)

# Try to import from the new package
try:
    from waterdynamics.waterdynamics import (
        WaterOrientationalRelaxation,
        AngularDistribution,
        MeanSquareDisplacement,
        SurvivalProbability,
    )
    _HAS_WATERDYNAMICS = True
except ImportError:
    _HAS_WATERDYNAMICS = False
    
    # Create dummy classes for type checking and documentation
    class _DummyAnalysis:
        """Dummy class for when waterdynamics is not installed."""
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            raise ImportError(
                "The waterdynamics MDAKit is not installed. "
                "Please install it with 'pip install waterdynamics'"
            )
    
    # Create dummy classes for each analysis class
    class WaterOrientationalRelaxation(_DummyAnalysis):  # type: ignore
        pass
        
    class AngularDistribution(_DummyAnalysis):  # type: ignore
        pass
        
    class MeanSquareDisplacement(_DummyAnalysis):  # type: ignore
        pass
        
    class SurvivalProbability(_DummyAnalysis):  # type: ignore
        pass
    
    # Show installation instructions
    _INSTALL_MSG = (
        "The waterdynamics MDAKit is required but not installed.\n"
        "Please install it with:\n"
        "    pip install waterdynamics\n\n"
        "For more information, visit: "
        "https://www.mdanalysis.org/waterdynamics/getting_started.html"
    )
    warnings.warn(_INSTALL_MSG, category=ImportWarning, stacklevel=2)
