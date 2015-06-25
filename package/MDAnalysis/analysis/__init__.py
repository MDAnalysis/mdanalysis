# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
:mod:`MDAnalysis.analysis` --- Analysis code based on MDAnalysis
================================================================

The :mod:`MDAnalysis.analysis` sub-package contains various recipes and
algorithms that can be used to analyze MD trajectories.

If you use them please check if the documentation mentions any specific caveats
and also if there are any published papers associated with these algorithms.

Available analysis modules
--------------------------

:mod:`~MDAnalysis.analysis.align`
    Fitting and aligning of coordinate frames, including the option to
    use a sequence alignment to define equivalent atoms to fit on.

:mod:`~MDAnalysis.analysis.contacts`
    Analyse the number of native contacts relative to a reference
    state, also known as a "q1-q2" analysis.

:mod:`~MDAnalysis.analysis.density`
    Creating and manipulating densities such as the density ow water
    molecules around a protein. Makes use of the external
    GridDataFormats_ package.

:mod:`~MDAnalysis.analysis.distances`
    Functions to calculate distances between atoms and selections; it
    contains the often-used
    :func:`~MDAnalysis.analysis.distances.distance_array` function.

:mod:`~MDAnalysis.analysis.hbonds`
    Analyze hydrogen bonds, including both the per frame results as well
    as the dynamic properties and lifetimes.

:mod:`~MDAnalysis.analysis.helanal`
    Analysis of helices with the HELANAL_ algorithm.

:mod:`~MDAnalysis.analysis.hole`
    Run and process output from the :program:`HOLE` program
    to analyze pores, tunnels and cavities in proteins.

:mod:`~MDAnalysis.analysis.gnm`
    Gaussian normal mode analysis of MD trajectories with the
    help of an elastic network.

:mod:`~MDAnalysis.analysis.leaflet`
    Find lipids in the upper and lower (or inner and outer) leaflet of
    a bilayer; the algorithm can deal with any deformations as long as
    the two leaflets are topologically distinct.

:mod:`~MDAnalysis.analysis.nuclinfo`
    Analyse the nucleic acid for the backbone dihedrals, chi, sugar
    pucker, and Watson-Crick distance (minor and major groove
    distances).

:mod:`~MDAnalysis.analysis.psa`
    Perform Path Similarity Analysis (PSA) on a set of trajectories to measure
    their mutual similarities, including the ability to perform hierarchical
    clustering and generate heat map-dendrogram plots.

:mod:`~MDAnalysis.analysis.rms`
    Calculation of RMSD and RMSF.

:mod:`~MDAnalysis.analysis.waterdynamics`
    Analysis of water.

:mod:`~MDAnalysis.analysis.x3dna`
    Analysis of helicoidal parameters driven by X3DNA_.

.. _GridDataFormats: https://github.com/orbeckst/GridDataFormats
.. _HELANAL: http://www.ccrnp.ncifcrf.gov/users/kumarsan/HELANAL/helanal.html
.. _X3DNA: http://x3dna.org/

.. versionchanged:: 0.10.0
   The analysis submodules are not automatically imported any more. Manually
   import any submodule that you need.

"""

__all__ = [
    'align', 'contacts', 'density', 'distances',
    'gnm', 'hbonds', 'helanal', 'hole', 'leaflet',
    'nuclinfo', 'psa', 'rms', 'waterdynamics',
    'x3dna',
]

