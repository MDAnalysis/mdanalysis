# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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

import mmtf


from . import Universe


def fetch_mmtf(pdb_id):
    """Create a Universe from the RCSB Protein Data Bank using mmtf format

    Parameters
    ----------
    pdb_id : string
        PDB code of the desired data, eg '4UCP'


    Returns
    -------
    MDAnalysis Universe of the corresponding PDB system


    See Also
    --------
    mmtf.fetch : Function for fetching raw mmtf data


    .. versionadded:: 0.16.0
    """
    top = mmtf.fetch(pdb_id)
    return Universe(top)
