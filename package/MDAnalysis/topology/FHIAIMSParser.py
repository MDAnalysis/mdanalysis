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

"""
FHI-AIMS Topology Parser --- :mod:`MDAnalysis.topolgy.FHIAIMSParser`
====================================================================

Reads an `FHI-AIMS`_ ``.in`` file and pulls the atom information from it.
Because an FHI-AIMS input file only has atom name information, any
information about residues and segments will not be populated.

.. _`FHI-AIMS`: https://aimsclub.fhi-berlin.mpg.de/


See Also
--------
:mod:`MDAnalysis.coordinates.FHIAIMS`


Classes
-------

.. autoclass:: FHIAIMSParser
   :members:
   :inherited-members:

"""
import numpy as np

from ..lib.util import openany
from .base import TopologyReaderBase
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomnames,
    Atomids,
    Resids,
    Resnums,
    Segids,
    Elements,
)


class FHIAIMSParser(TopologyReaderBase):
    """Parse a list of atoms from an FHI-AIMS file

    Creates the following attributes:
     - Atomnames


    .. note::

        By default, atomtypes and masses will be guessed on Universe creation.
        This may change in release 3.0.
        See :ref:`Guessers`_ for more information.

    .. versionchanged:: 2.8.0
        Removed type and mass guessing (attributes guessing takes place now
        through universe.guess_TopologyAttrs() API).
    """
    format = ['IN', 'FHIAIMS']

    def parse(self, **kwargs):
        """Read the file and return the structure.

        Returns
        -------
        MDAnalysis Topology object
        """
        # FHIAIMS geometry files are only single frames
        names = []
        skip_tags = ["#", "lattice_vector", "initial_moment", "velocity"]
        with openany(self.filename) as inf:
            for line in inf:
                line = line.strip()
                if line.startswith("atom"):
                    names.append(line.split()[-1])
                    continue
                if any([line.startswith(tag) for tag in skip_tags]):
                    continue
                # we are now seeing something that's neither atom nor lattice
                raise ValueError(
                    'Non-conforming line: ({0})in FHI-AIMS input file {0}'.format(line, self.filename))
            names = np.asarray(names)
            natoms = len(names)

        attrs = [Atomnames(names),
                 Atomids(np.arange(natoms) + 1),
                 Resids(np.array([1])),
                 Resnums(np.array([1])),
                 Segids(np.array(['SYSTEM'], dtype=object)),
                 Elements(names)]

        top = Topology(natoms, 1, 1,
                       attrs=attrs)

        return top
