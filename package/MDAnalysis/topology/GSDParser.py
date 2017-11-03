# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
GSD topology parser
=========================

.. versionadded:: 0.17.0

The :class:`GSDParser` generates a topology from files for the HOOMD_ code.

Parser for the HOOMD-blue GSD topology/trajectory format. The GSD file stores
information on both the topology and the trajectory in the same file, and allows
for varying atom numbers/identities and topologies during the course of the
simulation. At the moment MDAnalysis can deal only with the case in which there
is no variation.

.. _HOOMD: http://codeblue.umich.edu/hoomd-blue/index.html
.. _HOOMD GSD: https://bitbucket.org/glotzer/gsd

Classes
-------

.. autoclass:: GSDParser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import

try :
    import gsd.hoomd
except :
    pass
import numpy as np

from . import guessers
from ..lib.util import openany
from .base import TopologyReaderBase
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomtypes,
    Atomnames,
    Atomids,
    Angles,
    Bonds,
    Charges,
    Dihedrals,
    Impropers,
    Masses,
    Radii,
    Resids,
    Resnums,
    Segids,
)


class GSDParser(TopologyReaderBase):
    """Parses a Hoomd GSD file to create a Topology

    Reads the following Attributes:
     - Atomtypes
     - Bonds
     - Angles
     - Dihedrals
     - Impropers
     - Radii
     - Masses

    """
    format = 'GSD'

    def parse(self):
        """Parse Hoomd GSD file

        .. versionadded:: 0.17.0
        """
        attrs = {}

        with gsd.hoomd.open(self.filename,mode='rb') as t :
            # Here it is assumed that the particle data does not change in the
            # trajectory.
            snap = t[0]

            natoms = snap.particles.N


            ptypes = snap.particles.types
            atypes = [ptypes[idx] for idx in snap.particles.typeid]
            if len(atypes) != natoms:
                raise IOError("Number of types does not equal natoms.")
            attrs['types'] = Atomtypes(np.array(atypes, dtype=object))

            # set radii, masses, charges
            p = snap.particles
            attrs['diameter'] = Radii (np.array(p.diameter / 2.,dtype=np.float32))
            attrs['mass'] = Masses (np.array(p.mass,dtype=np.float64))
            attrs['charge'] = Charges (np.array(p.charge,dtype=np.float32))

            # set bonds, angles, dihedrals, impropers
            for attrname, attr, in (
                    ('bond', Bonds),
                    ('angle', Angles),
                    ('dihedral', Dihedrals),
                    ('improper', Impropers),
            ):
                try:
                    val = getattr(snap,attrname)
                    vals = val.group
                except:
                    pass
                else:
                    attrs[attrname] = attr(vals)

        attrs = list(attrs.values())
        attrs.append(Atomnames(np.array(atypes, dtype=object)))
        attrs.append(Atomids(np.arange(natoms) + 1))
        attrs.append(Resids(np.array([1])))
        attrs.append(Resnums(np.array([1])))
        attrs.append(Segids(np.array(['SYSTEM'], dtype=object)))

        top = Topology(natoms, 1, 1,
                       attrs=attrs)

        return top
