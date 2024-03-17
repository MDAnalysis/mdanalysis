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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
HOOMD XML topology parser
=========================

.. versionadded:: 0.11.0

The :class:`HoomdXMLParser` generates a topology from files for the HOOMD_ code.

Read a list of atoms from a `HOOMD XML`_ file to build a basic topology.
Masses and charges are set to zero if not present in the XML file.
Hoomd XML does not identify molecules or residues, so placeholder values
are used for residue numbers.
Bonds and angles are read if present.

.. _HOOMD: http://codeblue.umich.edu/hoomd-blue/index.html
.. _HOOMD XML: http://codeblue.umich.edu/hoomd-blue/doc/page_xml_file_format.html

Classes
-------

.. autoclass:: HoomdXMLParser
   :members:
   :inherited-members:

"""
import xml.etree.ElementTree as ET
import numpy as np

from ..lib.util import openany
from .base import TopologyReaderBase
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomtypes,
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


class HoomdXMLParser(TopologyReaderBase):
    """Parses a Hoomd XML file to create a Topology

    Reads the following Attributes:
     - Atomtypes
     - Bonds
     - Angles
     - Dihedrals
     - Impropers
     - Radii
     - Masses

    """
    format = 'XML'

    def parse(self, **kwargs):
        """Parse Hoomd XML file

        Hoomd XML format does not contain a node for names. The parser will
        look for a name node anyway, and if it doesn't find one, it will use
        the atom types as names. If the Hoomd XML file doesn't contain a type
        node (it should), then all atom types will be \'none\'. Similar to the
        names, the parser will try to read atom type, mass, and charge from the XML
        file, but it will use placeholder values if they are not present.

        Because Hoomd uses unitless mass, charge, etc., if they are not present
        they will not be guessed - they will be set to zero.


        .. versionadded:: 0.11.0
        """
        with openany(self.filename) as stream:
            tree = ET.parse(stream)
        root = tree.getroot()
        configuration = root.find('configuration')
        natoms = int(configuration.get('natoms'))

        attrs = {}

        atype = configuration.find('type')
        atypes = atype.text.strip().split('\n')
        if len(atypes) != natoms:
            raise IOError("Number of types does not equal natoms.")
        attrs['types'] = Atomtypes(np.array(atypes, dtype=object))

        for attrname, attr, mapper, dtype in (
                ('diameter', Radii, lambda x: float(x) / 2., np.float32),
                ('mass', Masses, float, np.float64),
                ('charge', Charges, float, np.float32),
        ):
            try:
                val = configuration.find(attrname)
                vals = [mapper(el) for el in val.text.strip().split()]
            except:
                pass
            else:
                attrs[attrname] = attr(np.array(vals, dtype=dtype))
        for attrname, attr, in (
                ('bond', Bonds),
                ('angle', Angles),
                ('dihedral', Dihedrals),
                ('improper', Impropers),
        ):
            try:
                val = configuration.find(attrname)
                vals = [tuple(int(el) for el in line.split()[1:])
                        for line in val.text.strip().split('\n')
                        if line.strip()]
            except:
                vals = []
            attrs[attrname] = attr(vals)

        if 'mass' not in attrs:
            attrs['mass'] = Masses(np.zeros(natoms))
        if 'charge' not in attrs:
            attrs['charge'] = Charges(np.zeros(natoms, dtype=np.float32))

        attrs = list(attrs.values())

        attrs.append(Atomids(np.arange(natoms) + 1))
        attrs.append(Resids(np.array([1])))
        attrs.append(Resnums(np.array([1])))
        attrs.append(Segids(np.array(['SYSTEM'], dtype=object)))

        top = Topology(natoms, 1, 1,
                       attrs=attrs)

        return top
