# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
HOOMD XML topology parser
========================
Read a list of atoms from a HOOMD XML file to build a basic topology.

Atom atypes, masses, and charges are guessed if not present in the XML file.

Classes
-------

.. autoclass:: HoomdXMLParser
   :members:
   :inherited-members:
"""
from __future__ import absolute_import
from ..core.AtomGroup import Atom
from .core import get_atom_mass, guess_atom_charge, guess_atom_element
from .base import TopologyReader
import xml.etree.ElementTree as ET

class HoomdXMLParser(TopologyReader):
    def parse(self):
        """Parse Hoomd XML file *filename* and return the dict `structure`.
            
            Hoomd XML format does not contain a node for names. The parser will
            look for a name node anyway, and if it doesn't find one, it will use
            the atom types as names. If the Hoomd XML file doesn't contain a type
            node (it should), then all atom types will be \'none\'. Similar to the
            names, the parser will try to read element, mass, and charge from the XML
            file, but it will guess if they are not present.

            :Returns: MDAnalysis internal *structure* dict

            .. SeeAlso:: The *structure* dict is defined in
            :func:`MDAnalysis.topology.base`.

        """
        tree = ET.parse(self.filename)
        root = tree.getroot()
        configuration = root.find('configuration')
        natoms = int(configuration.get('natoms'))

        atypes = []
        elems = []
        masses = []
        charges = []
        resname = "SYSTEM"
        resid = 0
        segid = "SYSTEM"

        try:
            atype = configuration.find('type')
            atypes = atype.text.strip().split('\n')
            if len(atypes) != 0 and len(atypes) != natoms:
                raise IOError("Number of types is neither zero nor natoms")
        except:
            atypes = ['none']*natoms
        try:
            name = configuration.find('name')
            names = name.text.strip().split('\n')
            if len(names) != 0 and len(names) != natoms:
                raise IOError("Number of types is neither zero nor natoms")
        except:
            names = [atype for atype in atypes]
        try:
            elem = configuration.find('element')
            elems = elem.text.strip().split('\n')
            if len(elems) != 0 and len(elems) != natoms:
                raise IOError("Number of elements is neither zero nor natoms")
        except:
            for atype in atypes:
                elems.append(guess_atom_element(atype))
        try:
            mass = configuration.find('mass')
            masses = [float(x) for x in mass.text.strip().split('\n')]
            if len(masses) != 0 and len(masses) != natoms:
                raise IOError("Number of masses is neither zero nor natoms")
        except:
            for elem in elems:
                masses.append(get_atom_mass(elem))

        try:
            charge = configuration.find('charge')
            charges = [float(x) for x in charge.text.strip().split('\n')]
            if len(charges) != 0 and len(charges) != natoms:
                raise IOError("Number of charges is neither zero nor natoms")
        except:
            for name in names:
                charges.append(guess_atom_charge(name))

        atoms = []
        for i in range(natoms):
            atoms.append(Atom(i, names[i], atypes[i], resname, resid, segid, masses[i], charges[i], universe=self._u))
        structure = {'atoms': atoms}

        return structure
