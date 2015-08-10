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
=========================

.. versionadded:: 0.11.0

The :class:`HoomdXMLParser` generates a topology from files for the HOOMD_ code.

Read a list of atoms from a `HOOMD XML`_ file to build a basic topology.
Atom names are set to atom type if not present (which they probably aren't).
Elements are guessed based on atom types.
Masses and charges are set to zero if not present in the XML file.
Hoomd XML does not identify molecules or residues, so placeholder values
are used for residue numbers and residue names.
Bonds and angles are read if present.

.. _HOOMD: http://codeblue.umich.edu/hoomd-blue/index.html
.. _HOOMD XML: http://codeblue.umich.edu/hoomd-blue/doc/page_xml_file_format.html

Classes
-------

.. autoclass:: HoomdXMLParser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import
from ..lib.util import openany
from ..core.AtomGroup import Atom
from .core import guess_atom_element
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
           file, but it will use placeholder values if they are not present.

           Because Hoomd uses unitless mass, charge, etc., if they are not present
           they will not be guessed - they will be set to zero.

           :Returns: MDAnalysis internal *structure* dict

           .. SeeAlso:: The *structure* dict is defined in
                        :func:`MDAnalysis.topology.base`.

           .. versionadded:: 0.11.0

        """
        with openany(self.filename) as stream:
            tree = ET.parse(stream)
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

        atype = configuration.find('type')
        atypes = atype.text.strip().split('\n')
        if len(atypes) != natoms:
            raise IOError("Number of types does not equal natoms.")

        try:
            name = configuration.find('name')
            names = name.text.strip().split('\n')
        except:
            names = [atype for atype in atypes]
        if len(names) != natoms:
            raise IOError("Number of names does not equal natoms.")

        try:
            elem = configuration.find('element')
            elems = elem.text.strip().split('\n')
        except:
            for atype in atypes:
                elems.append(guess_atom_element(atype))
        if len(elems) != natoms:
            raise IOError("Number of elements does not equal natoms.")

        try:
            mass = configuration.find('mass')
            masses = [float(x) for x in mass.text.strip().split('\n')]
        except:
            masses = [0]*natoms
        if len(masses) != natoms:
            raise IOError("Number of masses does not equal natoms.")

        try:
            charge = configuration.find('charge')
            charges = [float(x) for x in charge.text.strip().split('\n')]
        except:
            charges = [0]*natoms
        if len(charges) != natoms:
            raise IOError("Number of charges does not equal natoms.")

        atoms = []
        bonds = []
        angles = []
        dihedrals = []
        impropers = []

        for i in range(natoms):
            atoms.append(Atom(i, names[i], atypes[i], resname, resid, segid, masses[i], charges[i], universe=self._u))

        try:
            bond = configuration.find('bond')
            bondlines = bond.text.strip().split('\n')
            for bondline in bondlines:
                bondwords = bondline.split()
                bonds.append((int(bondwords[1]),int(bondwords[2])))
        except:
            bonds = []

        try:
            angle = configuration.find('angle')
            anglelines = angle.text.strip().split('\n')
            for angleline in anglelines:
                anglewords = angleline.split()
                angles.append((int(anglewords[1]),int(anglewords[2]),int(anglewords[3])))
        except:
            angles = []

        try:
            torsion = configuration.find('dihedral')
            torsionlines = torsion.text.strip().split('\n')
            for torsionline in torsionlines:
                torsionwords = torsionline.split()
                dihedrals.append((int(torsionwords[1]),int(torsionwords[2]),int(torsionwords[3]),int(torsionwords[4])))
        except:
            dihedrals = []

        try:
            improper = configuration.find('improper')
            improperlines = improper.text.strip().split('\n')
            for improperline in improperlines:
                improperwords = improperline.split()
                impropers.append((int(improperwords[1]),int(improperwords[2]),int(improperwords[3]),int(improperwords[4])))
        except:
            impropers = []


        structure = {'atoms': atoms, 'bonds': bonds, 'angles': angles, 'dihedrals': dihedrals, 'impropers': impropers}
        return structure
