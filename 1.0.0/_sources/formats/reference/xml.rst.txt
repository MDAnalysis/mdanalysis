.. -*- coding: utf-8 -*-
.. _XML-format:

=======================
XML (HOOMD)
=======================

    .. include:: classes/XML.txt

MDAnalysis can read topology informatin from a `HOOMD`_ `XML`_ file.
Masses and charges are set to zero if not present in the XML file.
Hoomd XML does not identify molecules or residues, so placeholder values
are used for residue numbers.
Bonds and angles are read if present.

.. _HOOMD: http://codeblue.umich.edu/hoomd-blue/index.html
.. _XML: http://codeblue.umich.edu/hoomd-blue/doc/page_xml_file_format.html

Hoomd XML format does not contain a node for names. The parser will
look for a name node anyway, and if it doesn't find one, it will use
the atom types as names. If the Hoomd XML file doesn't contain a type
node (it should), then all atom types will be \'none\'. 

Similar to the names, the parser will try to read atom type, mass, and charge from the XML
file. Therefore, they are not present, masses and charges will not be guessed. Instead, they will be set to zero, as Hoomd uses unitless mass, charge, etc.