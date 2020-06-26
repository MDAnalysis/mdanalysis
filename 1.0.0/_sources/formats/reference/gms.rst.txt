.. -*- coding: utf-8 -*-
.. _GMS-format:

=======================
GMS (Gamess trajectory)
=======================

.. include:: classes/GMS.txt

The GMS output format is a common output format for different
GAMESS distributions: `GAMESS-US`_, `Firefly`_ (PC-GAMESS) and `GAMESS-UK`_. 
The current version has been tested with US GAMESS and Firefly only.

.. _`GAMESS-US`: http://www.msg.ameslab.gov/gamess/
.. _Firefly: http://classic.chem.msu.su/gran/gamess/index.html
.. _`GAMESS-UK`: http://www.cfs.dl.ac.uk/

Reading in
==========

MDAnalysis can read a GAMESS output file and pull out atom information. Atom names are their elements. Information about residues and segments is not read.