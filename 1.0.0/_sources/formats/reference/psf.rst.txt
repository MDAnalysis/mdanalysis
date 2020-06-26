.. -*- coding: utf-8 -*-
.. _PSF-format:

===================================================
PSF (CHARMM, NAMD, or XPLOR protein structure file)
===================================================

.. include:: classes/PSF.txt

A `protein structure file (PSF) <https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node23.html>`_ contains topology information for CHARMM, NAMD, and XPLOR. The MDAnalysis PSFParser only reads information about atoms, bonds, angles, dihedrals, and impropers. While PSF files can include information on hydrogen-bond donor and acceptor groups, MDAnalysis does not read these in.

.. important:: **Atom ids**

    Although PSF files index atoms from 1 in the file, the MDAnalysis PSFParser subtracts 1 to create atom :code:`ids`. This means that if your atom is numbered 3 in your PSF file, it will have an :code:`Atom.id` of 2 in MDAnalysis.

    Atom :code:`indices` are MDAnalysis derived and always index from 0, no matter the file type.

Reading in
==========

PSF files can come in a number of "flavours": STANDARD, EXTENDED, and NAMD. If your file is not a standard file, it must have a NAMD or EXT flag to tell MDAnalysis to how to :ref:`parse the atom section <psf-spec>`. 

As a NAMD file is space-separated, files with missing columns can cause MDAnalysis to read information incorrectly. `This can cause issues for PSF files written from VMD.`_

.. _`This can cause issues for PSF files written from VMD.`: https://github.com/MDAnalysis/mdanalysis/issues/2061

PSF files can encode insertion codes. However, MDAnalysis `does not currently support reading PSF files with insertion codes`_. 

.. _`does not currently support reading PSF files with insertion codes`: https://github.com/MDAnalysis/mdanalysis/issues/2053



.. _psf-spec:

PSF specification
=================

**CHARMM**

Normal (standard) and extended (EXT) PSF format are
supported. CHEQ is supported in the sense that CHEQ data is simply
ignored.


CHARMM Format from ``source/psffres.src``:

CHEQ::

    II,LSEGID,LRESID,LRES,TYPE(I),IAC(I),CG(I),AMASS(I),IMOVE(I),ECH(I),EHA(I)

    standard format:
    (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8,2G14.6)
    (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8,2G14.6)  XPLOR
    expanded format EXT:
    (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8,2G14.6)
    (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A4,1X,2G14.6,I8,2G14.6) XPLOR

no CHEQ::

    II,LSEGID,LRESID,LRES,TYPE(I),IAC(I),CG(I),AMASS(I),IMOVE(I)

    standard format:
    (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8)
    (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8)  XPLOR
    expanded format EXT:
    (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8)
    (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A4,1X,2G14.6,I8) XPLOR

**NAMD**

This format is space separated (see the `release notes for VMD 1.9.1, psfplugin <http://www.ks.uiuc.edu/Research/vmd/current/devel.html>`_). 