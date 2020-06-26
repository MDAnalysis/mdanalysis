.. -*- coding: utf-8 -*-
.. _TPR-format:

====================================
TPR (GROMACS run topology files)
====================================

.. include:: classes/TPR.txt

A GROMACS TPR file is a portable binary run input file. It contains both topology and coordinate information. However, MDAnalysis currently only reads topology information about atoms, bonds, dihedrals, and impropers; it does not read the coordinate information.  

.. important:: **Atom ids, residue resids, and molnums**

    GROMACS indexes atom numbers and residue numbers from 1 in user-readable files. However, the MDAnalysis TPRParser creates atom :code:`ids` and residue :code:`resids` *from 0*. This means that if your atom is numbered 3 in your GRO, ITP, or TOP file, it will have an :code:`Atom.id` of 2 in MDAnalysis. Likewise, if your residue ALA has a residue number of 4 in your GRO file, it will have a :code:`Residue.resid` number of 3 in MDAnalysis. Finally, molecules are also numbered from 0, in the attribute :code:`molnums`.

    Atom :code:`indices` and residue :code:`resindices` are MDAnalysis derived and always index from 0, no matter the file type.

Supported versions
==================

.. table:: TPR format versions and generations read by :func:`MDAnalysis.topology.TPRParser.parse`.

   ========== ============== ==================== =====
   TPX format TPX generation Gromacs release      read
   ========== ============== ==================== =====
   ??         ??             3.3, 3.3.1           no

   58         17             4.0, 4.0.2, 4.0.3,   yes
                             4.0.4, 4.0.5, 4.0.6,
                             4.0.7

   73         23             4.5.0, 4.5.1, 4.5.2, yes
                             4.5.3, 4.5.4, 4.5.5

   83         24             4.6, 4.6.1           yes

   100        26             5.0, 5.0.1, 5.0.2,   yes
                             5.0.3,5.0.4, 5.0.5

   103        26             5.1                  yes

   110        26             2016                 yes
   112        26             2018                 yes
   116        26             2019                 yes
   ========== ============== ==================== =====

For further discussion and notes see `Issue 2`_. Please *open a new issue* in
the `Issue Tracker`_ when a new or different TPR file format version should be
supported.

TPR specification
=================

The TPR reader is a pure-python implementation of a basic TPR
parser. Currently the following sections of the topology are parsed:

* Atoms: number, name, type, resname, resid, segid, mass, charge,
  [residue, segment, radius, bfactor, resnum, moltype]
* Bonds
* Angles
* Dihedrals
* Impropers

Bonded interactions available in Gromacs are described in the
`Gromacs manual`_. The following ones are used to build the topology (see
`Issue 463`_):

.. _Gromacs: http://www.gromacs.org
.. _`Gromacs manual`: http://manual.gromacs.org/current/reference-manual/index.html
.. _TPR file: http://manual.gromacs.org/current/reference-manual/file-formats.html#tpr
.. _`Issue Tracker`: https://github.com/MDAnalysis/mdanalysis/issues
.. _`Issue 2`: https://github.com/MDAnalysis/mdanalysis/issues/2
.. _`Issue 463`: https://github.com/MDAnalysis/mdanalysis/pull/463

.. table:: GROMACS entries used to create bonds.

    =============	======	=====================================================	
    Directive   	 Type 	 Description                                         	
    =============	======	=====================================================	
    bonds       	 1    	 `regular bond`_                                     	
    bonds       	 2    	 `G96 bond`_                                         	
    bonds       	 3    	 `Morse bond`_                                       	
    bonds       	 4    	 `cubic bond`_                                       	
    bonds       	 5    	 `connections`_                                      	
    bonds       	 6    	 `harmonic potentials`_                              	
    bonds       	 7    	 `FENE bonds`_                                       	
    bonds       	 8    	 `tabulated potential with exclusion/connection`_    	
    bonds       	 9    	 `tabulated potential without exclusion/connection`_ 	
    bonds       	 10   	 `restraint potentials`_                             	
    constraints 	 1    	 `constraints with exclusion/connection`_            	
    constraints 	 2    	 `constraints without exclusion/connection`_         	
    settles     	 1    	 `SETTLE constraints`_                               	
    =============	======	=====================================================	

.. _`regular bond`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#harmonic-potential
.. _`G96 bond`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#fourth-power-potential
.. _`Morse bond`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#morse-potential-bond-stretching
.. _`cubic bond`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#cubic-bond-stretching-potential
.. _`connections`: http://manual.gromacs.org/current/reference-manual/topologies/molecule-definition.html#exclusions
.. _`harmonic potentials`: http://manual.gromacs.org/current/reference-manual/functions/restraints.html#distance-restraints
.. _`FENE bonds`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#fene-bond-stretching-potential
.. _`tabulated potential with exclusion/connection`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#tabulated-bonded-interaction-functions
.. _`tabulated potential without exclusion/connection`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#tabulated-bonded-interaction-functions
.. _`restraint potentials`: http://manual.gromacs.org/current/reference-manual/functions/restraints.html#distance-restraints
.. _`constraints with exclusion/connection`: http://manual.gromacs.org/current/reference-manual/functions/free-energy-interactions.html#constraints
.. _`constraints without exclusion/connection`: http://manual.gromacs.org/current/reference-manual/functions/free-energy-interactions.html#constraints
.. _`SETTLE constraints`: http://manual.gromacs.org/current/reference-manual/algorithms/constraint-algorithms.html#settle


.. table:: GROMACS entries used to create angles.

    ===========	======	=================================	
    Directive 	 Type 	 Description                     	
    ===========	======	=================================	
    angles    	 1    	 `regular angle`_                	
    angles    	 2    	 `G96 angle`_                    	
    angles    	 3    	 `Bond-bond cross term`_         	
    angles    	 4    	 `Bond-angle cross term`_        	
    angles    	 5    	 `Urey-Bradley`_                 	
    angles    	 6    	 `Quartic angles`_               	
    angles    	 8    	 `Tabulated angles`_             	
    angles    	 10   	 `restricted bending potential`_ 	
    ===========	======	=================================	


.. _`regular angle`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#harmonic-angle-potential
.. _`G96 angle`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#cosine-based-angle-potential
.. _`Bond-bond cross term`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#bond-bond-cross-term
.. _`Bond-angle cross term`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#bond-angle-cross-term
.. _`Urey-Bradley`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#urey-bradley-potential
.. _`Quartic angles`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#quartic-angle-potential
.. _`Tabulated angles`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#tabulated-bonded-interaction-functions
.. _`restricted bending potential`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#restricted-bending-potential

.. table:: GROMACS entries used to create dihedrals.

    ===========	======	=======================================	
    Directive 	 Type 	 Description                           	
    ===========	======	=======================================	
    dihedrals 	 1    	 `proper dihedral`_                    	
    dihedrals 	 3    	 `Ryckaert-Bellemans dihedral`_        	
    dihedrals 	 5    	 `Fourier dihedral`_                   	
    dihedrals 	 8    	 `Tabulated dihedral`_                 	
    dihedrals 	 9    	 `Periodic proper dihedral`_           	
    dihedrals 	 10   	 `Restricted dihedral`_                	
    dihedrals 	 11   	 `Combined bending-torsion potential`_ 	
    ===========	======	=======================================	

.. _`proper dihedral`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#proper-dihedrals
.. _`Ryckaert-Bellemans dihedral`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#proper-dihedrals-ryckaert-bellemans-function
.. _`Fourier dihedral`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#proper-dihedrals-fourier-function
.. _`Tabulated dihedral`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#tabulated-bonded-interaction-functions
.. _`Periodic proper dihedral`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#proper-dihedrals-periodic-type
.. _`Restricted dihedral`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#proper-dihedrals-restricted-torsion-potential
.. _`Combined bending-torsion potential`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#proper-dihedrals-combined-bending-torsion-potential

.. table:: GROMACS entries used to create improper dihedrals.

    ===========	======	===============================	
    Directive 	 Type 	 Description                   	
    ===========	======	===============================	
    dihedrals 	 2    	 `improper dihedral`_          	
    dihedrals 	 4    	 `periodic improper dihedral`_ 	
    ===========	======	===============================	

.. _`improper dihedral`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#improper-dihedrals-harmonic-type
.. _`periodic improper dihedral`: http://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html#improper-dihedrals-periodic-type


Developer notes
===============

This tpr parser is written according to the following files

- :file:`{gromacs_dir}/src/kernel/gmxdump.c`
- :file:`{gromacs_dir}/src/gmxlib/tpxio.c` (the most important one)
- :file:`{gromacs_dir}/src/gmxlib/gmxfio_rw.c`
- :file:`{gromacs_dir}/src/gmxlib/gmxfio_xdr.c`
- :file:`{gromacs_dir}/include/gmxfiofio.h`

or their equivalent in more recent versions of Gromacs.

The function :func:`read_tpxheader` is based on the
`TPRReaderDevelopment`_ notes.  Functions with names starting with
``read_`` or ``do_`` are trying to be similar to those in
:file:`gmxdump.c` or :file:`tpxio.c`, those with ``extract_`` are new.

Wherever ``fver_err(fver)`` is used, it means the tpx version problem
has not been solved. Versions prior to Gromacs 4.0.x are not supported.

.. _TPRReaderDevelopment: https://github.com/MDAnalysis/mdanalysis/wiki/TPRReaderDevelopment