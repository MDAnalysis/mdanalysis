==========================
 Overview over MDAnalysis
==========================

**MDAnalysis** is a Python package that provides classes to access
data in molecular dynamics trajectories. It is object oriented so it
treats atoms, groups of atoms, trajectories, etc as different
objects. Each object has a number of operations defined on itself
(also known as "methods") and also contains values describing the
object ("attributes"). For example, a
:class:`~MDAnalysis.core.AtomGroup.AtomGroup` object has a
:meth:`~MDAnalysis.core.AtomGroup.AtomGroup.center_of_mass` method that
returns the center of mass of the group of atoms. It also contains an
attribute called :attr:`~MDAnalysis.core.AtomGroup.AtomGroup.residues`
that lists all the residues that belong to the group. Using methods
such as :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.select_atoms`
(which uses `CHARMM-style`_ atom :ref:`selection-commands-label`) one
can create new objects (in this case, another
:class:`~MDAnalysis.core.AtomGroup.AtomGroup`).

A typical usage pattern is to iterate through a trajectory and analyze
coordinates for every frame. In the following example the end-to-end distance
of a protein and the radius of gyration of the backbone atoms are calculated::

  import MDAnalysis
  from MDAnalysis.tests.datafiles import PSF,DCD   # test trajectory
  import numpy.linalg
  u = MDAnalysis.Universe(PSF,DCD)                 # always start with a Universe
  nterm = u.s4AKE.N[0]   # can access structure via segid (s4AKE) and atom name
  cterm = u.s4AKE.C[-1]  # ... takes the last atom named 'C'
  bb = u.select_atoms('protein and backbone')  # a selection (a AtomGroup)
  for ts in u.trajectory:     # iterate through all frames
    r = cterm.pos - nterm.pos # end-to-end vector from atom positions
    d = numpy.linalg.norm(r)  # end-to-end distance
    rgyr = bb.radius_of_gyration()  # method of a AtomGroup; updates with each frame
    print "frame = %d: d = %f Angstroem, Rgyr = %f Angstroem" % (ts.frame, d, rgyr)


.. _NumPy:   http://numpy.scipy.org
.. _CHARMM:  http://www.charmm.org/
.. _LAMMPS:  http://lammps.sandia.gov/
.. _NAMD:    http://www.ks.uiuc.edu/Research/namd/
.. _Gromacs: http://www.gromacs.org/

.. _CHARMM-style:
   http://www.charmm.org/documentation/c37b1/select.html

.. TODO: more about philosophy etc... copy and paste from paper

Using MDAnalysis in python
==========================

If you've installed MDAnalysis in the standard python modules location, load
from within the interpreter::

 from MDAnalysis import *

or ::
 
 import MDAnalysis

The idea behind MDAnalysis is to get trajectory data into NumPy_
:class:`numpy.ndarray` arrays, where it can then be easily manipulated using
all the power in NumPy_ and SciPy_. 

MDAnalysis works well both in scripts and in interactive use. The developers
very much recommend using MDAnalysis from within the IPython_ Python shell.  It
allows one to interactively explore the objects (using TAB-completion and
online help), do analysis and immediately plot results. The examples in this manual
are typically run from an interactive :program:`ipython` session.

Invariably, a MDAnalysis session starts with loading data into the
:class:`~MDAnalysis.core.AtomGroup.Universe` class (which can be accessed
as :class:`MDAnalysis.Universe`)::

 from MDAnalysis import *
 universe = Universe(topology, trajectory)

- The *topology* file lists the atoms and residues (and also their
  connectivity). It can be a CHARMM/XPLOR/NAMD PSF file or a coordinate file
  such as a Protein Databank Brookhaven PDB file, a CHARMM card coordinate file
  (CRD), or a GROMOS/Gromacs GRO file.

- The *trajectory* contains a list of coordinates in the order defined in the
  *topology*. It can either be a single frame (PDB, CRD, and GRO are all read)
  or a time series of coordinate frames such as a CHARMM/NAMD/LAMMPS DCD
  binary file, a Gromacs XTC/TRR trajectory, or a XYZ trajectory (possibly
  compressed with gzip or bzip2).

For the remainder of this introduction we are using a short example trajectory
that is provided with MDAnalysis (as part of the `MDAnalysis test suite`_). The
trajectory is loaded with ::
 
  >>> from MDAnalysis import Universe
  >>> from MDAnalysis.tests.datafiles import PSF,DCD
  >>> u = Universe(PSF, DCD)

(The ``>>>`` signs are the Python input prompt and are not to be typed; they
just make clear in the examples what is input and what is output.)

The :class:`~MDAnalysis.core.AtomGroup.Universe` contains a number of important attributes,
the most important ones of which is
:attr:`~MDAnalysis.core.AtomGroup.Universe.atoms`::

  >>> print u.atoms
  <AtomGroup with 3341 atoms>

:attr:`Universe.atoms` is a
:class:`~MDAnalysis.core.AtomGroup.AtomGroup` and can be thought of as
list consisting of :class:`~MDAnalysis.core.AtomGroup.Atom`
objects. The :class:`~MDAnalysis.core.AtomGroup.Atom` is the
elementary and fundamental object in MDAnalysis.

The :attr:`MDAnalysis.Universe.trajectory` attribute gives access to the coordinates
over time::

  >>> print u.trajectory
  < DCDReader '/..../MDAnalysis/tests/data/adk_dims.dcd' with 98 frames of 3341 atoms (0 fixed) >

Finally, the :meth:`MDAnalysis.Universe.select_atoms` method generates a new
:class:`~MDAnalysis.core.AtomGroup.AtomGroup` according to a selection criterion

  >>> calphas = u.select_atoms("name CA")
  >>> print calphas
  <AtomGroup with 214 atoms>

as described in :ref:`selection-commands-label`.

.. _SciPy: http://www.scipy.org/
.. _IPython: http://ipython.scipy.org/
.. _MDAnalysis test suite: https://github.com/MDAnalysis/mdanalysis/wiki/UnitTests


Examples
========

The easiest way to get started with MDAnalysis is to read this
introduction and the chapter on :ref:`selection-commands-label` and then
explore the package interactively in IPython_ or another interactive
Python interpreter.

Included trajectories
---------------------

MDAnalysis comes with a number of real trajectories for testing. You
can also use them to explore the functionality and ensure that
everything is working properly::

  from MDAnalysis import *
  from MDAnalysis.tests.datafiles import PSF,DCD, PDB,XTC
  u_dims_adk = Universe(PSF,DCD)
  u_eq_adk = Universe(PDB, XTC)

The PSF and DCD file are a closed-form-to-open-form transition of
Adenylate Kinase (from [Beckstein2009]_) and the PDB+XTC file are ten
frames from a Gromacs simulation of AdK solvated in TIP4P water with
the OPLS/AA force field.

.. [Beckstein2009] O. Beckstein, E.J. Denning, J.R. Perilla, and
                   T.B. Woolf. Zipping and Unzipping of Adenylate
                   Kinase: Atomistic Insights into the Ensemble of
                   Open <--> Closed Transitions. *J Mol Biol* **394**
                   (2009), 160--176, doi:`10.1016/j.jmb.2009.09.009`_

.. _`10.1016/j.jmb.2009.09.009`: http://dx.doi.org/10.1016/j.jmb.2009.09.009

Code snippets
-------------

The source code distribution comes with a directory `examples`_ that
contains a number of code snippets that show how to use certain
aspects of MDAnalysis. 

For instance, there is code that shows how to

* fit a trajectory to a reference structure using the QCP
  RMSD-alignment code in :mod:`MDAnalysis.core.qcprot`
  (`rmsfit_qcp.py`_);

* do a block-averaging error analysis (`blocks.py`_);

* calculate a potential profile across a membrane (`potential_profile.py`_);

* do a native contact analysis using :mod:`MDAnalysis.analysis.contacts` (`nativecontacts.py`_)

* get the lipid composition of the individual leaflets of a bilayer
  using :mod:`MDAnalysis.analysis.leaflet` (`membrane-leaflets.py`_);

* define the multimeric states of a number of transmembrane peptides
  via clustering (`multimers-analysis.py`_);

* convert between trajectory formats (e.g. `dcd2xtc.py`_ or `amber2dcd.py`_)

* use MDAnalysis for simple model building (`make_MthK_tetramer.py`_);

and more.

.. Links to the stable git repository:

.. _examples:
   https://github.com/MDAnalysis/mdanalysis/blob/master/package/examples/

.. _`rmsfit_qcp.py`:
   https://github.com/MDAnalysis/mdanalysis/blob/master/package/examples/rmsfit_qcp.py
.. _`blocks.py`:
   https://github.com/MDAnalysis/mdanalysis/blob/master/package/examples/blocks.py
.. _`potential_profile.py`:
   https://github.com/MDAnalysis/mdanalysis/blob/master/package/examples/potential_profile.py
.. _`nativecontacts.py`:
   https://github.com/MDAnalysis/mdanalysis/blob/master/package/examples/nativecontacts.py
.. _`membrane-leaflets.py`:
   https://github.com/MDAnalysis/mdanalysis/blob/master/package/examples/membrane-leaflets.py
.. _`multimers-analysis.py`:
   https://github.com/MDAnalysis/mdanalysis/blob/master/package/examples/multimers-analysis.py
.. _`dcd2xtc.py`:
   https://github.com/MDAnalysis/mdanalysis/blob/master/package/examples/dcd2xtc.py
.. _`amber2dcd.py`:
   https://github.com/MDAnalysis/mdanalysis/blob/master/package/examples/amber2dcd.py
.. _`make_MthK_tetramer.py`:
   https://github.com/MDAnalysis/mdanalysis/blob/master/package/examples/make_MthK_tetramer.py
