=============================
 MDAnalysis Virtual Machines
=============================

This directory contains the set up for virtual machines (VM) using
Vagrant_ with VirtualBox_

.. _Vagrant: https://www.vagrantup.com/
.. _VirtualBox: https://www.virtualbox.org/
.. _ipython: http://ipython.org/

You first need to install VirtualBox_ and Vagrant_ --- either through
your package manager or by downloading the appropriate packages from
their websites.

Assuming that you have been successful, you can bring up a virtual
Linux machine that will automatically install itself with
MDAnalysis. For the following to work, you will need an internet
connection.

Start/stop VM
=============

For example, Ubuntu 14.04::
  
  cd Ubuntu/14.04
  vagrant up

The ``vagrant up`` command will take a while because additional
packages are installed and finally the latest MDAnalysis release is
installed from the internet. 

The next time you run ``vagrant up`` it will be much faster.

You can suspend the VM (saving its instantaneous state)::
  vagrant suspend
or you can shut it down ::
  vagrant halt
(Then the next ``up`` will boot up the VM.)

.. Warning:: 
   If you want to start completely fresh, run ``vagrant destroy``.


Working with the VM
===================

Once you have brought up your VM you then ``ssh`` into it::

   vagrant ssh

(This VM only supports commandline access, no graphical user
interface.)

When you are done, just exit the ssh session (``exit`` or Control-D).

You can access your real home directory in the virtual machine at the
path ``/myhome``. Anything that you do in this directory will be
reflected in your real home directory, including deletion of files!

.. Note:: There can be peformance issues when accessing shared
   directories in the VM so the virtual machine is not recommended for
   production analysis.)

Within the VM you can then use MDAnlysis. For instance, start ipython_::
  ipython

Within ``ipython`` (you can copy and paste line-by-line for
interactive exploration or copy and the use ``%paste`` in
``ipython``)::
  
  import MDAnalysis
  from MDAnalysis.tests.datafiles import PSF, DCD
  import numpy.linalg
  import matplotlib.pyplot as plt
  
  # load data
  u = MDAnalysis.Universe(PSF, DCD)
  
  # selections
  nterm = u.s4AKE.N[0]   # can access structure via segid (s4AKE) and atom name
  cterm = u.s4AKE.C[-1]  # ... takes the last atom named 'C'
  bb = u.selectAtoms('protein and backbone')  # a selection (a AtomGroup)

  # iterate through the trajectory and collect data
  data = []
  for ts in u.trajectory:     # iterate through all frames
     r = cterm.pos - nterm.pos # end-to-end vector from atom positions
     d = numpy.linalg.norm(r)  # end-to-end distance
     rgyr = bb.radius_of_gyration()  # method of a AtomGroup; updates with each frame
     data.append((ts.frame, u.trajectory.time, d, rgyr))
  data = numpy.array(data).transpose()
  frames, t, e2e, Rg = data  

  # plotting
  plt.close('all')
  fig = plt.figure(figsize=(5,5))
  ax = fig.add_subplot(111)
  ax.plot(t, e2e, 'k-', lw=2, label="end to end distance")
  ax.plot(t, Rg, 'r--', lw=2, label="radius of gyration")
  ax.set_xlabel(r"time $t$ (ps)")
  ax.set_ylabel(r"length ($\AA$)")
  ax.legend(loc="best")

  # you can view the figure in your *real* home directory:
  fig.savefig("/myhome/adk_e2e_rgyr.pdf")
  print("Figure was saved to your real home directory: ~/adk_e2e_rgyr.pdf")

The script will produce a figure in PDF format that will appear in
your real home directory so that you can view it with a PDF viewer.

It is generally recommended that you create a separate work directory
in your home directory to which you redirect output from working with
the tutorial VM.





