=======================
Examples for MDAnalysis
=======================

These python scripts show how MDAnalysis has been used to analyze MD
trajectories produced by LAMMPS or CHARMM. They are provided as
real-world examples but no effort has been made to make them
userfriendly; in fact, they will probably not work right away for you
as many assumptions have been simply hard coded.

Many of those scripts originated with older versions of MDAnalysis and
might be broken. Feedback is appreciated as are your own example
scripts.

Please discuss these scripts on the mailing list
http://groups.google.com/group/mdnalysis-discussion or file bugs
through the issue tracker
http://issues.mdanalysis.org


.. Note:: If you try to run example scripts from the top level
   directory (such as ``python ./examples/nativecontacts.py``) then
   you will likely get an error such as *ImportError: No module named
   _dcdmodule*. The problem is that python is trying to pick up
   MDAnalysis from the *source* directory (``./MDAnalysis``) but the
   compiled libraries were installed elsewhere. The solution is to not
   be in the source directory: simply ``cd examples/`` or copy the
   examples elsewhere into your work space.


Descriptions
============

For the time being, please see the scripts for details. They are all
fairly short.
