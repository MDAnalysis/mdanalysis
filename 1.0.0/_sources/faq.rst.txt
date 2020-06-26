.. -*- coding: utf-8 -*-

==========================
Frequently asked questions
==========================


Trajectories
==========================

Why do the atom positions change over trajectories?
---------------------------------------------------

A fundamental concept in MDAnalysis is that at any one time, 
only one time frame of the trajectory is being accessed. The 
:code:`trajectory` attribute of a :code:`Universe` is actually (usually) a file reader. 
Think of the trajectory as a function :math:`X(t)` of the frame index :math:`t` 
that makes the data from this specific frame available. This structure is important 
because it allows MDAnalysis to work with trajectory files too large to fit 
into the computer's memory. See :ref:`trajectories` for more information. 

