.. -*- coding: utf-8 -*-
.. _slicing-trajectories:

====================
Slicing trajectories
====================

MDAnalysis trajectories can be indexed to return a :class:`~MDAnalysis.coordinates.base.Timestep`, or sliced to give a :class:`~MDAnalysis.coordinates.base.FrameIterator`. 

.. ipython:: python

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import PSF, DCD

    u = mda.Universe(PSF, DCD)
    u.trajectory[4]


Indexing a trajectory shifts the :class:`~MDAnalysis.core.universe.Universe` to point towards that particular frame, updating dynamic data such as ``Universe.atoms.positions``. 

.. note::

    The trajectory frame is not read from the MD data. It is the internal index assigned by MDAnalysis.

.. ipython:: python

    u.trajectory.frame

*Creating* a :class:`~MDAnalysis.coordinates.base.FrameIterator` by slicing a trajectory does not shift the :class:`~MDAnalysis.core.universe.Universe` to a new frame, but *iterating* over the sliced trajectory will rewind the trajectory back to the first frame.

.. ipython::: python

    fiter = u.trajectory[10::2]
    u.trajectory.frame

.. ipython:: python

    fiter = u.trajectory[10::10]
    frames = [ts.frame for ts in fiter]
    print(frames, u.trajectory.frame)

You can also create a sliced trajectory with boolean indexing and fancy indexing. Boolean indexing allows you to select only frames that meet a certain condition, by passing a :class:`~numpy.ndarray` with the same length as the original trajectory. Only frames that have a boolean value of ``True`` will be in the resulting :class:`~MDAnalysis.coordinates.base.FrameIterator`. For example, to select only the frames of the trajectory with an RMSD under 2 angstrom:

.. ipython:: python

    from MDAnalysis.analysis import rms

    protein = u.select_atoms('protein')
    rmsd = rms.RMSD(protein, protein).run()
    bools = rmsd.rmsd.T[-1] < 2
    print(bools)

.. ipython:: python

    fiter = u.trajectory[bools]
    print([ts.frame for ts in fiter])

You can also use fancy indexing to control the order of specific frames.

.. ipython:: python

    indices = [10, 2, 3, 9, 4, 55, 2]
    print([ts.frame for ts in u.trajectory[indices]])

You can even slice a :class:`~MDAnalysis.coordinates.base.FrameIterator` to create a new :class:`~MDAnalysis.coordinates.base.FrameIterator`.

.. ipython:: python

    print([ts.frame for ts in fiter[::3]])
