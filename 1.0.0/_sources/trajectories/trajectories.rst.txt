.. -*- coding: utf-8 -*-
.. _trajectories:

============
Trajectories
============

In MDAnalysis, static data is contained in your universe Topology, while dynamic data is drawn from its trajectory at ``Universe.trajectory``. This is typically loaded from a trajectory file and includes information such as:

    * atom coordinates (``Universe.atoms.positions``)
    * box size (``Universe.dimensions``)
    * velocities and forces (if your file format contains the data) (``Universe.atoms.velocities``)

Although these properties look static, they are actually dynamic, and the data contained within can change.
In order to remain memory-efficient, MDAnalysis does not load every frame of your trajectory into memory at once. Instead, a Universe has a state: the particular timestep that it is currently associated with in the trajectory. When the timestep changes, the data in the properties above shifts accordingly.

The typical way to change a timestep is to index it. ``Universe.trajectory`` can be thought of as a list of :class:`~MDAnalysis.coordinates.base.Timestep`\ s, a data structure that holds information for the current time frame. For example, you can query its length.

.. ipython:: python

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import PSF, DCD

    u = mda.Universe(PSF, DCD)
    len(u.trajectory)

When a trajectory is first loaded from a file, it is set to the first frame (with index 0), by default.

.. ipython:: python

    print(u.trajectory.ts, u.trajectory.time)

Indexing the trajectory returns the timestep for that frame, and sets the Universe to point to that frame until the timestep next changes.

.. ipython:: python

    u.trajectory[3]

.. ipython:: python

    print('Time of fourth frame', u.trajectory.time)

Many tasks involve applying a function to each frame of a trajectory. For these, you need to iterate through the frames, *even if you don't directly use the timestep*. This is because the act of iterating moves the Universe onto the next frame, changing the dynamic atom coordinates. 

Trajectories can also be :ref:`sliced <slicing-trajectories>` if you only want to work on a subset of frames.

.. ipython:: python

    protein = u.select_atoms('protein')
    for ts in u.trajectory[:20:4]:
        # the radius of gyration depends on the changing positions
        rad = protein.radius_of_gyration()
        print('frame={}: radgyr={}'.format(ts.frame, rad))
    
Note that after iterating over the trajectory, the frame is always set back to the first frame, even if your loop stopped before the trajectory end.

.. ipython:: python

    u.trajectory.frame

Because MDAnalysis will pull trajectory data directly from the file it is reading from, changes to atom coordinates and box dimensions will not persist once the frame is changed. The only way to make these changes permanent is to load the trajectory into memory, or to write a new trajectory to file for every frame. For example, to set a cubic box size for every frame and write it out to a file::

    with mda.Writer('with_box.trr', 'w', n_atoms=u.atoms.n_atoms) as w:
        for ts in u.trajectory:
            ts.dimensions = [10, 10, 10, 90, 90, 90]
            w.write(u.atoms)
    
    u_with_box = mda.Universe(PSF, 'with_box.trr')


Sometimes you may wish to only transform part of the trajectory, or to not write a file out. In these cases, MDAnalysis supports :ref:`"on-the-fly" transformations <transformations>` that are performed on a frame when it is read. 


