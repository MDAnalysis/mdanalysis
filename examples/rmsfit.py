#!/usr/bin/env python
# Example script, part of MDAnalysis
"""
RMS-fit a trajectory to a reference structure
=============================================

Simple example implementation that shows how to fit by rms and access
the translated and rotated coordinates. For a more complicated example
look at ``rmsfit_alignment.py``.

"""

import numpy
import MDAnalysis

import sys

def rmsd(a,b):
    """Returns RMSD between two coordinate sets a and b."""
    return numpy.sqrt(numpy.sum(numpy.power(a-b,2))/a.shape[0])

def echo(s=''):
    """Simple string output that immediately prints to the console."""
    sys.stderr.write(s)
    sys.stderr.flush()

def rmsd_trj(traj,ref,select='backbone',):
    """Calculate the RMSD of a selection over a trajectory.

      rmsd_trj(traj, ref, 'backbone or name CB or name OT*')

    :Arguments:
      traj
         trajectory, :class:`MDAnalysis.Universe` object
      ref
         reference coordinates; :class:`MDAnalysis.Universe` object
         (uses the current time step of the object)
      select    
         any valid selection string for
         MDAnalysis.AtomGroup.selectAtom that produces identical
         selections in traj and ref or dictionary {'reference':sel1,
         'target':sel2}.  The fasta2select() function returns such a
         dictionary based on a ClustalW or STAMP sequence alignment.

    Both reference and trajectory must be MDAnalysis.Universe
    instances.
    """
    from MDAnalysis import SelectionError
    import MDAnalysis.core.rms_fitting

    # setup new dcd
    frames = traj.trajectory

    ref_atoms = ref.selectAtoms(select)
    traj_atoms = traj.selectAtoms(select)
    if len(ref_atoms) != len(traj_atoms):
        raise SelectionError("Reference and trajectory atom selections do not contain "+
                             "the same number of atoms: N_ref=%d, N_traj=%d" % \
                             (len(ref_atoms), len(traj_atoms)))

    # could check that ref and traj selection have same masses (see rmsfit_alignment.py)
    masses = ref_atoms.masses()

    # reference centre of mass system
    # (compatibility with pre 1.0 numpy: explicitly cast coords to float32)
    ref_com = ref_atoms.centerOfMass().astype(numpy.float32)
    ref_coordinates = ref_atoms.coordinates() - ref_com

    # allocate the array for selection atom coords
    # (allocating and re-filling is MUCH more efficient that re-allocating for 
    # every frame!)
    traj_coordinates = traj_atoms.coordinates().copy()

    # collect the time series
    rmsds = []

    # R: rotation matrix that aligns r-r_com, x~-x~com   
    #    (x~: selected coordinates, x: all coordinates)
    # Final transformed traj coordinates: x' = (x-x~_com)*R + ref_com
    for ts in frames:
        # shift coordinates for rotation fitting
        # selection is updated with the time frame
        x_com = traj_atoms.centerOfMass().astype(numpy.float32)
        # re-fill the coordinates array
        traj_coordinates[:] = traj_atoms.coordinates() - x_com
        R = numpy.matrix(MDAnalysis.core.rms_fitting.rms_rotation_matrix(
                traj_coordinates,ref_coordinates,masses),dtype=numpy.float32)
        # Transform each atom in the trajectory (use inplace ops to avoid copying arrays)
        # Note that this destroys the original coordinate information in the current timestep.
        ts._pos   -= x_com
        ts._pos[:] = ts._pos * R # R acts to the left (!) & is broadcasted N times (numpy magic!)
        ts._pos   += ref_com

        rmsd_old = rmsd(ref_atoms.coordinates(),traj_coordinates)
        rmsd_new = rmsd(ref_atoms.coordinates(),traj_atoms.coordinates())

        rmsds.append(rmsd_new)

        echo("Fitted frame %5d/%d  [%5.1f%%]  %5.2fA --> %5.2fA  |translation|=%.2fA\r" % \
            (ts.frame,frames.numframes,100.0*ts.frame/frames.numframes,
             rmsd_old, rmsd_new, rmsd(x_com,ref_com)))
        #if ts.frame % 10 == 0:
        #    echo("Fitted frame %5d/%d  [%5.1f%%]\r" % \
        #        (ts.frame,frames.numframes,100.0*ts.frame/frames.numframes))
    echo('\n')
    return numpy.array(rmsds)

if __name__ == '__main__':
   from MDAnalysis import *
   from MDAnalysis.tests.datafiles import PSF, DCD, PDB_small
   ref = Universe(PSF, PDB_small)   # reference structure 4AKE
   trj = Universe(PSF, DCD)         # trajectory of change 1AKE->4AKE
   
   print "CA RMSD for %(DCD)r versus %(PDB_small)r"
   rmsds1 = rmsd_trj(trj, ref, select='name CA')

   print "CA RMSD for %(DCD)r versus first frame"
   ref = Universe(PSF, DCD)
   ref.trajectory[0]  # go to first frame
   rmsds2 = rmsd_trj(trj, ref, select='name CA')

   print "CA RMSD for %(DCD)r versus last frame"
   ref = Universe(PSF, DCD)
   ref.trajectory[-1]  # go to last frame
   rmsds3 = rmsd_trj(trj, ref, select='name CA')

   try:
       from pylab import plot, xlabel, ylabel, legend, title
       plot(rmsds1, linewidth=2, color='black', label='4AKE')
       plot(rmsds2, linewidth=2, color='red', label='first frame')       
       plot(rmsds3, linewidth=2, color='blue', label='last frame')       
   
       xlabel('trajectory frame')
       ylabel(r'root mean square distance ($\mathrm{\AA}$)')
       title('DIMS trajectory of Adenylate Kinase from 1AKE to 4AKE')
       legend(loc='best')
   except ImportError:
       pass
       

