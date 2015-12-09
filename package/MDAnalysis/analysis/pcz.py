"""
A PCA module for MDAnalysis
==========================================================================

:Authors: The pyPcazip Developers - Ardita Shkurti, Ramon Goni, Pau Andrio, Elena Breitmoser, Iain
Bethune, Modesto Orozco, Charles Laughton
:Year: 2015
:Copyright: BSD


The module contains code to perform Principal Component Analysis of trajectory
data. For examples of the application of PCA to MD analysis see Wlodek et al.
[Wlodek1997]_ and Sherer et al. [Sherer1999]_.

Much of the code has been adapted from that in the `pyPcazip package`_.

.. Rubric:: References

..[Wlodek1997] Wlodek ST, Clark TW, Scott LR, McCammon JA (1997)
               J Am Chem Soc 119, 9513-9522.

..[Sherer1999] Sherer EC, Harris SA, Soliva R, Orozco M, Laughton CA (1999)
               J Am Chem Soc 121, 5981-5991.

.._`pyPcazip package`: https://bitbucket.org/ramonbsc/pypcazip



.. examples

Examples
--------

    import MDAnalysis as mda
    from MDAnalysis.analysis.pca import Pca

    u = mda.Universe(pdb, trajectory)
    p = Pca(u, atom_selection='name CA')
    # print the eigenvalues:
    print p.evals
    # print the projection of the trajectory along the second eigenvector:
    print p.proj(1)  # note zero-based indexing!
    # print the scores for trajectory frame 7:
    print p.scores(7)
    # print the scores for a different coordinate set:
    u2 = mda.Universe(pdb, trajectory2)
    x = u.atomSelection('name CA').ts._pos
    print p.transform(x)
    # generate a coordinate set from a scores vector defining
    # a point in (2D) PC space:
    s = [5.0, 3.0]
    x = p.transform(s, inverse=True)


Helper functions
----------------
The following are convenience classes and functions primarily for internal use.

.. autoclass:: Fasu
.. autofunction:: lsqfit
.. autofunction:: rmsd

Classes, methods, and functions
-------------------------------

.. autoclass:: Pca
   :members:

   .. attribute:: natoms

      int, number of atoms selected for PC analysis

   .. attribute:: numframes

      int, number of trajectory frames selected for PC analysis

   .. attribute:: avg

      :class:`numpy.ndarray` containing the coordinates of the average structure

   .. attribute:: nvecs

      int, the number of eigenvectors/eigenvalues calculated

   .. attribute:: captured_variance

      float, the fraction (between 0.0 and 1.0) of the total variance captured 
      by the *nvecs* eigenvectors

   .. attribute:: evals

      :class:`numpy.ndarray` containing the eigenvalues (largest first)

   .. attribute:: evecs

      :class:`numpy.ndarray` containing the eigenvectors

"""
import MDAnalysis as mda
import numpy as np

from MDAnalysis.analysis.align import rotation_matrix
from scipy.linalg import eigh


def lsqfit(mobile, reference, weights = None):
    """
    Utility function to fit the coordinates in the [n, 3] numpy array
    *mobile* to those in the [n, 3] array *reference*, with optional mass
    weighting.

    :Arguments:
      *mobile*
        :class:`numpy.ndarray` containing the coordinates to be fitted

      *reference*
        :class:`numpy.ndarray` containing the target coordinates

      *weights*
        :class: numpy.ndarray`  (optional) for mass-weighted fitting

    :Returns:
       :class:`numpy.ndarray` containing the transformed coordinates of *mobile* 
    """
    lmob = mobile - mobile.mean(axis=0)
    targ_com = reference.mean(axis=0)
    ltarg = reference - targ_com
    R, rms = rotation_matrix(ltarg, lmob, weights = weights)
    fitted = np.dot(lmob, np.array(R))
    fitted += targ_com
    return fitted

def rmsd(mobile, reference, weights = None):
    """
    Utility function to calculate the rmsd between the coordinates in 
    the [n, 3] numpy array *mobile* and those in the [n, 3] array
    *reference*, with optional mass weighting.

    :Arguments:
      *mobile*
        :class:`numpy.ndarray` containing the coordinates to be fitted

      *reference*
        :class:`numpy.ndarray` containing the target coordinates

      *weights*
        :class: numpy.ndarray`  (optional) for mass-weighted fitting

    :Returns:
       float, the RMS deviation between *mobile* and *reference* after
       least-squares fiting (with optional mass-weighting) 
    """
    lmob = mobile - mobile.mean(axis=0)
    targ_com = reference.mean(axis=0)
    ltarg = reference - targ_com
    R, rms = rotation_matrix(ltarg, lmob, weights = weights)
    return rms

class Fasu:
    def __init__(self, universe, atom_selection='name *', start=None, stop=None, step=None):
        """
        A "Fasu" is a convenience wrapper round an MDAnalysis universe,
        providing slicing and pseudo-direct access methods for any 
        trajectory format.

        :Arguments:
          *universe*
            :class:`MDAnalysis.Universe` object

          *atom_selection* 
            a valid selection string for 
            :meth:`MDAnalysis.core.AtomGroup.AtomGroup.select_atoms`

          *start*
            int, the index (zero-based) of the first frame in the trajectory
            to include in the PC analysis

          *stop*
            int, the index (zero-based) of the last frame in the trajectory 
            to include in the PC analysis

          *step*
            int, the increment between frames to be included in the analysis
        """
        self._u = universe
        self._sel = self._u.select_atoms(atom_selection)
        self.natoms = len(self._sel)
        if start is None:
            self._start = 0
        else:
            self._start = start
        if stop is None:
            self._stop = self._u.trajectory.n_frames
        else: 
            self._stop = stop
        if step is None:
            self._step = 1
        else:
            self._step = step

        self.numframes = abs(self._stop - self._start)/self._step

    def coords(self, frame_number, reference_structure = None):
        """
        Returns a [natoms, 3] numpy array of the selected coordinates at the
        chosen frame (zero-based). The coordinates are least-squares
        fitted to reference_structure, if this is specified.

        :Arguments:
          *frame_number*
            int, the index (zero-based) of the desired trajectory frame

          *reference_structure*
            :class:`numpy.ndarray` containing the coordinates of the 
            reference structure

        :Returns:
          :class:`numpy.ndarray` containing the coordinates of the selected 
          snapshot
        """
        f_real = self._start + frame_number * self._step

        if f_real < (self._u.trajectory.frame - 1):
            self._u.trajectory.rewind()
        while (self._u.trajectory.frame ) < f_real:
            try:
                self._u.trajectory.next()
            except EOFError:
                print 'End of file error trying to access frame ', frame
                raise EOFError

        crds = self._sel.ts._pos
        if reference_structure is not None:
            crds = lsqfit(crds, reference_structure)
        return crds

class Pca:
    def __init__(self, universe, atom_selection=None, start=None, stop=None,
    step=None, reference_coordinates=None, fitting_tolerance = 0.01,
    max_iterations = 10, captured_variance = None):
        """
        Initialises a new Pca object with the data from the trajectory 
        associated with the supplied universe, subject to the (optional) 
        atom and frame selections. The selected atoms are least-squares 
        fitted using an iterative procedure to remove global rotation and 
        translation before the covariance matrix is calculated. 
        If reference_coordinates, in the form of a [N, 3] numpy array (where 
        N is the number of selected atoms) is provided, this will be used to 
        initialise the least-squares fitting process, otherwise the first 
        snapshot in the trajectory will be used.

        :Arguments:
          *universe*
            :class:`MDAnalysis.Universe`

          *atom_selection* 
            a valid selection string for 
            :meth:`MDAnalysis.core.AtomGroup.AtomGroup.select_atoms`

          *start*
            int, the index (zero-based) of the first frame in the trajectory
            to include in the PC analysis

          *stop*
            int, the index (zero-based) of the last frame in the trajectory 
            to include in the PC analysis

          *step*
            int, the increment between frames to be included in the analysis

          *fitting_tolerance*
            float, the iterative least-squares fitting process will terminate
            when the average structure has an RMSD less than this from the
            structure generated in the previous cycle (or if *max_iterations*
            is reached first).

          *max_iterations*
            int, the maximum number of least-squares fitting cycles

          *captured_variance*
            int or float, if < 1 then enough eigenvectors will be kept to
            explain *captured_variance* fraction of the total variance, if
            >=1 (and integral), this is the number of eigenvectors that will
            be kept. Providing this argument permits the precalculation of 
            projections, which can improve performance.
        """

        # Step 1: Create the "Fasu" object from the selected trajectory data.
        self._f = Fasu(universe, atom_selection=atom_selection, start=start, stop=stop, step=step)
        self.natoms = self._f.natoms
        self.numframes = self._f.numframes

        # Step 2: Calculate the mean structure, later used to remove global
        #         translation and rotation from the trajectory frames before
        #         the covariance matrix calculation. This iterative method
        #         is based on the "Procrustes" approach (Dryden and Mardia,
        #         "Statistical Shape Analysis", ISBN 978-0-471-95816-1).

        if reference_coordinates is None:
            reference_coordinates = self._f.coords(0)

        rms = fitting_tolerance + 1.0
        old_avg = reference_coordinates
        new_avg = reference_coordinates
        iteration = 0
        while rms > fitting_tolerance and iteration < max_iterations:
            old_avg = new_avg
            new_avg = 0.0
            for i in range(self.numframes):
                new_avg += self._f.coords(i, reference_structure=old_avg)

            new_avg = new_avg/self.numframes
            rms  = rmsd(new_avg, old_avg)
            iteration += 1

        self.avg = new_avg

        # Step 3: Calculation of the covariance matrix. Each snapshot is
        #         first least-squares fitted to the mean structure, then 
        #         the mean structure is subtracted from it. The covariance
        #         matrix is built up in as big blocks as memory will permit.
        covar = np.zeros((3*self._f.natoms, 3*self._f.natoms))
        blocksize = self.numframes * 2
        toobig = True
        while toobig:
            blocksize = blocksize/2
            try:
                trj = np.zeros((blocksize, self.natoms * 3))
                toobig = False
            except MemoryError:
                blocksize = blocksize/2
                
        num_full_blocks = self.numframes/blocksize
        remainder = self.numframes%blocksize

        frames_so_far = 0
        for i in range(num_full_blocks):
            for j in range(blocksize):
                mob = self._f.coords(j + frames_so_far, reference_structure = self.avg) - self.avg
                trj[j,:] = mob.flatten()
            covar += np.dot(trj.T, trj.conj())
            frames_so_far += blocksize

        trj = np.zeros((remainder, self.natoms * 3))
        for j in range(remainder):
            mob = self._f.coords(j + frames_so_far, reference_structure = self.avg) - self.avg
            trj[j,:] = mob.flatten()
        covar += np.dot(trj.T, trj.conj())
        covar = covar/self.numframes

        # Step 4: Diagonalise the covariance matrix; reorder the eigenvectors
        #         and eigenvalues into descending importance; reduce to the
        #         number required to explain the desired fraction of the total
        #         variance (if specified).
        self.evals, self.evecs = eigh(covar)
        self.evals = self.evals[-1::-1]
        self.evecs = self.evecs[:, -1::-1].T
        self.nvecs = len(self.evals)

        cs = np.cumsum(self.evals)
        cs = cs/cs[-1]

        if captured_variance is not None:
            if captured_variance < 1.0:
                i = 0
                while cs[i] < captured_variance:
                    i += 1
                self.nvecs = i +1 
            else:
                self.nvecs = int(captured_variance)

            self.evals = self.evals[:self.nvecs]
            self.evecs = self.evecs[:self.nvecs]

        self.captured_variance = cs[self.nvecs - 1]

        # Step 5: If only a subset of the 3N eigenvectors have been
        #         saved, precalculate and store the projection data.
        if self.captured_variance < 1.0:
            self._proj = np.zeros((self.nvecs, self.numframes))
            for i in range(self.numframes):
                x = self._f.coords(i, reference_structure = self.avg) - self.avg
                self._proj[:,i] = np.dot(x.flatten(), self.evecs.T)
            
            
    def proj(self, ev):
        """
        Returns the projection of each frame in the trajectory along
        eigenvector *ev* (zero-based), as a numpy array.

        :Arguments:
          *ev*
            int, the index (zero-based) of the desired eigenvector

        :Returns:
          :class:`numpy.ndarray` containing the projection data (one value 
          per frame in the trajectory)
        """
        if ev >= len(self.evals):
            print 'Error - there are only ', len(self.evals), 'eigenvectors.'
            return None
        else:
            # Use precalculated projection data, if available
            if self.captured_variance < 1.0:
                return self._proj[ev]
            else:
                p = np.zeros(self.numframes)
                blocksize = self.numframes * 2
                toobig = True
                while toobig:
                    blocksize = blocksize/2
                    try:
                        trj = np.zeros((blocksize, self.natoms * 3))
                        toobig = False
                    except MemoryError:
                        blocksize = blocksize/2
                
                num_full_blocks = self.numframes/blocksize
                remainder = self.numframes%blocksize

                frames_so_far = 0
                for i in range(num_full_blocks):
                    for j in range(blocksize):
                        mob = self._f.coords(j + frames_so_far, reference_structure = self.avg) - self.avg
                        trj[j,:] = mob.flatten()
                    p[frames_so_far:frames_so_far+blocksize] = np.dot(trj, self.evecs[ev].T).T
                    frames_so_far += blocksize

                trj = np.zeros((remainder, self.natoms * 3))
                for j in range(remainder):
                    mob = self._f.coords(j + frames_so_far, reference_structure = self.avg) - self.avg
                    trj[j,:] = mob.flatten()
                p[frames_so_far:frames_so_far+remainder] = np.dot(trj, self.evecs[ev].T).T
                return p

    def scores(self, frame_number):
        """
        Returns the scores (coordinates in PC space) corresponding to frame
        *frame_number* (zero-based) in the trajectory.

        :Arguments:
          *frame_number*
            int, the (zero-based) index of the frame whose scores are required

        :Returns:
            :class:`numpy.ndarray` containing the vector of scores (one 
            value per eigenvector)
        """
        # Use precalculated projection/score data if available
        if self.captured_variance < 1.0:
            return self._proj[:, frame_number]
        else:
            x = self._f.coords(frame_number, reference_structure = self.avg) - self.avg
            return np.dot(x.flatten(), self.evecs.T).T

    def transform(self, coordinates, inverse = False):
        """
        Transforms a [N, 3] array of Cartesian coordinates into an [3*N]
        array of PC coordinates, or (if inverse = True) the reverse.

        :Arguments:
          *coordinates*
            :class:`numpy.ndarray` of either Cartesian coordinates (if
            *inverse* = False), or of scores (if *inverse* = True). In the
            former case the array must be of size [N, 3] exactly (where N is
            the number of atoms), but in the latter case the scores vector
            may contain less than [3*N] elements, in which case the Cartesian
            coordinate transformation will use just the most significant 
            eigenvectors

          *inverse*
            logical, controls the direction of the transformation
        """
        if not inverse:
            x = lsqfit(coordinates, self.avg) - self.avg
            return np.dot(x.flatten(), self.evecs.T).T
        else:
            x = self.avg
            for i in range(len(coordinates)):
                x += (self.evecs[i]*coordinates[i]).reshape((self.natoms, 3))
            return x
