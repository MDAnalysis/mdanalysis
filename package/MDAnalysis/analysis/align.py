# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""Coordinate fitting and alignment --- :mod:`MDAnalysis.analysis.align`
=====================================================================

:Author: Oliver Beckstein, Joshua Adelman
:Year: 2010--2013
:Copyright: GNU Public License v3

The module contains functions to fit a target structure to a reference
structure. They use the fast QCP algorithm to calculate the root mean
square distance (RMSD) between two coordinate sets [Theobald2005]_ and
the rotation matrix *R* that minimizes the RMSD [Liu2010]_. (Please
cite these references when using this module.).

Typically, one selects a group of atoms (such as the C-alphas),
calculates the RMSD and transformation matrix, and applys the
transformation to the current frame of a trajectory to obtain the
rotated structure. The :func:`alignto` and :func:`rms_fit_trj`
functions can be used to do this for individual frames and
trajectories respectively.

The :ref:`RMS-fitting-tutorial` shows how to do the individual steps
manually and explains the intermediate steps.

.. SeeAlso::

   :mod:`MDAnalysis.analysis.rms`
        contains functions to compute RMSD (when structural alignment is not
        required)
   :mod:`MDAnalysis.lib.qcprot`
        implements the fast RMSD algorithm.


.. _RMS-fitting-tutorial:

RMS-fitting tutorial
--------------------

The example uses files provided as part of the MDAnalysis test suite
(in the variables :data:`~MDAnalysis.tests.datafiles.PSF`,
:data:`~MDAnalysis.tests.datafiles.DCD`, and
:data:`~MDAnalysis.tests.datafiles.PDB_small`). For all further
examples execute first ::

   >>> from MDAnalysis import *
   >>> from MDAnalysis.analysis.align import *
   >>> from MDAnalysis.analysis.rms import rmsd
   >>> from MDAnalysis.tests.datafiles import PSF, DCD, PDB_small


In the simplest case, we can simply calculate the C-alpha RMSD between
two structures, using :func:`rmsd`::

   >>> ref = Universe(PDB_small)
   >>> mobile = Universe(PSF,DCD)
   >>> rmsd(mobile.atoms.CA.positions, ref.atoms.CA.positions)
   16.282308620224068

Note that in this example translations have not been removed. In order
to look at the pure rotation one needs to superimpose the centres of
mass (or geometry) first:

   >>> rmsd(mobile.atoms.CA.positions, ref.atoms.CA.positions, center=True)
   12.639693690256898

This has only done a translational superposition. If you want to also do a
rotational superposition use the superposition keyword. This will calculate a
minimized RMSD between the reference and mobile structure.

   >>> rmsd(mobile.atoms.CA.positions, ref.atoms.CA.positions,
   >>>      superposition=True)
   6.8093965864717951

The rotation matrix that superimposes *mobile* on *ref* while
minimizing the CA-RMSD is obtained with the :func:`rotation_matrix`
function ::

   >>> mobile0 = mobile.atoms.CA.positions - mobile.atoms.center_of_mass()
   >>> ref0 = ref.atoms.CA.positions - ref.atoms.center_of_mass()
   >>> R, rmsd = rotation_matrix(mobile0, ref0)
   >>> print rmsd
   6.8093965864717951
   >>> print R
   [[ 0.14514539 -0.27259113  0.95111876]
    [ 0.88652593  0.46267112 -0.00268642]
    [-0.43932289  0.84358136  0.30881368]]

Putting all this together one can superimpose all of *mobile* onto *ref*::

   >>> mobile.atoms.translate(-mobile.atoms.CA.center_of_mass())
   >>> mobile.atoms.rotate(R)
   >>> mobile.atoms.translate(ref.atoms.CA.center_of_mass())
   >>> mobile.atoms.write("mobile_on_ref.pdb")


Common usage
------------

To **fit a single structure** with :func:`alignto`::

   >>> ref = Universe(PSF, PDB_small)
   >>> mobile = Universe(PSF, DCD)     # we use the first frame
   >>> alignto(mobile, ref, select="protein and name CA", mass_weighted=True)

This will change *all* coordinates in *mobile* so that the protein
C-alpha atoms are optimally superimposed (translation and rotation).

To **fit a whole trajectory** to a reference structure with the
:func:`rms_fit_trj` function::

   >>> ref = Universe(PSF, PDB_small)   # reference structure 1AKE
   >>> trj = Universe(PSF, DCD)         # trajectory of change 1AKE->4AKE
   >>> rms_fit_trj(trj, ref, filename='rmsfit.dcd')

It is also possible to align two arbitrary structures by providing a
mapping between atoms based on a sequence alignment. This allows
fitting of structural homologs or wild type and mutant.

If a alignment was provided as "sequences.aln" one would first produce
the appropriate MDAnalysis selections with the :func:`fasta2select`
function and then feed the resulting dictionary to :func:`rms_fit_trj`::

   >>> seldict = fasta2select('sequences.aln')
   >>> rms_fit_trj(trj, ref, filename='rmsfit.dcd', select=seldict)

(See the documentation of the functions for this advanced usage.)


Functions
---------

.. autofunction:: alignto
.. autofunction:: rms_fit_trj
.. autofunction:: rotation_matrix

.. versionchanged:: 0.10.0
   Function :func:`~MDAnalysis.analysis.rms.rmsd` was removed from
   this module and is now exclusively accessible as
   :func:`~MDAnalysis.analysis.rms.rmsd`.

Helper functions
----------------

The following functions are used by the other functions in this
module. They are probably of more interest to developers than to
normal users.

.. autofunction:: fasta2select
.. autofunction:: get_matching_atoms

"""

import os.path
from six.moves import range, zip, zip_longest
import numpy as np
import warnings
import logging

import MDAnalysis.lib.qcprot as qcp
from MDAnalysis.exceptions import SelectionError, SelectionWarning
from MDAnalysis.lib.log import ProgressMeter
import MDAnalysis.analysis.rms as rms


logger = logging.getLogger('MDAnalysis.analysis.align')


def rotation_matrix(a, b, weights=None):
    r"""Returns the 3x3 rotation matrix for RMSD fitting coordinate sets *a* and
    *b*.

    The rotation matrix *R* transforms *a* to overlap with *b* (i.e. *b* is the
    reference structure):

    .. math::
        \vec{b} = \bold{R} \dot \vec{a}

    Parameters
    ----------
    a : array-like
          coordinates that are to be rotated ("mobile set"); array of N atoms
          of shape N*3 as generated by, e.g.,
          :meth:`MDAnalysis.core.AtomGroup.AtomGroup.coordinates`.
    b : array-like
          reference coordinates; array of N atoms of shape N*3 as generated by,
          e.g., :meth:`MDAnalysis.core.AtomGroup.AtomGroup.coordinates`.
    weights : array-like (optional)
          array of floats of size N for doing weighted RMSD fitting (e.g. the
          masses of the atoms)

    Returns
    -------
    R : ndarray
        rotation matrix
    rmsd : float
        RMSD between *a* and *b* before rotation
    ``(R, rmsd)`` rmsd and rotation matrix *R*

    Example
    -------
    *R* can be used as an argument for
    :meth:`MDAnalysis.core.AtomGroup.AtomGroup.rotate` to generate a rotated
    selection, e.g. ::

    >>> R = rotation_matrix(A.select_atoms('backbone').positions,
    >>>                     B.select_atoms('backbone').positions)[0]
    >>> A.atoms.rotate(R)
    >>> A.atoms.write("rotated.pdb")

    Notes
    -----
    The function does *not* shift the centers of mass or geometry;
    this needs to be done by the user.

    See Also
    --------
    :func:`rmsd` calculates the RMSD between *a* and *b*; for fitting a whole
    trajectory it is more efficient to use :func:`rms_fit_trj`. A complete fit
    of two structures can be done with :func:`alignto`. """
    if weights is not None:
        # weights are constructed as relative to the mean
        weights = np.asarray(weights, dtype=np.float64) / np.mean(weights)

    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    N = b.shape[0]

    rot = np.zeros(9, dtype=np.float64)
    rmsd = qcp.CalcRMSDRotationalMatrix(a.T, b.T, N, rot, weights)
    return np.matrix(rot.reshape(3, 3)), rmsd


def alignto(mobile, reference, select="all", mass_weighted=False,
            subselection=None, tol_mass=0.1, strict=False):
    """Spatially align *mobile* to *reference* by doing a RMSD fit on *select* atoms.

    The superposition is done in the following way:

    1. A rotation matrix is computed that minimizes the RMSD between
       the coordinates of `mobile.select_atoms(sel1)` and
       `reference.select_atoms(sel2)`; before the rotation, *mobile* is
       translated so that its center of geometry (or center of mass)
       coincides with the one of *reference*. (See below for explanation of
       how *sel1* and *sel2* are derived from *select*.)

    2. All atoms in :class:`~MDAnalysis.core.AtomGroup.Universe` that
       contains *mobile* are shifted and rotated. (See below for how
       to change this behavior through the *subselection* keyword.)

    The *mobile* and *reference* atom groups can be constructed so that they
    already match atom by atom. In this case, *select* should be set to "all"
    (or ``None``) so that no further selections are applied to *mobile* and
    *reference*, therefore preserving the exact atom ordering (see
    :ref:`ordered-selections-label`).

    .. Warning:: The atom order for *mobile* and *reference* is *only*
       preserved when *select* is either "all" or ``None``. In any other case,
       a new selection will be made that will sort the resulting AtomGroup by
       index and therefore destroy the correspondence between the two groups. **It
       is safest not to mix ordered AtomGroups with selection strings.**

    :Arguments:
      *mobile*
         structure to be aligned, a :class:`~MDAnalysis.core.AtomGroup.AtomGroup`
         or a whole :class:`~MDAnalysis.core.AtomGroup.Universe`
      *reference*
         reference structure, a :class:`~MDAnalysis.core.AtomGroup.AtomGroup`
         or a whole :class:`~MDAnalysis.core.AtomGroup.Universe`
      *select*
         1. any valid selection string for
            :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.select_atoms` that produces identical
            selections in *mobile* and *reference*; or
         2. dictionary ``{'mobile':sel1, 'reference':sel2}``.
            (the :func:`fasta2select` function returns such a
            dictionary based on a ClustalW_ or STAMP_ sequence alignment); or
         3.  tuple ``(sel1, sel2)``

         When using 2. or 3. with *sel1* and *sel2* then these selections can also each be
         a list of selection strings (to generate a AtomGroup with defined atom order as
         described under :ref:`ordered-selections-label`).
      *mass_weighted* : boolean
         ``True`` uses the masses :meth:`reference.masses` as weights for the
         RMSD fit.
      *tol_mass*
         Reject match if the atomic masses for matched atoms differ by more than
         *tol_mass* [0.1]
      *strict*
         ``True``
             Will raise :exc:`SelectioError` if a single atom does not
             match between the two selections.
         ``False`` [default]
             Will try to prepare a matching selection by dropping
             residues with non-matching atoms. See :func:`get_matching_atoms`
             for details.
      *subselection*
         Apply the transformation only to this selection.

         ``None`` [default]
             Apply to `mobile.universe.atoms` (i.e. all atoms in the
             context of the selection from *mobile* such as the rest of a
             protein, ligands and the surrounding water)
         *selection-string*
             Apply to `mobile.select_atoms(selection-string)`
         :class:`~MDAnalysis.core.AtomGroup.AtomGroup`
             Apply to the arbitrary group of atoms

    :Returns: RMSD before and after spatial alignment.

    .. SeeAlso:: For RMSD-fitting trajectories it is more efficient to
                 use :func:`rms_fit_trj`.

    .. versionchanged:: 0.8
       Added check that the two groups describe the same atoms including
       the new *tol_mass* keyword.

    .. versionchanged:: 0.10.0
       Uses :func:`get_matching_atoms` to work with incomplete selections
       and new *strict* keyword. The new default is to be lenient whereas
       the old behavior was the equivalent of *strict* = ``True``.
    """
    if select in ('all', None):
        # keep the EXACT order in the input AtomGroups; select_atoms('all')
        # orders them by index, which can lead to wrong results if the user
        # has crafted mobile and reference to match atom by atom
        mobile_atoms = mobile.atoms
        ref_atoms = reference.atoms
    else:
        select = rms._process_selection(select)
        mobile_atoms = mobile.select_atoms(*select['mobile'])
        ref_atoms = reference.select_atoms(*select['reference'])

    ref_atoms, mobile_atoms = get_matching_atoms(ref_atoms, mobile_atoms,
                                                 tol_mass=tol_mass, strict=strict)

    if mass_weighted:
        weights = ref_atoms.masses / np.mean(ref_atoms.masses)
        ref_com = ref_atoms.center_of_mass()
        mobile_com = mobile_atoms.center_of_mass()
    else:
        weights = None
        ref_com = ref_atoms.center_of_geometry()
        mobile_com = mobile_atoms.center_of_geometry()

    ref_coordinates = ref_atoms.positions - ref_com
    mobile_coordinates = mobile_atoms.positions - mobile_com

    old_rmsd = rms.rmsd(mobile_coordinates, ref_coordinates)

    R, new_rmsd = rotation_matrix(mobile_coordinates, ref_coordinates, weights=weights)

    if subselection is None:
        atoms = mobile.universe.atoms
    elif type(subselection) is str:
        atoms = mobile.select_atoms(subselection)
    else:
        try:
            atoms = subselection.atoms
        except AttributeError:
            raise TypeError("subselection must be a selection string, a AtomGroup or Universe or None")

    atoms.translate(-mobile_com)
    atoms.rotate(R)
    atoms.translate(ref_com)

    return old_rmsd, new_rmsd


def rms_fit_trj(traj, reference, select='all', filename=None, rmsdfile=None, prefix='rmsfit_',
                mass_weighted=False, tol_mass=0.1, strict=False, force=True, quiet=False, **kwargs):
    """RMS-fit trajectory to a reference structure using a selection.

    Both reference *ref* and trajectory *traj* must be
    :class:`MDAnalysis.Universe` instances. If they contain a
    trajectory then it is used. The output file format is determined
    by the file extension of *filename*. One can also use the same
    universe if one wants to fit to the current frame.

    :Arguments:
      *traj*
         trajectory, :class:`MDAnalysis.Universe` object
      *reference*
         reference coordinates; :class:`MDAnalysis.Universe` object
         (uses the current time step of the object)
      *select*
         1. any valid selection string for
            :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.select_atoms` that produces identical
            selections in *mobile* and *reference*; or
         2. a dictionary ``{'mobile':sel1, 'reference':sel2}`` (the
            :func:`fasta2select` function returns such a
            dictionary based on a ClustalW_ or STAMP_ sequence alignment); or
         3. a tuple ``(sel1, sel2)``

         When using 2. or 3. with *sel1* and *sel2* then these selections can also each be
         a list of selection strings (to generate a AtomGroup with defined atom order as
         described under :ref:`ordered-selections-label`).
      *filename*
         file name for the RMS-fitted trajectory or pdb; defaults to the
         original trajectory filename (from *traj*) with *prefix* prepended
      *rmsdfile*
         file name for writing the RMSD timeseries [``None``]
      *prefix*
         prefix for autogenerating the new output filename
      *mass_weighted*
         do a mass-weighted RMSD fit
      *tol_mass*
         Reject match if the atomic masses for matched atoms differ by more than
         *tol_mass* [0.1]
      *strict*
         Default: ``False``
         - ``True``: Will raise :exc:`SelectioError` if a single atom does not
           match between the two selections.
         - ``False``: Will try to prepare a matching selection by dropping
           residues with non-matching atoms. See :func:`get_matching_atoms`
           for details.
      *force*
         - ``True``: Overwrite an existing output trajectory (default)
         - ``False``: simply return if the file already exists
      *quiet*
         - ``True``: suppress progress and logging for levels INFO and below.
         - ``False``: show all status messages and do not change the the logging
           level (default)

         .. Note:: If


      *kwargs*
         All other keyword arguments are passed on the trajectory
         :class:`~MDAnalysis.coordinates.base.Writer`; this allows manipulating/fixing
         trajectories on the fly (e.g. change the output format by changing the extension of *filename*
         and setting different parameters as described for the corresponding writer).

    :Returns: *filename* (either provided or auto-generated)

    .. _ClustalW: http://www.clustal.org/
    .. _STAMP: http://www.compbio.dundee.ac.uk/manuals/stamp.4.2/

    .. versionchanged:: 0.8
       Added *kwargs* to be passed to the trajectory :class:`~MDAnalysis.coordinates.base.Writer` and
       *filename* is returned.

    .. versionchanged:: 0.10.0
       Uses :func:`get_matching_atoms` to work with incomplete selections
       and new *strict* keyword. The new default is to be lenient whereas
       the old behavior was the equivalent of *strict* = ``True``.

    """
    frames = traj.trajectory
    if quiet:
        # should be part of a try ... finally to guarantee restoring the log level
        logging.disable(logging.WARN)

    kwargs.setdefault('remarks', 'RMS fitted trajectory to reference')
    if filename is None:
        path, fn = os.path.split(frames.filename)
        filename = os.path.join(path, prefix + fn)
        _Writer = frames.Writer
    else:
        _Writer = frames.OtherWriter
    if os.path.exists(filename) and not force:
        logger.warn("{0} already exists and will NOT be overwritten; use force=True if you want this".format(filename))
        return filename
    writer = _Writer(filename, **kwargs)
    del _Writer

    select = rms._process_selection(select)
    ref_atoms = reference.select_atoms(*select['reference'])
    traj_atoms = traj.select_atoms(*select['mobile'])
    natoms = traj_atoms.n_atoms

    ref_atoms, traj_atoms = get_matching_atoms(ref_atoms, traj_atoms,
                                                 tol_mass=tol_mass, strict=strict)

    logger.info("RMS-fitting on {0:d} atoms.".format(len(ref_atoms)))
    if mass_weighted:
        # if performing a mass-weighted alignment/rmsd calculation
        weight = ref_atoms.masses / ref_atoms.masses.mean()
    else:
        weight = None

    # reference centre of mass system
    ref_com = ref_atoms.center_of_mass()
    ref_coordinates = ref_atoms.positions - ref_com

    # allocate the array for selection atom coords
    traj_coordinates = traj_atoms.positions.copy()

    # RMSD timeseries
    nframes = len(frames)
    rmsd = np.zeros((nframes,))

    # R: rotation matrix that aligns r-r_com, x~-x~com
    #    (x~: selected coordinates, x: all coordinates)
    # Final transformed traj coordinates: x' = (x-x~_com)*R + ref_com
    rot = np.zeros(9, dtype=np.float64)  # allocate space for calculation
    R = np.matrix(rot.reshape(3, 3))

    percentage = ProgressMeter(nframes, interval=10, quiet=quiet,
                               format="Fitted frame %(step)5d/%(numsteps)d  [%(percentage)5.1f%%]\r")

    for k, ts in enumerate(frames):
        # shift coordinates for rotation fitting
        # selection is updated with the time frame
        x_com = traj_atoms.center_of_mass().astype(np.float32)
        traj_coordinates[:] = traj_atoms.positions - x_com

        # Need to transpose coordinates such that the coordinate array is
        # 3xN instead of Nx3. Also qcp requires that the dtype be float64
        # (I think we swapped the position of ref and traj in CalcRMSDRotationalMatrix
        # so that R acts **to the left** and can be broadcasted; we're saving
        # one transpose. [orbeckst])
        rmsd[k] = qcp.CalcRMSDRotationalMatrix(ref_coordinates.T.astype(np.float64),
                                               traj_coordinates.T.astype(np.float64),
                                               natoms, rot, weight)
        R[:, :] = rot.reshape(3, 3)

        # Transform each atom in the trajectory (use inplace ops to avoid copying arrays)
        # (Marginally (~3%) faster than "ts.positions[:] = (ts.positions - x_com) * R + ref_com".)
        ts.positions -= x_com
        ts.positions[:] = ts.positions * R  # R acts to the left & is broadcasted N times.
        ts.positions += ref_com

        writer.write(traj.atoms)  # write whole input trajectory system
        percentage.echo(ts.frame)
    logger.info("Wrote %d RMS-fitted coordinate frames to file %r",
                frames.n_frames, filename)
    if rmsdfile is not None:
        np.savetxt(rmsdfile, rmsd)
        logger.info("Wrote RMSD timeseries  to file %r", rmsdfile)

    if quiet:
        # should be part of a try ... finally to guarantee restoring the log level
        logging.disable(logging.NOTSET)

    return filename


def sequence_alignment(mobile, reference, **kwargs):
    """Generate a global sequence alignment between residues in *reference* and *mobile*.

    The global alignment uses the Needleman-Wunsch algorith as
    implemented in :mod:`Bio.pairwise2`. The parameters of the dynamic
    programming algorithm can be tuned with the keywords. The defaults
    should be suitable for two similar sequences. For sequences with
    low sequence identity, more specialized tools such as clustalw,
    muscle, tcoffee, or similar should be used.

    :Arguments:
       *mobile*
          protein atom group
       *reference*
          protein atom group

    :Keywords:
      *match_score*
         score for matching residues [2]
      *mismatch_penalty*
         penalty for residues that do not match [-1]
      *gap_penalty*
         penalty for opening a gap; the high default value creates compact
         alignments for highly identical sequences but might not be suitable
         for sequences with low identity [-2]
      *gapextension_penalty*
         penalty for extending a gap [-0.1]

    .. versionadded:: 0.10.0
    """
    import Bio.pairwise2
    kwargs.setdefault('match_score', 2)
    kwargs.setdefault('mismatch_penalty', -1)
    kwargs.setdefault('gap_penalty', -2)
    kwargs.setdefault('gapextension_penalty', -0.1)

    aln = Bio.pairwise2.align.globalms(
        reference.sequence(format="string"), mobile.sequence(format="string"),
        kwargs['match_score'], kwargs['mismatch_penalty'],
        kwargs['gap_penalty'], kwargs['gapextension_penalty'])
    # choose top alignment
    return aln[0]


def fasta2select(fastafilename, is_aligned=False,
                 ref_resids=None, target_resids=None,
                 ref_offset=0, target_offset=0, verbosity=3,
                 alnfilename=None, treefilename=None, clustalw="clustalw2"):
    """Return selection strings that will select equivalent residues.

    The function aligns two sequences provided in a FASTA file and
    constructs MDAnalysis selection strings of the common atoms. When
    these two strings are applied to the two different proteins they
    will generate AtomGroups of the aligned residues.

    *fastafilename* contains the two un-aligned sequences in FASTA
    format. The reference is assumed to be the first sequence, the
    target the second. ClustalW_ produces a pairwise
    alignment (which is written to a file with suffix .aln).  The
    output contains atom selection strings that select the same atoms
    in the two structures.

    Unless *ref_offset* and/or *target_offset* are specified, the resids
    in the structure are assumed to correspond to the positions in the
    un-aligned sequence, namely the first residue has resid == 1.

    In more complicated cases (e.g. when the resid numbering in the
    structure/psf has gaps due to missing parts), simply provide the
    sequence of resids as they appear in the psf in *ref_resids* or
    *target_resids*, e.g. ::

       target_resids = [a.resid for a in trj.select_atoms('name CA')]

    (This translation table *is* combined with any value for *xxx_offset*!)

    :Arguments:
      *fastafilename*
         FASTA file with first sequence as reference and
         second the one to be aligned (ORDER IS IMPORTANT!)
      *is_aligned*
         False: run clustalw for sequence alignment; True: use
         the alignment in the file (e.g. from STAMP) [``False``]
      *ref_offset*
         add this number to the column number in the FASTA file
         to get the original residue number
      *target_offset*
         same for the target
      *ref_resids*
         sequence of resids as they appear in the reference structure
      *target_resids*
         sequence of resids as they appear in the target
      *alnfilename*
         filename of ClustalW alignment (clustal format) that is
         produced by *clustalw* when *is_aligned* = ``False``.
         ``None`` uses the name and path of *fastafilename* and
         subsititutes the suffix with '.aln'.[``None``]
      *treefilename*
         filename of ClustalW guide tree (Newick format);
         if ``None``  the the filename is generated from *alnfilename*
         with the suffix '.dnd' instead of '.aln' [``None``]
      *clustalw*
         path to the ClustalW (or ClustalW2) binary; only
         needed for *is_aligned* = ``False`` ["clustalw2"]

    :Returns:
      *select_dict*
          dictionary with 'reference' and 'mobile' selection string
          that can be used immediately in :func:`rms_fit_trj` as
          ``select=select_dict``.
    """
    import Bio.SeqIO
    import Bio.AlignIO
    import Bio.Alphabet
    import numpy as np

    protein_gapped = Bio.Alphabet.Gapped(Bio.Alphabet.IUPAC.protein)
    if is_aligned:
        logger.info("Using provided alignment %r", fastafilename)
        with open(fastafilename) as fasta:
            alignment = Bio.AlignIO.read(fasta, "fasta", alphabet=protein_gapped)
    else:
        from Bio.Align.Applications import ClustalwCommandline
        import os.path

        if alnfilename is None:
            filepath, ext = os.path.splitext(fastafilename)
            alnfilename = filepath + '.aln'
        if treefilename is None:
            filepath, ext = os.path.splitext(alnfilename)
            treefilename = filepath + '.dnd'
        run_clustalw = ClustalwCommandline(clustalw, infile=fastafilename, type="protein",
                                           align=True, outfile=alnfilename, newtree=treefilename)
        logger.debug("Aligning sequences in %(fastafilename)r with %(clustalw)r.", vars())
        logger.debug("ClustalW commandline: %r", str(run_clustalw))
        try:
            stdout, stderr = run_clustalw()
        except:
            logger.exception("ClustalW %(clustalw)r failed", vars())
            logger.info("(You can get clustalw2 from http://www.clustal.org/clustal2/)")
            raise
        with open(alnfilename) as aln:
            alignment = Bio.AlignIO.read(aln, "clustal", alphabet=protein_gapped)
        logger.info("Using clustalw sequence alignment {0!r}".format(alnfilename))
        logger.info("ClustalW Newick guide tree was also produced: {0!r}".format(treefilename))

    nseq = len(alignment)
    if nseq != 2:
        raise ValueError("Only two sequences in the alignment can be processed.")

    orig_resids = [ref_resids, target_resids]  # implict assertion that
    # we only have two sequences in the alignment
    offsets = [ref_offset, target_offset]
    for iseq, a in enumerate(alignment):  # need iseq index to change orig_resids
        if orig_resids[iseq] is None:
            # build default: assume consecutive numbering of all
            # residues in the alignment
            GAP = a.seq.alphabet.gap_char
            length = len(a.seq) - a.seq.count(GAP)
            orig_resids[iseq] = np.arange(1, length + 1)
        else:
            orig_resids[iseq] = np.asarray(orig_resids[iseq])
    # add offsets to the sequence <--> resid translation table
    seq2resids = [resids + offset for resids, offset in zip(orig_resids, offsets)]
    del orig_resids
    del offsets

    def resid_factory(alignment, seq2resids):
        """Return a function that gives the resid for a position ipos in
        the nseq'th alignment.

        resid = resid_factory(alignment,seq2resids)
        r = resid(nseq,ipos)

        It is based on a look up table that translates position in the
        alignment to the residue number in the original
        sequence/structure.

        The first index of resid() is the alignmment number, the
        second the position in the alignment.

        seq2resids translates the residues in the sequence to resid
        numbers in the psf. In the simplest case this is a linear map
        but if whole parts such as loops are ommitted from the protein
        the seq2resids may have big gaps.

        Format: a tuple of two numpy arrays; the first array is for
        the reference, the second for the target, The index in each
        array gives the consecutive number of the amino acid in the
        sequence, the value the resid in the structure/psf.

        Note: assumes that alignments have same length and are padded if
        necessary.
        """
        # could maybe use Bio.PDB.StructureAlignment instead?
        nseq = len(alignment)
        t = np.zeros((nseq, alignment.get_alignment_length()), dtype=int)
        for iseq, a in enumerate(alignment):
            GAP = a.seq.alphabet.gap_char
            t[iseq, :] = seq2resids[iseq][np.cumsum(np.where(
                np.array(list(a.seq)) == GAP, 0, 1)) - 1]
            # -1 because seq2resid is index-1 based (resids start at 1)

        def resid(nseq, ipos, t=t):
            return t[nseq, ipos]

        return resid

    resid = resid_factory(alignment, seq2resids)

    res_list = []  # collect individual selection string
    # could collect just resid and type (with/without CB) and
    # then post-process and use ranges for continuous stretches, eg
    # ( resid 1:35 and ( backbone or name CB ) ) or ( resid 36 and backbone ) ...

    GAP = alignment[0].seq.alphabet.gap_char  # should be the same for both seqs
    if GAP != alignment[1].seq.alphabet.gap_char:
        raise ValueError("Different gap characters in sequence 'target' and 'mobile'.")
    for ipos in range(alignment.get_alignment_length()):
        aligned = list(alignment[:, ipos])
        if GAP in aligned:
            continue  # skip residue
        template = "resid %i"
        if 'G' not in aligned:
            # can use CB
            template += " and ( backbone or name CB )"
        else:
            template += " and backbone"
        template = "( " + template + " )"

        res_list.append([template % resid(iseq, ipos) for iseq in range(nseq)])

    sel = np.array(res_list).transpose()

    ref_selection = " or ".join(sel[0])
    target_selection = " or ".join(sel[1])
    return {'reference': ref_selection, 'mobile': target_selection}


def get_matching_atoms(ag1, ag2, tol_mass=0.1, strict=False):
    """Return two atom groups with one-to-one matched atoms.

    The function takes two :class:`~MDAnalysis.core.AtomGroup.AtomGroup`
    instances *ag1* and *ag2* and returns two atom groups *g1* and *g2* that
    consist of atoms so that the mass of atom ``g1[0]`` is the same as the mass
    of atom ``g2[0]``, ``g1[1]`` and ``g2[1]`` etc.

    The current implementation is very simplistic and works on a per-residue basis:

    1. The two groups must contain the same number of residues.
    2. Any residues in each group that have differing number of atoms are discarded.
    3. The masses of corresponding atoms are compared. and if any masses differ
       by more than *tol_mass* the test is considered failed and a
       :exc:`SelectionError` is raised.

    The log file (see :func:`MDAnalysis.start_logging`) will contain detailed
    information about mismatches.

    :Arguments:
      *ag1*, *ag2*
         :class:`~MDAnalysis.core.AtomGroup.AtomGroup` instances that are compared
    :Keywords:
      *tol_mass*
         Reject if the atomic masses for matched atoms differ by more than
         *tol_mass* [0.1]
      *strict*
         ``True``
             Will raise :exc:`SelectioError` if a single atom does not
             match between the two selections.
         ``False`` [default]
             Will try to prepare a matching selection by dropping
             residues with non-matching atoms. See :func:`get_matching_atoms`
             for details.

    :Returns: Tuple ``(g1, g2)`` with :class:`~MDAnalysis.core.AtomGroup.AtomGroup` instances
              that match, atom by atom. The groups are either the original groups if all matches
              or slices of the original groups.

    :Raises: :exc:`SelectionError` if the number of residues does not match or if in the final
             matching masses differ by more than *tol*.

    The algorithm could be improved by using e.g. the Needleman-Wunsch
    algorithm in :mod:`Bio.profile2` to align atoms in each residue (doing a
    global alignment is too expensive).

    .. versionadded:: 0.8

    .. versionchanged:: 0.10.0
       Renamed from :func:`check_same_atoms` to :func:`get_matching_atoms` and now returns
       matching atomgroups (possibly with residues removed)

    """

    if ag1.n_atoms != ag2.n_atoms:
        if ag1.n_residues != ag2.n_residues:
            errmsg = "Reference and trajectory atom selections do not contain "
            "the same number of atoms: \n"
            "atoms:    N_ref={0}, N_traj={1}\n"
            "and also not the same number of residues:\n"
            "residues: N_ref={2}, N_traj={3}\n"
            "\n"
            "(More details can be found in the log file "
            "which can be enabled with 'MDAnalysis.start_logging()')".format(
                ag1.n_atoms, ag2.n_atoms,
                ag1.n_residues, ag2.n_residues)
            dbgmsg = "mismatched residue numbers\n" + \
                "\n".join(["{0} | {1}"  for r1, r2 in
                           zip_longest(ag1.resids, ag2.resids)])
            logger.error(errmsg)
            logger.debug(dbgmsg)
            raise SelectionError(errmsg)
        else:
            msg = ("Reference and trajectory atom selections do not contain "
                   "the same number of atoms: \n"
                   "atoms:    N_ref={0}, N_traj={1}").format(
                ag1.n_atoms, ag2.n_atoms)
            if strict:
                raise SelectionError(msg)

            # continue with trying to creating a valid selection
            warnings.warn(msg + "\nbut we attempt to create a valid selection.",
                          category=SelectionWarning)

        # continue with trying to salvage the selection:
        # - number of atoms is different
        # - number of residues is the same
        # We will remove residues with mismatching number of atoms (e.g. not resolved
        # in an X-ray structure)
        assert ag1.n_residues == ag2.n_residues

        # Alternatively, we could align all atoms but Needleman-Wunsch
        # pairwise2 consumes too much memory for thousands of characters in
        # each sequence. Perhaps a solution would be pairwise alignment per residue.
        #
        # aln_elem = Bio.pairwise2.align.globalms("".join([MDAnalysis.topology.core.guess_atom_element(n) for n in gref.atoms.names]),
        #    "".join([MDAnalysis.topology.core.guess_atom_element(n) for n in models[0].atoms.names]),
        #                               2, -1, -1, -0.1,
        #                               one_alignment_only=True)

        # For now, just remove the residues that don't have matching numbers
        rsize1 = np.array([r.n_atoms for r in ag1.residues])
        rsize2 = np.array([r.n_atoms for r in ag2.residues])
        rsize_mismatches = np.absolute(rsize1 - rsize2)
        mismatch_mask = (rsize_mismatches > 0)
        if np.any(mismatch_mask):
            if strict:
                # diagnostics
                mismatch_resindex = np.arange(ag1.n_residues)[mismatch_mask]
                def log_mismatch(number, ag, rsize, mismatch_resindex=mismatch_resindex):
                    logger.error("Offending residues: group {0}: {1}".format(
                            number,
                            ", ".join(["{0[0]}{0[1]} ({0[2]})".format(r) for r in
                                       zip(ag.resnames[mismatch_resindex],
                                           ag.resids[mismatch_resindex],
                                           rsize[mismatch_resindex]
                                       )])))
                logger.error("Found {0} residues with non-matching numbers of atoms (#)".format(
                        mismatch_mask.sum()))
                log_mismatch(1, ag1, rsize1)
                log_mismatch(2, ag2, rsize2)

                raise SelectionError("Different number of atoms in some residues. "
                                     "(Use strict=False to attempt using matching atoms only.)")

            def get_atoms_byres(g, match_mask=np.logical_not(mismatch_mask)):
                # not pretty... but need to do things on a per-atom basis in order
                # to preserve original selection
                ag = g.atoms
                good = ag.resids[match_mask]
                resids = np.array([a.resid for a in ag])  # resid for each atom
                ix_good = np.in1d(resids, good)   # boolean array for all matching atoms
                return ag[np.arange(len(ag))[ix_good]]   # workaround for missing boolean indexing
            _ag1 = get_atoms_byres(ag1)
            _ag2 = get_atoms_byres(ag2)

            # diagnostics
            # (ugly workaround for missing boolean indexing of AtomGroup)
            # note: ag[arange(len(ag))[boolean]] is ~2x faster than ag[where[boolean]]
            mismatch_resindex = np.arange(ag1.n_residues)[mismatch_mask]
            logger.warn("Removed {0} residues with non-matching numbers of atoms".format(
                    mismatch_mask.sum()))
            logger.debug("Removed residue ids: group 1: {0}".format(ag1.resids[mismatch_resindex]))
            logger.debug("Removed residue ids: group 2: {0}".format(ag2.resids[mismatch_resindex]))
            # replace after logging (still need old ag1 and ag2 for diagnostics)
            ag1 = _ag1
            ag2 = _ag2
            del _ag1, _ag2

    mass_mismatches = (np.absolute(ag1.masses - ag2.masses) > tol_mass)
    if np.any(mass_mismatches):
        # Test 2 failed.
        # diagnostic output:
        # (ugly workaround because boolean indexing is not yet working for atomgroups)
        assert ag1.n_atoms == ag2.n_atoms
        mismatch_atomindex = np.arange(ag1.n_atoms)[mass_mismatches]

        logger.error("Atoms: reference | trajectory")
        for ar, at in zip(ag1[mismatch_atomindex], ag2[mismatch_atomindex]):
            logger.error("{0!s:>4} {1:3d} {2!s:>3} {3!s:>3} {4:6.3f}  |  {5!s:>4} {6:3d} {7!s:>3} {8!s:>3} {9:6.3f}".format(ar.segid, ar.resid, ar.resname, ar.name, ar.mass,
                          at.segid, at.resid, at.resname, at.name, at.mass))
        errmsg = ("Inconsistent selections, masses differ by more than {0}; " + \
            "mis-matching atoms are shown above.").format(tol_mass)
        logger.error(errmsg)
        raise SelectionError(errmsg)
    return ag1, ag2
