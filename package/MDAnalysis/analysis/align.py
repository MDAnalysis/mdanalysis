# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

"""
Coordinate fitting and alignment --- :mod:`MDAnalysis.analysis.align`
=====================================================================

:Author: Oliver Beckstein, Joshua Adelman
:Year: 2010--2011
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
   :mod:`MDAnalysis.core.qcprot`
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
   >>> from MDAnalysis.tests.datafiles import PSF, DCD, PDB_small


In the simplest case, we can simply calculate the C-alpha RMSD between
two structures, using :func:`rmsd`::

   >>> ref = Universe(PDB_small)
   >>> mobile = Universe(PSF,DCD)
   >>> rmsd(mobile.atoms.CA.coordinates(), ref.atoms.CA.coordinates())
   18.858259026820352

Note that in this example translations have not been removed. In order
to look at the pure rotation one needs to superimpose the centres of
mass (or geometry) first:

   >>> ref0 =  ref.atoms.CA.coordinates() - ref.atoms.CA.centerOfMass()
   >>> mobile0 =  mobile.atoms.CA.coordinates() - mobile.atoms.CA.centerOfMass()
   >>> rmsd(mobile0, ref0)
    6.8093965864717951

The rotation matrix that superimposes *mobile* on *ref* while
minimizing the CA-RMSD is obtained with the :func:`rotation_matrix`
function ::

   >>> R, rmsd = rotation_matrix(mobile0, ref0)
   >>> print rmsd
   6.8093965864717951
   >>> print R
   [[ 0.14514539 -0.27259113  0.95111876]
    [ 0.88652593  0.46267112 -0.00268642]
    [-0.43932289  0.84358136  0.30881368]]

Putting all this together one can superimpose all of *mobile* onto *ref*::

   >>> mobile.atoms.translate(-mobile.atoms.CA.centerOfMass())
   >>> mobile.atoms.rotate(R)
   >>> mobile.atoms.translate(ref.atoms.CA.centerOfMass())
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

.. autofunction:: rmsd
.. autofunction:: alignto
.. autofunction:: rms_fit_trj


Helper functions
----------------

The following functions are used by the other functions in this
module. They are probably of more interest to developers than to
normal users.

.. autofunction:: fasta2select

"""

import numpy

import MDAnalysis
import MDAnalysis.core.qcprot as qcp
from MDAnalysis import SelectionError
from MDAnalysis.core.log import ProgressMeter
from MDAnalysis.core.util import asiterable
from MDAnalysis.analysis.rms import rmsd, _process_selection

import os.path

import logging
logger = logging.getLogger('MDAnalysis.analysis.align')


def rotation_matrix(a,b, weights=None):
    """Returns the 3x3 rotation matrix for RMSD fitting coordinate sets *a* and *b*.

    The rotation matrix *R* transforms *a* to overlap with *b* (i.e. *b* is the
    reference structure):

       *b* = *R* . *a*

    :Arguments:
       *a*
          coordinates that are to be rotated ("mobile set"); array of N atoms
          of shape N*3 as generated by, e.g.,
          :meth:`MDAnalysis.core.AtomGroup.AtomGroup.coordinates`.
       *b*
          reference coordinates; array of N atoms of shape N*3 as generated by,
          e.g., :meth:`MDAnalysis.core.AtomGroup.AtomGroup.coordinates`.
       *weights*
          array of floats of size N for doing weighted RMSD fitting (e.g. the
          masses of the atoms)

    :Returns: ``(R, rmsd)`` rmsd and rotation matrix *R*


    *R* can be used as an argument for
    :meth:`MDAnalysis.core.AtomGroup.AtomGroup.rotate` to generate a rotated
    selection, e.g. ::

      >>> R = rotation_matrix(A.selectAtoms('backbone').coordinates(), B.selectAtoms('backbone').coordinates())
      >>> A.atoms.rotate(R)
      >>> A.atoms.write("rotated.pdb")

    Note that the function does *not* shift the centers of mass or geometry;
    this needs to be done by the user.

    .. SeeAlso:: :func:`rmsd` calculates the RMSD between *a* and *b*; for
                 fitting a whole trajectory it is more efficient to use
                 :func:`rms_fit_trj`. A complete fit of two structures can be
                 done with :func:`alignto`.
    """
    if not weights is None:
        # weights are constructed as relative to the mean
        weights = numpy.asarray(weights)/numpy.mean(weights)
    rot = numpy.zeros(9, dtype=numpy.float64)
    rmsd = qcp.CalcRMSDRotationalMatrix(a.T.astype(numpy.float64), b.T.astype(numpy.float64),
                                        b.shape[0], rot, weights)
    return numpy.matrix(rot.reshape(3,3)), rmsd


def alignto(mobile, reference, select="all", mass_weighted=False,
            subselection=None):
    """Spatially align *mobile* to *reference* by doing a RMSD fit on *select* atoms.

    The superposition is done in the following way:

    1. A rotation matrix is computed that minimizes the RMSD between
       the coordinates of `mobile.selectAtoms(sel1)` and
       `reference.selectAtoms(sel2)`; before the rotation, *mobile* is
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
            :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.selectAtoms` that produces identical
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
      *subselection*
         Apply the transformation only to this selection.

         ``None`` [default]
             Apply to `mobile.universe.atoms` (i.e. all atoms in the
             context of the selection from *mobile* such as the rest of a
             protein, ligands and the surrounding water)
         *selection-string*
             Apply to `mobile.selectAtoms(selection-string)`
         :class:`~MDAnalysis.core.AtomGroup.AtomGroup`
             Apply to the arbitrary group of atoms

    :Returns: RMSD before and after spatial alignment.

    .. SeeAlso:: For RMSD-fitting trajectories it is more efficient to
                 use :func:`rms_fit_trj`.
    """
    if select in ('all', None):
        # keep the EXACT order in the input AtomGroups; selectAtoms('all')
        # orders them by index, which can lead to wrong results if the user
        # has crafted mobile and reference to match atom by atom
        mobile_atoms = mobile.atoms
        ref_atoms = reference.atoms
    else:
        select = _process_selection(select)
        mobile_atoms = mobile.selectAtoms(*select['mobile'])
        ref_atoms = reference.selectAtoms(*select['reference'])
    if mass_weighted:
        weights = ref_atoms.masses()/numpy.mean(ref_atoms.masses())
        ref_com = ref_atoms.centerOfMass()
        mobile_com = mobile_atoms.centerOfMass()
    else:
        weights = None
        ref_com = ref_atoms.centerOfGeometry()
        mobile_com = mobile_atoms.centerOfGeometry()

    ref_coordinates = ref_atoms.coordinates() - ref_com
    mobile_coordinates = mobile_atoms.coordinates() - mobile_com

    old_rmsd = rmsd(mobile_atoms.coordinates(), ref_atoms.coordinates())

    R, new_rmsd = rotation_matrix(mobile_coordinates, ref_coordinates, weights=weights)

    if subselection is None:
        atoms = mobile.universe.atoms
    elif type(subselection) is str:
        atoms = mobile.selectAtoms(subselection)
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
                mass_weighted=False, tol_mass=0.1):
    """RMS-fit trajectory to a reference structure using a selection.

    :Arguments:
      *traj*
         trajectory, :class:`MDAnalysis.Universe` object
      *reference*
         reference coordinates; :class:`MDAnalysis.Universe` object
         (uses the current time step of the object)
      *select*
         1. any valid selection string for
            :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.selectAtoms` that produces identical
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

    Both reference and trajectory must be :class:`MDAnalysis.Universe`
    instances. If they contain a trajectory then it is used. The
    output file format is the same as the input *traj*.

    .. _ClustalW: http://www.clustal.org/
    .. _STAMP: http://www.compbio.dundee.ac.uk/manuals/stamp.4.2/
    """

    frames = traj.trajectory

    if filename is None:
        path,fn = os.path.split(frames.filename)
        filename = os.path.join(path,prefix+fn)
        writer = frames.Writer(filename, remarks='RMS fitted trajectory to reference')
    else:
        writer = frames.OtherWriter(filename)

    select = _process_selection(select)

    ref_atoms = reference.selectAtoms(*select['reference'])
    traj_atoms = traj.selectAtoms(*select['mobile'])
    natoms = traj_atoms.numberOfAtoms()
    if len(ref_atoms) != len(traj_atoms):
        raise SelectionError("Reference and trajectory atom selections do not contain "+
                             "the same number of atoms: N_ref=%d, N_traj=%d" % \
                             (len(ref_atoms), len(traj_atoms)))
    logger.info("RMS-fitting on %d atoms." % len(ref_atoms))
    mass_mismatches = (numpy.absolute(ref_atoms.masses() - traj_atoms.masses()) > tol_mass)
    if numpy.any(mass_mismatches):
        # diagnostic output:
        logger.error("Atoms: reference | trajectory")
        for ar,at in zip(ref_atoms,traj_atoms):
            if ar.name != at.name:
                logger.error("%4s %3d %3s %3s %6.3f  |  %4s %3d %3s %3s %6.3f" %  \
                      (ar.segid, ar.resid, ar.resname, ar.name, ar.mass,
                       at.segid, at.resid, at.resname, at.name, at.mass,))
        errmsg = "Inconsistent selections, masses differ by more than %f; mis-matching atoms are shown above." % tol_mass
        logger.error(errmsg)
        raise SelectionError(errmsg)
    del mass_mismatches

    if mass_weighted:
        # if performing a mass-weighted alignment/rmsd calculation
        weight = ref_atoms.masses()/ref_atoms.masses().mean()
    else:
        weight = None

    # reference centre of mass system
    # (compatibility with pre 1.0 numpy: explicitly cast coords to float32)
    ref_com = ref_atoms.centerOfMass().astype(numpy.float32)
    ref_coordinates = ref_atoms.coordinates() - ref_com

    # allocate the array for selection atom coords
    traj_coordinates = traj_atoms.coordinates().copy()

    # RMSD timeseries
    nframes = len(frames)
    rmsd = numpy.zeros((nframes,))

    # R: rotation matrix that aligns r-r_com, x~-x~com
    #    (x~: selected coordinates, x: all coordinates)
    # Final transformed traj coordinates: x' = (x-x~_com)*R + ref_com
    rot = numpy.zeros(9,dtype=numpy.float64)      # allocate space for calculation
    R = numpy.matrix(rot.reshape(3,3))

    percentage = ProgressMeter(nframes, interval=10,
                               format="Fitted frame %(step)5d/%(numsteps)d  [%(percentage)5.1f%%]\r")

    for k,ts in enumerate(frames):
        # shift coordinates for rotation fitting
        # selection is updated with the time frame
        x_com = traj_atoms.centerOfMass().astype(numpy.float32)
        traj_coordinates[:] = traj_atoms.coordinates() - x_com

        # Need to transpose coordinates such that the coordinate array is
        # 3xN instead of Nx3. Also qcp requires that the dtype be float64
        # (I think we swapped the position of ref and traj in CalcRMSDRotationalMatrix
        # so that R acts **to the left** and can be broadcasted; we're saving
        # one transpose. [orbeckst])
        rmsd[k] = qcp.CalcRMSDRotationalMatrix(ref_coordinates.T.astype(numpy.float64),
                                               traj_coordinates.T.astype(numpy.float64),
                                               natoms, rot, weight)
        R[:,:] = rot.reshape(3,3)

        # Transform each atom in the trajectory (use inplace ops to avoid copying arrays)
        # (Marginally (~3%) faster than "ts._pos[:] = (ts._pos - x_com) * R + ref_com".)
        ts._pos   -= x_com
        ts._pos[:] = ts._pos * R # R acts to the left & is broadcasted N times.
        ts._pos   += ref_com

        writer.write(traj.atoms) # write whole input trajectory system
        percentage.echo(ts.frame)
    logger.info("Wrote %d RMS-fitted coordinate frames to file %r",
                frames.numframes, filename)
    if not rmsdfile is None:
        numpy.savetxt(rmsdfile,rmsd)
        logger.info("Wrote RMSD timeseries  to file %r", rmsdfile)

def fasta2select(fastafilename,is_aligned=False,
                 ref_resids=None, target_resids=None,
                 ref_offset=0,target_offset=0,verbosity=3,
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

       target_resids = [a.resid for a in trj.selectAtoms('name CA')]

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
    import Bio.SeqIO, Bio.AlignIO, Bio.Alphabet
    import numpy

    protein_gapped = Bio.Alphabet.Gapped(Bio.Alphabet.IUPAC.protein)
    if is_aligned:
        logger.info("Using provided alignment %r", fastafilename)
        with open(fastafilename) as fasta:
            alignment = Bio.AlignIO.read(fasta, "fasta", alphabet=protein_gapped)
    else:
        from Bio.Align.Applications import ClustalwCommandline
        import os.path
        if alnfilename is None:
            filepath,ext = os.path.splitext(fastafilename)
            alnfilename = filepath + '.aln'
        if treefilename is None:
            filepath,ext = os.path.splitext(alnfilename)
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
        logger.info("Using clustalw sequence alignment %r" % alnfilename)
        logger.info("ClustalW Newick guide tree was also produced: %r" % treefilename)

    nseq = len(alignment)
    if nseq != 2:
        raise ValueError("Only two sequences in the alignment can be processed.")

    orig_resids = [ref_resids,target_resids] # implict assertion that
                                             # we only have two sequences in the alignment
    offsets = [ref_offset,target_offset]
    for iseq,a in enumerate(alignment):      # need iseq index to change orig_resids
        if orig_resids[iseq] is None:
            # build default: assume consecutive numbering of all
            # residues in the alignment
            GAP = a.seq.alphabet.gap_char
            length = len(a.seq) - a.seq.count(GAP)
            orig_resids[iseq] = numpy.arange(1,length+1)
        else:
            orig_resids[iseq] = numpy.asarray(orig_resids[iseq])
    # add offsets to the sequence <--> resid translation table
    seq2resids = [resids + offset for resids,offset in zip(orig_resids,offsets)]
    del orig_resids
    del offsets

    def resid_factory(alignment,seq2resids):
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
        t = numpy.zeros((nseq,alignment.get_alignment_length()),dtype=int)
        for iseq,a in enumerate(alignment):
            GAP = a.seq.alphabet.gap_char
            t[iseq,:] = seq2resids[iseq][numpy.cumsum(numpy.where(
                        numpy.array(list(a.seq))==GAP,0,1)) - 1]
            # -1 because seq2resid is index-1 based (resids start at 1)

        def resid(nseq,ipos,t=t):
            return t[nseq,ipos]
        return resid
    resid = resid_factory(alignment,seq2resids)

    res_list = []     # collect individual selection string
    # could collect just resid and type (with/without CB) and
    # then post-process and use ranges for continuous stretches, eg
    # ( resid 1:35 and ( backbone or name CB ) ) or ( resid 36 and backbone ) ...

    GAP = alignment[0].seq.alphabet.gap_char        # should be the same for both seqs
    if GAP != alignment[1].seq.alphabet.gap_char:
        raise ValueError("Different gap characters in sequence 'target' and 'mobile'.")
    for ipos in xrange(alignment.get_alignment_length()):
        aligned = list(alignment[:, ipos])
        if GAP in aligned:
            continue       # skip residue
        template = "resid %i"
        if 'G' not in aligned:
            # can use CB
            template += " and ( backbone or name CB )"
        else:
            template += " and backbone"
        template = "( "+template+" )"

        res_list.append([template % resid(iseq,ipos) for iseq in xrange(nseq)])

    sel = numpy.array(res_list).transpose()

    ref_selection =  " or ".join(sel[0])
    target_selection =  " or ".join(sel[1])
    return {'reference':ref_selection, 'mobile':target_selection}

