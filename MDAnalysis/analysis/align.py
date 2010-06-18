# aligning of coordinate frames
# Part of MDAnalysis http://mdanalysis.googlecode.com
# Copyright (c) 2006-2010 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# Released under the GNU Public Licence, v2+
"""
Coordinate fitting and alignment
================================

:Author: Oliver Beckstein
:Year: 2010
:Copyright: GNU Public License v3

The script show how to select a group of atoms (such as the C-alphas),
calculate the RMSD and transformation matrix, and apply the
transformation to the current frame of a trajectory to obtain the
rotated structure.

In this example we use files provided as part of the MDAnalysis test
suite (in the variables PSF, DCD, and PDB_small)::

   >>> from MDAnalysis import *
   >>> from MDAnalysis.analysis.align import *
   >>> from MDAnalysis.tests.datafiles import PSF, DCD, PDB_small
   >>> ref = Universe(PSF, PDB_small)   # reference structure 1AKE
   >>> trj = Universe(PSF, DCD)         # trajectory of change 1AKE->4AKE
   >>> rms_fit_trj(trj, ref, filename='rmsfit.dcd')


It is also possible to align two arbitrary structures by providing a
mapping between atoms based on a sequence alignment. This allows
fitting of structural homologs or wild type and mutant.

If a alignment was provided as "sequences.aln" one would first produce
the appropriate MDAnalysis selections with the :func:`fasta2select`
function and the feed the resulting dictionary to :func:`rms_fit_trj`::

   >>> seldict = fasta2select('sequences.aln')
   >>> rms_fit_trj(trj, ref, filename='rmsfit.dcd', select=seldict)

(See the documentation of the functions for this advanced usage.)
"""

import numpy
import MDAnalysis
import MDAnalysis.core.rms_fitting
from MDAnalysis import SelectionError

import os.path
import sys

import logging
logger = logging.getLogger('MDAnalysis.analysis.align')

def rmsd(a,b):
    """Returns RMSD between two coordinate sets a and b."""
    return numpy.sqrt(numpy.sum(numpy.power(a-b,2))/a.shape[0])

def echo(s=''):
    """Simple string output that immediately prints to the console."""
    sys.stderr.write(s)
    sys.stderr.flush()

def rms_fit_trj(traj,ref,select='backbone',filename=None,prefix='rmsfit_'):
    """RMS-fit trajectory to a reference structure using a selection.

      rms_fit_trj(traj, ref, 'backbone or name CB or name OT*')

    :Arguments:
      *traj*
         trajectory, :class:`MDAnalysis.Universe` object
      *ref*
         reference coordinates; :class:`MDAnalysis.Universe` object
         (uses the current time step of the object)
      *filename*
         file name for the RMS-fitted trajectory or pdb; defaults to the 
         original trajectory filename (from *traj*) with *prefix* prepended
      *prefix*
         prefix for autogenerating the new output filename
      *select*
         any valid selection string for
         :meth:`MDAnalysis.AtomGroup.selectAtoms` that produces identical
         selections in *traj* and *ref* or dictionary {'reference':sel1,
         'target':sel2}.  The :func:`fasta2select` function returns such a
         dictionary based on a ClustalW_ or STAMP_ sequence alignment.

    Both reference and trajectory must be :class:`MDAnalysis.Universe`
    instances. If they contain a trajectory then it is used. The
    output file format is the same as the input *traj*.

    .. Note:: This function is flexible but slow. In order to fit long
       trajectory it is suggested to only align one frame and then use
       that frame with a more efficient program (that otherwise could
       not easily deal with aligning different molecules).

    .. _ClustalW: http://www.clustal.org/
    .. _STAMP: http://www.compbio.dundee.ac.uk/manuals/stamp.4.2/
    """

    frames = traj.trajectory

    if filename is None:
        path,fn = os.path.split(frames.filename)
        filename = os.path.join(path,prefix+fn)
    if type(select) is str:
        select = {'reference':select,'target':select}

    # TODO: dealing with different formats is not nice; this should be
    # rewritten with the trajectory API in mind...
    if isinstance(frames,MDAnalysis.coordinates.DCD.DCDReader):
        writer = MDAnalysis.coordinates.DCD.DCDWriter(
            filename,frames.numatoms,
            frames.start_timestep,
            frames.skip_timestep,
            frames.delta,
            remarks='RMS fitted trajectory to ref')
    if isinstance(frames,MDAnalysis.coordinates.XTC.XTCReader):
        writer = MDAnalysis.coordinates.XTC.XTCWriter(  # untested!
            filename,frames.numatoms,
            frames.start_timestep,
            frames.skip_timestep,
            frames.delta,
            remarks='RMS fitted trajectory to ref')
    if isinstance(frames,MDAnalysis.coordinates.TRR.TRRReader):
        writer = MDAnalysis.coordinates.TRR.TRRWriter(  # untested!
            filename,frames.numatoms,
            frames.start_timestep,
            frames.skip_timestep,
            frames.delta,
            remarks='RMS fitted trajectory to ref')
    elif isinstance(frames,MDAnalysis.coordinates.PDB.PDBReader):
        writer = MDAnalysis.coordinates.PDB.PDBWriter( # still working like this?
            filename, universe=traj,
            remarks='RMS fitted pdb frame to ref')
    else:
        raise TypeError('traj must contain a DCD, XTC, TRR, or a PDB.')

    ref_atoms = ref.selectAtoms(select['reference'])
    traj_atoms = traj.selectAtoms(select['target'])
    if len(ref_atoms) != len(traj_atoms):
        raise SelectionError("Reference and trajectory atom selections do not contain "+
                             "the same number of atoms: N_ref=%d, N_traj=%d" % \
                             (len(ref_atoms), len(traj_atoms)))
    logger.info("RMS-fitting on %d atoms." % len(ref_atoms))
    mass_mismatches = (numpy.absolute(ref_atoms.masses() - traj_atoms.masses()) > 0)
    if numpy.any(mass_mismatches):
        # diagnostic output:
        logger.error("Atoms: reference | trajectory")
        for ar,at in zip(ref_atoms,traj_atoms):
            if ar.name != at.name:
                logger.error("%4s %3d %3s %3s %6.3f  |  %4s %3d %3s %3s %6.3f" %  \
                      (ar.segid, ar.resid, ar.resname, ar.name, ar.mass,
                       at.segid, at.resid, at.resname, at.name, at.mass,))
        errmsg = "Inconsistent selections, masses don't match; mis-matching atoms are shown above."
        logger.error(errmsg)
        raise SelectionError(errmsg)
    del mass_mismatches
    masses = ref_atoms.masses()

    # reference centre of mass system
    # (compatibility with pre 1.0 numpy: explicitly cast coords to float32)
    ref_com = ref_atoms.centerOfMass().astype(numpy.float32)
    ref_coordinates = ref_atoms.coordinates() - ref_com

    # allocate the array for selection atom coords
    traj_coordinates = traj_atoms.coordinates().copy()

    # R: rotation matrix that aligns r-r_com, x~-x~com   
    #    (x~: selected coordinates, x: all coordinates)
    # Final transformed traj coordinates: x' = (x-x~_com)*R + ref_com
    for ts in frames:
        # shift coordinates for rotation fitting
        # selection is updated with the time frame
        x_com = traj_atoms.centerOfMass().astype(numpy.float32)
        traj_coordinates[:] = traj_atoms.coordinates() - x_com
        R = numpy.matrix(MDAnalysis.core.rms_fitting.rms_rotation_matrix(
                traj_coordinates,ref_coordinates,masses),dtype=numpy.float32)
        # Transform each atom in the trajectory (use inplace ops to avoid copying arrays)
        ts._pos   -= x_com
        ts._pos[:] = ts._pos * R # R acts to the left & is broadcasted N times.
        ts._pos   += ref_com
        writer.write_next_timestep(ts)
        # for debugging:
        # rmsd_old = rmsd(ref_atoms.coordinates(),traj_coordinates)
        # rmsd_new = rmsd(ref_atoms.coordinates(),traj_atoms.coordinates())
        # logger.debug("Fitted frame %5d/%d  [%5.1f%%]  %5.2fA --> %5.2fA  |translation|=%.2fA\r" % \
        #            (ts.frame,frames.numframes,100.0*ts.frame/frames.numframes,
        #             rmsd_old, rmsd_new, rmsd(x_com,ref_com)) )
        if ts.frame % 10 == 0:
            echo("Fitted frame %5d/%d  [%5.1f%%]\r" % \
                (ts.frame,frames.numframes,100.0*ts.frame/frames.numframes))
    # done
    echo("Fitted frame %5d/%d  [%5.1f%%]\r\n" % \
             (ts.frame,frames.numframes,100.0*ts.frame/frames.numframes))
    logger.info("Wrote %d RMS-fitted coordinate frames to file %r" % (frames.numframes, filename))

def fasta2select(fastafilename,is_aligned=False,
                 ref_resids=None, target_resids=None,
                 ref_offset=0,target_offset=0,verbosity=3):
    """Align two sequences from a fasta file and construct a MDAnalysis
    selection string of the common atoms.

      fasta2select(fastafilename) -> selection_dict

    *fastafilename* contains the the two un-aligned sequences in FASTA
    format. The reference is assumed to be be the first sequence, the
    target the second. :program:`ClustalW` produces a pairwise
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
      fastafilename    
         FASTA file with first sequence as reference and
         second the one to be aligned (ORDER IS IMPORTANT!)
      is_aligned
         False: run clustalw for sequence alignment; True: use
         the alignment in the file (e.g. from STAMP)
      ref_offset
         add this number to the column number in the FASTA
         to get the original residue number
      target_offset
         same for the target
      ref_resids
         sequence of resids as they appear in the reference structure
      target_resids
         sequence of resids as they appear in the target
 
    :Returns:
      select_dict      
          dictionary with 'reference' and 'target' selection string
          that can be used immediately in :func:`rms_fit_trj` as
          ``select=select_dict``.
    """
    import Bio
    import numpy
    
    if is_aligned:
        import Bio.SeqIO, Bio.Alphabet
        protein_gapped = Bio.Alphabet.Gapped(Bio.Alphabet.IUPAC.protein)
        fasta = open(fastafilename)
        try:
            alignment = Bio.SeqIO.to_alignment(
                Bio.SeqIO.FastaIO.FastaIterator(fasta,alphabet=protein_gapped),
                alphabet=protein_gapped)
        finally:
            fasta.close()
        logger.info("Using provided alignment, %s", fastafilename)
    else:
        import Bio.Clustalw
        import os.path
        filepath,ext = os.path.splitext(fastafilename)
        alnfilename = filepath + '.aln'
        cline = Bio.Clustalw.MultipleAlignCL(fastafilename)
        cline.set_output(alnfilename)
        cline.set_type('protein')
        alignment = Bio.Clustalw.do_alignment(cline)
        logger.info("Using clustalw sequence alignment, %s.\n" % alnfilename)

    nseq = len(alignment._records)    # the stupid class should provide __len__ !
    if nseq != 2:
        raise ValueError("Only two sequences in the alignment can be processed.")

    orig_resids = [ref_resids,target_resids] # implict assertion that
                                             # we only have to sequenceses in the alignment
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
        nseq = len(alignment._records)
        t = numpy.zeros((nseq,alignment.get_alignment_length()),dtype=int)
        for iseq,a in enumerate(alignment):
            GAP = a.seq.alphabet.gap_char
            t[iseq,:] = seq2resids[iseq][numpy.cumsum(numpy.where(
                        numpy.array(list(a.seq))==GAP,0,1)) - 1]  
            # -1 because seq2resid is index-1 based (resids start at 1)

        def resid(nseq,ipos):
            return t[nseq,ipos]
        return resid
    resid = resid_factory(alignment,seq2resids)

    res_list = []     # collect individual selection string
    # could collect just resid and type (with/without CB) and
    # then post-process and use ranges for continuous stretches, eg
    # ( resid 1:35 and ( backbone or name CB ) ) or ( resid 36 and backbone ) ...

    GAP = alignment.get_seq_by_num(0).alphabet.gap_char # should be same for all seqs
    for ipos in xrange(alignment.get_alignment_length()):
        aligned = list(alignment.get_column(ipos))
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
    return {'reference':ref_selection, 'target':target_selection}

