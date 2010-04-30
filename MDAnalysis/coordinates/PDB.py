# $Id: PDB.py 101 2008-05-18 13:19:06Z orbeckst $
"""
PDB structure files in MDAnalysis
=================================

MDAnalysis reads coordinates from PDB files and additional optional
data such as B-factors. It is also possible to substitute a PDB file
instead of PSF file in order to define the list of atoms (but no
connectivity information will be  available in this case).

The :mod:`PDB` module makes heavy use of Biopython's :mod:`Bio.PDB`:

  Hamelryck, T., Manderick, B. (2003) PDB parser and structure class 
  implemented in Python. Bioinformatics, 19, 2308-2310.

  http://biopython.org

but replaces the standard PDB file parser with one that uses the
:class:`MDAnalysis.coordinates.pdb.extensions.SloppyStructureBuilder`
to cope with very large pdb files as commonly encountered in MD
simulations.
"""

try:
    # BioPython is overkill but potentially extensible (altLoc etc)
    import Bio.PDB
except ImportError:
    raise ImportError("No PDB I/O functionality. Install biopython.")

import numpy
from DCD import Timestep
import MDAnalysis.core.util as util
import pdb.extensions

class PDBReader:
    """Read a pdb file into a BioPython pdb structure.

    The coordinates are also supplied as one numpy array and wrapped
    into a Timestep object; attributes are set so that the PDBReader
    object superficially resembles the DCDReader object.
    """
    def __init__(self,pdbfilename):
        pdb_id = "0UNK"
        self.pdb = pdb.extensions.get_structure(pdbfilename, pdb_id)
        pos = numpy.array([atom.coord for atom in self.pdb.get_atoms()])
        self.pdbfilename = pdbfilename
        self.filename = self.pdbfilename
        self.numatoms = pos.shape[0]
        self.numframes = 1
        self.fixed = 0          # parse B field for fixed atoms?
        self.skip = 1
        self.periodic = False
        self.ts = Timestep(pos)
        del pos

    def __len__(self):
        return self.numframes
    def __iter__(self):
        def iterPDB():
            yield self.ts   # just a single frame available
            raise StopIteration
        return iterPDB()
    def __getitem__(self, frame):
        if frame != 0:
            raise IndexError('PDBReader only contains a single frame at index 0')
        return self.ts
    def __repr__(self):
            return "<MDAnalysis.PDB.PDBReader '"+ self.filename + "' with " + repr(self.numframes) + " frames of " + repr(self.numatoms) + " atoms (" + repr(self.fixed) + " fixed)>"


class PDBWriter:
    """Write out the current time step as a pdb file.

    This is not cleanly implemented at the moment. One must supply a
    universe, even though this is nominally an optional argument. The
    class behaves slightly differently depending on if the structure
    was loaded from a PDB (then the full-fledged Bio.PDB writer is
    used) or if this is really only an atom selection (then a less
    sophistiocated writer is employed).
    """

    # PDBWriter is a bit more complicated than the DCDWriter in the
    # sense that a DCD frame only contains coordinate information. The
    # PDB contains atom data as well and hence it MUST access the
    # universe. In order to present a unified (and backwards
    # compatible) interface we must keep the universe argument an
    # optional keyword argument even though it really is required.
    
    def __init__(self,pdbfilename,universe=None,multi=False,**kwargs):
        """pdbwriter = PDBWriter(<pdbfilename>,PDBstructure==universe.pdb.pdb,**kwargs)
        :Arguments:
        pdbfilename     filename; if multi=True, embed a %%d formatstring
                        so that write_next_timestep() can insert the frame number
        universe        supply a universe [really REQUIRED; optional only for compatibility]
        multi           False: write a single structure to a single pdb
                        True: write all frames to multiple pdb files
        """
        import Bio.PDB.Structure
        self.universe = universe
        try:
            self.PDBstructure = universe.pdb.pdb
        except AttributeError:
            self.PDBstructure = None
        self.filename = pdbfilename
        self.multi = multi
        if self.multi:
            raise NotImplementedError('Sorry, multi=True does not work yet.')
        if self.PDBstructure is not None and not isinstance(self.PDBstructure,Bio.PDB.Structure.Structure):
            raise TypeError('If defined, PDBstructure must be a Bio.PDB.Structure.Structure, eg '
                            'Universe.pdb.pdb.')
    def write_next_timestep(self,ts=None):
        self.write(ts)
    def write(self,ts=None):
        """Write timestep as a pdb file.

        If ts=None then we try to get the current one from the universe.
        """        
        if self.PDBstructure is None:
            if self.universe is None:
                warnings.warn("PDBWriter: Not writing frame as no universe supplied.")
                return            
            # primitive PDB writing (ignores timestep argument)
            ppw = PrimitivePDBWriter(self.filename)
            ppw.write(self.universe.selectAtoms('all'))
            ppw.close()
        else:
            # full fledged PDB writer
            # Let's cheat and use universe.pdb.pdb: modify coordinates
            # and save...
            if ts is None:
                try:
                    ts = self.universe.trajectory.ts
                except AttributeError:
                    warnings.warn("PDBWriter: Not writing frame as neither universe nor timestep supplied.")
                    return            
            for a,pos in zip(self.PDBstructure.get_atoms(), ts._pos):
                a.set_coord(pos)
            io = pdb.extensions.SloppyPDBIO()
            io.set_structure(self.PDBstructure)
            io.save(self.filename)
    

class PrimitivePDBWriter(object):
    """PDB writer that implements a subset of the PDB 3.2 standard.
    http://www.wwpdb.org/documentation/format32/v3.2.html
    """
    #          1         2         3         4         5         6         7         8
    # 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.
    # ATOM__seria nameAres CressI   xxxxxxxxyyyyyyyyzzzzzzzzOCCUPAtempft          elCH
    # ATOM  %5d   %-4s %-3s %4d %1s %8.3f   %8.3f   %8.3f   %6.2f %6.2f           %2s
    #                 %1s  %1s                                                      %2d
    #            =        =      ===                                    ==========
    # ATOM  %5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2d
    # ATOM  %(serial)5d %(name)-4s%(altLoc)1s%(resName)-3s %(chainID)1s%(resSeq)4d%(iCode)1s   %(x)8.3f%(y)8.3f%(z)8.3f%(occupancy)6.2f%(tempFactor)6.2f          %(element)2s%(charge)2d

    fmt = {'ATOM':   "ATOM  %(serial)5d %(name)-4s%(altLoc)1s%(resName)-3s %(chainID)1s%(resSeq)4d%(iCode)1s   %(x)8.3f%(y)8.3f%(z)8.3f%(occupancy)6.2f%(tempFactor)6.2f          %(element)2s%(charge)2d\n",
           'REMARK': "REMARK     %s\n",
           'TITLE':  "TITLE    %s\n",
           }

    def __init__(self,filename):
        self.filename = util.filename(filename,ext='pdb')
        self.pdb = open(self.filename,'w')

    def close(self):
        self.pdb.close()

    def write(self,selection,frame=None):
        """Write selection at current trajectory frame to file.

        write(selection,frame=FRAME)

        selection         MDAnalysis AtomGroup
        frame             optionally move to frame FRAME
        """
        u = selection.universe
        if frame is not None:            
            u.trajectory[frame]  # advance to frame
        else:
            try:
                frame = u.trajectory.ts.frame
            except AttributeError:
                frame = 1   # should catch cases when we are analyzing a single PDB (?)
        coor = selection.coordinates()
        
        self.TITLE("FRAME "+str(frame)+" FROM "+str(u.trajectory.filename))
        for i, atom in enumerate(selection.atoms):
            self.ATOM(serial=i+1, name=atom.name.strip(), resName=atom.resname.strip(), resSeq=atom.resid,
                      x=coor[i,0], y=coor[i,1], z=coor[i,2])
        # get bfactor, too?
        self.close()


    def TITLE(self,*title):
        """Write TITLE record.
        http://www.wwpdb.org/documentation/format32/sect2.html
        """        
        line = " ".join(title)    # should do continuation automatically
        self.pdb.write(self.fmt['TITLE'] % line)

    def REMARK(self,*remark):
        """Write generic REMARK record (without number).
        http://www.wwpdb.org/documentation/format32/remarks1.html
        http://www.wwpdb.org/documentation/format32/remarks2.html
        """
        line = " ".join(remark)
        self.pdb.write(self.fmt['REMARK'] % line)
        
    def ATOM(self,serial=None,name=None,altLoc=None,resName=None,chainID=None,
             resSeq=None,iCode=None,x=None,y=None,z=None,occupancy=1.0,tempFactor=0.0,
             element=None,charge=0):
        """Write ATOM record. 
        http://www.wwpdb.org/documentation/format32/sect9.html
        Only some keword args are optional (altLoc, iCode, chainID), for some defaults are set.

        All inputs are cut to the maximum allowed length. For integer
        numbers the highest-value digits are chopped (so that the
        serial and reSeq wrap); for strings the trailing characters
        are chopped.

        Note: Floats are not checked and can potentially screw up the format.
        """
        for arg in ('serial','name','resName','resSeq','x','y','z',
                    'occupancy','tempFactor','charge'):
            if locals()[arg] is None:
                raise ValueError('parameter '+arg+' must be defined.')
        serial = int(str(serial)[-5:])  # check for overflow here?
        name = name[:4]
        if len(name) < 4:
            name = " "+name   # customary to start in column 14
        altLoc = altLoc or " "
        altLoc= altLoc[:1]
        resName = resName[:3]
        chainID = chainID or ""   # or should we provide a chainID such as 'A'?
        chainId = chainID[:1]
        resSeq = int(str(resSeq)[-4:]) # check for overflow here?
        iCode = iCode or ""
        iCode = iCode[:1]
        element = element or name.strip()[0]  # could have a proper dict here...
        element = element[:2]
        self.pdb.write(self.fmt['ATOM'] % vars())        
        
    def __del__(self):
        self.close()
