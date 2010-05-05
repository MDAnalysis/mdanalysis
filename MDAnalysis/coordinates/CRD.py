# $Id: CRD.py 101 2009-08-11 13:19:06Z denniej0 $
"""CRD structure files in MDAnalysis

Only coordinates are written at the moment.

VERY Primative CRD generator (may still have to be debugged!)

 It may need some debugging (i.e.-it might not work for large systems
 as they usually need the extended version of crd writing).

"""

import MDAnalysis.core.util as util
import base

class CRDWriter(base.Writer):
    """CRD writer that implements the standard CRD coordinate format.
    """
    format = 'CRD'
    units = {'time': None, 'length': 'Angstrom'}

    #          1         2         3         4         5         6         7         8
    # 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.
    # ATOM__seria nameAres CressI   xxxxxxxxyyyyyyyyzzzzzzzzOCCUPAtempft          elCH
    # %5d   %-4s %-3s %4d %1s %8.3f   %8.3f   %8.3f   %6.2f %6.2f           %2s
    #                 %1s  %1s                                                      %2d
    #            =        =      ===                                    ==========
    # %5d %4d %-4s %-4s %9.5f %9.5f %9.5f %-4s %-4d  %9.5f\n
    # %(serial)5d %(TotRes)4d %(resName)-4s %(name)-4s %(x)9.5f %(y)9.5f %(z)9.5f %(chainID)-4s %(resSeq)-4d %(tempFactor)9.5f
    fmt = {'ATOM':"%(serial)5d %(TotRes)4d %(resName)-4s %(name)-4s %(x)9.5f %(y)9.5f %(z)9.5f %(chainID)-4s %(resSeq)-4d %(tempFactor)9.5f\n",
           'TITLE':  "*%s\n",
	   'NUMATOMS':"%5d\n",
           }

    def __init__(self,filename):
        self.filename = util.filename(filename,ext='crd')
        self.crd = None

    def close_trajectory(self):
        pass

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
       
        self.crd = open(self.filename,'w')
        try:
            self._TITLE("FRAME "+str(frame)+" FROM "+str(u.trajectory.filename))
            self._TITLE("")
            self._NUMATOMS(len(u.atoms))
            inst_resid = 0
            for i, atom in enumerate(selection.atoms):
          	#print selection[i].resname, selection[i-1].resname
	        #print selection[i].resid, selection[i-1].resid
                if selection[i].resid != selection[i-1].resid:
                    inst_resid += 1
                    totres = inst_resid
                else:
                    inst_resid = inst_resid 
                    totres = inst_resid
                self._ATOM(serial=i+1, resSeq=atom.resid, resName=atom.resname, name=atom.name,
                          x=coor[i,0], y=coor[i,1], z=coor[i,2], chainID=atom.segid,tempFactor=0.0,TotRes=totres)
                # get bfactor, too?
        finally:
            self.crd.close()

    def _TITLE(self,*title):
        """Write TITLE record.
        """        
        line = " ".join(title)    # should do continuation automatically
        line = line.strip()
        if len(line) > 0:
            line = " "+line
        self.crd.write(self.fmt['TITLE'] % line)

    def _NUMATOMS(self,numatoms):
	"""Write generic total number of atoms in system)
	"""
	self.crd.write(self.fmt['NUMATOMS'] % numatoms)
    
    def _ATOM(self,serial=None,resSeq=None,resName=None,name=None,x=None,y=None,z=None,chainID=None,tempFactor=0.0,TotRes=None):
        """Write ATOM record. 

        All inputs are cut to the maximum allowed length. For integer
        numbers the highest-value digits are chopped (so that the
        serial and reSeq wrap); for strings the trailing characters
        are chopped.

        Note: Floats are not checked and can potentially screw up the format.
        """
        for arg in ('serial','name','resName','resSeq','x','y','z','tempFactor'):
            if locals()[arg] is None:
                raise ValueError('parameter '+arg+' must be defined.')
        serial = int(str(serial)[-5:])  # check for overflow here?
        name = name[:4]
        if len(name) < 4:
            name = name   # customary to start in column 14
        resName = resName[:4]
        chainID = chainID or ""   # or should we provide a chainID such as 'A'?
        chainID = chainID[:4]
        resSeq = int(str(resSeq)[-4:]) # check for overflow here?
	totres = int(str(TotRes)[-4:])
	self.crd.write(self.fmt['ATOM'] % vars())
        

