# $Id: CRD.py 101 2009-08-11 13:19:06Z denniej0 $
"""CRD structure files in MDAnalysis

Only coordinates are written at the moment.

VERY Primative CRD generator (may still have to be debugged!)

 It may need some debugging (i.e.-it might not work for large systems
 as they usually need the extended version of crd writing).

"""

import MDAnalysis.core.util as util

class CRDWriter(object):
    """CRD writer that implements the standard CRD coordinate format.
    """

    fmt = {'ATOM':"%(serial)5d %(TotRes)4d %(resName)-4s %(name)-4s %(x)9.5f %(y)9.5f %(z)9.5f %(chainID)-4s %(resSeq)-4d %(tempFactor)9.5f\n",
           'TITLE':  "* %s\n",
	   'NUMATOMS':"%5d\n",
           }

    def __init__(self,filename):
        self.filename = util.filename(filename,ext='crd')
        self.crd = open(self.filename,'w')

    def close(self):
        self.crd.close()

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
       
        self.TITLE("  FRAME "+str(frame)+" FROM "+str(u.trajectory.filename))
        self.TITLE("")
	self.NUMATOMS(len(u.atoms))
	inst_resid = 0
	for i, atom in enumerate(selection.atoms):
	    #print selection[i].resname, selection[i-1].resname
	    if atom.resid < inst_resid:
		if selection[i].resname != selection[i-1].resname:
			totres = inst_resid + atom.resid
			inst_resid = totres 
	    else:
		inst_resid = atom.resid
		totres = inst_resid
            self.ATOM(serial=i+1, resSeq=atom.resid, resName=atom.resname, name=atom.name,
                      x=coor[i,0], y=coor[i,1], z=coor[i,2], chainID=atom.segid,tempFactor=0.0,TotRes=totres)
        # get bfactor, too?
        self.close()


    def TITLE(self,*title):
        """Write TITLE record.
        """        
        line = " ".join(title)    # should do continuation automatically
        self.crd.write(self.fmt['TITLE'] % line)

    def NUMATOMS(self,numatoms):
	"""Write genertic total number of atoms in system)
	"""
	self.crd.write(self.fmt['NUMATOMS'] % numatoms)
    
    def ATOM(self,serial=None,resSeq=None,resName=None,name=None,x=None,y=None,z=None,chainID=None,tempFactor=0.0,TotRes=None):
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
        resName = resName[:3]
        chainID = chainID or ""   # or should we provide a chainID such as 'A'?
        chainID = chainID[:4]
        resSeq = int(str(resSeq)[-4:]) # check for overflow here?
	totres = int(str(TotRes)[-4:])
	self.crd.write(self.fmt['ATOM'] % vars())
        
        
    def __del__(self):
        self.close()

