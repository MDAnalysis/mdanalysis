"""PDB structure files in MDAnalysis

Only coordinates are read; the topology must still be provided by a
psf file.

The PRB module makes heavy use of Biopython's Bio.PDB:

  Hamelryck, T., Manderick, B. (2003) PDB parser and structure class 
  implemented in Python. Bioinformatics, 19, 2308-2310.

  http://biopython.org
"""

        
# BioPython is overkill but potentially extensible (altLoc etc)
import Bio.PDB
import numpy
from DCD import Timestep

class PDBReader:
    """Read a pdb file into a BioPython pdb structure.

    The coordinates are also supplied as one numpy array and wrapped
    into a Timestep object; attributes are set so that the PDBReader
    object superficially resembles the DCDReader object.
    """
    def __init__(self,pdbfilename):
        p=Bio.PDB.PDBParser(PERMISSIVE=1)
        pdb_id = "0UNK"
        self.pdb = p.get_structure(pdb_id, pdbfilename)
        coord_list = [atom.coord for atom in self.pdb.get_atoms()]
        pos = numpy.array(coord_list)
        del coord_list
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

    Currently, this only works when the structure coordinates were
    loaded from a pdb because we are only modifying all the
    coordinates of this structure and writing it out again.
    """
    def __init__(self,PDBstructure,pdbfilename,multi=False,**kwargs):
        """pdbwriter = PDBWriter(universe,<pdbfilename>,**kwargs)
        :Arguments:
        PDBstructure    Bio.PDB.Structure.Structure (universe.pdb.pdb)
        pdbfilename     filename; if multi=True, embed a %%d formatstring
                        so that write_next_timestep() can insert the frame number
        multi           False: write a single structure to a single pdb
                        True: write all frames to multiple pdb files
        """
        import Bio.PDB.Structure
        self.PDBstructure = PDBstructure
        self.filename = pdbfilename
        self.multi = multi
        if not isinstance(PDBstructure,Bio.PDB.Structure.Structure):
            raise TypeError('PDBstructure must be a Bio.PDB.Structure.Structure, eg '
                            'Universe.pdb.pdb.')
    def write_next_timestep(self,ts):
        self.write(ts)
    def write(self,ts):
        # Let's cheat and use universe.pdb.pdb: modify coordinates
        # and save...
        for a,pos in zip(self.PDBstructure.get_atoms(), ts._pos):
            a.set_coord(pos)
        io = Bio.PDB.PDBIO()
        io.set_structure(self.PDBstructure)
        io.save(self.filename)
    
