"""PDB structure files in MDAnalysis

(very thin...)
"""

        
# BioPython is overkill but potentially extensible (altLoc etc)
import Bio.PDB
import numpy
from DCD import Timestep

class PDBReader:
    """Read a pdb file into a BioPython pdb structure.

    The coordinates are also supplied as one numpy array and wrapped
    into a Timestep object; attributes are set so that the PDBReader
    obejct superficially resembles the DCDReader object.
    """
    def __init__(self,pdbfilename):
        p=Bio.PDB.PDBParser(PERMISSIVE=1)
        pdb_id = "0UNK"
        self.pdb = p.get_structure(pdb_id, pdbfilename)
        coord_list = [atom.coord for atom in self.pdb.get_atoms()]
        pos = numpy.array(coord_list)
        del coord_list
        self.pdbfilename = pdbfilename
        self.numatoms = pos.shape[0]
        self.numframes = 1
        self.fixed = 0          # parse B field for fixed atoms?
        self.skip = 1
        self.periodic = False
        self.ts = Timestep(pos)
        del pos
