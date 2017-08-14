import numpy as np
import MDAnalysis as mda
from MDAnalysis.topology.base import TopologyReaderBase
from MDAnalysis.topology import guessers
from MDAnalysis.lib.util import openany
from MDAnalysis.core.topology import Topology
from MDAnalysis.core.topologyattrs import (
    Atomnames,
    Atomids,
    Atomtypes,
    Bonds,
    Masses,
    Resids,
    Resnums,
    Segids,
)


class TXYZParser(TopologyReaderBase):
    """Parse a list of atoms from a Tinker XYZ file.
    Creates the following attributes:
     - Atomnames
     - Atomtypes
    .. versionadded:: 0.17.0
    """
    format = ['TXYZ', 'ARC']

    def parse(self):
        """Read the file and return the structure.
        Returns
        -------
        MDAnalysis Topology object
        """
        with openany(self.filename) as inf:
            #header
            natoms = int(inf.readline().split()[0])

            atomids = np.zeros(natoms, dtype=np.int)
            names = np.zeros(natoms, dtype=object)
            types = np.zeros(natoms, dtype=np.int)
            bonds = []
            # Can't infinitely read as XYZ files can be multiframe
            for i in range(natoms):
                line = inf.readline().split()
                atomids[i]= line[0]
                names[i] = line[1]
                types[i] = line[5]
                bonded_atoms = line[6:]
                for other_atom in bonded_atoms:
                    other_atom = int(other_atom) - 1
                    if i < other_atom: 
                         bonds.append((i, other_atom))
                    

        # Guessing time
        masses = guessers.guess_masses(names)

        attrs = [Atomnames(names),
                 Atomids(atomids),
                 Atomtypes(types),
                 Bonds(tuple(bonds)),
                 Masses(masses, guessed=True),
                 Resids(np.array([1])),
                 Resnums(np.array([1])),
                 Segids(np.array(['SYSTEM'], dtype=object)),
                 ]

        top = Topology(natoms, 1, 1,
                       attrs=attrs)

        return top

        
