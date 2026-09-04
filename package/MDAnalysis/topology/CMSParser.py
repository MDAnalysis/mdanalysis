import numpy as np
import re

from ..lib.util import openany
from ..core.topologyattrs import (
    Atomnames,
    Atomids,
    Resids,
    Resnames,
    Resnums,
    Segids,
    Resindices,
)
from ..core.topology import Topology
from .base import TopologyReaderBase, change_squash


class CMSParser(TopologyReaderBase):
    """Parser for CMS file format."""
    format = 'CMS'
    
    def __init__(self, filename):
        super().__init__(filename=filename)

    def parse(self, **kwargs):
        """Parse the CMS file and debug data extraction."""
        
        with open(self.filename, 'rt') as inf:
            # Read all lines
            lines = inf.readlines()

        # Extract the number of atoms (n_atoms) from the m_atom[] line
        n_atoms = 0
        for line in lines:
            line = line.strip()
            if line.startswith("m_atom["):
                # Extract the number of atoms from m_atom[] (e.g., m_atom[100])
                parts = line.split('[')
                n_atoms = int(parts[1].split(']')[0])
                break

        # If the number of atoms wasn't found, raise an error
        if n_atoms == 0:
            raise ValueError("Number of atoms (n_atoms) could not be found in the file.")
        
        # Create arrays
        resids = np.zeros(n_atoms, dtype=np.int32)
        resnames = np.zeros(n_atoms, dtype=object)
        segids = np.zeros(n_atoms, dtype=object)
        atomnames = np.zeros(n_atoms, dtype=object)
        atom_ids = np.zeros(n_atoms, dtype=np.int32)

        #Parse the atom data after the third occurrence of ":::"
        colon_count = 0  # Counter for ":::" markers
        atom_block = False

        atom_idx = 0

        # Regex pattern to split while keeping quoted strings together
        split_pattern = r'".*?"|\S+'

        for line in lines:
            line = line.strip()

            # Count occurrences of ":::" and start atom block after the third occurrence
            if ":::" in line:
                colon_count += 1
                if colon_count < 3:
                    continue
                else:
                    atom_block = True
                    continue  # Skip the ":::" marker line

            if atom_block:
                # Stop processing when encountering a line starting with '}'
                if line.startswith("}") or line.startswith('{'):
                    break

                details = re.findall(split_pattern, line)
                details = [item.strip('"') for item in details]

                if len(details) >= 7:
                    try:
                        # Extract atom data
                        atom_id = int(details[0])  # Atom ID
                        resid = int(details[5])  # Residue ID
                        segid = details[7].strip()  # Segment ID

                        resname = details[11].strip()  # Residue name
                        atomname = details[12].strip()  # Atom name

                        if atomname == '': 
                            atomname = details[14].strip()

                        # Fill the allocated arrays with parsed data
                        resids[atom_idx] = resid
                        resnames[atom_idx] = resname
                        segids[atom_idx] = segid
                        atomnames[atom_idx] = atomname
                        atom_ids[atom_idx] = atom_id

                        atom_idx += 1

                    except (ValueError, IndexError):
                        # Skip malformed lines
                        print(f"Skipping invalid line: {line}")

        print(segids)

        attrs = [
            Atomnames(atomnames),
            Atomids(atom_ids),
            Resids(resids),
            Resnames(resnames),
            Segids(segids),
        ]

        topology = Topology(n_atoms=n_atoms, n_res=len(resids), n_seg=len(segids), attrs=attrs)

        return topology
