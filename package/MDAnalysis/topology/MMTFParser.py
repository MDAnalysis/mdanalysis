"""MMTF Parser

More docs here
"""
import mmtf
from six.moves import zip

from . import base
from ..core.AtomGroup import Atom

class MMTFParser(base.TopologyReader):
    format = 'MMTF'

    def parse(self):
        top = mmtf.parse(self.filename)

        atoms = []
        atom_id = 0
        # Loop over groups in mmtf top object
        for i, group_id in enumerate(top.group_type_list):
            group = top.group_list[group_id]

            resid = i + 1
            resname = group['groupName']

            for name, charge, element in zip(
                    group['atomNameList'],
                    group['formalChargeList'],
                    group['elementList']):
                altLoc = top.alt_loc_list[atom_id]
                occupancy = top.occupancy_list[atom_id]
                bfactor = top.b_factor_list[atom_id]

                # TODO: Fix mass!
                mass = 1.0
                # TODO: Fix segid!
                segid = 'FIXME'

                atoms.append(Atom(
                    atom_id,
                    name,
                    element,
                    resname,
                    resid,
                    segid,
                    mass,
                    charge,
                    radius=None,
                    bfactor=bfactor,
                    resnum=None,
                    serial=None,
                    altLoc=altLoc,
                    universe=self._u,
                ))
                atom_id += 1

        struc = {'atoms':atoms}

        # TODO: Parse bonds!

        return struc

        
