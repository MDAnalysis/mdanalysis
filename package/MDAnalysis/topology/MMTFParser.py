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
            # Select the corresponding group
            group = top.group_list[group_id]
            # Grab per Group data
            resid = i + 1
            resnum = top.ins_code_list[group_id]
            # this field seems to be filled with '\x00' ie null
            if resnum == '\x00':
                resnum = None
            resname = group['groupName']

            # Unused group properties
            sec_struc = top.sec_struc_list[group_id]
            sequence = top.sequence_index_list[group_id]
            single_letter_code = group['singleLetterCode']
            chem_comp_type = group['chemCompType']
            

            for name, charge, element in zip(
                    group['atomNameList'],
                    group['formalChargeList'],
                    group['elementList']):
                altLoc = top.alt_loc_list[atom_id]
                occupancy = top.occupancy_list[atom_id]
                bfactor = top.b_factor_list[atom_id]

                # atom id as defined in mmtf
                mmtf_atom_id = top.atom_id_list[atom_id]

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
                    resnum=resnum,
                    serial=None,
                    altLoc=altLoc,
                    universe=self._u,
                ))
                atom_id += 1

        struc = {'atoms':atoms}

        # TODO: Parse bonds!

        return struc

        
