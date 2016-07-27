"""MMTF Parser

More docs here
"""
import mmtf
from six.moves import zip, range

from . import base
from ..core.AtomGroup import Atom


class MMTFParser(base.TopologyReader):
    """
    Mapping of MMTF terms to MDAnalysis terms:
      models -> frames
      chains -> segments
      groups -> residues
      atoms  -> atoms
    """
    format = 'MMTF'

    def parse(self):
        top = mmtf.parse(self.filename)

        # TODO: What is the difference between these two?
        #       Which is closest to segid?
        #       Gets resolved in 363 anyway as we can use both fields
        chain_id = 0
        segid = top.chain_name_list[chain_id]
        #chain_id = top.chain_id_list[chain_id]
        chain_cum_sum = [sum(top.groups_per_chain[:i+1])
                         for i in range(top.num_chains)]

        atoms = []
        bonds = []
        bondorders = {}

        atom_offset = 0
        # Loop over groups in mmtf top object
        for i, group_id in enumerate(top.group_type_list):
            # Check if we've jumped into the next chain
            if i > chain_cum_sum[0]:
                # Move on to next tally
                chain_cum_sum.pop(0)
                chain_id += 1
                segid = top.chain_name_list[chain_id]

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
            sec_struc = top.sec_struct_list[group_id]
            sequence = top.sequence_index_list[group_id]
            single_letter_code = group['singleLetterCode']
            chem_comp_type = group['chemCompType']

            # Atom properties
            for j, (name, charge, element) in enumerate(zip(
                    group['atomNameList'],
                    group['formalChargeList'],
                    group['elementList'])):
                altLoc = top.alt_loc_list[atom_offset + j]
                occupancy = top.occupancy_list[atom_offset + j]
                bfactor = top.b_factor_list[atom_offset + j]

                # atom id as defined in mmtf
                mmtf_atom_id = top.atom_id_list[atom_offset + j]

                # TODO: Fix mass!
                mass = 1.0

                atoms.append(Atom(
                    atom_offset + j,
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
            # Intra group bonds
            bondlist = group['bondAtomList']
            for ibond, jbond, order in zip(
                    map(lambda x: x + atom_offset, bondlist[::2]),
                    map(lambda x: x + atom_offset, bondlist[1::2]),
                    group['bondOrderList'],
            ):
                bondtuple = tuple(sorted([ibond, jbond]))
                bonds.append(bondtuple)
                bondorders[bondtuple] = order

            # Jump offset by size of group
            atom_offset += len(group['atomNameList'])

        # Inter group bonds
        for ibond, jbond, order in zip(
            top.bond_atom_list[::2],
            top.bond_atom_list[1::2],
            top.bond_order_list
        ):
            bondtuple = tuple(sorted([ibond, jbond]))
            bonds.append(bondtuple)
            bondorders[bondtuple] = order

        struc = {
            'atoms': atoms,
            'bonds': bonds,
            'bondorder': bondorders,
        }

        return struc

        
