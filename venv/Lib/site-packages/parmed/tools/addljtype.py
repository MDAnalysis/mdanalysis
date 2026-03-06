" This command adds a new Lennard Jones atom type from the selected atoms. "
from math import sqrt
from ..constants import PrmtopPointers

def AddLJType(parm, sel_atms, radius, epsilon, radius14, epsilon14):
    """ Adds a new Lennard Jones type to a topology file """
    # Set the new atom type
    old_type = None
    for i in range(len(sel_atms)):
        if sel_atms[i] == 1:
            if old_type is None:
                old_type = parm.parm_data['ATOM_TYPE_INDEX'][i]
            parm.parm_data['ATOM_TYPE_INDEX'][i] = parm.ptr('ntypes') + 1

    # Now increment NTYPES
    parm.parm_data['POINTERS'][PrmtopPointers.NTYPES] += 1
    parm.pointers['NTYPES'] += 1
   
    # Now create a whole new array for NONBONDED_PARM_INDEX
    start_idx = max(parm.parm_data['NONBONDED_PARM_INDEX']) + 1
    current_idx = max(parm.parm_data['NONBONDED_PARM_INDEX']) + 1
    old_ntypes = parm.ptr('ntypes') - 1
    for i in range(old_ntypes):
        # Copy over the first ntypes terms
        parm.parm_data['NONBONDED_PARM_INDEX'].insert(parm.ptr('ntypes')*(i+1)-1, current_idx)
        current_idx += 1

    # Now add the interaction of the last type with every other
    for i in range(old_ntypes):
        parm.parm_data['NONBONDED_PARM_INDEX'].append(start_idx)
        start_idx += 1

    # New type interacting with itself
    parm.parm_data['NONBONDED_PARM_INDEX'].append(current_idx)

    # Now we need to add onto the ACOEF and BCOEF arrays
    for i in range(old_ntypes):
        rad = parm.LJ_radius[i] + radius
        depth = sqrt(parm.LJ_depth[i] * epsilon)
        parm.parm_data['LENNARD_JONES_ACOEF'].append(depth * rad**12)
        parm.parm_data['LENNARD_JONES_BCOEF'].append(2 * depth * rad**6)
    # Do the same for the 1-4 interactions if we're doing a chamber prmtop
    if parm.chamber:
        for i in range(old_ntypes):
            rad = parm.LJ_14_radius[i] + radius14
            depth = sqrt(parm.LJ_14_depth[i] * epsilon14)
            parm.parm_data['LENNARD_JONES_14_ACOEF'].append(depth * rad**12)
            parm.parm_data['LENNARD_JONES_14_BCOEF'].append(2 * depth * rad**6)
   
    # Now add the last type interacting with itself, and add it to the
    # LJ_radius/depth arrays
    depth = epsilon
    rad = 2 * radius
    parm.parm_data['LENNARD_JONES_ACOEF'].append(depth * rad**12)
    parm.parm_data['LENNARD_JONES_BCOEF'].append(2 * depth * rad**6)
    parm.LJ_radius.append(radius)
    parm.LJ_depth.append(depth)

    # Do the same for chamber prmtops
    if parm.chamber:
        depth = epsilon14
        rad = 2 * radius14
        parm.parm_data['LENNARD_JONES_14_ACOEF'].append(depth * rad**12)
        parm.parm_data['LENNARD_JONES_14_BCOEF'].append(2 * depth * rad**6)
        parm.LJ_14_radius.append(radius14)
        parm.LJ_14_depth.append(epsilon14)

    # If our prmtop had a C-coefficient array, just copy the terms from the old
    # atom type to the new atom type, as defined by the type of the first atom
    # that got assigned the "new" type.
    if 'LENNARD_JONES_CCOEF' in parm.parm_data:
        ccoeffs = parm.parm_data['LENNARD_JONES_CCOEF']
        for i in range(old_ntypes):
            nbi = parm.ptr('ntypes')*(old_type-1) + i
            idx = parm.parm_data['NONBONDED_PARM_INDEX'][nbi] - 1
            ccoeffs.append(ccoeffs[idx])

        # Now the last type interacting with itself
        nbi = parm.ptr('ntypes')*(old_type-1) + old_type - 1
        idx = parm.parm_data['NONBONDED_PARM_INDEX'][nbi] - 1
        ccoeffs.append(ccoeffs[idx])
