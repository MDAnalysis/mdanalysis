"""
Contains the (very long) routine that does in-depth sanity checks on the
topology file
"""
from math import sqrt

from .exceptions import (AmberIncompatibleWarning, BadParmWarning,
        FixableParmWarning, NonfatalWarning, NonUniversalWarning,
        MissingDisulfide, LongBondWarning)
from ..constants import TINY
from ..amber.mask import AmberMask
from ..exceptions import MoleculeError

def check_validity(parm, warnings):

    # local shortcut for _check_exist_nvals
    def checkme(*args, **kwargs):
        kwargs['warnings'] = warnings
        return _check_exist_nvals(parm, *args, **kwargs)

    # Make sure all of our sections that we expect are present and have the
    # right number of variables in it
    nttyp = parm.ptr('ntypes') * (parm.ptr('ntypes') + 1) // 2 # for LJ stuff
    checkme('ATOM_NAME', parm.ptr('natom'), True, False, None, str)
    checkme('CHARGE', parm.ptr('natom'), True, False, None, float)
    checkme('MASS', parm.ptr('natom'), True, False, None, float)
    checkme('ATOM_TYPE_INDEX', parm.ptr('natom'), True, False, None, int)
    checkme('NUMBER_EXCLUDED_ATOMS', parm.ptr('natom'), True, False, None, int)
    checkme('NONBONDED_PARM_INDEX', parm.ptr('ntypes')**2, True, False, None, int)
    checkme('RESIDUE_LABEL', parm.ptr('nres'), True, False, None, str)
    checkme('RESIDUE_POINTER', parm.ptr('nres'), True, False, None, int)
    checkme('BOND_FORCE_CONSTANT', parm.ptr('numbnd'), True, False, None, float)
    checkme('BOND_EQUIL_VALUE', parm.ptr('numbnd'), True, False, None, float)
    checkme('ANGLE_FORCE_CONSTANT', parm.ptr('numang'), True, False, None, float)
    checkme('ANGLE_EQUIL_VALUE', parm.ptr('numang'), True, False, None, float)
    checkme('DIHEDRAL_FORCE_CONSTANT', parm.ptr('nptra'), True, False, None, float)
    checkme('DIHEDRAL_PERIODICITY', parm.ptr('nptra'), True, False, None, float)
    checkme('DIHEDRAL_PHASE', parm.ptr('nptra'), True, False, None, float)
    checkme('SCEE_SCALE_FACTOR', parm.ptr('nptra'), False, True, 'scee 1.2', float)
    checkme('SCNB_SCALE_FACTOR', parm.ptr('nptra'), False, True, 'scnb 2.0', float)
    checkme('SOLTY', parm.ptr('natyp'), True, False, None, float)
    checkme('LENNARD_JONES_ACOEF', nttyp, True, False, None, float)
    checkme('LENNARD_JONES_BCOEF', nttyp, True, False, None, float)
    checkme('BONDS_INC_HYDROGEN', parm.ptr('nbonh')*3, True, False, None, int)
    checkme('BONDS_WITHOUT_HYDROGEN', parm.ptr('nbona')*3, True, False, None, int)
    checkme('ANGLES_INC_HYDROGEN', parm.ptr('ntheth')*4, True, False, None, int)
    checkme('ANGLES_WITHOUT_HYDROGEN', parm.ptr('ntheta')*4, True, False, None, int)
    checkme('DIHEDRALS_INC_HYDROGEN', parm.ptr('nphih')*5, True, False, None, int)
    checkme('DIHEDRALS_WITHOUT_HYDROGEN', parm.ptr('nphia')*5, True, False, None, int)
    checkme('EXCLUDED_ATOMS_LIST', parm.ptr('next'), True, False, None, int)
    checkme('HBOND_ACOEF', parm.ptr('nphb'), True, False, None, float)
    checkme('HBOND_BCOEF', parm.ptr('nphb'), True, False, None, float)
    checkme('HBCUT', parm.ptr('nphb'), True, False, None, float)
    checkme('AMBER_ATOM_TYPE', parm.ptr('natom'), True, False, None, str)
    checkme('TREE_CHAIN_CLASSIFICATION', parm.ptr('natom'), True, False, None, str)
    checkme('JOIN_ARRAY', parm.ptr('natom'), True, False, None, int)
    checkme('IROTAT', parm.ptr('natom'), True, False, None, int)
    # PBC tests
    if parm.ptr('ifbox') > 0:
        if checkme('SOLVENT_POINTERS', 3, True, True, 'setMolecules', int):
            checkme('ATOMS_PER_MOLECULE', parm.parm_data['SOLVENT_POINTERS'][1],
                    True, True, 'setMolecules', int)
        checkme('BOX_DIMENSIONS', 4, True, False, None, float)
        # Check for contiguous molecules. Back up ATOMS_PER_MOLECULE and
        # SOLVENT_POINTERS so we can use rediscover_molecules without clobbering
        # any changes users may have made
        apm_cpy = parm.parm_data['ATOMS_PER_MOLECULE'][:]
        sp_cpy = parm.parm_data['SOLVENT_POINTERS'][:]
        try:
            parm.rediscover_molecules(fix_broken=False)
            parm.parm_data['ATOMS_PER_MOLECULE'] = apm_cpy
            parm.parm_data['SOLVENT_POINTERS'] = sp_cpy
        except MoleculeError:
            warnings.warn('ATOMS_PER_MOLECULE section corrupt! Molecules are not contiguous!')
    # Cap info
    if parm.ptr('ifcap') > 0:
        checkme('CAP_INFO', 1, True, False, None, int)
        checkme('CAP_INFO2', 4, True, False, None, float)
    # GB stuff
    if parm.ptr('ifbox') == 0:
        checkme('RADII', parm.ptr('natom'), True, True, 'changeRadii', float)
        checkme('SCREEN', parm.ptr('natom'), True, True, 'changeRadii', float)
        checkme('RADIUS_SET', 1, True, True, 'changeRadii', str)
        checkme('ATOMIC_NUMBER', parm.ptr('natom'), False, True, 'addAtomicNumber', int)
    # Some chamber checks
    if parm.chamber:
        # See if we have any CMAP terms
        has_cmap = 'CHARMM_CMAP_COUNT' in parm.flag_list
        # Ugh, don't want short-circuiting, so have to do this...
        try:
            nvals = parm.parm_data['FORCE_FIELD_TYPE'][0]
            nvals = int(nvals)
        except IndexError:
            warnings.warn('%FLAG FORCE_FIELD_TYPE is empty', BadParmWarning)
        except ValueError:
            warnings.warn('%FLAG FORCE_FIELD_TYPE does not have integer first', BadParmWarning)
        if len(parm.parm_data['FORCE_FIELD_TYPE']) - 1 < nvals:
            warnings.warn('Not enough lines in %FLAG FORCE_FIELD_TYPE', BadParmWarning)
        else:
            for i, val in enumerate(parm.parm_data['FORCE_FIELD_TYPE']):
                if i % 2 == 0 and not isinstance(val, int):
                    warnings.warn('int type mismatch in FORCE_FIELD_TYPE', BadParmWarning)
                elif i % 2 == 1 and not isinstance(val, str):
                    warnings.warn('str type mismatch in FORCE_FIELD_TYPE', BadParmWarning)
        hasallkeys = checkme('CHARMM_UREY_BRADLEY_COUNT', 2, False, False, None, int)
        hk = checkme('CHARMM_NUM_IMPROPERS', 1, True, False, None, int)
        hasallkeys = hasallkeys and hk
        hk = checkme('CHARMM_NUM_IMPR_TYPES', 1, True, False, None, int)
        hasallkeys = hasallkeys and hk
        if has_cmap:
            hk = checkme('CHARMM_CMAP_COUNT', 2, True, False, None, int)
            hasallkeys = hasallkeys and hk
     
        if hasallkeys:
            # Get the number of terms
            nub = parm.parm_data['CHARMM_UREY_BRADLEY_COUNT'][0]
            nubtypes = parm.parm_data['CHARMM_UREY_BRADLEY_COUNT'][1]
            nimphi = parm.parm_data['CHARMM_NUM_IMPROPERS'][0]
            nimprtyp = parm.parm_data['CHARMM_NUM_IMPR_TYPES'][0]
            cmap_term_cnt = parm.parm_data['CHARMM_CMAP_COUNT'][0]
            cmap_type_cnt = parm.parm_data['CHARMM_CMAP_COUNT'][1]
            # Look for terms
            checkme('LENNARD_JONES_14_ACOEF', nttyp, True, False, None, float)
            checkme('LENNARD_JONES_14_BCOEF', nttyp, True, False, None, float)
            checkme('CHARMM_UREY_BRADLEY', nub, True, False, None, int)
            checkme('CHARMM_UREY_BRADLEY_FORCE_CONSTANT', nubtypes, True, False, None, float)
            checkme('CHARMM_UREY_BRADLEY_EQUIL_CONSTANT', nubtypes, True, False, None, float)
            checkme('CHARMM_IMPROPERS', nimphi*5, True, False, None, int)
            checkme('CHARMM_IMPROPER_FORCE_CONSTANT', nimprtyp, True, False, None, float)
            checkme('CHARMM_IMPROPER_PHASE', nimprtyp, True, False, None, float)
            if has_cmap:
                checkme('CHARMM_CMAP_INDEX', cmap_term_cnt*6, True, False, None,int)
                if checkme('CHARMM_CMAP_RESOLUTION', cmap_type_cnt, True, False, True, int):
                    for i in range(cmap_type_cnt):
                        checkme('CHARMM_CMAP_PARAMETER_%s' % (str(i).zfill(2)),
                                parm.parm_data['CHARMM_CMAP_RESOLUTION'][i],
                                True, False, True, float)

    # Check that ATOMS_PER_MOLECULE == NATOM
    if 'ATOMS_PER_MOLECULE' in parm.parm_data:
        if sum(parm.parm_data['ATOMS_PER_MOLECULE']) != parm.ptr('natom'):
            warnings.warn('sum(ATOMS_PER_MOLECULE) != NATOM. Use "setMolecules" to fix this',
                          FixableParmWarning)

    # Duplicate pmemd's checks
    if parm.ptr('ifpert') != 0:
        warnings.warn('IFPERT must be 0! Parm will not work in Amber', AmberIncompatibleWarning)
    if (parm.ptr('mbona') != parm.ptr('nbona') or parm.ptr('mtheta') != parm.ptr('ntheta') or
                parm.ptr('mphia') != parm.ptr('nphia')):
        warnings.warn('Constraints can no longer be put in the prmtop in Amber programs!',
                      AmberIncompatibleWarning)

    # Check that we didn't change off-diagonal LJ pairs, since this will not
    # work for all programs necessarily...  Check relative error, though, since
    # Lennard Jones coefficients are very large
    parm.fill_LJ()
    ntypes = parm.ptr('ntypes')
    try:
        for i in range(ntypes):
            for j in range(ntypes):
                idx = parm.parm_data['NONBONDED_PARM_INDEX'][ntypes*i+j] - 1
                rij = parm.LJ_radius[i] + parm.LJ_radius[j]
                wdij = sqrt(parm.LJ_depth[i] * parm.LJ_depth[j])
                acoef = parm.parm_data['LENNARD_JONES_ACOEF'][idx]
                bcoef = parm.parm_data['LENNARD_JONES_BCOEF'][idx]
                if acoef == 0 or bcoef == 0:
                    if acoef != 0 or bcoef != 0 or (wdij != 0 and rij != 0):
                        warnings.warn('Modified off-diagonal LJ parameters detected',
                                      NonUniversalWarning)
                        raise StopIteration
                elif (abs((acoef - (wdij * rij**12)) / acoef) > TINY or
                      abs((bcoef - (2 * wdij * rij**6)) / bcoef) > TINY):
                    warnings.warn('Modified off-diagonal LJ parameters detected',
                                  NonUniversalWarning)
                    raise StopIteration
    except StopIteration:
        pass

    # Check to see if we forgot to add disulfides. Unless we have coordinates,
    # though, we can only check that CYX sulfurs are bonded to another Sulfur
    mask = AmberMask(parm, ':CYX@SG')
    for i, sel in enumerate(mask.Selection()):
        if not sel: continue
        atm = parm.atoms[i]
        # We expect 2 bonds
        bondedatms = [a.name for a in atm.bond_partners]
        if len(bondedatms) != 2 or 'SG' not in bondedatms:
            warnings.warn('Detected CYX residue with a Sulfur atom not bonded '
                          'to another CYX Sulfur! Did you forget to make the '
                          'disulfide bond?', MissingDisulfide)
            break

    if parm.coordinates is not None:
        # Check if we think any disulfide bonds might be missing.
        mask = AmberMask(parm, ':CYS,CYM@SG')
        s_atms = []
        for i, sel in enumerate(mask.Selection()):
            if not sel: continue
            s_atms.append(parm.atoms[i])
        try:
            for i in range(len(s_atms)-1):
                for j in range(i, len(s_atms)):
                    dx = s_atms[i].xx - s_atms[j].xx
                    dy = s_atms[i].xy - s_atms[j].xy
                    dz = s_atms[i].xz - s_atms[j].xz
                    if (dx * dx + dy * dy + dz * dz) < 9.0:
                        warnings.warn(
                            "Detected two cysteine residues whose sulfur atoms are within 3 "
                            "Angstroms. Rename CYS to CYX in the PDB file and use the 'bond' "
                            "command in tleap to create the disulfide bond",
                            MissingDisulfide,
                        )
                        raise StopIteration
        except StopIteration:
            pass
        # Now check if we have any bonds that seem unreasonably large. This
        # would indicate gaps in the structure.
        for bnd in parm.bonds:
            dx = bnd.atom1.xx - bnd.atom2.xx
            dy = bnd.atom1.xy - bnd.atom2.xy
            dz = bnd.atom1.xz - bnd.atom2.xz
            d2 = dx*dx + dy*dy + dz*dz
            req = bnd.type.req
            # Warn if any bond starts at > 3 times its equilibrium length
            if d2 > 9 * req*req:
                warnings.warn(
                    f'Atoms {bnd.atom1.idx + 1} ({bnd.atom1.residue.name} {bnd.atom1.residue.idx + 1} [{bnd.atom1.name}]) '
                    f'and {bnd.atom2.idx + 1} ({bnd.atom2.residue.name} {bnd.atom2.residue.idx + 1} [{bnd.atom2.name}]) '
                    f'are bonded (equilibrium length {req:.3f} A) but are {sqrt(d2):.3f} A apart. This '
                    'often indicates gaps in the original sequence and should be checked carefully.',
                    LongBondWarning
                )

def _check_exist_nvals(parm, key, nvals, required=True, addable=False,
                       addaction=None, typ=str, warnings=None):
    """
    Checks if a prmtop has a given key with the given number of values. We
    issue warnings based on whether it's a required flag and print advice
    if parmed can add it (i.e., addable)
    """
    if warnings is None:
        import warnings

    if addable and addaction is None:
        warnings.warn('Implementation Error: addable/addaction must be specified together!')
        addable = False

    if not key in parm.parm_data:
        if not required:
            if addable:
                warnings.warn(
                    f'%FLAG {key} not found, but it is not required. It can be added to the prmtop '
                    f'using the [[ {addaction} ]] action in ParmEd',
                    FixableParmWarning
                )
            else:
                warnings.warn(f'%FLAG {key} not found, but it is not required.', NonfatalWarning)
        else:
            if addable:
                warnings.warn(
                    f'%FLAG {key} not found! Use the [[ {addaction} ]] command in ParmEd to add it.',
                    FixableParmWarning,
                )
            else:
                warnings.warn(f'%FLAG {key} not found!', BadParmWarning)
        # Bail out here or we get a KeyError
        return

    if len(parm.parm_data[key]) != nvals:
        warnings.warn(f'%FLAG {key} has {len(parm.parm_data[key])} values, expected {nvals}!',
                      BadParmWarning)

    for val in parm.parm_data[key]:
        if not isinstance(val, typ):
            warnings.warn(
                f'Unexpected type for %FLAG {key}: expected {typ.__name__}, but got '
                f'{type(val).__name__}',
                BadParmWarning,
            )
            break
