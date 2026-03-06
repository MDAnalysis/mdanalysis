from .exceptions import ChangeRadiiError, ParmWarning
import warnings

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def bondi(parm):
    """ Sets the bondi radii """
    for i, atom in enumerate(parm.atoms):
        # Radius of C atom depends on what type it is
        if atom.atomic_number == 6:
            if isinstance(atom.type, int):
                atom.solvent_radius = 1.7
            elif atom.type.startswith('C1') and atom.mass > 13.0:
                atom.solvent_radius = 2.2
            elif atom.type.startswith('C2') and atom.mass > 14.0:
                atom.solvent_radius = 2.2
            elif atom.type.startswith('C3') and atom.mass > 15.0:
                atom.solvent_radius = 2.2
            else:
                atom.solvent_radius = 1.7
        # All other elements have fixed radii for all types/partners
        elif atom.atomic_number == 1:
            atom.solvent_radius = 1.2
        elif atom.atomic_number == 7:
            atom.solvent_radius = 1.55
        elif atom.atomic_number == 8:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 9:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 14:
            atom.solvent_radius = 2.1
        elif atom.atomic_number == 15:
            atom.solvent_radius = 1.85
        elif atom.atomic_number == 16:
            atom.solvent_radius = 1.8
        elif atom.atomic_number == 17:
            atom.solvent_radius = 1.7
        else:
            atom.solvent_radius = 1.5

    try:
        parm.parm_data['RADIUS_SET'][0] = 'Bondi radii (bondi)'
    except AttributeError:
        pass
    _screen1(parm)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def amber6(parm):
    """ Sets the amber6 radii """
    for i, atom in enumerate(parm.atoms):
        # Radius of H atom depends on element it is bonded to
        bondeds = list(atom.bond_partners)
        if atom.atomic_number == 1:
            if not bondeds:
                warnings.warn(f"hydrogen atom {i} has no bonds; this is unusual", ParmWarning)
                atom.solvent_radius = 1.2
            elif bondeds[0].atomic_number == 6: # carbon
                atom.solvent_radius = 1.3
            elif bondeds[0].atomic_number in (8, 16): # oxygen or sulfur
                atom.solvent_radius = 0.8
            else: # anything else
                atom.solvent_radius = 1.2
        # Radius of C atom depends on what type it is
        elif atom.atomic_number == 6:
            if isinstance(atom.type, int):
                atom.solvent_radius = 1.7
            elif atom.type.startswith('C1') and atom.mass > 13.0:
                atom.solvent_radius = 2.2
            elif atom.type.startswith('C2') and atom.mass > 14.0:
                atom.solvent_radius = 2.2
            elif atom.type.startswith('C3') and atom.mass > 15.0:
                atom.solvent_radius = 2.2
            else:
                atom.solvent_radius = 1.7
        # All other elements have fixed solvent_radius for all types/partners
        elif atom.atomic_number == 7:
            atom.solvent_radius = 1.55
        elif atom.atomic_number == 8:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 9:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 14:
            atom.solvent_radius = 2.1
        elif atom.atomic_number == 15:
            atom.solvent_radius = 1.85
        elif atom.atomic_number == 16:
            atom.solvent_radius = 1.8
        elif atom.atomic_number == 17:
            atom.solvent_radius = 1.7
        else:
            atom.solvent_radius = 1.5
    try:
        parm.parm_data['RADIUS_SET'][0] = 'amber6 modified Bondi radii (amber6)'
    except AttributeError:
        pass
    _screen1(parm)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def mbondi(parm):
    """ Sets the mbondi radii """
    for i, atom in enumerate(parm.atoms):
        # Radius of H atom depends on element it is bonded to
        if atom.atomic_number == 1:
            bondeds = list(atom.bond_partners)
            if not bondeds:
                warnings.warn(f"hydrogen atom {i} has no bonds; this is unusual", ParmWarning)
                atom.solvent_radius = 0.8
            elif bondeds[0].atomic_number in (6, 7): # C or N
                atom.solvent_radius = 1.3
            elif bondeds[0].atomic_number in (8, 16): # O or S
                atom.solvent_radius = 0.8
            else:
                atom.solvent_radius = 1.2
        # Radius of C atom depends on what type it is
        elif atom.atomic_number == 6:
            if isinstance(atom.type, int):
                atom.solvent_radius = 1.7
            elif atom.type.startswith('C1') and atom.mass > 13.0:
                atom.solvent_radius = 2.2
            elif atom.type.startswith('C2') and atom.mass > 14.0:
                atom.solvent_radius = 2.2
            elif atom.type.startswith('C3') and atom.mass > 15.0:
                atom.solvent_radius = 2.2
            else:
                atom.solvent_radius = 1.7
        # All other elements have fixed solvent_radius for all types/partners
        elif atom.atomic_number == 7:
            atom.solvent_radius = 1.55
        elif atom.atomic_number == 8:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 9:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 14:
            atom.solvent_radius = 2.1
        elif atom.atomic_number == 15:
            atom.solvent_radius = 1.85
        elif atom.atomic_number == 16:
            atom.solvent_radius = 1.8
        elif atom.atomic_number == 17:
            atom.solvent_radius = 1.7
        else:
            atom.solvent_radius = 1.5
    try:
        parm.parm_data['RADIUS_SET'][0] = 'modified Bondi radii (mbondi)'
    except AttributeError:
        pass
    _screen1(parm)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def mbondi_pb2(parm):
    """
    Sets the mbondi_pb2 radii

    Use ropt values for halogens in Table 3 (https://pubs.acs.org/doi/full/10.1021/acs.jctc.9b00106)

                     pb2
    halogen	 rstd	 MAEstd	 MAEopt (ropt)
    Cl	     1.70    1.544	 1.529  (1.76)
    Br	     1.85    0.792	 0.743  (1.97)
    I	     1.98    0.905	 0.858  (2.09)

    These ropt values have been obtained using the "pb2" routine described
    in https://pubs.acs.org/doi/full/10.1021/acs.jctc.9b00106

    setup	radii set(s)	ΔGnonpolar	      γ	        β        inp value in gmx_MMPBSA
    pb2	    mbondi	        γ * SASA + β	  0.00720	0.0000   1

    """
    for i, atom in enumerate(parm.atoms):
        # Radius of H atom depends on element it is bonded to
        if atom.atomic_number == 1:
            bondeds = list(atom.bond_partners)
            if bondeds[0].atomic_number in (6, 7): # C or N
                atom.solvent_radius = 1.3
            elif bondeds[0].atomic_number in (8, 16): # O or S
                atom.solvent_radius = 0.8
            else:
                atom.solvent_radius = 1.2
        # Radius of C atom depends on what type it is
        elif atom.atomic_number == 6:
            if atom.type.startswith('C1') and atom.mass > 13.0:
                atom.solvent_radius = 2.2
            if atom.type.startswith('C2') and atom.mass > 14.0:
                atom.solvent_radius = 2.2
            if atom.type.startswith('C3') and atom.mass > 15.0:
                atom.solvent_radius = 2.2
            else:
                atom.solvent_radius = 1.7
        # All other elements have fixed solvent_radius for all types/partners
        elif atom.atomic_number == 7:
            atom.solvent_radius = 1.55
        elif atom.atomic_number == 8:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 9:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 14:
            atom.solvent_radius = 2.1
        elif atom.atomic_number == 15:
            atom.solvent_radius = 1.85
        elif atom.atomic_number == 16:
            atom.solvent_radius = 1.8
        # Cl
        elif atom.atomic_number == 17:
            atom.solvent_radius = 1.76
        # Br
        elif atom.atomic_number == 35:
            atom.solvent_radius = 1.97
        # I
        elif atom.atomic_number == 53:
            atom.solvent_radius = 2.09
        else:
            atom.solvent_radius = 1.5
    try:
        parm.parm_data['RADIUS_SET'][0] = 'modified Bondi radii with ropt for halogens (mbondi_pb2)'
    except AttributeError:
        pass
    _screen1(parm)


def mbondi_pb3(parm):
    """
    Sets the mbondi_pb3 radii

    Use ropt values for halogens in Table 3 (https://pubs.acs.org/doi/full/10.1021/acs.jctc.9b00106)

                     pb3
    halogen	 rstd	 MAEstd	 MAEopt (ropt)
    Cl	     1.70    1.358	 1.236 (2.20)
    Br	     1.85    0.770	 0.681 (2.04)
    I	     1.98    0.769	 0.697 (2.19)

    These ropt values have been obtained using the "pb3" routine described
    in https://pubs.acs.org/doi/full/10.1021/acs.jctc.9b00106

    setup	radii set(s)	ΔGnonpolar	                γ	     β        inp value in gmx_MMPBSA
    pb3	    mbondi + Rmin	γ * SAV + β + ΔGdispersion	0.03780	 –0.5692  2

    The mbondi radii set is used in the ΔGpolar calculation, while the value corresponding to Rmin in GAFF is
    used for the calculation of ΔGnonpolar

    """
    for i, atom in enumerate(parm.atoms):
        # Radius of H atom depends on element it is bonded to
        if atom.atomic_number == 1:
            bondeds = list(atom.bond_partners)
            if bondeds[0].atomic_number in (6, 7): # C or N
                atom.solvent_radius = 1.3
            elif bondeds[0].atomic_number in (8, 16): # O or S
                atom.solvent_radius = 0.8
            else:
                atom.solvent_radius = 1.2
        # Radius of C atom depends on what type it is
        elif atom.atomic_number == 6:
            if atom.type.startswith('C1') and atom.mass > 13.0:
                atom.solvent_radius = 2.2
            if atom.type.startswith('C2') and atom.mass > 14.0:
                atom.solvent_radius = 2.2
            if atom.type.startswith('C3') and atom.mass > 15.0:
                atom.solvent_radius = 2.2
            else:
                atom.solvent_radius = 1.7
        # All other elements have fixed solvent_radius for all types/partners
        elif atom.atomic_number == 7:
            atom.solvent_radius = 1.55
        elif atom.atomic_number == 8:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 9:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 14:
            atom.solvent_radius = 2.1
        elif atom.atomic_number == 15:
            atom.solvent_radius = 1.85
        elif atom.atomic_number == 16:
            atom.solvent_radius = 1.8
        # Cl
        elif atom.atomic_number == 17:
            atom.solvent_radius = 2.20
        # Br
        elif atom.atomic_number == 35:
            atom.solvent_radius = 2.04
        # I
        elif atom.atomic_number == 53:
            atom.solvent_radius = 2.19
        else:
            atom.solvent_radius = 1.5
    try:
        parm.parm_data['RADIUS_SET'][0] = 'modified Bondi radii with ropt for halogens (mbondi_pb3)'
    except AttributeError:
        pass
    _screen1(parm)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def mbondi2(parm):
    """ Sets the mbondi2 radii """
    for i, atom in enumerate(parm.atoms):
        # Radius of H atom depends on element it is bonded to
        if atom.atomic_number == 1:
            if not atom.bond_partners:
                warnings.warn(f"hydrogen atom {i} has no bonds; this is unusual", ParmWarning)
                atom.solvent_radius = 1.2
            if atom.bond_partners[0].atomic_number == 7:
                atom.solvent_radius = 1.3
            else:
                atom.solvent_radius = 1.2
        # Radius of C atom depends on what type it is
        elif atom.atomic_number == 6:
            if isinstance(atom.type, int):
                atom.solvent_radius = 1.7
            elif atom.type.startswith('C1') and atom.mass > 13.0:
                atom.solvent_radius = 2.2
            elif atom.type.startswith('C2') and atom.mass > 14.0:
                atom.solvent_radius = 2.2
            elif atom.type.startswith('C3') and atom.mass > 15.0:
                atom.solvent_radius = 2.2
            else:
                atom.solvent_radius = 1.7
        # All other elements have fixed solvent_radius for all types/partners
        elif atom.atomic_number == 7:
            atom.solvent_radius = 1.55
        elif atom.atomic_number == 8:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 9:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 14:
            atom.solvent_radius = 2.1
        elif atom.atomic_number == 15:
            atom.solvent_radius = 1.85
        elif atom.atomic_number == 16:
            atom.solvent_radius = 1.8
        elif atom.atomic_number == 17:
            atom.solvent_radius = 1.7
        else:
            atom.solvent_radius = 1.5
    try:
        parm.parm_data['RADIUS_SET'][0] = 'H(N)-modified Bondi radii (mbondi2)'
    except AttributeError:
        pass
    _screen1(parm)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def mbondi3(parm):
    """ Sets mbondi3 radii """
    mbondi2(parm) # start from mbondi2 radii
    for i, atom in enumerate(parm.atoms):
        # Adjust OE (GLU), OD (ASP), and HH/HE (ARG)
        if atom.residue.name in ('GLU', 'ASP', 'GL4', 'AS4'):
            if atom.name.startswith('OE') or atom.name.startswith('OD'):
                atom.solvent_radius = 1.4
        elif atom.residue.name == 'ARG':
            if atom.name.startswith('HH') or atom.name.startswith('HE'):
                atom.solvent_radius = 1.17
        # Adjust carboxylate O radii on C-Termini. Don't just do the end
        # residue, since we can have C-termini in the middle as well
        # (i.e., 2-chain dimers)
        if atom.name == 'OXT':
            atom.solvent_radius = 1.4
            parm.atoms[i-1].solvent_radius = 1.4

    try:
        parm.parm_data['RADIUS_SET'][0] = 'ArgH and AspGluO modified Bondi2 radii (mbondi3)'
    except AttributeError:
        pass
    _screen1(parm)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def charmm_radii(parm):
    """
    Sets charmm radii

    * Atomic radii for Poisson-Boltzmann calculations derived from average
    * solvent electrostatic charge distribution with explicit solvent.
    * Accuracy tested with free energy perturbation with explicit solvent

    * Most of the values were taken from a *_radii.str file used in PBEQ Solver in charmm-gui

    ! Radii for protein atoms in 20 standard amino acids
    ! Authors:  Mafalda Nina, Dmitrii Belogv, and Benoit Roux
    ! University of Montreal, June 1996.
    ! M. Nina and B. Roux. "Atomic Radii for Continuum Electrostatics
    ! Calculations Based on Molecular Dynamics Free Energy Simulations",
    ! J. Phys. Chem. B 101: 5239-5248 (1997).

    ! Radii for nucleic acid atoms (RNA and DNA)
    ! Authors:  Nilesh Banavali and Benoit Roux
    ! Cornell University, June 2002.
    ! N.K. Banavali and B. Roux, "Atomic Radii for Continuum Electrostatics
    ! Calculations on Nucleic Acids", J. Phys.  Chem. B 106:11026-11035 (2002).

    ! UPDATES:
    ! --------
    !
    ! December 1998:  GLU and ASP modified by Mafalda Nina.
    ! January  1999:  Protonated histidine HSP has been included by Mafalda Nina,
    !                 dG_elec = -68.15 kcal/mol (MD/FES) or -68.10 kcal/mol (PBEQ).
    ! January  1999:  TEA and ions added by Benoit Roux.
    ! November 2000:  Na+ added by Benoit Roux.
    ! June     2001:  All nucleic acid atoms added by Nilesh Banavali.
    ! February 2022:  Other atoms (https://pubs.acs.org/doi/10.1021/acs.jcim.1c00177)

    """
    for i, atom in enumerate(parm.atoms):
        # PROTEIN ATOMS
        if atom.residue.name in ['ALA', 'ARG', 'ASN', 'ASP', 'ASPP', 'CYS', 'GLN', 'GLU', 'GLUP', 'GLY', 'HSD',
                                 'HSE', 'HSP', 'ILE', 'LEU', 'LYS', 'LSN', 'MET', 'PHE', 'PRO', 'SER', 'SERD',
                                 'THR', 'TRP', 'TYR', 'VAL', 'ALAD', 'CYM', 'CYSD']:
            # !Patches CT3 N-Methylamide C-terminus
            # !        ACE acetylated N-terminus (ACP for PRO)
            # scalar wmain set 2.06 sele (type CAY .or. type CAT) end
            if atom.name == 'CAY' or atom.name == 'CAT':
                atom.solvent_radius = 2.06
            # scalar wmain set 2.04 sele type CY end
            if atom.name == 'CY':
                atom.solvent_radius = 2.04
            # scalar wmain set 1.52 sele type OY end
            if atom.name == 'OY':
                atom.solvent_radius = 1.52
            # scalar wmain set 2.23 sele type NT end
            if atom.name == 'NT':
                atom.solvent_radius = 2.23
            # ! for COO- terminus
            # scalar wmain set 1.40 sele type OT* end
            if atom.name.startswith('OT'):
                atom.solvent_radius = 1.40

            # !Backbone
            # scalar wmain set 2.04 sele type C  end
            #        ! for peptide bond
            if atom.name == 'C':
                atom.solvent_radius = 2.04
            # scalar wmain set 1.52 sele type O  end
            #        ! for peptide bond
            if atom.name == 'O':
                atom.solvent_radius = 1.52
            # scalar wmain set 2.23 sele type N  end
            #        ! for peptide bond
            if atom.name == 'N':
                atom.solvent_radius = 2.23
            # scalar wmain set 2.86 sele type CA  end
            #        ! for all CA except GLY
            # scalar wmain set 2.38 sele (resnam GLY .and. type CA) end
            # scalar wmain set 2.38 sele (resnam BGLY .and. type CA) end
            #        ! for GLY only
            if atom.name == 'CA':
                if atom.residue.name == 'GLY':
                    atom.solvent_radius = 2.38
                else:
                    atom.solvent_radius = 2.86
            # scalar wmain set 2.38 sele type CB2 end
            if atom.name == 'CB2':
                atom.solvent_radius = 2.38

            # !Hydrogens
            # scalar wmain set 0.00 sele type H* end
            #        ! for all hydrogens
            if atom.name.startswith('H'):
                atom.solvent_radius = 0.00

            # !Carbons
            # scalar wmain set 2.67 sele type CB end
            #        ! for all residues
            if atom.name == 'CB':
                atom.solvent_radius = 2.67
            # scalar wmain set 2.46 sele type CG* end
            #        ! for ARG, GLN, ILE, LYS, MET, PHE, THR, TRP, VAL,
            #              HSP, HSD, HSE
            if atom.residue.name in ['ARG', 'GLN', 'ILE', 'LYS', 'MET', 'PHE', 'THR', 'TRP', 'VAL',
                                     'HSP', 'HSD', 'HSE', 'LEU']:
                if atom.name.startswith('CG'):
                    atom.solvent_radius = 2.46
            # scalar wmain set 2.77 sele resnam GLU .and. type CG end
            # scalar wmain set 2.77 sele resnam BGLU .and. type CG end
            #        ! for GLU only
            if atom.residue.name == 'GLU':
                if atom.name == 'CG':
                    atom.solvent_radius = 2.77
            # scalar wmain set 2.44 sele type CD* end
            #        ! for ARG, ILE, LEU, LYS
            if atom.residue.name in ['ARG', 'ILE', 'LEU', 'LYS']:
                if atom.name.startswith('CD'):
                    atom.solvent_radius = 2.44
            # scalar wmain set 1.98 sele (resnam GLN .and. type CD) .or. (resnam ASN .and. type CG) .or. -
            #                            (resnam GLU .and. type CD) .or. (resnam ASP .and. type CG) end
            # scalar wmain set 1.98 sele (resnam BGLN .and. type CD) .or. (resnam BASN .and. type CG) .or. -
            #                            (resnam BGLU .and. type CD) .or. (resnam BASP .and. type CG) end
            #        ! for ASP, GLU, ASN, GLN
            if atom.residue.name in ['GLU', 'GLN']:
                if atom.name == 'CD':
                    atom.solvent_radius = 1.98
            if atom.residue.name in ['ASP', 'ASN']:
                if atom.name == 'CG':
                    atom.solvent_radius = 1.98
            # scalar wmain set 1.98 sele (resnam PRO .and. (type CB .or. type CG .or. type CD)) end
            # scalar wmain set 1.98 sele (resnam BPRO .and. (type CB .or. type CG .or. type CD)) end
            #        ! for PRO only
            if atom.residue.name == 'PRO':
                if atom.name in ['CB', 'CG', 'CD']:
                    atom.solvent_radius = 1.98
            # scalar wmain set 2.00 sele (resnam TYR .and. (type CE* .or. type CD* .or. -
            #                           type CZ)) .or. (resnam PHE .and. (type CE* .or. -
            #                           type CD* .or. type CZ))  end
            # scalar wmain set 2.00 sele (resnam BTYR .and. (type CE* .or. type CD* .or. -
            #                           type CZ)) .or. (resnam BPHE .and. (type CE* .or. -
            #                           type CD* .or. type CZ))  end
            #        ! for TYR, PHE rings
            if atom.residue.name in ['TYR', 'PHE']:
                if atom.name.startswith('CE') or atom.name.startswith('CD') or atom.name.startswith('CG') or \
                        atom.name == 'CZ':
                    atom.solvent_radius = 2.00
            # scalar wmain set 1.78 sele (resnam TRP .and. (type CE* .or. type CD* .or. -
            #                           type CZ* .or. type CH2)) end
            # scalar wmain set 1.78 sele (resnam BTRP .and. (type CE* .or. type CD* .or. -
            #                           type CZ* .or. type CH2)) end
            #        ! for TRP ring only
            if atom.residue.name == 'TRP':
                if atom.name.startswith('CE') or atom.name.startswith('CD') or atom.name.startswith('CZ') or \
                        atom.name == 'CH2':
                    atom.solvent_radius = 1.78
            # scalar wmain set 2.10 sele type CE end
            #        ! for MET only
            if atom.residue.name == 'MET':
                if atom.name == 'CE':
                    atom.solvent_radius = 2.10
            # scalar wmain set 2.80 sele (resnam ARG .and. type CZ) .or. (resnam LYS .and. type CE) end
            # scalar wmain set 2.80 sele (resnam BARG .and. type CZ) .or. (resnam BLYS .and. type CE) end
            # scalar wmain set 2.80 sele (resnam BORN .and. type CE) end
            #        ! for ARG, LYS
            if atom.residue.name in ['ARG', 'LYS', 'ORN']:
                if atom.name == 'CE' or atom.name == 'CZ':
                    atom.solvent_radius = 2.80
            # scalar wmain set 1.98 select (( resnam HSD .or. resnam HSE .or. resnam HSP) .and. type CE1) -
            #        .or. (( resnam HSD .or. resnam HSE .or. resnam HSP)  .and. type CD2) end
            # scalar wmain set 1.98 select (( resnam BHSD .or. resnam BHSE .or. resnam BHSP) .and. type CE1) -
            #        .or. (( resnam BHSD .or. resnam BHSE .or. resnam BHSP)  .and. type CD2) end
            #        ! for neutral HSD and protonated HSP
            if atom.residue.name in ['HSD', 'HSE', 'HSP']:
                if atom.name == 'CE1' or atom.name == 'CD2':
                    atom.solvent_radius = 1.98

            # !Oxygens
            # scalar wmain set 1.40 sele (resnam GLU .or. resnam ASP) .and. (type OE* .or. type OD*) end
            # scalar wmain set 1.40 sele (resnam BGLU .or. resnam BASP) .and. (type OE* .or. type OD*) end
            #        ! for GLU, ASP
            if atom.residue.name in ['GLU', 'ASP']:
                if atom.name.startswith('OE') or atom.name.startswith('OD'):
                    atom.solvent_radius = 1.40
            # scalar wmain set 1.42 sele (resnam ASN .or. resnam GLN) .and. (type OE* .or. type OD*) end
            # scalar wmain set 1.42 sele (resnam BASN .or. resnam BGLN) .and. (type OE* .or. type OD*) end
            #        ! for ASN, GLN
            if atom.residue.name in ['ASN', 'GLN']:
                if atom.name.startswith('OE') or atom.name.startswith('OD'):
                    atom.solvent_radius = 1.42
            # scalar wmain set 1.64 sele type OG* end
            #        ! for SER, THR
            if atom.residue.name in ['SER', 'THR']:
                if atom.name.startswith('OG'):
                    atom.solvent_radius = 1.64
            # scalar wmain set 1.85 sele (resnam TYR .and. type OH) end
            # scalar wmain set 1.85 sele (resnam BTYR .and. type OH) end
            #        ! for TYR only
            if atom.residue.name == 'TYR':
                if atom.name == 'OH':
                    atom.solvent_radius = 1.85
            # scalar wmain set 2.2 select resname TIP3 .and. type OH2 end
            #        ! for explicit water molecules
            if atom.residue.name == 'TIP3':
                if atom.name == 'OH2':
                    atom.solvent_radius = 2.2

            # !Nitrogens
            # scalar wmain set 1.80 sele resnam HSD  .and. (type NE2 .or. type ND1) end
            # scalar wmain set 1.80 sele resnam BHSD  .and. (type NE2 .or. type ND1) end
            #        ! for neutral HSD
            # scalar wmain set 1.80 sele resnam HSE  .and. (type NE2 .or. type ND1) end
            # scalar wmain set 1.80 sele resnam BHSE  .and. (type NE2 .or. type ND1) end
            #        ! for neutral HSE
            if atom.residue.name in ['HSD', 'HSE']:
                if atom.name == 'NE2' or atom.name == 'ND1':
                    atom.solvent_radius = 1.80
            # scalar wmain set 2.30 sele resnam HSP  .and. (type NE2 .or. type ND1) end
            # scalar wmain set 2.30 sele resnam BHSP  .and. (type NE2 .or. type ND1) end
            #        ! for protonated HSP
            if atom.residue.name == 'HSP':
                if atom.name == 'NE2' or atom.name == 'ND1':
                    atom.solvent_radius = 2.30
            # scalar wmain set 2.13 sele resnam ARG .and. (type NH* .or. type NE) .or. -
            #                       (resnam LYS .and. type NZ) end
            # scalar wmain set 2.13 sele resnam BARG .and. (type NH* .or. type NE) .or. -
            #                       (resnam BLYS .and. type NZ) end
            # scalar wmain set 2.13 sele (resnam BORN .and. type NZ) end
            #        ! for ARG, LYS
            if atom.residue.name in ['ARG', 'LYS', 'ORN']:
                if atom.name.startswith('NH') or atom.name == 'NE' or atom.name == 'NZ':
                    atom.solvent_radius = 2.13
            # scalar wmain set 2.15 sele (resnam GLN .and. type NE2) .or. (resnam ASN .and. type ND2) end
            # scalar wmain set 2.15 sele (resnam BGLN .and. type NE2) .or. (resnam BASN .and. type ND2) end
            #        ! for GLN, ASN
            if atom.residue.name in ['GLN', 'ASN']:
                if atom.name == 'NE2' or atom.name == 'ND2':
                    atom.solvent_radius = 2.15
            # scalar wmain set 2.40 sele resnam TRP .and. type NE1 end
            # scalar wmain set 2.40 sele resnam BTRP .and. type NE1 end
            #        ! for TRP
            if atom.residue.name == 'TRP':
                if atom.name == 'NE1':
                    atom.solvent_radius = 2.40

            # !Sulphur
            # scalar wmain set 2.00 sele type S* end
            #        ! for MET, CYS
            if atom.residue.name in ['MET', 'CYS']:
                if atom.name.startswith('S'):
                    atom.solvent_radius = 2.00

            # ! Phosphate atoms in phosphotyrosine, phosphoserine and phosphothreonine
            # calc radius = 2.35*@factor
            # scalar wmain set @radius sele (resn tyr .and. type p1) end
            # scalar wmain set @radius sele (resn ser .and. type p ) end
            if atom.residue.name in ['SER', 'TYR']:
                if atom.name == 'p1' or atom.name == 'p':
                    atom.solvent_radius = 2.35
            # calc radius = 1.49*@factor
            # scalar wmain set @radius sele (resn tyr .and. (type o2  .or. type o3  .or. type o4)) end
            # scalar wmain set @radius sele (resn ser .and. (type o1p .or. type o2p .or. type o3p .or. type ot )) end
            if atom.residue.name == 'TYR':
                if atom.name == 'o2' or atom.name == 'o3' or atom.name == 'o4':
                    atom.solvent_radius = 1.49
            if atom.residue.name in ['SER', 'THR']:
                if atom.name == 'o1p' or atom.name == 'o2p' or atom.name == 'o3p' or atom.name == 'ot':
                    atom.solvent_radius = 1.49

        # NUCLEIC ACID ATOMS
        elif atom.residue.name in ['ADE', 'GUA' , 'CYT', 'THY', 'URA', 'DMPA', 'ATP', 'ADP']:
            factor = 1
            # if @?factor eq 0 set factor = 1
            #
            # ! Set all H radii to 0.0
            # calc radius = 0.0*@factor
            # scalar wmain set @radius sele chem H* end
            if atom.name.startswith('H'):
                atom.solvent_radius = 0.00

            # ! Purine base atoms
            # calc radius = 1.75*@factor
            # scalar wmain set @radius sele resn ade .and. type n1 end
            if atom.residue.name == 'ADE':
                if atom.name == 'N1':
                    atom.solvent_radius = 1.75 * factor
            # calc radius = 2.17*@factor
            # scalar wmain set @radius sele resn ade .and. type n6 end
                elif atom.name == 'N6':
                    atom.solvent_radius = 2.17 * factor
            # calc radius = 2.15*@factor
            # scalar wmain set @radius sele resn gua .and. type n1 end
            if atom.residue.name == 'GUA':
                if atom.name == 'N1':
                    atom.solvent_radius = 2.15 * factor
            # calc radius = 1.55*@factor
            # scalar wmain set @radius sele resn gua .and. type o6 end
                elif atom.name == 'O6':
                    atom.solvent_radius = 1.55 * factor
            # calc radius = 2.12*@factor
            # scalar wmain set @radius sele resn gua .and. type n2 end
                elif atom.name == 'N2':
                    atom.solvent_radius = 2.12 * factor
            # calc radius = 2.15*@factor
            # scalar wmain set @radius sele (resn gua .or. resn ade) .and. type c2 end
            if atom.residue.name in ['GUA', 'ADE']:
                if atom.name == 'C2':
                    atom.solvent_radius = 2.15 * factor
            # calc radius = 1.69*@factor
            # scalar wmain set @radius sele (resn gua .or. resn ade) .and. type n3 end
                elif atom.name == 'N3':
                    atom.solvent_radius = 1.69 * factor
            # calc radius = 2.12*@factor
            # scalar wmain set @radius sele (resn gua .or. resn ade) .and. type c4 end
                elif atom.name == 'C4':
                    atom.solvent_radius = 2.12 * factor
            # calc radius = 2.12*@factor
            # scalar wmain set @radius sele (resn gua .or. resn ade) .and. type c5 end
                elif atom.name == 'C5':
                    atom.solvent_radius = 2.12 * factor
            # calc radius = 2.12*@factor
            # scalar wmain set @radius sele (resn gua .or. resn ade) .and. type c6 end
                elif atom.name == 'C6':
                    atom.solvent_radius = 2.12 * factor
            # calc radius = 1.69*@factor
            # scalar wmain set @radius sele (resn gua .or. resn ade) .and. type n7 end
                elif atom.name == 'N7':
                    atom.solvent_radius = 1.69 * factor
            # calc radius = 2.12*@factor
            # scalar wmain set @radius sele (resn gua .or. resn ade) .and. type c8 end
                elif atom.name == 'C8':
                    atom.solvent_radius = 2.12 * factor
            # calc radius = 2.13*@factor
            # scalar wmain set @radius sele (resn gua .or. resn ade) .and. type n9 end
                elif atom.name == 'N9':
                    atom.solvent_radius = 2.13 * factor
            # calc radius = 2.30*@factor
            # scalar wmain set @radius sele (resn gua .or. resn ade) .and. type c9 end
                elif atom.name == 'C9':
                    atom.solvent_radius = 2.30 * factor

            # ! Pyrimidine base atoms
            # calc radius = 2.20*@factor
            # scalar wmain set @radius sele (resn cyt .or. resn thy .or. resn ura) .and. type n1 end
            if atom.residue.name in ['CYT', 'THY', 'URA']:
                if atom.name == 'N1':
                    atom.solvent_radius = 2.20 * factor
            # calc radius = 2.04*@factor
            # scalar wmain set @radius sele (resn cyt .or. resn thy .or. resn ura) .and. type c2 end
                elif atom.name == 'C2':
                    atom.solvent_radius = 2.04 * factor
            # calc radius = 2.12*@factor
            # scalar wmain set @radius sele (resn cyt .or. resn thy .or. resn ura) .and. type c4 end
                elif atom.name == 'C4':
                    atom.solvent_radius = 2.12 * factor
            # calc radius = 2.25*@factor
            # scalar wmain set @radius sele (resn cyt .or. resn thy .or. resn ura) .and. type c5 end
                elif atom.name == 'C5':
                    atom.solvent_radius = 2.25 * factor
            # calc radius = 2.25*@factor
            # scalar wmain set @radius sele (resn cyt .or. resn thy .or. resn ura) .and. type c6 end
                elif atom.name == 'C6':
                    atom.solvent_radius = 2.25 * factor
            # calc radius = 1.60*@factor
            # scalar wmain set @radius sele (resn cyt .or. resn thy .or. resn ura) .and. type o2 end
                elif atom.name == 'O2':
                    atom.solvent_radius = 1.60 * factor
            # calc radius = 2.30*@factor
            # scalar wmain set @radius sele (resn cyt .or. resn thy .or. resn ura) .and. type c1 end
                elif atom.name == 'C1':
                    atom.solvent_radius = 2.30 * factor
            # calc radius = 1.68*@factor
            # scalar wmain set @radius sele resn cyt .and. type n3 end
            if atom.residue.name == 'CYT':
                if atom.name == 'N3':
                    atom.solvent_radius = 1.68 * factor
            # calc radius = 2.08*@factor
            # scalar wmain set @radius sele resn cyt .and. type n4 end
                elif atom.name == 'N4':
                    atom.solvent_radius = 2.08 * factor
            # calc radius = 2.20*@factor
            # scalar wmain set @radius sele (resn thy .or. resn ura) .and. type n3 end
            if atom.residue.name in ['THY', 'URA']:
                if atom.name == 'N3':
                    atom.solvent_radius = 2.20 * factor
            # calc radius = 1.60*@factor
            # scalar wmain set @radius sele (resn thy .or. resn ura) .and. type o4 end
                elif atom.name == 'O4':
                    atom.solvent_radius = 1.60 * factor
            # calc radius = 2.30*@factor
            # scalar wmain set @radius sele resn thy .and. type c5m end
            if atom.residue.name == 'THY':
                if atom.name == 'C5M':
                    atom.solvent_radius = 2.30 * factor

            # ! sugar atoms
            # calc radius = 2.57*@factor
            # scalar wmain set @radius sele type c1' end
            if atom.name == "C1'":
                atom.solvent_radius = 2.57 * factor
            # calc radius = 2.70*@factor
            # scalar wmain set @radius sele type c2' end
            if atom.name == "C2'":
                atom.solvent_radius = 2.70 * factor
            # calc radius = 2.73*@factor
            # scalar wmain set @radius sele type c3' end
            if atom.name == "C3'":
                atom.solvent_radius = 2.73 * factor
            # calc radius = 2.50*@factor
            # scalar wmain set @radius sele type c4' end
            if atom.name == "C4'":
                atom.solvent_radius = 2.50 * factor
            # calc radius = 1.55*@factor
            # scalar wmain set @radius sele type o4' end
            if atom.name == "O4'":
                atom.solvent_radius = 1.55 * factor
            # calc radius = 2.57*@factor
            # scalar wmain set @radius sele type c5' end
            if atom.name == "C5'":
                atom.solvent_radius = 2.57 * factor
            # calc radius = 1.75*@factor
            # scalar wmain set @radius sele type o2' end
            if atom.name == "O2'":
                atom.solvent_radius = 1.75 * factor

            if atom.name in ["O3'", "O5'"]:
                bondeds = list(atom.bond_partners)
                # ! add radii for blocking group hydroxyl oxygens
                # define oter  sele (type o3' .or. type o5') .and. .bonded. (type h3t .or. type h5t) end
                # calc radius = 1.72*@factor
                # scalar wmain set @radius sele oter end
                if "H3T" in [at.name for at in bondeds] or "H5T" in [at.name for at in bondeds]:
                    atom.solvent_radius = 1.72 * factor
                else:
                    # calc radius = 1.65*@factor
                    # scalar wmain set @radius sele type o3' end
                    # calc radius = 1.65*@factor
                    # scalar wmain set @radius sele type o5' end
                    atom.solvent_radius = 1.65 * factor

            # ! atoms for sugar2phos
            # calc radius = 1.65*@factor
            # scalar wmain set @radius sele resn ade .and. type o5t  end
            if atom.residue.name in ['ADE', 'GUA', 'CYT', 'THY', 'URA']:
                if atom.name == 'O5T':
                    atom.solvent_radius = 1.65 * factor
            # calc radius = 2.30*@factor
            # scalar wmain set @radius sele resn ade .and. type c5t  end
                elif atom.name == 'C5T':
                    atom.solvent_radius = 2.30 * factor
            # calc radius = 2.35*@factor
            # scalar wmain set @radius sele resn ade .and. type p3   end
                elif atom.name == 'P3':
                    atom.solvent_radius = 2.35 * factor
            # calc radius = 1.49*@factor
            # scalar wmain set @radius sele resn ade .and. type o1p3 end
                elif atom.name == 'O1P3':
                    atom.solvent_radius = 1.49 * factor
            # calc radius = 1.49*@factor
            # scalar wmain set @radius sele resn ade .and. type o2p3 end
                elif atom.name == 'O2P3':
                    atom.solvent_radius = 1.49 * factor
            # calc radius = 1.65*@factor
            # scalar wmain set @radius sele resn ade .and. type o3t  end
                elif atom.name == 'O3T':
                    atom.solvent_radius = 1.65 * factor
            # calc radius = 2.30*@factor
            # scalar wmain set @radius sele resn ade .and. type c3t  end
                elif atom.name == 'C3T':
                    atom.solvent_radius = 2.30 * factor

            # ! phosphate atoms
            # calc radius = 2.35*@factor
            # scalar wmain set @radius sele (resn gua .or. resn ade .or. resn cyt -
            #                     .or. resn thy .or. resn ura) .and. type p end
            if atom.residue.name in ['GUA', 'ADE', 'CYT', 'THY', 'URA']:
                if atom.name == 'P':
                    atom.solvent_radius = 2.35 * factor
            # calc radius = 1.49*@factor
            # scalar wmain set @radius sele (resn gua .or. resn ade .or. resn cyt -
            #                     .or. resn thy .or. resn ura) .and. type o1p end
                elif atom.name == 'O1P':
                    atom.solvent_radius = 1.49 * factor
            # calc radius = 1.49*@factor
            # scalar wmain set @radius sele (resn gua .or. resn ade .or. resn cyt -
            #                     .or. resn thy .or. resn ura) .and. type o2p end
                elif atom.name == 'O2P':
                    atom.solvent_radius = 1.49 * factor

            # ! phosphate atoms in ADP or ATP
            # calc radius = 2.35*@factor
            # scalar wmain set @radius sele type pa .or. type pb .or. type pg end
            if atom.name == 'PA' or atom.name == 'PB' or atom.name == 'PG':
                atom.solvent_radius = 2.35 * factor
            # calc radius = 1.65*@factor
            # scalar wmain set @radius sele (resn atp .and. ( type o3a .or. type o3b)) .or. -
            #                               (resn adp .and. ( type o3a)) end
            # calc radius = 1.49*@factor
            # scalar wmain set @radius sele type o1a .or. type o2a .or. type o1b .or. type o2b -
            #                  .or. type o1g .or. type o2g .or. type o3g .or. (resn adp .and. type o3b) end
            if atom.residue.name == 'ATP':
                if atom.name == 'O3A' or atom.name == 'O3B':
                    atom.solvent_radius = 1.65 * factor

            if atom.residue.name == 'ADP':
                if atom.name == 'O3A':
                    atom.solvent_radius = 1.65 * factor
                elif atom.name == 'O3B':
                    atom.solvent_radius = 1.49 * factor

            if atom.name in ['O1A', 'O2A', 'O1B', 'O2B', 'O1G', 'O2G', 'O3G']:
                atom.solvent_radius = 1.49 * factor

        # OTHER ATOMS
        else:
            # ! Set to zero all H radii
            # scalar wmain set 0.0 sele chem H* end
            if atom.name.startswith('H'):
                atom.solvent_radius = 0.00

            # !Carbons
            elif atom.name.startswith('C'):
                # from Table S5 (https://pubs.acs.org/doi/10.1021/acs.jcim.1c00177)
                if atom.type in [x[:4] for x in
                                 ['CG1N1', 'CG2D1', 'CG2D1O', 'CG2D2', 'CG2DC1', 'CG2O1', 'CG2O2',
                                  'CG2O4', 'CG2R51', 'CG2R53', 'CG2R61', 'CG2R62', 'CG2R63', 'CG2R67',
                                  'CG301', 'CG302', 'CG311', 'CG312', 'CG321', 'CG322', 'CG331',
                                  'CG3C50', 'CG3C51', 'CG3C52', 'CG3RC1']]:
                    atom.solvent_radius = 2.46
                # from Table S5 (https://pubs.acs.org/doi/10.1021/acs.jcim.1c00177)
                elif atom.type == 'CG2RC0'[:4]:
                    atom.solvent_radius = 2.12
                else:
                    # scalar wmain set 2.3 sele chem C* end
                    atom.solvent_radius = 2.30

            # !Oxygens
            elif atom.name.startswith('O'):
                # from Table S5 (https://pubs.acs.org/doi/10.1021/acs.jcim.1c00177)
                if atom.type in [x[:4] for x in ['OG2D1', 'OG2N1', 'OG301', 'OG302', 'OG303', 'OG311', 'OG312',
                                                 'OG3R60']]:
                    atom.solvent_radius = 1.64
                # from Table S5 (https://pubs.acs.org/doi/10.1021/acs.jcim.1c00177)
                elif atom.type == 'OG2D4'[:4]:
                    atom.solvent_radius = 1.60
                # from Table S5 (https://pubs.acs.org/doi/10.1021/acs.jcim.1c00177)
                elif atom.type == 'OG2P1'[:4]:
                    atom.solvent_radius = 1.49
                else:
                    # scalar wmain set 1.8 sele chem O* end
                    atom.solvent_radius = 1.8

            # !Nitrogens
            elif atom.name.startswith('N'):
                # from Table S5 (https://pubs.acs.org/doi/10.1021/acs.jcim.1c00177)
                if atom.type == 'NG2S1'[:4]:
                    atom.solvent_radius = 2.23
                else:
                    # scalar wmain set 2.3 sele chem N* end
                    atom.solvent_radius = 2.30

            # !Sulphur
            # scalar wmain set 2.3 sele chem S* end
            elif atom.name.startswith('S'):
                atom.solvent_radius = 2.30

            # !Phosphorus
            elif atom.name.startswith('P'):
                # from Table S5 (https://pubs.acs.org/doi/10.1021/acs.jcim.1c00177)
                atom.solvent_radius = 2.35

            # !Halogens rstd from Table S8 (https://pubs.acs.org/doi/10.1021/acs.jcim.1c00177)
            elif atom.name.startswith('CL'):
                atom.solvent_radius = 1.86
            elif atom.name.startswith('BR'):
                atom.solvent_radius = 1.98
            elif atom.name.startswith('I'):
                atom.solvent_radius = 2.24

            else:
                atom.solvent_radius = 2.00

    try:
        parm.parm_data['RADIUS_SET'][0] = \
                'charmm radii (charmm_radii)'
    except AttributeError:
        pass
    _screen1(parm)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _screen1(parm):
    """ Applies the first set of screening parameters found in tleap source """
    for i, atom in enumerate(parm.atoms):
        if atom.atomic_number == 1:
            atom.screen = 0.85
        elif atom.atomic_number == 6:
            atom.screen = 0.72
        elif atom.atomic_number == 7:
            atom.screen = 0.79
        elif atom.atomic_number == 8:
            atom.screen = 0.85
        elif atom.atomic_number == 9:
            atom.screen = 0.88
        elif atom.atomic_number == 15:
            atom.screen = 0.86
        elif atom.atomic_number == 16:
            atom.screen = 0.96
        else:
            atom.screen = 0.8

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _screen2(parm):
    """
    Applies the second set of screening parameters found in tleap source
    (unused as far as I can tell)
    """
    for i, atom in enumerate(parm.atoms):
        if atom.atomic_number == 1:
            atom.screen = 0.8461
        elif atom.atomic_number == 6:
            atom.screen = 0.9615
        elif atom.atomic_number == 7:
            atom.screen = 0.9343
        elif atom.atomic_number == 8:
            atom.screen = 1.0088
        elif atom.atomic_number == 11:
            atom.screen = 1.0000
        elif atom.atomic_number == 12:
            atom.screen = 1.0000
        elif atom.atomic_number == 15:
            atom.screen = 1.0700
        elif atom.atomic_number == 16:
            atom.screen = 1.1733
        else:
            atom.screen = 0.8

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _screen3(parm):
    """
    Applies the third and final set of screening parameters found in tleap
    source (unused as far as I can tell)
    """
    for i, atom in enumerate(parm.atoms):
        if atom.atomic_number == 1:
            atom.screen = 0.8846
        elif atom.atomic_number == 6:
            atom.screen = 0.9186
        elif atom.atomic_number == 7:
            atom.screen = 0.8733
        elif atom.atomic_number == 8:
            atom.screen = 0.8836
        elif atom.atomic_number == 11:
            atom.screen = 1.0000
        elif atom.atomic_number == 12:
            atom.screen = 1.0000
        elif atom.atomic_number == 15:
            atom.screen = 0.9604
        elif atom.atomic_number == 16:
            atom.screen = 0.9323
        else:
            atom.screen = 0.8

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

_call_method = dict(bondi=bondi, mbondi=mbondi, mbondi_pb2=mbondi_pb2, mbondi_pb3=mbondi_pb3, mbondi2=mbondi2,
                    mbondi3=mbondi3, amber6=amber6, charmm_radii=charmm_radii)

def ChRad(parm, radii_set):
    if radii_set not in _call_method:
        raise ChangeRadiiError(f"You must choose from {list(_call_method.keys())} radii sets")
    _call_method[radii_set](parm)
