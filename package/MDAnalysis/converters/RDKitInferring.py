import warnings

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    pass
else:
    RDBONDORDER = {
        1: Chem.BondType.SINGLE,
        1.5: Chem.BondType.AROMATIC,
        "ar": Chem.BondType.AROMATIC,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
    }
    # add string version of the key for each bond
    RDBONDORDER.update({str(key): value for key, value in RDBONDORDER.items()})
    PERIODIC_TABLE = Chem.GetPeriodicTable()


# charges that should be assigned to monatomic cations
# structure --> atomic number : formal charge
# anion charges are directly handled by the code using the typical valence
# of the atom
MONATOMIC_CATION_CHARGES = {
    3: 1, 11: 1, 19: 1, 37: 1, 47: 1, 55: 1,
    12: 2, 20: 2, 29: 2, 30: 2, 38: 2, 56: 2,
    26: 2,  # Fe could also be +3
    13: 3,
}
# reactions uses by _standardize_patterns to fix challenging cases
# must have single reactant and product, and cannot add any atom
STANDARDIZATION_REACTIONS = [
    "[C-;X2;H0:1]=[O:2]>>[C+0:1]=[O:2]",  # Cterm
    "[N-;X2;H1;$(N-[*^3]):1]>>[N+0:1]",  # Nterm
    "[#6-:1]-[#6:2]=[O:3]>>[#6+0:1]=[#6:2]-[O-:3]",  # keto-enolate
    "[C-;v3:1]-[#7+0;v3;H2:2]>>[#6+0:1]=[#7+:2]",  # ARG
    "[#6+0;H0:1]=[#6+0:2]-[#7;X3:3]-[#6-;X3:4]"
    ">>[#6:1]=[#6:2]-[#7+:3]=[#6+0:4]",  # HIP
    "[S;D4;!v6:1]-[*-:2]>>[S;v6:1]=[*+0:2]",  # sulfone
    "[#7+0;X3:1]-[*-:2]>>[#7+:1]=[*+0:2]",  # charged-N
]


def _atom_sorter(atom):
    """Sorts atoms in the molecule in a way that makes it easy for the bond
    order and charge infering code to get the correct state on the first
    try. Currently sorts by number of unpaired electrons, then by number of
    heavy atom neighbors (i.e. atoms at the edge first)."""
    num_heavy_neighbors = len([
        neighbor for neighbor in atom.GetNeighbors()
        if neighbor.GetAtomicNum() > 1]
    )
    return (-_get_nb_unpaired_electrons(atom)[0], num_heavy_neighbors)


def _infer_bo_and_charges(mol):
    """Infer bond orders and formal charges from a molecule.

    Since most MD topology files don't explicitly retain information on bond
    orders or charges, it has to be guessed from the topology. This is done by
    looping over each atom and comparing its expected valence to the current
    valence to get the Number of Unpaired Electrons (NUE).
    If an atom has a negative NUE, it needs a positive formal charge (-NUE).
    If two neighbouring atoms have UEs, the bond between them most
    likely has to be increased by the value of the smallest NUE.
    If after this process, an atom still has UEs, it needs a negative formal
    charge of -NUE.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.RWMol
        The molecule is modified inplace and must have all hydrogens added

    Notes
    -----
    This algorithm is order dependant. For example, for a carboxylate group
    R-C(-O)-O the first oxygen read will receive a double bond and the other
    one will be charged. It will also affect more complex conjugated systems.
    """
    # heavy atoms sorted by number of heavy atom neighbors (lower first) then
    # NUE (higher first)
    atoms = sorted([atom for atom in mol.GetAtoms()
                    if atom.GetAtomicNum() > 1],
                   key=_atom_sorter)

    for atom in atoms:
        # monatomic ions
        if atom.GetDegree() == 0:
            atom.SetFormalCharge(MONATOMIC_CATION_CHARGES.get(
                                 atom.GetAtomicNum(),
                                 -_get_nb_unpaired_electrons(atom)[0]))
            mol.UpdatePropertyCache(strict=False)
            continue
        # get NUE for each possible valence
        nue = _get_nb_unpaired_electrons(atom)
        # if there's only one possible valence state and the corresponding
        # NUE is negative, it means we can only add a positive charge to
        # the atom
        if (len(nue) == 1) and (nue[0] < 0):
            atom.SetFormalCharge(-nue[0])
            mol.UpdatePropertyCache(strict=False)
        # go to next atom if above case or atom has no unpaired electron
        if (len(nue) == 1) and (nue[0] <= 0):
            continue
        else:
            neighbors = sorted(atom.GetNeighbors(), reverse=True,
                               key=lambda a: _get_nb_unpaired_electrons(a)[0])
            # check if one of the neighbors has a common NUE
            for na in neighbors:
                # get NUE for the neighbor
                na_nue = _get_nb_unpaired_electrons(na)
                # smallest common NUE
                common_nue = min(
                    min([i for i in nue if i >= 0], default=0),
                    min([i for i in na_nue if i >= 0], default=0)
                )
                # a common NUE of 0 means we don't need to do anything
                if common_nue != 0:
                    # increase bond order
                    bond = mol.GetBondBetweenAtoms(
                        atom.GetIdx(), na.GetIdx())
                    order = common_nue + 1
                    bond.SetBondType(RDBONDORDER[order])
                    mol.UpdatePropertyCache(strict=False)
                    # go to next atom if one of the valences is complete
                    nue = _get_nb_unpaired_electrons(atom)
                    if any([n == 0 for n in nue]):
                        break

            # if atom valence is still not filled
            nue = _get_nb_unpaired_electrons(atom)
            if not any([n == 0 for n in nue]):
                # transform nue to charge
                atom.SetFormalCharge(-nue[0])
                atom.SetNumRadicalElectrons(0)
                mol.UpdatePropertyCache(strict=False)


def _get_nb_unpaired_electrons(atom):
    """Calculate the number of unpaired electrons (NUE) of an atom

    Parameters
    ----------
    atom: rdkit.Chem.rdchem.Atom
        The atom for which the NUE will be computed

    Returns
    -------
    nue : list
        The NUE for each possible valence of the atom
    """
    expected_vs = PERIODIC_TABLE.GetValenceList(atom.GetAtomicNum())
    current_v = atom.GetTotalValence() - atom.GetFormalCharge()
    return [v - current_v for v in expected_vs]


def _standardize_patterns(mol, max_iter=200):
    """Standardizes functional groups

    Uses :func:`_rebuild_conjugated_bonds` to standardize conjugated systems,
    and SMARTS reactions for other functional groups.
    Due to the way reactions work, we first have to split the molecule by
    fragments. Then, for each fragment, we apply the standardization reactions.
    Finally, the fragments are recombined.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.RWMol
        The molecule to standardize
    max_iter : int
        Maximum number of iterations to standardize conjugated systems

    Returns
    -------
    mol : rdkit.Chem.rdchem.Mol
        The standardized molecule

    Notes
    -----
    The following functional groups are transformed in this order:

    +---------------+------------------------------------------------------------------------------+
    | Name          | Reaction                                                                     |
    +===============+==============================================================================+
    | conjugated    | ``[*-:1]-[*:2]=[*:3]-[*-:4]>>[*+0:1]=[*:2]-[*:3]=[*+0:4]``                   |
    +---------------+------------------------------------------------------------------------------+
    | conjugated N+ | ``[N;X3;v3:1]-[*:2]=[*:3]-[*-:4]>>[N+:1]=[*:2]-[*:3]=[*+0:4]``               |
    +---------------+------------------------------------------------------------------------------+
    | conjugated O- | ``[O:1]=[#6+0,#7+:2]-[*:3]=[*:4]-[*-:5]>>[O-:1]-[*:2]=[*:3]-[*:4]=[*+0:5]``  |
    +---------------+------------------------------------------------------------------------------+
    | conjug. S=O   | ``[O-:1]-[S;D4;v4:2]-[*:3]=[*:4]-[*-:5]>>[O+0:1]=[*:2]=[*:3]-[*:4]=[*+0:5]`` |
    +---------------+------------------------------------------------------------------------------+
    | Cterm         | ``[C-;X2;H0:1]=[O:2]>>[C+0:1]=[O:2]``                                        |
    +---------------+------------------------------------------------------------------------------+
    | Nterm         | ``[N-;X2;H1;$(N-[*^3]):1]>>[N+0:1]``                                         |
    +---------------+------------------------------------------------------------------------------+
    | keto-enolate  | ``[#6-:1]-[#6:2]=[O:3]>>[#6+0:1]=[#6:2]-[O-:3]``                             |
    +---------------+------------------------------------------------------------------------------+
    | arginine      | ``[C-;v3:1]-[#7+0;v3;H2:2]>>[#6+0:1]=[#7+:2]``                               |
    +---------------+------------------------------------------------------------------------------+
    | histidine     | ``[#6+0;H0:1]=[#6+0:2]-[#7;X3:3]-[#6-;X3:4]>>[#6:1]=[#6:2]-[#7+:3]=[#6+0:4]``|
    +---------------+------------------------------------------------------------------------------+
    | sulfone       | ``[S;D4;!v6:1]-[*-:2]>>[S;v6:1]=[*+0:2]``                                    |
    +---------------+------------------------------------------------------------------------------+
    | charged N     | ``[#7+0;X3:1]-[*-:2]>>[#7+:1]=[*+0:2]``                                      |
    +---------------+------------------------------------------------------------------------------+

    """

    # standardize conjugated systems
    _rebuild_conjugated_bonds(mol, max_iter)

    # list of sanitized reactions
    reactions = []
    for rxn in STANDARDIZATION_REACTIONS:
        reaction = AllChem.ReactionFromSmarts(rxn)
        reactions.append(reaction)

    # fragment mol (reactions must have single reactant and product)
    fragments = list(Chem.GetMolFrags(mol, asMols=True))
    for reactant in fragments:
        _apply_reactions(reactions, reactant)

    # recombine fragments
    newmol = fragments.pop(0)
    for fragment in fragments:
        newmol = Chem.CombineMols(newmol, fragment)

    return newmol


def _apply_reactions(reactions, reactant):
    """Applies a series of unimolecular reactions to a molecule. The reactions
    must not add any atom to the product. The molecule is modified inplace.

    Parameters
    ----------
    reactions : List[rdkit.Chem.rdChemReactions.ChemicalReaction]
        Reactions from SMARTS. Each reaction is applied until no more
        transformations can be made.
    reactant : rdkit.Chem.rdchem.Mol
        The molecule to transform inplace

    """
    reactant.UpdatePropertyCache(strict=False)
    Chem.Kekulize(reactant)
    for reaction in reactions:
        while reaction.RunReactantInPlace(reactant):
            reactant.UpdatePropertyCache(strict=False)
        reactant.UpdatePropertyCache(strict=False)
        Chem.Kekulize(reactant)


def _rebuild_conjugated_bonds(mol, max_iter=200):
    """Rebuild conjugated bonds without negatively charged atoms at the
    beginning and end of the conjugated system

    Depending on the order in which atoms are read during the conversion, the
    :func:`_infer_bo_and_charges` function might write conjugated systems with
    a double bond less and both edges of the system as anions instead of the
    usual alternating single and double bonds. This function corrects this
    behaviour by using an iterative procedure.
    The problematic molecules always follow the same pattern:
    ``anion[-*=*]n-anion`` instead of ``*=[*-*=]n*``, where ``n`` is the number
    of successive single and double bonds. The goal of the iterative procedure
    is to make ``n`` as small as possible by consecutively transforming
    ``anion-*=*`` into ``*=*-anion`` until it reaches the smallest pattern with
    ``n=1``. This last pattern is then transformed from ``anion-*=*-anion`` to
    ``*=*-*=*``.
    Since ``anion-*=*`` is the same as ``*=*-anion`` in terms of SMARTS, we can
    control that we don't transform the same triplet of atoms back and forth by
    adding their indices to a list.
    This function also handles conjugated systems with triple bonds.
    The molecule needs to be kekulized first to also cover systems
    with aromatic rings.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.RWMol
        The molecule to transform, modified inplace
    max_iter : int
        Maximum number of iterations performed by the function

    Notes
    -----
    The molecule is modified inplace
    """
    mol.UpdatePropertyCache(strict=False)
    Chem.Kekulize(mol)
    # pattern used to find problematic conjugated bonds
    # there's usually an even number of matches for this
    pattern = Chem.MolFromSmarts("[*-{1-2}]-,=[*+0]=,#[*+0]")
    # pattern used to finish fixing a series of conjugated bonds
    base_end_pattern = Chem.MolFromSmarts(
        "[*-{1-2}]-,=[*+0]=,#[*+0]-,=[*-{1-2}]")
    # used when there's an odd number of matches for `pattern`
    odd_end_pattern = Chem.MolFromSmarts(
        "[*-]-[*+0]=[*+0]-[*-,$([#7;X3;v3]),$([#6+0,#7+1]=O),"
        "$([S;D4;v4]-[O-])]")
    # number of unique matches with the pattern
    n_matches = len(set([match[0]
                         for match in mol.GetSubstructMatches(pattern)]))
    # nothing to standardize
    if n_matches == 0:
        return
    # single match (unusual)
    elif n_matches == 1:
        # as a last resort, the only way to standardize is to find a nitrogen
        # that can accept a double bond and a positive charge
        # or a carbonyl that will become an enolate
        end_pattern = odd_end_pattern
    # at least 2 matches
    else:
        # give priority to base end pattern and then deal with the odd one if
        # necessary
        end_pattern = base_end_pattern
    backtrack = []
    backtrack_cycles = 0
    for _ in range(max_iter):
        # check for ending pattern
        end_match = mol.GetSubstructMatch(end_pattern)
        if end_match:
            # index of each atom
            anion1, a1, a2, anion2 = end_match
            term_atom = mol.GetAtomWithIdx(anion2)
            # edge-case 1: C-[O-] or [N+]-[O-]
            # [*-]-*=*-[C,N+]=O --> *=*-*=[C,N+]-[O-]
            # transform the =O to -[O-]
            if (
                term_atom.GetAtomicNum() == 6
                and term_atom.GetFormalCharge() == 0
            ) or (
                term_atom.GetAtomicNum() == 7
                and term_atom.GetFormalCharge() == 1
            ):
                for neighbor in term_atom.GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(anion2, neighbor.GetIdx())
                    if (neighbor.GetAtomicNum() == 8 and
                            bond.GetBondTypeAsDouble() == 2):
                        bond.SetBondType(Chem.BondType.SINGLE)
                        neighbor.SetFormalCharge(-1)
                        break
            # edge-case 2: S=O
            # [*-]-*=*-[Sv4]-[O-] --> *=*-*=[Sv6]=O
            # transform -[O-] to =O
            elif (term_atom.GetAtomicNum() == 16 and
                  term_atom.GetFormalCharge() == 0):
                for neighbor in term_atom.GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(anion2, neighbor.GetIdx())
                    if (neighbor.GetAtomicNum() == 8 and
                        neighbor.GetFormalCharge() == -1 and
                            bond.GetBondTypeAsDouble() == 1):
                        bond.SetBondType(Chem.BondType.DOUBLE)
                        neighbor.SetFormalCharge(0)
                        break
            # not an edge case:
            # increment the charge
            else:
                term_atom.SetFormalCharge(term_atom.GetFormalCharge() + 1)
            # common to all cases:
            # [*-]-*=*-[R] --> *=*-*=[R]
            # increment the charge and switch single and double bonds
            a = mol.GetAtomWithIdx(anion1)
            a.SetFormalCharge(a.GetFormalCharge() + 1)
            b = mol.GetBondBetweenAtoms(anion1, a1)
            b.SetBondType(RDBONDORDER[b.GetBondTypeAsDouble() + 1])
            b = mol.GetBondBetweenAtoms(a1, a2)
            b.SetBondType(RDBONDORDER[b.GetBondTypeAsDouble() - 1])
            b = mol.GetBondBetweenAtoms(a2, anion2)
            b.SetBondType(RDBONDORDER[b.GetBondTypeAsDouble() + 1])
            mol.UpdatePropertyCache(strict=False)
            continue

        # switch the position of the charge and the double bond
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            # check if we haven't already transformed this triplet
            for match in matches:
                # store order-independent atom indices
                g = set(match)
                # already transformed --> try the next one
                if g in backtrack:
                    continue
                # add to backtracking and start the switch
                else:
                    anion, a1, a2 = match
                    backtrack.append(g)
                    break
            # already performed all changes
            else:
                if backtrack_cycles == 1:
                    # might be stuck because there were 2 "odd" end patterns
                    # misqualified as a single base one
                    end_pattern = odd_end_pattern
                elif backtrack_cycles > 1:
                    # likely stuck changing the same bonds over and over so
                    # let's stop here
                    mol.UpdatePropertyCache(strict=False)
                    return
                match = matches[0]
                anion, a1, a2 = match
                backtrack = [set(match)]
                backtrack_cycles += 1

            # switch charges
            a = mol.GetAtomWithIdx(anion)
            a.SetFormalCharge(a.GetFormalCharge() + 1)
            a = mol.GetAtomWithIdx(a2)
            a.SetFormalCharge(a.GetFormalCharge() - 1)
            # switch bond orders
            b = mol.GetBondBetweenAtoms(anion, a1)
            b.SetBondType(RDBONDORDER[b.GetBondTypeAsDouble() + 1])
            b = mol.GetBondBetweenAtoms(a1, a2)
            b.SetBondType(RDBONDORDER[b.GetBondTypeAsDouble() - 1])
            mol.UpdatePropertyCache(strict=False)
            # update number of matches for the end pattern
            n_matches = len(set([match[0] for match in matches]))
            if n_matches == 1:
                end_pattern = odd_end_pattern
            # start new iteration
            continue

        # no more changes to apply
        mol.UpdatePropertyCache(strict=False)
        return

    # reached max_iter
    warnings.warn("The standardization could not be completed within a "
                  "reasonable number of iterations")
