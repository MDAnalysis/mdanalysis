"""
This module contains all of the information for the titratable residues,
including reference energies, model compound pKas, and charge vectors for every
titratable residue treated.
"""
from math import log
from io import StringIO
from ..exceptions import AmberWarning, AmberError
import warnings

titratable_residues = ['AS4', 'GL4', 'CYS', 'TYR', 'HIP', 'LYS', 'DAP', 'DCP',
                       'DG', 'DT', 'AP', 'CP', 'G', 'U', 'HEH', 'PRN', 'TYX']

class _State(object):
    """ A protonation state """
    sort_by_resnum = True

    def __init__(self, charges, refene, refene_old=None, protcnt=None, pka_corr=None, eleccnt=None, eo_corr=None):
        self.charges = charges
        self.refene = refene
        self.refene_old = refene_old
        self.protcnt = protcnt
        self.pka_corr = pka_corr
        self.eleccnt = eleccnt
        self.eo_corr = eo_corr

class _ReferenceEnergy(object):
    """ Reference energies for various solvent models """
    LN_TO_LOG = log(10.0)
    KB = 0.00199
    TEMP = 300.0
    def __init__(self, igb1=None, igb2=None, igb5=None, igb7=None, igb8=None):
        self.pKa_is_set = False
        self.igb1 = igb1
        self.igb2 = igb2
        self.igb5 = igb5
        self.igb7 = igb7
        self.igb8 = igb8

    def solvent_energies(self, igb1=None, igb2=None, igb5=None, igb7=None, igb8=None):
        """
        Add solvent reference energies, copying the GB reference energies if
        none are explicitly given
        """
        self.solvent = _ReferenceEnergy(igb1, igb2, igb5, igb7, igb8)

    def dielc2_energies(self, igb1=None, igb2=None, igb5=None, igb7=None, igb8=None):
        """
        Add reference energies for a dielectric constant of 2.0. Since this has
        no reason to be near the original reference energies, do not use those
        in place of energies that aren't provided
        """
        self.dielc2 = _ReferenceEnergy(igb1, igb2, igb5, igb7, igb8)

    def set_pKa(self, pKa, deprotonated=False):
        """
        Adjusts the reference energies based on the pKa. If the reference energy
        is given for a deprotonated state, the pKa adjustment is subtracted from
        the reference energy. Otherwise, this state is a 'protonated' state, so
        the pKa adjustment is added to the reference energy
        """
        if self.pKa_is_set: return
        self.pKa_is_set = True
        factor = self.KB * self.LN_TO_LOG * self.TEMP * pKa

        if deprotonated:
            if self.igb1 is not None: self.igb1 -= factor
            if self.igb2 is not None: self.igb2 -= factor
            if self.igb5 is not None: self.igb5 -= factor
            if self.igb7 is not None: self.igb7 -= factor
            if self.igb8 is not None: self.igb8 -= factor
        else:
            if self.igb1 is not None: self.igb1 += factor
            if self.igb2 is not None: self.igb2 += factor
            if self.igb5 is not None: self.igb5 += factor
            if self.igb7 is not None: self.igb7 += factor
            if self.igb8 is not None: self.igb8 += factor
        # Now apply the correction to our solvent reference energies if we have
        # them (and dielectric-2 energies, if we have them)
        if hasattr(self, 'solvent'):
            self.solvent.set_pKa(pKa, deprotonated)
        if hasattr(self, 'dielc2'):
            self.dielc2.set_pKa(pKa, deprotonated)

class _LineBuffer(object):
    """ Buffer to add lines to the cpin file """

    CHARS_PER_LINE = 80

    def __init__(self, file):
        self.file = file
        self.linebuffer = ''

    def add_word(self, word):
        if len(self.linebuffer) + len(word) > self.CHARS_PER_LINE:
            self.file.write(f"{self.linebuffer}\n")
            self.linebuffer = f' {word}'
        else:
            self.linebuffer += word

    def add_words(self, words, space_delimited=False):
        """ Adds multiple words """
        extra = ''
        if space_delimited:
            extra = ' '
        for word in words:
            self.add_word(word + extra)

    def flush(self):
        """ Flushes this buffer to the file """
        if len(self.linebuffer) == 0:
            return
        self.file.write(self.linebuffer + '\n')
        self.linebuffer = ''

class TitratableResidue(object):
    """
    A residue with different protonation states defined for Amber for use in
    the Constant pH MD method implemented in sander
    """

    def __init__(self, resname, atom_list, typ, pka=None, eo=None):
        self.resname = resname
        self.atom_list = list(atom_list) # list of atom names
        self.states = []
        self.first_state = -1
        self.first_charge = -1
        self.typ = typ
        self.pKa = pka
        self.Eo = eo

    def _str_refenes(self, solvent=False, igb=2, dielc=1.0):
        """
        Converts all reference energies into a formatted string with a message
        saying if the energy is not set
        """
        ret_str = ''
        igb_str = f'igb{igb}'
        for state in self.states:
            if dielc == 2:
                if solvent:
                    refene = getattr(state.refene.dielc2.solvent, igb_str)
                else:
                    refene = getattr(state.refene.dielc2, igb_str)
            else:
                if solvent:
                    refene = getattr(state.refene.solvent, igb_str)
                else:
                    refene = getattr(state.refene, igb_str)

            if refene is None:
                ret_str += f"{'Not Set':>12s}"
            else:
                ret_str += f"{refene:12.5f}"

        return ret_str

    def __str__(self):
        delineator = '-' * (8 + 12 * len(self.states))
        if self.typ == "ph":
            ret_strs = [f"{self.resname:<4s}\tpKa = {self.pKa:5.1f}"]
        elif self.typ == "redox":
            ret_strs = [f"{self.resname:<4s}\tEo = {self.Eo:7.3f} V"]
        else:
            ret_strs = [f"{self.resname:<4s}"]
        ret_strs.append(
            f"{'ATOM':>8s}" + ''.join([f"STATE {i}".rjust(12) for i in range(len(self.states))])
        )
        for i, atom in enumerate(self.atom_list):
            ret_strs.append(f"{atom:>8s}" + "".join([f"{state.charges[i]:12.4f}" for state in self.states]))
        ret_strs.append(delineator)
        if self.typ == "ph" or self.typ == "phredox":
            ret_strs.append("Prot Cnt".rjust(8) + "".join([f"{state.protcnt:12d}" for state in self.states]))
            ret_strs.append(delineator)
            ret_strs.append("pKa Corr".rjust(8) + "".join([f"{state.pka_corr:12.4f}" for state in self.states]))
        if (self.typ == "phredox"):
            ret_strs.append(delineator)
        if (self.typ == "redox" or self.typ == "phredox"):
            ret_strs.append("Elec Cnt".rjust(8) + "".join([f"{state.eleccnt:12d}" for state in self.states]))
            ret_strs.append(delineator)
            ret_strs.append("Eo Corr".rjust(8) + "".join([f"{state.eo_corr:12.4f}" for state in self.states]))
        ret_strs.append(delineator)
        ret_strs.extend(["Reference Energies (ES = Explicit solvent, IS = Implicit solvent)", ""])
        ret_strs.append(f"{'igb=1 IS':8s}{self._str_refenes(False, 1)}")
        ret_strs.append(f"{'igb=2 IS':8s}{self._str_refenes(False, 2)}")
        ret_strs.append(f"{'igb=5 IS':8s}{self._str_refenes(False, 5)}")
        ret_strs.append(f"{'igb=7 IS':8s}{self._str_refenes(False, 7)}")
        ret_strs.append(f"{'igb=8 IS':8s}{self._str_refenes(False, 8)}")
        ret_strs.append(f"{'igb=1 ES':8s}{self._str_refenes(True, 1)}")
        ret_strs.append(f"{'igb=2 ES':8s}{self._str_refenes(True, 2)}")
        ret_strs.append(f"{'igb=5 ES':8s}{self._str_refenes(True, 5)}")
        ret_strs.append(f"{'igb=7 ES':8s}{self._str_refenes(True, 7)}")
        ret_strs.append(f"{'igb=8 ES':8s}{self._str_refenes(True, 8)}")
        ret_strs.append(delineator)
        ret_strs.extend(["Reference Energies for Internal Dielectric of 2.0", ""])
        ret_strs.append(f"{'igb=1 IS':8s}{self._str_refenes(False, 1, 2)}")
        ret_strs.append(f"{'igb=2 IS':8s}{self._str_refenes(False, 2, 2)}")
        ret_strs.append(f"{'igb=5 IS':8s}{self._str_refenes(False, 5, 2)}")
        ret_strs.append(f"{'igb=7 IS':8s}{self._str_refenes(False, 7, 2)}")
        ret_strs.append(f"{'igb=8 IS':8s}{self._str_refenes(False, 8, 2)}")
        ret_strs.append(f"{'igb=1 ES':8s}{self._str_refenes(True, 1, 2)}")
        ret_strs.append(f"{'igb=2 ES':8s}{self._str_refenes(True, 2, 2)}")
        ret_strs.append(f"{'igb=5 ES':8s}{self._str_refenes(True, 5, 2)}")
        ret_strs.append(f"{'igb=7 ES':8s}{self._str_refenes(True, 7, 2)}")
        ret_strs.append(f"{'igb=8 ES':8s}{self._str_refenes(True, 8, 2)}")
        return "\n".join(ret_strs) + "\n"

    def add_state(self, charges, refene, refene_old=None, protcnt=None, pka_corr=None, eleccnt=None, eo_corr=None):
        """ Add a single titratable state for this titratable residue """
        new_state = _State(charges, refene, refene_old=refene_old, protcnt=protcnt, pka_corr=pka_corr, eleccnt=eleccnt, eo_corr=eo_corr)
        if len(new_state.charges) != len(self.atom_list):
            raise AmberError('Wrong number of charges for new state')
        self.states.append(new_state)

    def add_states(self, charges, refenes, refenes_old=None, protcnts=None, pka_corrs=None, eleccnts=None, eo_corrs=None):
        """ Add multiple titratable states for this titratable residue """
        if len(charges) != len(refenes) or (len(charges) != len(refenes_old) and refenes_old) or (len(charges) != len(protcnts) and protcnts) \
           or (len(charges) != len(pka_corrs) and pka_corrs) or (len(charges) != len(eleccnts) and eleccnts) \
           or (len(charges) != len(eo_corrs) and eo_corrs):
            raise AmberError('Inconsistent list of parameters for TitratableResidue.add_states')
        for i in range(len(charges)):
            self.add_state(charges[i], refenes[i], refenes_old[i], protcnts[i], pka_corrs[i], eleccnts[i], eo_corrs[i])

    def cpin_pointers(self, first_atom):
        """ Sets and returns the cpin info """
        if self.first_state == -1 or self.first_charge == -1:
            raise AmberError('Must set residue pointers before writing cpin info!')
        return {'FIRST_ATOM': first_atom,
                'FIRST_CHARGE': self.first_charge,
                'FIRST_STATE': self.first_state,
                'NUM_ATOMS': len(self.atom_list),
                'NUM_STATES': len(self.states)}

    def set_first_state(self, index):
        """ Sets the first state index """
        # Has the first state already been set?
        if self.first_state != -1:
            if index != self.first_state:
                raise AmberError('First state already set differently')
        self.first_state = index

    def set_first_charge(self, index):
        """ Sets the first charge index """
        # Has it already been set?
        if self.first_charge != -1:
            if index != self.first_charge:
                raise AmberError('First charge already set differently')
        self.first_charge = index

    def reset(self):
        """ Resets the pointers """
        self.first_state = -1
        self.first_charge = -1

    def check(self):
        """ Checks that the charges are consistent w/ the protonation states """
        sum_charges = [sum(state.charges) for state in self.states]
        if (self.typ == "ph" or self.typ == "phredox"):
            protcnts = [state.protcnt for state in self.states]
        if (self.typ == "redox" or self.typ == "phredox"):
            eleccnts = [state.eleccnt for state in self.states]
        # All we have to do is make sure that the charges/proton counts are
        # consistent between the first state and every other state
        for i in range(1, len(sum_charges)):
            charge_diff = sum_charges[i] - sum_charges[0]
            if (self.typ == "ph"):
                diff = protcnts[i] - protcnts[0]
            elif (self.typ == "redox"):
                diff = eleccnts[0] - eleccnts[i]
            elif (self.typ == "phredox"):
                diff = protcnts[i] - protcnts[0] + eleccnts[0] - eleccnts[i]
            if abs(charge_diff - diff) >= 0.0001:
                warnings.warn(f'Charge definitions inconsistent in {self.resname}', AmberWarning)
        # Check all of the reference energies to make sure that the pKa was set
        # for all but one of them
        notset = 0
        valid = True
        for state in self.states:
            if not state.refene_old:
                valid = False
                break
            notset += int(state.refene_old.pKa_is_set)
        if (notset != len(self.states) - 1 and valid):
            warnings.warn(f"Not enough states are pKa-adjusted in {self.resname}")

class TitratableResidueList(list):
    """ List of all titratable residues """
    def __init__(self, system_name='Unknown', solvated=False, first_solvent=0):
        list.__init__(self)
        self.first_atoms = []
        self.residue_nums = []
        self.resstates = []
        self.system_name = system_name
        self.solvated = solvated
        self.first_sol = first_solvent

    def add_residue(self, residue, resnum, first_atom, state=0):
        """ Adds a residue to the list """
        list.append(self, residue)
        self.first_atoms.append(first_atom)
        self.residue_nums.append(resnum)
        if state < 0 or state >= len(residue.states):
            raise AmberError(
                f"Residue {residue.resname} only has states 0-{len(residue.states)} ({state} chosen)"
            )
        self.resstates.append(state)

    def set_states(self, statelist):
        """
        Sets the initial protonation states from a list -- make sure there are
        enough states in the list to set every residue, or emit a warning
        """
        if len(statelist) != len(self):
            warnings.warn(
                f"Number of states ({len(statelist)}) does not equal number of "
                f"residues ({len(self)}). Using default initial states.", AmberWarning
            )
            return
        # Check that all states are allowable
        for i, state in enumerate(statelist):
            if state < 0 or state >= len(self[i].states):
                raise AmberError(f"Bad state choice ({state}). Minimum is 0, maximum is {len(self[i].states)}")
        # If we got here, then we are OK
        self.resstates = statelist

    def sort(self):
        """ Sorts by residue number """
        atoms_residues = sorted(list(zip(self.first_atoms, self.residue_nums)))
        self.first_atoms = [atom_residue[0] for atom_residue in atoms_residues]
        self.residue_nums = [atom_residue[1] for atom_residue in atoms_residues]

    def write_cpin(self, output, igb=2, intdiel=1.0, oldfmt=False, typ="ph", coions=False):
        """ Writes the CPIN file based on the titrated residues """
        end = StringIO()
        # Reset all residues
        for res in self: res.reset()
        # Sort our residue list
        self.sort()
        limit_buf = _LineBuffer(output)
        buf = _LineBuffer(end)
        limit_buf.add_word('&CNSTPHE_LIMITS')
        limit_buf.flush()
        limit_buf.add_word(f' ntres={len(self)},')
        limit_buf.add_word(f' maxh={max(len(r.states) for r in self)},')
        if (typ == "ph"):
            buf.add_word('&CNSTPH')
        elif (typ == "redox"):
            buf.add_word('&CNSTE')
        elif (typ == "phredox"):
            buf.add_word('&CNSTPHE')
        buf.flush()
        buf.add_word(' CHRGDAT=')
        charges, energies, protcnts, pka_corrs, eleccnts, eo_corrs, pointers = [], [], [], [], [], [], []
        first_charge = 0
        first_state = 0
        for i, res in enumerate(self):
            if res.first_charge == -1:
                res.set_first_charge(first_charge)
                res.set_first_state(first_state)

                for state in res.states:
                    # See which dielectric reference energies we want
                    if intdiel == 2:
                        if (not oldfmt):
                            refene = state.refene.dielc2
                        else:
                            refene = state.refene_old.dielc2
                    else:
                        if (not oldfmt):
                            refene = state.refene
                        else:
                            refene = state.refene_old
                    # See if we want the explicit solvent refene or not
                    if self.solvated:
                        energies.append(getattr(refene.solvent, f"igb{igb}"))
                    else:
                        energies.append(getattr(refene, f"igb{igb}"))
                    if (typ == "ph" or typ == "phredox"):
                        # Add protonation count of this state
                        if (state.protcnt):
                            protcnts.append(state.protcnt)
                        else:
                            protcnts.append(0)
                        # Add pka reference of this state
                        if (state.pka_corr):
                            pka_corrs.append(state.pka_corr)
                        else:
                            pka_corrs.append(0.0)
                    if (typ == "redox" or typ == "phredox"):
                        # Add electron count of this state
                        if (state.eleccnt):
                            eleccnts.append(state.eleccnt)
                        else:
                            eleccnts.append(0)
                        # Add Eo reference of this state
                        if (state.eo_corr):
                            eo_corrs.append(state.eo_corr)
                        else:
                            eo_corrs.append(0.0)

                first_state += len(res.states)
                new_charges = []
                for state in res.states:
                    new_charges.extend(state.charges)
                charges.extend(new_charges)
                first_charge += len(new_charges)
            pointers.append(res.cpin_pointers(self.first_atoms[i]))

        limit_buf.add_word(f' natchrg={len(charges)},')
        limit_buf.add_word(f' ntstates={max(len(protcnts), len(eleccnts))},')
        limit_buf.flush()
        limit_buf.add_word('/')
        limit_buf.flush()
        # Print the charges
        for charge in charges:
            buf.add_word(f'{charge},')
        buf.flush()
        # Print the protcnts
        if (typ == "ph" or typ == "phredox"):
            buf.add_word(' PROTCNT=')
            for protcnt in protcnts:
                buf.add_word(f'{protcnt:d},')
        buf.flush()
        if (typ == "redox" or typ == "phredox"):
            buf.add_word(' ELECCNT=')
            for eleccnt in eleccnts:
                buf.add_word('%d,' % eleccnt)
        buf.flush()
        # Print the residue names
        buf.add_word(f" RESNAME='System: {self.system_name}',")
        for i, res in enumerate(self):
            buf.add_word(f"'Residue: {res.resname} {self.residue_nums[i]:d}',")
        buf.flush()
        # Print the residue states
        buf.add_word(" RESSTATE=")
        for state in self.resstates:
            buf.add_word(f'{state:d},')
        buf.flush()
        # Print the residue pointers
        buf.add_word(' ') # get a leading space
        for i, p in enumerate(pointers):
            buf.add_word(f"STATEINF({i})%FIRST_ATOM={p['FIRST_ATOM']}, ")
            buf.add_word(f"STATEINF({i})%FIRST_CHARGE={p['FIRST_CHARGE']}, ")
            buf.add_word(f"STATEINF({i})%FIRST_STATE={p['FIRST_STATE']}, ")
            buf.add_word(f"STATEINF({i})%NUM_ATOMS={p['NUM_ATOMS']}, ")
            buf.add_word(f"STATEINF({i})%NUM_STATES={p['NUM_STATES']}, ")
        buf.flush()
        # Print the reference energies
        buf.add_word(' STATENE=')
        for i, energy in enumerate(energies):
            if energy is None:
                raise AmberError(f"{i}'th reference energy not known for igb = {igb}")
            buf.add_word(f'{energy:.6f},')
        buf.flush()
        # Print the pKa or Eo reference
        if not oldfmt:
            if typ == "ph" or typ == "phredox":
                buf.add_word(' PKA_CORR=')
                for pka_corr in pka_corrs:
                    buf.add_word(f'{pka_corr:.4f},')
            buf.flush()
            if (typ == "redox" or typ == "phredox"):
                buf.add_word(' EO_CORR=')
                for eo_corr in eo_corrs:
                    buf.add_word(f'{eo_corr:.4f},')
            buf.flush()
        # Print the # of residues and explicit solvent info if required
        buf.add_word(' TRESCNT=%d,' % len(self))
        if self.solvated:
            if (typ == "ph"):
                buf.add_word(f'CPHFIRST_SOL={self.first_sol}, CPH_IGB={igb}, CPH_INTDIEL={intdiel}, ')
            elif (typ == "redox"):
                buf.add_word(f'CEFIRST_SOL={self.first_sol}, CE_IGB={igb}, CE_INTDIEL={intdiel}, ')
            elif (typ == "phredox"):
                buf.add_word(f'CPHEFIRST_SOL={self.first_sol}, CPHE_IGB={igb}, CPHE_INTDIEL={intdiel}, ')
            buf.flush()
            # Now scan through all of the waters
        buf.flush()
        buf.add_word('/')
        buf.flush()
        end.seek(0)
        output.write(end.read())

# Now define all of the titratable residues

# Aspartate
refene1 = _ReferenceEnergy(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.dielc2_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.dielc2.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene2 = _ReferenceEnergy(igb1=21.4298008, igb2=26.8894581, igb5=26.5980488,
                 igb7=23.4181107, igb8=26.3448911)
refene2.solvent_energies(igb2=33.2613028, igb5=32.064349, igb7=28.262350, igb8=31.286037)
refene2.dielc2_energies(igb2=12.676908, igb5=13.084913)
refene2.dielc2.solvent_energies()
# Copying the reference energy to be printted on the old CPIN format
refene2_old = _ReferenceEnergy(igb1=21.4298008, igb2=26.8894581, igb5=26.5980488,
                 igb7=23.4181107, igb8=26.3448911)
refene2_old.solvent_energies(igb2=33.2613028, igb5=32.064349, igb7=28.262350, igb8=31.286037)
refene2_old.dielc2_energies(igb2=12.676908, igb5=13.084913)
refene2_old.dielc2.solvent_energies()
refene2_old.set_pKa(4.0, deprotonated=False)

AS4 = TitratableResidue('AS4', ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG',
                        'OD1', 'OD2', 'HD21', 'C', 'O', 'HD22', 'HD11', 'HD12'],
                        pka=4.0, typ="ph")
AS4.add_state(protcnt=0, refene=refene1, refene_old=refene1, pka_corr=0.0, # deprotonated
              charges=[-0.4157, 0.2719, 0.0341, 0.0864, -0.1783, -0.0122,
              -0.0122, 0.7994, -0.8014, -0.8014, 0.0, 0.5973, -0.5679, 0.0, 0.0,
              0.0])
AS4.add_state(protcnt=1, refene=refene2, refene_old=refene2_old, pka_corr=4.0, # protonated syn-O2
              charges=[-0.4157, 0.2719, 0.0341, 0.0864, -0.0316, 0.0488, 0.0488,
              0.6462, -0.5554, -0.6376, 0.4747, 0.5973, -0.5679, 0.0, 0.0, 0.0])
AS4.add_state(protcnt=1, refene=refene2, refene_old=refene2_old, pka_corr=4.0, # protonated anti-O2
              charges=[-0.4157, 0.2719, 0.0341, 0.0864, -0.0316, 0.0488, 0.0488,
              0.6462, -0.5554, -0.6376, 0.0, 0.5973, -0.5679, 0.4747, 0.0, 0.0])
AS4.add_state(protcnt=1, refene=refene2, refene_old=refene2_old, pka_corr=4.0, # protonated syn-O1
              charges=[-0.4157, 0.2719, 0.0341, 0.0864, -0.0316, 0.0488, 0.0488,
              0.6462, -0.6376, -0.5554, 0.0, 0.5973, -0.5679, 0.0, 0.4747, 0.0])
AS4.add_state(protcnt=1, refene=refene2, refene_old=refene2_old, pka_corr=4.0, # protonated anti-O1
              charges=[-0.4157, 0.2719, 0.0341, 0.0864, -0.0316, 0.0488, 0.0488,
              0.6462, -0.6376, -0.5554, 0.0, 0.5973, -0.5679, 0.0, 0.0, 0.4747])
AS4.check() # check that everything is consistent

# Glutamate
refene1 = _ReferenceEnergy(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.dielc2_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.dielc2.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene2 = _ReferenceEnergy(igb1=3.89691326, igb2=8.4057785, igb5=8.0855764,
               igb7=5.305949, igb8=8.3591335)
refene2.solvent_energies(igb2=15.20019319, igb5=13.533409, igb7=10.682953, igb8=13.242410)
refene2.dielc2_energies(igb2=3.455596, igb5=3.957270)
refene2.dielc2.solvent_energies()
# Copying the reference energy to be printted on the old CPIN format
refene2_old = _ReferenceEnergy(igb1=3.89691326, igb2=8.4057785, igb5=8.0855764,
               igb7=5.305949, igb8=8.3591335)
refene2_old.solvent_energies(igb2=15.20019319, igb5=13.533409, igb7=10.682953, igb8=13.242410)
refene2_old.dielc2_energies(igb2=3.455596, igb5=3.957270)
refene2_old.dielc2.solvent_energies()
refene2_old.set_pKa(4.4, deprotonated=False)

GL4 = TitratableResidue('GL4', ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG',
                        'HG2', 'HG3', 'CD', 'OE1', 'OE2', 'HE21', 'C', 'O',
                        'HE22', 'HE11', 'HE12'], pka=4.4, typ="ph")
GL4.add_state(protcnt=0, refene=refene1, refene_old=refene1, pka_corr=0.0, # deprotonated
              charges=[-0.4157, 0.2719, 0.0145, 0.0779, -0.0398, -0.0173,
              -0.0173, 0.0136, -0.0425, -0.0425, 0.8054, -0.8188, -0.8188, 0.0,
              0.5973, -0.5679, 0.0, 0.0, 0.0])
GL4.add_state(protcnt=1, refene=refene2, refene_old=refene2_old, pka_corr=4.4, # protonated syn-O2
              charges=[-0.4157, 0.2719, 0.0145, 0.0779, -0.0071, 0.0256, 0.0256,
              -0.0174, 0.0430, 0.0430, 0.6801, -0.5838, -0.6511, 0.4641, 0.5973,
              -0.5679, 0.0, 0.0, 0.0])
GL4.add_state(protcnt=1, refene=refene2, refene_old=refene2_old, pka_corr=4.4, # protonated anti-O2
              charges=[-0.4157, 0.2719, 0.0145, 0.0779, -0.0071, 0.0256, 0.0256,
              -0.0174, 0.0430, 0.0430, 0.6801, -0.5838, -0.6511, 0.0, 0.5973,
              -0.5679, 0.4641, 0.0, 0.0])
GL4.add_state(protcnt=1, refene=refene2, refene_old=refene2_old, pka_corr=4.4, # protonated syn-O1
              charges=[-0.4157, 0.2719, 0.0145, 0.0779, -0.0071, 0.0256, 0.0256,
              -0.0174, 0.0430, 0.0430, 0.6801, -0.6511, -0.5838, 0.0, 0.5973,
              -0.5679, 0.0, 0.4641, 0.0])
GL4.add_state(protcnt=1, refene=refene2, refene_old=refene2_old, pka_corr=4.4, # protonated syn-O2
              charges=[-0.4157, 0.2719, 0.0145, 0.0779, -0.0071, 0.0256, 0.0256,
              -0.0174, 0.0430, 0.0430, 0.6801, -0.6511, -0.5838, 0.0, 0.5973,
              -0.5679, 0.0, 0.0, 0.4641])
GL4.check()

# Tyrosine
refene1 = _ReferenceEnergy(igb2=0, igb5=0, igb8=0)
refene1.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.dielc2_energies(igb2=0, igb5=0, igb8=0)
refene1.dielc2.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene2 = _ReferenceEnergy(igb2=-65.113428, igb5=-64.166385, igb8=-61.3305355)
refene2.solvent_energies(igb2=-65.003415, igb5=-64.047229)
refene2.dielc2_energies(igb2=-32.167520, igb5=-31.751177)
refene2.dielc2.solvent_energies()
# Copying the reference energy to be printted on the old CPIN format
refene2_old = _ReferenceEnergy(igb2=-65.113428, igb5=-64.166385, igb8=-61.3305355)
refene2_old.solvent_energies(igb2=-65.003415, igb5=-64.047229)
refene2_old.dielc2_energies(igb2=-32.167520, igb5=-31.751177)
refene2_old.dielc2.solvent_energies()
refene2_old.set_pKa(9.6, deprotonated=True)

TYR = TitratableResidue('TYR', ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG',
                        'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'OH', 'HH', 'CE2',
                        'HE2', 'CD2', 'HD2', 'C', 'O'], pka=9.6, typ="ph")
TYR.add_state(protcnt=1, refene=refene1, refene_old=refene1, pka_corr=9.6, # protonated
              charges=[-0.4157, 0.2719, -0.0014, 0.0876, -0.0152, 0.0295,
              0.0295, -0.0011, -0.1906, 0.1699, -0.2341, 0.1656, 0.3226,
              -0.5579, 0.3992, -0.2341, 0.1656, -0.1906, 0.1699, 0.5973,
              -0.5679])
TYR.add_state(protcnt=0, refene=refene2, refene_old=refene2_old, pka_corr=0.0, # deprotonated
              charges=[-0.4157, 0.2719, -0.0014, 0.0876, -0.0858, 0.0190,
              0.0190, -0.2130, -0.1030, 0.1320, -0.4980, 0.1320, 0.7770,
              -0.8140, 0.0, -0.4980, 0.1320, -0.1030, 0.1320, 0.5973, -0.5679])
TYR.check()

# Histidine
refene1 = _ReferenceEnergy(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.dielc2_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.dielc2.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene2 = _ReferenceEnergy(igb1=-4.208863, igb2=-2.84183, igb5=-2.86001,
               igb7=-1.741947, igb8=-3.4000)
refene2.solvent_energies(igb2=-2.77641, igb5=-2.90517)
refene2.dielc2_energies(igb2=-1.628110, igb5=-1.691093)
refene2.dielc2.solvent_energies()
# Copying the reference energy to be printted on the old CPIN format
refene2_old = _ReferenceEnergy(igb1=-4.208863, igb2=-2.84183, igb5=-2.86001,
               igb7=-1.741947, igb8=-3.4000)
refene2_old.solvent_energies(igb2=-2.77641, igb5=-2.90517)
refene2_old.dielc2_energies(igb2=-1.628110, igb5=-1.691093)
refene2_old.dielc2.solvent_energies()
refene2_old.set_pKa(6.5, deprotonated=True)
refene3 = _ReferenceEnergy(igb1=-8.230643, igb2=-6.58793, igb5=-6.70726,
               igb7=-5.118453, igb8=-6.3190)
refene3.solvent_energies(igb2=-6.483630, igb5=-6.82684)
refene3.dielc2_energies(igb2=-3.444200, igb5=-3.070113)
refene3.dielc2.solvent_energies()
# Copying the reference energy to be printted on the old CPIN format
refene3_old = _ReferenceEnergy(igb1=-8.230643, igb2=-6.58793, igb5=-6.70726,
               igb7=-5.118453, igb8=-6.3190)
refene3_old.solvent_energies(igb2=-6.483630, igb5=-6.82684)
refene3_old.dielc2_energies(igb2=-3.444200, igb5=-3.070113)
refene3_old.dielc2.solvent_energies()
refene3_old.set_pKa(7.1, deprotonated=True)

HIP = TitratableResidue('HIP', ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG',
                        'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2',
                        'C', 'O'], pka=6.6, typ="ph")
HIP.add_state(protcnt=2, refene=refene1, refene_old=refene1, pka_corr=7.1, # HIP
              charges=[-0.3479, 0.2747, -0.1354, 0.1212, -0.0414, 0.0810,
              0.0810, -0.0012, -0.1513, 0.3866, -0.0170, 0.2681, -0.1718,
              0.3911, -0.1141, 0.2317, 0.7341, -0.5894])
HIP.add_state(protcnt=1, refene=refene2, refene_old=refene2_old, pka_corr=0.6, # HID
              charges=[-0.3479, 0.2747, -0.1354, 0.1212, -0.1110, 0.0402,
              0.0402, -0.0266, -0.3811, 0.3649, 0.2057, 0.1392, -0.5727, 0.0,
              0.1292, 0.1147, 0.7341, -0.5894])
HIP.add_state(protcnt=1, refene=refene3, refene_old=refene3_old, pka_corr=0.0, # HIE
              charges=[-0.3479, 0.2747, -0.1354, 0.1212, -0.1012, 0.0367,
              0.0367, 0.1868, -0.5432, 0.0, 0.1635, 0.1435, -0.2795, 0.3339,
              -0.2207, 0.1862, 0.7341, -0.5894])
HIP.check()

# Lysine
refene1 = _ReferenceEnergy(igb2=-15.2423959, igb5=-14.5392838, igb8=-18.393654)
refene1.solvent_energies(igb2=-15.1417977, igb5=-14.3152107)
refene1.dielc2_energies(igb2=-7.239587, igb5=-6.825997)
refene1.dielc2.solvent_energies()
# Copying the reference energy to be printted on the old CPIN format
refene1_old = _ReferenceEnergy(igb2=-15.2423959, igb5=-14.5392838, igb8=-18.393654)
refene1_old.solvent_energies(igb2=-15.1417977, igb5=-14.3152107)
refene1_old.dielc2_energies(igb2=-7.239587, igb5=-6.825997)
refene1_old.dielc2.solvent_energies()
refene1_old.set_pKa(10.4, deprotonated=False)
refene2 = _ReferenceEnergy(igb2=0, igb5=0, igb8=0)
refene2.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene2.dielc2_energies(igb2=0, igb5=0, igb8=0)
refene2.dielc2.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)

LYS = TitratableResidue('LYS', ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG',
                        'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2' ,'HE3',
                        'NZ', 'HZ1', 'HZ2', 'HZ3', 'C', 'O'], pka=10.4, typ="ph")
LYS.add_state(protcnt=3, refene=refene1, refene_old=refene1_old, pka_corr=10.4, # protonated
              charges=[-0.3479, 0.2747, -0.2400, 0.1426, -0.0094, 0.0362,
              0.0362, 0.0187, 0.0103, 0.0103, -0.0479, 0.0621, 0.0621, -0.0143,
              0.1135, 0.1135, -0.3854, 0.3400, 0.3400, 0.3400, 0.7341, -0.5894])
LYS.add_state(protcnt=2, refene=refene2, refene_old=refene2, pka_corr=0.0, # deprotonated
              charges=[-0.3479, 0.2747, -0.2400, 0.1426, -0.10961, 0.0340,
              0.0340, 0.06612, 0.01041, 0.01041, -0.03768, 0.01155, 0.01155,
              0.32604, -0.03358, -0.03358, -1.03581, 0.0, 0.38604, 0.38604,
              0.7341, -0.5894])
LYS.check()

# Cysteine
refene1 = _ReferenceEnergy(igb2=77.4666763, igb5=76.2588331, igb8=71.5804519)
refene1.solvent_energies(igb2=77.6041407, igb5=76.2827217)
refene1.dielc2_energies(igb2=38.090523, igb5=37.454637)
refene1.dielc2.solvent_energies(igb2=38.489170)
# Copying the reference energy to be printted on the old CPIN format
refene1_old = _ReferenceEnergy(igb2=77.4666763, igb5=76.2588331, igb8=71.5804519)
refene1_old.solvent_energies(igb2=77.6041407, igb5=76.2827217)
refene1_old.dielc2_energies(igb2=38.090523, igb5=37.454637)
refene1_old.dielc2.solvent_energies(igb2=38.489170)
refene1_old.set_pKa(8.5, deprotonated=False)
refene2 = _ReferenceEnergy(igb2=0, igb5=0, igb8=0)
refene2.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene2.dielc2_energies(igb2=0, igb5=0, igb8=0)
refene2.dielc2.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)

CYS = TitratableResidue('CYS', ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'SG',
                        'HG', 'C', 'O'], pka=8.5, typ="ph")
CYS.add_state(protcnt=1, refene=refene1, refene_old=refene1_old, pka_corr=8.5, # protonated
              charges=[-0.4157, 0.2719, 0.0213, 0.1124, -0.1231, 0.1112, 0.1112,
                       -0.3119, 0.1933, 0.5973, -0.5679])
CYS.add_state(protcnt=0, refene=refene2, refene_old=refene2, pka_corr=0.0, # deprotonated
              charges=[-0.4157, 0.2719, 0.0213, 0.1124, -0.3593, 0.1122, 0.1122,
                       -0.8844, 0.0, 0.5973, -0.5679])
CYS.check()

# Deoxy-adenine
refene1 = _ReferenceEnergy(igb2=-19.8442, igb5=-19.8442)
refene1.solvent_energies()
refene1.dielc2_energies(igb2=-9.106013, igb5=-9.404867)
refene1.dielc2.solvent_energies(igb2=-9.779586)
# Copying the reference energy to be printted on the old CPIN format
refene1_old = _ReferenceEnergy(igb2=-19.8442, igb5=-19.8442)
refene1_old.solvent_energies()
refene1_old.dielc2_energies(igb2=-9.106013, igb5=-9.404867)
refene1_old.dielc2.solvent_energies(igb2=-9.779586)
refene1_old.set_pKa(3.9, deprotonated=True)
refene2 = _ReferenceEnergy(igb2=0, igb5=0, igb8=0)
refene2.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene2.dielc2_energies(igb2=0, igb5=0, igb8=0)
refene2.dielc2.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)

DAP = TitratableResidue('DAP', ['P', 'O1P', 'O2P', "O5'", "C5'", "H5'1", "H5'2",
                        "C4'", "H4'", "O4'", "C1'", "H1'", 'N9', 'C8', 'H8',
                        'N7', 'C5', 'C6', 'N6', 'H61', 'H62', 'N1', 'C2', 'H2',
                        'N3', 'C4', "C3'", "H3'", "C2'", "H2'1", "H2'2", "O3'",
                        'H1'], pka=3.9, typ="ph")
DAP.add_state(protcnt=1, refene=refene1, refene_old=refene1_old, pka_corr=0.0, # deprotonated
              charges=[1.1659, -0.7761, -0.7761, -0.4954, -0.0069, 0.0754,
              0.0754, 0.1629, 0.1176, -0.3691, 0.0431, 0.1838, -0.0268, 0.1607,
              0.1877, -0.6175, 0.0725, 0.6897, -0.9123, 0.4167, 0.4167, -0.7624,
              0.5716, 0.0598, -0.7417, 0.38, 0.0713, 0.0985, -0.0854, 0.0718,
              0.0718, -0.5232, 0.0])
DAP.add_state(protcnt=2, refene=refene2, refene_old=refene2, pka_corr=3.9, # protonated
              charges=[1.1659, -0.7761, -0.7761, -0.4954, -0.0069, 0.0754,
              0.0754, 0.1629, 0.1176, -0.3691, 0.0431, 0.1838, 0.0944, 0.1617,
              0.2281, -0.5674, 0.1358, 0.5711, -0.8251, 0.4456, 0.4456, -0.575,
              0.4251, 0.1437, -0.5611, 0.3421, 0.0713, 0.0985, -0.0854, 0.0718,
              0.0718, -0.5232, 0.4301])
DAP.check()

# Deoxy-cytosine
refene1 = _ReferenceEnergy(igb2=-40.526, igb5=-40.526)
refene1.solvent_energies()
refene1.dielc2_energies(igb2=-19.447553, igb5=-19.842087)
refene1.dielc2.solvent_energies(igb2=-20.121129)
# Copying the reference energy to be printted on the old CPIN format
refene1_old = _ReferenceEnergy(igb2=-40.526, igb5=-40.526)
refene1_old.solvent_energies()
refene1_old.dielc2_energies(igb2=-19.447553, igb5=-19.842087)
refene1_old.dielc2.solvent_energies(igb2=-20.121129)
refene1_old.set_pKa(4.3, deprotonated=True)
refene2 = _ReferenceEnergy(igb2=0, igb5=0)
refene2.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene2.dielc2_energies(igb2=0, igb5=0, igb8=0)
refene2.dielc2.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)

DCP = TitratableResidue('DCP', ['P', 'O1P', 'O2P', "O5'", "C5'", "H5'1", "H5'2",
                        "C4'", "H4'", "O4'", "C1'", "H1'", 'N1', 'C6', 'H6',
                        'C5', 'H5', 'C4', 'N4', 'H41', 'H42', 'N3', 'C2', 'O2',
                        "C3'", "H3'", "C2'", "H2'1", "H2'2", "O3'", 'H3'],
                        pka=4.3, typ="ph")
DCP.add_state(protcnt=1, refene=refene1, refene_old=refene1_old, pka_corr=0.0, # deprotonated
              charges=[1.1659, -0.7761, -0.7761, -0.4954, -0.0069, 0.0754,
              0.0754, 0.1629, 0.1176, -0.3691, -0.0116, 0.1963, -0.0339,
              -0.0183, 0.2293, -0.5222, 0.1863, 0.8439, -0.9773, 0.4314, 0.4314,
              -0.7748, 0.7959, -0.6548, 0.0713, 0.0985, -0.0854, 0.0718, 0.0718,
              -0.5232, 0.0])
DCP.add_state(protcnt=2, refene=refene2, refene_old=refene2, pka_corr=4.3, # protonated
              charges=[1.1659, -0.7761, -0.7761, -0.4954, -0.0069, 0.0754,
              0.0754, 0.1629, 0.1176, -0.3691, -0.0116, 0.1963, 0.2167, -0.0282,
              0.2713, -0.4162, 0.2179, 0.6653, -0.859, 0.4598, 0.4598, -0.4956,
              0.5371, -0.5028, 0.0713, 0.0985, -0.0854, 0.0718, 0.0718, -0.5232,
              0.4108])
DCP.check()

# Deoxy-guanine
refene1 = _ReferenceEnergy(igb2=0, igb5=0, igb8=0)
refene1.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.dielc2_energies(igb2=0, igb5=0, igb8=0)
refene1.dielc2.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene2 = _ReferenceEnergy(igb2=-90.0011, igb5=-90.0011)
refene2.solvent_energies()
refene2.dielc2_energies(igb2=-44.031593, igb5=-43.588343)
refene2.dielc2.solvent_energies(igb2=-45.090067)
# Copying the reference energy to be printted on the old CPIN format
refene2_old = _ReferenceEnergy(igb2=-90.0011, igb5=-90.0011)
refene2_old.solvent_energies()
refene2_old.dielc2_energies(igb2=-44.031593, igb5=-43.588343)
refene2_old.dielc2.solvent_energies(igb2=-45.090067)
refene2_old.set_pKa(9.2, deprotonated=True)

DG = TitratableResidue('DG', ['P', 'O1P', 'O2P', "O5'", "C5'", "H5'1", "H5'2",
                       "C4'", "H4'", "O4'", "C1'", "H1'", 'N9', 'C8', 'H8',
                       'N7', 'C5', 'C6', 'O6', 'N1', 'H1', 'C2', 'N2', 'H21',
                       'H22', 'N3', 'C4', "C3'", "H3'", "C2'", "H2'1", "H2'2",
                       "O3'"], pka=9.2, typ="ph")
DG.add_state(protcnt=1, refene=refene1, refene_old=refene1, pka_corr=9.2, # protonated
             charges=[1.1659, -0.7761, -0.7761, -0.4954, -0.0069, 0.0754,
             0.0754, 0.1629, 0.1176, -0.3691, 0.0358, 0.1746, 0.0577, 0.0736,
             0.1997, -0.5725, 0.1991, 0.4918, -0.5699, -0.5053, 0.352, 0.7432,
             -0.923, 0.4235, 0.4235, -0.6636, 0.1814, 0.0713, 0.0985, -0.0854,
             0.0718, 0.0718, -0.5232])
DG.add_state(protcnt=0, refene=refene2, refene_old=refene2_old, pka_corr=0.0, # deprotonated
             charges=[1.1659, -0.7761, -0.7761, -0.4954, -0.0069, 0.0754,
             0.0754, 0.1629, 0.1176, -0.3691, 0.0358, 0.1746, -0.0507, 0.0779,
             0.1516, -0.6122, 0.0806, 0.7105, -0.7253, -0.8527, 0.0, 0.9561,
             -0.9903, 0.3837, 0.3837, -0.8545, 0.2528, 0.0713, 0.0985, -0.0854,
             0.0718, 0.0718, -0.5232])
DG.check()

# Deoxy-thymine
refene1 = _ReferenceEnergy(igb2=0, igb5=0)
refene1.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.dielc2_energies(igb2=0, igb5=0, igb8=0)
refene1.dielc2.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene2 = _ReferenceEnergy(igb2=-56.7729, igb5=-56.7729)
refene2.solvent_energies(igb2=-28.429391)
refene2.dielc2_energies(igb2=-28.085730, igb5=-27.298290)
refene2.dielc2.solvent_energies()
# Copying the reference energy to be printted on the old CPIN format
refene2_old = _ReferenceEnergy(igb2=-56.7729, igb5=-56.7729)
refene2_old.solvent_energies(igb2=-28.429391)
refene2_old.dielc2_energies(igb2=-28.085730, igb5=-27.298290)
refene2_old.dielc2.solvent_energies()
refene2_old.set_pKa(9.7, deprotonated=True)

DT = TitratableResidue('DT', ['P', 'O1P', 'O2P', "O5'", "C5'", "H5'1", "H5'2",
                       "C4'", "H4'", "O4'", "C1'", "H1'", 'N1', 'C6', 'H6',
                       'C5', 'C7', 'H71', 'H72', 'H73', 'C4', 'O4', 'N3', 'H3',
                       'C2', 'O2', "C3'", "H3'", "C2'", "H2'1", "H2'2", "O3'"],
                       pka=9.7, typ="ph")
DT.add_state(protcnt=1, refene=refene1, refene_old=refene1, pka_corr=9.7, # protonated
             charges=[1.1659, -0.7761, -0.7761, -0.4954, -0.0069, 0.0754,
             0.0754, 0.1629, 0.1176, -0.3691, 0.068, 0.1804, -0.0239, -0.2209,
             0.2607, 0.0025, -0.2269, 0.077, 0.077, 0.077, 0.5194, -0.5563,
             -0.434, 0.342, 0.5677, -0.5881, 0.0713, 0.0985, -0.0854, 0.0718,
             0.0718, -0.5232])
DT.add_state(protcnt=0, refene=refene2, refene_old=refene2_old, pka_corr=0.0, # deprotonated
             charges=[1.1659, -0.7761, -0.7761, -0.4954, -0.0069, 0.0754,
             0.0754, 0.1629, 0.1176, -0.3691, 0.068, 0.1804, -0.2861, -0.1874,
             0.2251, -0.1092, -0.2602, 0.0589, 0.0589, 0.0589, 0.8263, -0.7396,
             -0.9169, 0.0, 0.9167, -0.7722, 0.0713, 0.0985, -0.0854, 0.0718,
             0.0718, -0.5232])
DT.check()

# Adenine
refene1 = _ReferenceEnergy(igb2=0, igb5=0)
refene1.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.dielc2_energies(igb2=0, igb5=0, igb8=0)
refene1.dielc2.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene2 = _ReferenceEnergy(igb2=14.851840012, igb5=15.1166001033)
refene2.solvent_energies(igb2=15.026098360790293,igb5=15.143651997474915)
refene2.solvent_energies()
refene2.dielc2_energies(igb2=6.953887, igb5=7.092043)
refene2.dielc2.solvent_energies(igb2=7.544988)
# Copying the reference energy to be printted on the old CPIN format
refene2_old = _ReferenceEnergy(igb2=14.851840012, igb5=15.1166001033)
refene2_old.solvent_energies(igb2=15.026098360790293,igb5=15.143651997474915)
refene2_old.solvent_energies()
refene2_old.dielc2_energies(igb2=6.953887, igb5=7.092043)
refene2_old.dielc2.solvent_energies(igb2=7.544988)
refene2_old.set_pKa(3.5, deprotonated=False)

AP = TitratableResidue('AP', ['P', 'O1P', 'O2P', "O5'", "C5'", "H5'1", "H5'2",
                       "C4'", "H4'", "O4'", "C1'", "H1'", 'N9', 'C8', 'H8',
                       'N7', 'C5', 'C6', 'N6', 'H61', 'H62', 'N1', 'C2', 'H2',
                       'N3', 'C4', "C3'", "H3'", "C2'", "H2'1", "O2'", "HO'2",
                       "O3'", 'H1'], pka=3.9, typ="ph")
AP.add_state(protcnt=0, refene=refene1, refene_old=refene1, pka_corr=0.0, # deprotonated
             charges=[1.1662, -0.776, -0.776, -0.4989, 0.0558, 0.0679, 0.0679,
             0.1065, 0.1174, -0.3548, 0.0394, 0.2007, -0.0251, 0.2006, 0.1553,
             -0.6073, 0.0515, 0.7009, -0.9019, 0.4115, 0.4115, -0.7615, 0.5875,
             0.0473, -0.6997, 0.3053, 0.2022, 0.0615, 0.067, 0.0972, -0.6139,
             0.4186, -0.5246, 0.0])
AP.add_state(protcnt=1, refene=refene2, refene_old=refene2_old, pka_corr=3.5, # protonated
             charges=[1.1662, -0.776, -0.776, -0.4989, 0.0558, 0.0679, 0.0679,
             0.1065, 0.1174, -0.3548, 0.0394, 0.2007, 0.0961, 0.2011, 0.1965,
             -0.5569, 0.1136, 0.5845, -0.8152, 0.4403, 0.4403, -0.5776, 0.4435,
             0.1307, -0.5201, 0.2681, 0.2022, 0.0615, 0.067, 0.0972, -0.6139,
             0.4186, -0.5246, 0.431])
AP.check()

# Cytosine
refene1 = _ReferenceEnergy(igb2=0, igb5=0)
refene1.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.dielc2_energies(igb2=0, igb5=0, igb8=0)
refene1.dielc2.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene2 = _ReferenceEnergy(igb2=37.501800178, igb5=38.0081251132)
refene2.solvent_energies(igb2=37.378544257354164, igb5=37.90444570773976)
refene2.dielc2_energies(igb2=18.483513, igb5=19.016390)
refene2.dielc2.solvent_energies()
# Copying the reference energy to be printted on the old CPIN format
refene2_old = _ReferenceEnergy(igb2=37.501800178, igb5=38.0081251132)
refene2_old.solvent_energies(igb2=37.378544257354164, igb5=37.90444570773976)
refene2_old.dielc2_energies(igb2=18.483513, igb5=19.016390)
refene2_old.dielc2.solvent_energies()
refene2_old.set_pKa(4.2, deprotonated=False)

CP = TitratableResidue('CP', ['P', 'O1P', 'O2P', "O5'", "C5'", "H5'1", "H5'2",
                       "C4'", "H4'", "O4'", "C1'", "H1'", 'N1', 'C6', 'H6',
                       'C5', 'H5', 'C4', 'N4', 'H41', 'H42', 'N3', 'C2', 'O2',
                       "C3'", "H3'", "C2'", "H2'1", "O2'", "HO'2", "O3'", 'H3'],
                       pka=4.3, typ="ph")
CP.add_state(protcnt=1, refene=refene1, refene_old=refene1, pka_corr=0.0, # deprotonated
             charges=[1.1662, -0.776, -0.776, -0.4989, 0.0558, 0.0679, 0.0679,
             0.1065, 0.1174, -0.3548, 0.0066, 0.2029, -0.0484, 0.0053, 0.1958,
             -0.5215, 0.1928, 0.8185, -0.953, 0.4234, 0.4234, -0.7584, 0.7538,
             -0.6252, 0.2022, 0.0615, 0.067, 0.0972, -0.6139, 0.4186, -0.5246,
             0.0])
CP.add_state(protcnt=2, refene=refene2, refene_old=refene2_old, pka_corr=4.2, # protonated
             charges=[1.1662, -0.776, -0.776, -0.4989, 0.0558, 0.0679, 0.0679,
             0.1065, 0.1174, -0.3548, 0.0066, 0.2029, 0.1954, 0.0028, 0.2366,
             -0.4218, 0.2253, 0.6466, -0.8363, 0.4518, 0.4518, -0.4871, 0.5039,
             -0.4753, 0.2022, 0.0615, 0.067, 0.0972, -0.6139, 0.4186, -0.5246,
             0.4128])
CP.check()

# Guanine
refene1 = _ReferenceEnergy(igb2=0, igb5=0)
refene1.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.dielc2_energies(igb2=0, igb5=0, igb8=0)
refene1.dielc2.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene2 = _ReferenceEnergy(igb2=-97.094725165, igb5=-96.0365352027)
refene2.solvent_energies(igb2=-97.31657849010276, igb5=-95.95654436492156)
refene2.dielc2_energies(igb2=-47.410980, igb5=-47.008233)
refene2.dielc2.solvent_energies(igb2=-48.222021)
# Copying the reference energy to be printted on the old CPIN format
refene2_old = _ReferenceEnergy(igb2=-97.094725165, igb5=-96.0365352027)
refene2_old.solvent_energies(igb2=-97.31657849010276, igb5=-95.95654436492156)
refene2_old.dielc2_energies(igb2=-47.410980, igb5=-47.008233)
refene2_old.dielc2.solvent_energies(igb2=-48.222021)
refene2_old.set_pKa(9.2, deprotonated=True)

G = TitratableResidue('G', ['P', 'O1P', 'O2P', "O5'", "C5'", "H5'1", "H5'2",
                      "C4'", "H4'", "O4'", "C1'", "H1'", 'N9', 'C8', 'H8', 'N7',
                      'C5', 'C6', 'O6', 'N1', 'H1', 'C2', 'N2', 'H21', 'H22',
                      'N3', 'C4', "C3'", "H3'", "C2'", "H2'1", "O2'", "HO'2",
                      "O3'"], pka=9.2, typ="ph")
G.add_state(protcnt=1, refene=refene1, refene_old=refene1, pka_corr=9.2, # protonated
            charges=[1.1662, -0.776, -0.776, -0.4989, 0.0558, 0.0679, 0.0679,
            0.1065, 0.1174, -0.3548, 0.0191, 0.2006, 0.0492, 0.1374, 0.164,
            -0.5709, 0.1744, 0.477, -0.5597, -0.4787, 0.3424, 0.7657, -0.9672,
            0.4364, 0.4364, -0.6323, 0.1222, 0.2022, 0.0615, 0.067, 0.0972,
            -0.6139, 0.4186, -0.5246])
G.add_state(protcnt=0, refene=refene2, refene_old=refene2_old, pka_corr=0.0, # deprotonated
            charges=[1.1662, -0.776, -0.776, -0.4989, 0.0558, 0.0679, 0.0679,
            0.1065, 0.1174, -0.3548, 0.0191, 0.2006, -0.0623, 0.1479, 0.1137,
            -0.6127, 0.0488, 0.7137, -0.7191, -0.8557, 0.0, 0.9976, -1.0387,
            0.3969, 0.3969, -0.8299, 0.1992, 0.2022, 0.0615, 0.067, 0.0972,
            -0.6139, 0.4186, -0.5246])
G.check()

# Uracil
refene1 = _ReferenceEnergy(igb2=0, igb5=0)
refene1.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.dielc2_energies(igb2=0, igb5=0, igb8=0)
refene1.dielc2.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene2 = _ReferenceEnergy(igb2=-136.326020191, igb5=-134.938275039)
refene2.solvent_energies(igb2=-136.5653533428478, igb5=-135.06973320905044)
refene2.dielc2_energies(igb2=-67.270690, igb5=-66.605330)
refene2.dielc2.solvent_energies()
# Copying the reference energy to be printted on the old CPIN format
refene2_old = _ReferenceEnergy(igb2=-136.326020191, igb5=-134.938275039)
refene2_old.solvent_energies(igb2=-136.5653533428478, igb5=-135.06973320905044)
refene2_old.dielc2_energies(igb2=-67.270690, igb5=-66.605330)
refene2_old.dielc2.solvent_energies()
refene2_old.set_pKa(9.2, deprotonated=True)

U = TitratableResidue('U', ['P', 'O1P', 'O2P', "O5'", "C5'", "H5'1", "H5'2",
                      "C4'", "H4'", "O4'", "C1'", "H1'", 'N1', 'C6', 'H6', 'C5',
                      'H5', 'C4', 'O4', 'N3', 'H3', 'C2', 'O2', "C3'", "H3'",
                      "C2'", "H2'1", "O2'", "HO'2", "O3'"], pka=9.3, typ="ph")
U.add_state(protcnt=1, refene=refene1, refene_old=refene1, pka_corr=9.2, # protonated
            charges=[1.1662, -0.776, -0.776, -0.4989, 0.0558, 0.0679, 0.0679,
            0.1065, 0.1174, -0.3548, 0.0674, 0.1824, 0.0418, -0.1126, 0.2188,
            -0.3635, 0.1811, 0.5952, -0.5761, -0.3549, 0.3154, 0.4687, -0.5477,
            0.2022, 0.0615, 0.067, 0.0972, -0.6139, 0.4186, -0.5246])
U.add_state(protcnt=0, refene=refene2, refene_old=refene2_old, pka_corr=0.0, # deprotonated
            charges=[1.1662, -0.776, -0.776, -0.4989, 0.0558, 0.0679, 0.0679,
            0.1065, 0.1174, -0.3548, 0.0674, 0.1824, -0.2733, 0.0264, 0.1501,
            -0.582, 0.156, 0.9762, -0.7808, -0.9327, 0.0, 0.8698, -0.7435,
            0.2022, 0.0615, 0.067, 0.0972, -0.6139, 0.4186, -0.5246])
U.check()

# HEH: HEME ring + parts of 2 HIS + 2 CYS
refene1 = _ReferenceEnergy(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.dielc2_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.dielc2.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene2 = _ReferenceEnergy(igb2=-15.493731, igb5=-16.349152, igb7=-16.509509, igb8=-22.025653) # Implicit
refene2.solvent_energies(igb2=-15.209270, igb5=-15.840853, igb7=-15.495868) # Explicit
refene2.dielc2_energies()
refene2.dielc2.solvent_energies()

HEH = TitratableResidue('HEH', ['FE', 'NA', 'C1A', 'C2A', 'C3A', 'CMA', 'HMA1', 'HMA2', 'HMA3', 'C4A',
                                'CHB', 'HHB', 'C1B', 'NB', 'C2B', 'CMB', 'HMB1', 'HMB2', 'HMB3', 'C3B',
                                'CAB', 'HAB', 'CBB', 'HBB1', 'HBB2', 'HBB3', 'C4B', 'CHC', 'HHC', 'C1C',
                                'NC', 'C2C', 'CMC', 'HMC1', 'HMC2', 'HMC3', 'C3C', 'CAC', 'HAC', 'CBC',
                                'HBC1', 'HBC2', 'HBC3', 'C4C', 'CHD', 'HHD', 'C1D', 'ND', 'C2D', 'CMD',
                                'HMD1', 'HMD2', 'HMD3', 'C3D', 'C4D', 'CHA', 'HHA', 'CBC1', 'HB2C',
                                'HB3C', 'SGC1', 'CB1', 'HB21', 'HB31', 'CG1', 'ND11', 'HD11', 'CE11',
                                'HE11', 'NE21', 'CD21', 'HD21', 'CBB2', 'HB2B', 'HB3B', 'SGB2', 'CB2',
                                'HB22', 'HB32', 'CG2', 'ND12', 'HD12', 'CE12', 'HE12', 'NE22', 'CD22',
                                'HD22'],
                        eo=-0.203, typ="redox")
HEH.add_state(eleccnt=2, refene=refene1, eo_corr=0.0, # FE3+ (oxidized state)
              charges=[ 0.6660, -0.1530, -0.0956,  0.1274,  0.1624, -0.2600,  0.0743,  0.0743,  0.0743,
                       -0.0766, -0.0586,  0.1300, -0.0206, -0.2560,  0.1394, -0.2240,  0.0663,  0.0663,
                        0.0663,  0.0734, -0.0935,  0.1765, -0.4035,  0.1375,  0.1375,  0.1375,  0.0504,
                       -0.1806,  0.1430,  0.0364, -0.2390,  0.1784, -0.2325,  0.0635,  0.0635,  0.0635,
                        0.0084, -0.0570,  0.2130, -0.4090,  0.1357,  0.1357,  0.1357, -0.0266, -0.0576,
                        0.1360, -0.1156, -0.1130,  0.1604, -0.2575,  0.0752,  0.0752,  0.0752,  0.1444,
                       -0.1126,  0.0104,  0.1090, -0.3349,  0.1297,  0.1297, -0.1760, -0.0803,  0.0269,
                        0.0269,  0.1990, -0.2930,  0.3650,  0.0120,  0.1180, -0.0400, -0.2020,  0.1790,
                       -0.3349,  0.1297,  0.1297, -0.1760, -0.0803,  0.0269,  0.0269,  0.1990, -0.2930,
                        0.3650,  0.0120,  0.1180, -0.0400, -0.2020,  0.1790]
              )
HEH.add_state(eleccnt=3, refene=refene2, eo_corr=-0.203, # FE2+ (reduced state)
              charges=[ 0.4800, -0.1337, -0.1455,  0.1285,  0.1325, -0.2545,  0.0608,  0.0608,  0.0608,
                       -0.0865, -0.0815,  0.1220, -0.0425, -0.2490,  0.1605, -0.1935,  0.0422,  0.0422,
                        0.0422, -0.0025, -0.0390,  0.1720, -0.4255,  0.1348,  0.1348,  0.1348,  0.0405,
                       -0.1725,  0.1320, -0.0525, -0.1980,  0.2125, -0.1985,  0.0405,  0.0405,  0.0405,
                       -0.0355, -0.0195,  0.1915, -0.4340,  0.1320,  0.1320,  0.1320, -0.0275, -0.0805,
                        0.1290, -0.1275, -0.1020,  0.1335, -0.2525,  0.0615,  0.0615,  0.0615,  0.1415,
                       -0.1575,  0.0275,  0.0950, -0.3343,  0.1300,  0.1300, -0.2315, -0.0967,  0.0187,
                        0.0187,  0.2030, -0.3450,  0.3550,  0.0110,  0.1090, -0.0270, -0.2170,  0.1750,
                       -0.3343,  0.1300,  0.1300, -0.2315, -0.0967,  0.0187,  0.0187,  0.2030, -0.3450,
                        0.3550,  0.0110,  0.1090, -0.0270, -0.2170,  0.1750]
              )
HEH.check()

# Propionate
refene1 = _ReferenceEnergy(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.dielc2_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.dielc2.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene2 = _ReferenceEnergy(igb2=10.356928, igb5=9.943308, igb7=7.020632, igb8=5.259028) # Implicit
refene2.solvent_energies(igb2=16.751825, igb5=15.661934, igb7=12.906876, igb8=15.024553) # Explicit
refene2.dielc2_energies()
refene2.dielc2.solvent_energies()
# Copying the reference energy to be printted on the old CPIN format
refene2_old = _ReferenceEnergy(igb2=10.356928, igb5=9.943308, igb7=7.020632, igb8=5.259028) # Implicit
refene2_old.solvent_energies(igb2=16.751825, igb5=15.661934, igb7=12.906876, igb8=15.024553) # Explicit
refene2_old.dielc2_energies()
refene2_old.dielc2.solvent_energies()
refene2_old.set_pKa(4.85, deprotonated=False)

PRN = TitratableResidue('PRN',
                        ['CA', 'HA1', 'HA2', 'CB', 'HB1', 'HB2', 'CG',
                        'O1', 'O2', 'H11', 'H12', 'H21', 'H22'], pka=4.85, typ="ph")

PRN.add_state(protcnt=0, refene=refene1, refene_old=refene1, pka_corr=0.0, # deprotonated
              charges=[-0.0508, -0.0173,
              -0.0173, 0.0026, -0.0425, -0.0425, 0.8054, -0.8188, -0.8188, 0.0,
              0.0, 0.0, 0.0])

PRN.add_state(protcnt=1, refene=refene2, refene_old=refene2_old, pka_corr=4.85, # protonated syn-O1
              charges=[-0.0181, 0.0256, 0.0256,
              -0.0284, 0.0430, 0.0430, 0.6801, -0.6511, -0.5838, 0.4641,
              0.0, 0.0, 0.0])

PRN.add_state(protcnt=1, refene=refene2, refene_old=refene2_old, pka_corr=4.85, # protonated anti-O1
              charges=[-0.0181, 0.0256, 0.0256,
              -0.0284, 0.0430, 0.0430, 0.6801, -0.6511, -0.5838, 0.0,
              0.4641, 0.0, 0.0])

PRN.add_state(protcnt=1, refene=refene2, refene_old=refene2_old, pka_corr=4.85, # protonated syn-O2
              charges=[-0.0181, 0.0256, 0.0256,
              -0.0284, 0.0430, 0.0430, 0.6801, -0.5838, -0.6511, 0.0,
              0.0, 0.4641, 0.0])

PRN.add_state(protcnt=1, refene=refene2, refene_old=refene2_old, pka_corr=4.85, # protonated anti-O2
              charges=[-0.0181, 0.0256, 0.0256,
              -0.0284, 0.0430, 0.0430, 0.6801, -0.5838, -0.6511, 0.0,
              0.0, 0.0, 0.4641])

PRN.check()

# Tyrosine, pH and redox active
refene1 = _ReferenceEnergy(igb2=0, igb5=0, igb8=0)
refene1.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene1.dielc2_energies(igb2=0, igb5=0, igb8=0)
refene1.dielc2.solvent_energies(igb1=0, igb2=0, igb5=0, igb7=0, igb8=0)
refene2 = _ReferenceEnergy(igb2=5.409084) # Implicit
refene2.solvent_energies(igb2=5.210075) # Explicit
refene2.dielc2_energies()
refene2.dielc2.solvent_energies()
refene3 = _ReferenceEnergy(igb2=-52.293859) # Implicit
refene3.solvent_energies(igb2=-52.185005) # Explicit
refene3.dielc2_energies()
refene3.dielc2.solvent_energies()
refene4 = _ReferenceEnergy(igb2=24.166408) # Implicit
refene4.solvent_energies(igb2=24.224913) # Explicit
refene4.dielc2_energies()
refene4.dielc2.solvent_energies()

TYX = TitratableResidue('TYX', ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG',
                        'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'OH', 'HH', 'CE2',
                        'HE2', 'CD2', 'HD2', 'C', 'O'], typ="phredox")
# Note: pka_corr differences only make sense between states with the same number of electrons
# Note: eo_corr differences only make sense between states with the same number of protons
TYX.add_state(protcnt=1, eleccnt=1, refene=refene1, pka_corr=9.6, eo_corr=1.4, # tyrOH
              charges=[-0.4157, 0.2719, -0.0014, 0.0876, -0.1163, 0.0548,
                       0.0548, -0.0139, -0.1142, 0.1615, -0.3410, 0.1911,
                       0.4198, -0.5278, 0.3621, -0.3410, 0.1911, -0.1142,
                       0.1615, 0.5973, -0.5679])
TYX.add_state(protcnt=1, eleccnt=0, refene=refene2, pka_corr=-2.0, eo_corr=0.0, # tyrOH+
              charges=[-0.4157, 0.2719, -0.0014, 0.0876, -0.2911, 0.1123,
                       0.1123, 0.4479, -0.1923, 0.2177, -0.1736, 0.2109,
                       0.5560, -0.4376, 0.4031, -0.1736, 0.2109, -0.1923,
                       0.2177, 0.5973, -0.5679])
TYX.add_state(protcnt=0, eleccnt=1, refene=refene3, pka_corr=0.0, eo_corr=0.71, # tyrO-
              charges=[-0.4157, 0.2719, -0.0014, 0.0876, -0.0500, 0.0300,
                       0.0300, -0.1786, -0.1545, 0.1418, -0.4793, 0.1318,
                       0.7197, -0.8026, 0.0000, -0.4793, 0.1318, -0.1545,
                       0.1418, 0.5974, -0.5678])
TYX.add_state(protcnt=0, eleccnt=0, refene=refene4, pka_corr=0.0, eo_corr=0.0, # tyrO
              charges=[-0.4157, 0.2719, -0.0014, 0.0876, -0.1941, 0.0931,
                       0.0931, 0.0538, -0.1216, 0.1629, -0.3224, 0.1643,
                       0.6679, -0.4519, 0.0000, -0.3224, 0.1643, -0.1216,
                       0.1629, 0.5973, -0.5679])
TYX.check()
