"""
This module contains the following classes:

   TinkerAnalout - Reads the output from the "analyze" command

"""
from ..exceptions import TinkerError
from .topologyobjects import (AtomList, BondList, AngleList,
      StretchBendList, UreyBradleyList, OutOfPlaneBendList, OutOfPlaneDistList,
      TorsionAngleList, PiTorsionList, TorsionTorsionList, AtomicMultipoleList,
      DipolePolarizabilityList, TorsionTorsionGrid)

class TinkerAnalout(object):
    """ Reads the output of "analyze" to determine system parameters """

    # Flags paired with the attributes in the pointers
    atom_inter_flags = {'Atoms in System' : 'natom', 
                        'Pisystem Atoms' : 'norbit',
                        'Bond Stretches' : 'nbond',
                        'Conjugated Pi-Bonds' : 'nbpi',
                        'Angle Bends' : 'nangle',
                        'Stretch-Bends' : 'nstrbnd',
                        'Urey-Bradley' : 'nurey',
                        'Angle-Angles' : 'nangang',
                        'Out-of-Plane Bends' : 'nopbend',
                        'Out-of-Plane Distances' : 'nopdist',
                        'Improper Dihedrals' : 'niprop',
                        'Improper Torsions' : 'nitors',
                        'Torsional Angles' : 'ntors',
                        'Pi-Orbital Torsions' : 'npitors',
                        'Stretch-Torsions' : 'nstrtor',
                        'Torsion-Torsions' : 'ntortor',
                        'Atomic Partial Charges' : 'nion',
                        'Bond Dipole Moments' : 'ndipole',
                        'Polarizable Multipoles' : 'npole',
                        'Number of 1-2 Pairs' : 'pair12',
                        'Number of 1-3 Pairs' : 'pair13',
                        'Number of 1-4 Pairs' : 'pair14',
                        'Number of 1-5 Pairs' : 'pair15',
    }

    def __init__(self, fname=None):
        if fname is not None:
            self.read(fname)

    def read(self, fname):
        """ Reads the analout file """
        self.pointers = dict(natom=0, norbit=0, nbond=0, nbpi=0, nangle=0,
                nstrbnd=0, nurey=0, nangang=0, nopbend=0, nopdist=0, niprop=0,
                nitors=0, ntors=0, npitors=0, nstrtor=0, ntortor=0, nion=0,
                ndipole=0, npole=0, pair12=0, pair13=0, pair14=0, pair15=0)
        self.fname = fname
        f = open(self.fname, 'r')
        try:
            line = f.readline()
            # Look for the TINKER watermark
            while True:
                if not line:
                    raise TinkerError(f'Could not find the TINKER watermark in {fname}')
                if line.lstrip().startswith('###            TINKER'):
                    break
                line = f.readline()
            # Get the numbers of all parameters
            while not line.lstrip().startswith(
                        'Total Numbers of Atoms and Interactions'):
                if not line:
                    raise TinkerError('Could not find the atom/ixn count')
                line = f.readline()
            # Eat the next line
            f.readline()
            line = f.readline()
            # Get all of the pointers
            while line.strip():
                try:
                    key = TinkerAnalout.atom_inter_flags[line[:27].strip()]
                except KeyError:
                    raise TinkerError(f'Unrecognized pointer keyword {key}')
                try:
                    self.pointers[key] = int(line[27:].strip())
                except ValueError:
                    raise TinkerError(f'Could not convert pointer {key} to int [{line.rstrip()}]')
                except KeyError:
                    raise Exception('Should not be here -- internal error')
                line = f.readline()
            # Check that we read in some pointers
            s = 0
            for key in self.pointers:
                s += self.pointers[key]
            if s <= 0:
                raise TinkerError('All pointers are 0')
            # Get the atoms in the next section
            while line.strip() != 'Atom Type Definition Parameters :':
                if not line:
                    raise TinkerError('Unexpected EOF when looking for atom definitions')
                line = f.readline()
            TinkerAnalout._read_atom_definitions(self, f)
            # Get all of the sections defined in _functionmap -- see bottom of
            # this file
            while True:
                line = f.readline()
                try:
                    self._functionmap[line.strip()](self, f)
                except KeyError:
                    break
        finally:
            f.close()

    # These are functions that take an open file and load a section of the file
    # into the data structures

    def _read_atom_definitions(self, f):
        # Make an atom list
        self.atom_list = AtomList()
        # Eat 3 lines, then begin
        f.readline(); f.readline(); f.readline()
        line = f.readline()
        for i in range(self.pointers['natom']):
            try:
                atnum = int(line[:6])
                self.atom_list.add(
                    line[11:14], line[14:21], line[21:28], line[28:34],
                    line[34:44], line[44:49], line[54:].strip()
                )
            except ValueError:
                raise TinkerError(f'Error parsing atomic properties\n\t[{line.rstrip()}]')
            line = f.readline()
            if atnum != i + 1:
                raise TinkerError(f'Atom number mismatch [{i + 1} vs {atnum}]')

    def _read_vdw(self, f):
        """ Reads the van der Waals parameters """
        if not hasattr(self, 'atom_list'):
            raise AttributeError('Atom definitions must be loaded prior to vdW!')
        # Eat the next 3 lines
        f.readline(); f.readline(); f.readline()
        line = f.readline()
        for i in range(self.pointers['natom']):
            try:
                int(line[0:6])
                atnum = int(line[9:15])
                self.atom_list[i].add_vdw(
                    line[22:32], line[32:42], line[43:53], line[53:63], line[64:74],
                )
            except ValueError:
                raise TinkerError('Error parsing van der Waals term')
            line = f.readline()
            if atnum != i + 1:
                raise TinkerError(f'Atom number mismatch in vdW [{i + 1} vs {atnum}]')

    def _read_bonds(self, f):
        """ Reads the bond stretching terms """
        self.bond_list = BondList()
        # Eat the next 3 lines
        f.readline(); f.readline(); f.readline()
        line = f.readline()
        for i in range(self.pointers['nbond']):
            try:
                int(line[0:6])
                at1 = int(line[9:15]) - 1
                at2 = int(line[15:21]) - 1
                self.bond_list.add(self.atom_list[at1], self.atom_list[at2], line[40:50], line[50:60])
            except ValueError:
                raise TinkerError('Error parsing bonded term')
            line = f.readline()

    def _read_angs(self, f):
        """ Reads the angle stretching terms """
        self.angle_list = AngleList()
        # Eat the next 3 lines
        f.readline(); f.readline(); f.readline()
        line = f.readline()
        for i in range(self.pointers['nangle']):
            try:
                int(line[0:6])
                at1 = int(line[9:15]) - 1
                at2 = int(line[15:21]) - 1
                at3 = int(line[21:27]) - 1
                self.angle_list.add(
                    self.atom_list[at1], self.atom_list[at2], self.atom_list[at3],
                    line[40:50], line[50:60], line[60:67], line[69:]
                )
            except ValueError as err:
                raise TinkerError('Error parsing angle term') from err
            line = f.readline()

    def _read_strbnd(self, f):
        """ Reads the stretch-bend terms """
        self.stretchbend_list = StretchBendList()
        # Eat the next 3 lines
        f.readline(); f.readline(); f.readline()
        line = f.readline()
        for i in range(self.pointers['nstrbnd']):
            try:
                int(line[0:6])
                at1 = int(line[9:15]) - 1
                at2 = int(line[15:21]) - 1
                at3 = int(line[21:27]) - 1
                self.stretchbend_list.add(
                    self.atom_list[at1], self.atom_list[at2], self.atom_list[at3],
                    line[27:40], line[40:50], line[50:60], line[60:70]
                )
            except ValueError:
                raise TinkerError('Error parsing stretch-bend term')
            line = f.readline()

    def _read_ureybrad(self, f):
        """ Reads the urey-bradley terms """
        self.ureybrad_list = UreyBradleyList()
        # Eat the next 3 lines
        f.readline(); f.readline(); f.readline()
        line = f.readline()
        for i in range(self.pointers['nurey']):
            try:
                int(line[0:6])
                at1 = int(line[9:15]) - 1
                at2 = int(line[15:21]) - 1
                # Support older-style analyze output that had all 3 atoms
                # involved in the angle. Newer analyze output contains only the
                # first and 3rd atoms, so that is all we pass into the
                # UreyBradleyList function
                try:
                    at3 = int(line[21:27]) - 1
                except ValueError:
                    at3 = at2
                self.ureybrad_list.add(
                    self.atom_list[at1], self.atom_list[at3], line[34:50], line[50:60]
                )
            except ValueError as err:
                raise TinkerError('Error parsing Urey-Bradley term') from err
            line = f.readline()

    def _read_opbend(self, f):
        """ Reads the out-of-plane bending terms """
        self.oopbend_list = OutOfPlaneBendList()
        # Eat the next 3 lines
        f.readline(); f.readline(); f.readline()
        line = f.readline()
        for i in range(self.pointers['nopbend']):
            try:
                int(line[0:6])
                at1 = int(line[9:15]) - 1
                at2 = int(line[15:21]) - 1
                at3 = int(line[21:27]) - 1
                at4 = int(line[27:33]) - 1
                self.oopbend_list.add(
                    self.atom_list[at1], self.atom_list[at2], self.atom_list[at3],
                    self.atom_list[at4], line[42:52]
                )
            except ValueError as err:
                raise TinkerError('Error parsing out-of-plane bending term') from err
            line = f.readline()

    def _read_opdist(self, f):
        """ Read the out-of-plane distance parameters """
        self.oopdist_list = OutOfPlaneDistList()
        # Eat the next 3 lines
        f.readline(); f.readline(); f.readline()
        line = f.readline()
        for i in range(self.pointers['nopdist']):
            try:
                int(line[0:6])
                at1 = int(line[9:15]) - 1
                at2 = int(line[15:21]) - 1
                at3 = int(line[21:27]) - 1
                at4 = int(line[27:33]) - 1
                self.oopdist_list.add(
                    self.atom_list[at1], self.atom_list[at2], self.atom_list[at3],
                    self.atom_list[at4], line[42:52]
                )
            except ValueError as err:
                raise TinkerError('Error parsing out-of-plane distance term') from err
            line = f.readline()

    def _read_torang(self, f):
        """ Read the torsion-angle parameters """
        self.torangle_list = TorsionAngleList()
        # Eat the next 3 lines
        f.readline(); f.readline(); f.readline()
        line = f.readline()
        for i in range(self.pointers['ntors']):
            try:
                int(line[0:6])
                at1 = int(line[9:15]) - 1
                at2 = int(line[15:21]) - 1
                at3 = int(line[21:27]) - 1
                at4 = int(line[27:33]) - 1
                # Get the rest of the terms (replace / with ' ' so we can do a
                # simple string split on whitespace)
                terms = line[33:].replace('/', ' ').split()
                self.torangle_list.add(
                    self.atom_list[at1], self.atom_list[at2], self.atom_list[at3],
                    self.atom_list[at4], terms
                )
            except ValueError as err:
                raise TinkerError('Error parsing torsion angle term') from err
            line = f.readline()

    def _read_pitors(self, f):
        """ Read the Pi-Orbital Torsion parameters """
        self.pitors_list = PiTorsionList()
        # Eat the next 3 lines
        f.readline(); f.readline(); f.readline()
        line = f.readline()
        for i in range(self.pointers['npitors']):
            try:
                int(line[0:6])
                at1 = int(line[9:15]) - 1
                at2 = int(line[15:21]) - 1
                self.pitors_list.add(self.atom_list[at1], self.atom_list[at2], line[40:50])
            except ValueError as err:
                raise TinkerError('Error parsing pi-torsion term') from err
            line = f.readline()

    def _read_tortors(self, f):
        """ Read the Torsion-Torsion parameters """
        self.tortor_list = TorsionTorsionList()
        # Eat the next 3 lines
        f.readline(); f.readline(); f.readline()
        line = f.readline()
        for i in range(self.pointers['ntortor']):
            try:
                int(line[0:6])
                at1 = int(line[9:15]) - 1
                at2 = int(line[15:21]) - 1
                at3 = int(line[21:27]) - 1
                at4 = int(line[27:33]) - 1
                at5 = int(line[33:39]) - 1
                dim1 = int(line[49:55])
                dim2 = int(line[55:61])
                self.tortor_list.add(
                    self.atom_list[at1], self.atom_list[at2], self.atom_list[at3],
                    self.atom_list[at4], self.atom_list[at5], dim1, dim2
                )
            except ValueError as err:
                raise TinkerError('Error parsing torsion-torsion term') from err
            line = f.readline()
            # The CMAP section was adjusted to print out the entire torsion
            # grid under each tor-tor definition. If this line has 3 words, we
            # need to eat the next dim1*dim2 lines
            if len(line.split()) == 3:
                gridvals = []
                for j in range(dim1*dim2):
                    gridvals.append(tuple([float(x) for x in line.split()]))
                    line = f.readline()
                self.tortor_list[-1].type = TorsionTorsionGrid.new(gridvals)

    def _read_multipoles(self, f):
        """ Read the atomic multipoles """
        self.multipole_list = AtomicMultipoleList()
        # Eat the next 3 lines
        f.readline(); f.readline(); f.readline()
        line = f.readline()
        for i in range(self.pointers['npole']):
            try:
                int(line[0:6])
                at = int(line[9:15]) - 1
                # Collect the rest of the arguments, kept as strings
                frame = [line[16:23], line[23:30], line[30:37]]
                typestr = line[40:48]
                moments = [line[50:59]]
                # Next numbers are on the next line
                line = f.readline()
                moments.extend([line[50:59], line[59:68], line[68:77]])
                line = f.readline()
                moments.append(line[50:59])
                line = f.readline()
                moments.extend([line[50:59], line[59:68]])
                line = f.readline()
                moments.extend([line[50:59], line[59:68], line[68:77]])
                self.multipole_list.add(self.atom_list[at], frame, typestr, moments)
            except ValueError as err:
                raise TinkerError('Error parsing multipole term') from err
            line = f.readline()

    def _read_dipoles(self, f):
        """ Read atomic dipole polarizabilities """
        self.dipole_list = DipolePolarizabilityList()
        # Eat the next 3 lines
        f.readline(); f.readline(); f.readline()
        line = f.readline()
        for i in range(self.pointers['npole']):
            try:
                int(line[0:6])
                at = int(line[9:15]) - 1
                self.dipole_list.add(self.atom_list[at], line[25:35], line[40:].split())
            except ValueError as err:
                raise TinkerError('Error parsing dipole polarizabilities') from err
            line = f.readline()

    def _read_interactions(self, f):
        """ Reads the number of interactions present in the system """
        # Eat the next line (1 only)
        f.readline()
        line = f.readline()
        while line.strip():
            try:
                key = TinkerAnalout.atom_inter_flags[line[1:20]]
            except KeyError:
                raise TinkerError(f'Unrecognized token in interaction count: [{line[1:20]}]')
            self.pointers[key] = int(line[21:])
            line = f.readline()

    def _read_12pairs(self, f):
        # Eat the next line
        self.pair12_list = set()
        f.readline()
        line = f.readline()
        for i in range(self.pointers['pair12']):
            at1, at2 = int(line[0:8]) - 1, int(line[8:16]) - 1
            self.pair12_list.add( (self.atom_list[at1], self.atom_list[at2]) )
            line = f.readline()

    def _read_13pairs(self, f):
        # Eat the next line
        self.pair13_list = set()
        f.readline()
        line = f.readline()
        for i in range(self.pointers['pair13']):
            at1, at2 = int(line[0:8]) - 1, int(line[8:16]) - 1
            self.pair13_list.add( (self.atom_list[at1], self.atom_list[at2]) )
            line = f.readline()

    def _read_14pairs(self, f):
        # Eat the next line
        self.pair14_list = set()
        f.readline()
        line = f.readline()
        for i in range(self.pointers['pair14']):
            at1, at2 = int(line[0:8]) - 1, int(line[8:16]) - 1
            self.pair14_list.add( (self.atom_list[at1], self.atom_list[at2]) )
            line = f.readline()

    def _read_15pairs(self, f):
        # Eat the next line
        self.pair15_list = set()
        f.readline()
        line = f.readline()
        for i in range(self.pointers['pair15']):
            at1, at2 = int(line[0:8]) - 1, int(line[8:16]) - 1
            self.pair15_list.add( (self.atom_list[at1], self.atom_list[at2]) )
            line = f.readline()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def __str__(self):
        return self.fname

    def __repr__(self):
        return f"<TinkerAnalout; {self.pointers['natom']} atoms>"


# Add the functionmap onto TinkerAnalout after the class has been built. Only
# then will the function references be complete

# Function map is a hash with the title string of an analout section
# mapping to the function responsible for reading that section and the
# key of the pointers flag (see self.pointers below) that declares how
# many numbers must be parsed from that section
TinkerAnalout._functionmap = {
    'Van der Waals Parameters :' : TinkerAnalout._read_vdw,
    'Bond Stretching Parameters :' : TinkerAnalout._read_bonds,
    'Angle Bending Parameters :' : TinkerAnalout._read_angs,
    'Stretch-Bend Parameters :' : TinkerAnalout._read_strbnd,
    'Urey-Bradley Parameters :' : TinkerAnalout._read_ureybrad,
    'Out-of-Plane Bending Parameters :' : TinkerAnalout._read_opbend,
    'Out-of-Plane Distance Parameters :' : TinkerAnalout._read_opdist,
    'Torsional Angle Parameters :' : TinkerAnalout._read_torang,
    'Pi-Orbital Torsion Parameters :' : TinkerAnalout._read_pitors,
    'Torsion-Torsion Parameters :' : TinkerAnalout._read_tortors,
    'Atomic Multipole Parameters :' : TinkerAnalout._read_multipoles,
    'Dipole Polarizability Parameters :' : TinkerAnalout._read_dipoles,
    'List of 1-2 Connected Atomic Interactions :' : TinkerAnalout._read_12pairs,
    'List of 1-3 Connected Atomic Interactions :' : TinkerAnalout._read_13pairs,
    'List of 1-4 Connected Atomic Interactions :' : TinkerAnalout._read_14pairs,
    'List of 1-5 Connected Atomic Interactions :' : TinkerAnalout._read_15pairs,
    'Total Number of Pairwise Atomic Interactions :' : TinkerAnalout._read_interactions,
}
