
class TimeseriesCollection:
    def __init__(self):
        self.timeseries = []
    def addTimeseries(self, ts):
        # if ts is not Timeseries do something
        self.timeseries.append(ts)
    def getAtomList(self):
        self._atomlist = []
        for ts in self.timeseries:
            self._atomlist += list(ts.getAtomList())
        return self._atomlist
    def getFormat(self):
        return ''.join([ts.getFormatCode() for ts in self.timeseries])
    def getBounds(self):
        if not hasattr(self, '_atomlist'):
            self.getAtomList()
        lowerb = min(self._atomlist)
        upperb = max(self._atomlist)
        return lowerb, upperb
    def getDataSize(self):
        return sum([ts.getDataSize() for ts in self.timeseries])
    def getAtomCounts(self):
        return [ts.getNumAtoms() for ts in self.timeseries]
    def getAuxData(self):
        return [ts.getAuxData() for ts in self.timeseries]
    def clear(self):
        self.timeseries = []

class Timeseries:
    def __init__(self, code, atoms, dsize):
        if type(atoms) == list:
            self.atoms = atoms
        else: # This is probably an AtomGroup
            self.atoms = atoms._atoms
        self.code = code
        self.numatoms = len(atoms)
        self.dsize = dsize
    def getAtomList(self):
        return [a.number for a in self.atoms]
    def getFormatCode(self):
        return self.code
    def getDataSize(self):
        return self.dsize
    def getNumAtoms(self):
        return self.numatoms
    def getAuxData(self):
        return [0.]*self.numatoms

class Atom(Timeseries):
    def __init__(self, code, atom):
        if code not in ('x', 'y', 'z', 'v'):
            raise Exception("Bad code")
        if code == 'v': size = 3
        else: size = 1
        numatoms = len(atoms)
        Timeseries.__init__(self, code*numatoms, atoms, size*numatoms)

class Angle(Timeseries):
    def __init__(self, atoms):
        if not len(atoms) == 3:
            raise Exception("Angle timeseries requires a 3 atom selection")
        Timeseries.__init__(self, 'a', atoms, 1)

class Dihedral(Timeseries):
    def __init__(self, atoms):
        if not len(atoms) == 4:
            raise Exception("Dihedral timeseries requires a 4 atom selection")
        Timeseries.__init__(self, 'h', atoms, 1)

class Distance(Timeseries):
    def __init__(self, code, atoms):
        if code not in ('d', 'r'):
            raise Exception("Bad code")
        if code == 'd': size = 3
        else: size = 1
        if not len(atoms) == 2:
            raise Exception("Distance timeseries requires a 2 atom selection")
        Timeseries.__init__(self, code, atoms, size)

class CenterOfGeometry(Timeseries):
    def __init__(self, atoms):
        Timeseries.__init__(self, 'g', atoms, 3)

class CenterOfMass(Timeseries):
    def __init__(self, atoms):
        Timeseries.__init__(self, 'm', atoms, 3)
    def getAuxData(self):
        return [a.mass for a in self.atoms]
