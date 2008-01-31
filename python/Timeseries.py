
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
        fmt = []
        for ts in self.timeseries:
            fmt.append(ts.getFormatCode())
        return ''.join(fmt)
    def getBounds(self):
        if not hasattr(self, '_atomlist'):
            self.getAtomList()
        lowerb = min(self._atomlist)
        upperb = max(self._atomlist)
        return lowerb, upperb
    def getDataSize(self):
        size = 0
        for ts in self.timeseries:
            size += ts.getDataSize()
        return size

class Timeseries:
    def __init__(self, code, numatoms, size):
        self.code = code
        self.numatoms = numatoms
        self.size = size
    def getAtomList(self):
        return [a.number for a in self.atoms]
    def getFormatCode(self):
        return self.code
    def getDataSize(self):
        return self.size

class Atom(Timeseries):
    def __init__(self, code, atom):
        if code not in ('x', 'y', 'z', 'v'):
            raise Exception("Bad code")
        if code == 'v': size = 3
        else: size = 1
        self.atoms = [atom]
        Timeseries.__init__(self, code, 1, size)

class Angle(Timeseries):
    def __init__(self, atoms):
        self.atoms = atoms
        Timeseries.__init__(self, 'a', 3, 1)

class Dihedral(Timeseries):
    def __init__(self, atoms):
        self.atoms = atoms
        Timeseries.__init__(self, 'h', 4, 1)

class Distance(Timeseries):
    def __init__(self, code, atoms):
        if code not in ('d', 'r'):
            raise Exception("Bad code")
        if code == 'd': size = 3
        else: size = 1
        self.atoms = atoms
        Timeseries.__init__(self, code, 2, size)
