
from Scientific.Geometry import Transformation

class Transformation:
    def __init__(self):
        pass
    def transform(self):
        raise NotImplemented()

class Recenter(Transformation):
    def __init__(self, system, asel):
        self.system = system
        self.asel = asel
    def transform(self):
        com = self.asel.centerOfMass()
        self.system.coord.x -= com[0]
        self.system.coord.y -= com[1]
        self.system.coord.z -= com[2]

class RMSOrient(Transformation):
    def __init__(self, system, asel):
        self.system = system
        self.asel = asel
    def transform(self):
        # XXX Not complete yet
        com = self.asel.centerOfMass()

