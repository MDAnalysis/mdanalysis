"""The all singing all dancing new Topology system"""


class SizedDict(dict):
    """Dict that only allows things of a certain size to be added"""
    def __init__(self, size):
        self.size = size

    def __setitem__(self, key, item):
        if not len(item) == self.size:
            raise ValueError("Setting element with incorrect size")
        dict.__setitem__(self, key, item)


class Topology(object):
    def __init__(self, natoms, nresidues, nsegments)
        self.atoms = SizedDict(natoms)
        self.residues = SizedDict(nresidues)
        self.segments = SizedDict(nsegments)
