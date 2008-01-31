"""
"""

import AtomGroup

class Analysis:
    def __init__(self,psffilename, trjfilename, skip):
        self.system = AtomGroup.System(psffilename, trjfilename)
        self.skip = skip
    def analyze_trj(self):
        raise NotImplementedError
    def analyze_timestep(self):
        raise NotImplementedError
    def run(self, verbose=False):
        # Get the name of the analysis function
        if not (hasattr(self, "analyze")):
            if not (hasattr(self, self.analysisfn)):
                raise SystemError("No analysis function defined")
            else:
                afunc = getattr(self, self.analysisfn)
        else:
            afunc = self.analyze
        for i, ts in enumerate(self.system._dcd):
            if (verbose): print "Analyzing timestep",i
            afunc(i)
