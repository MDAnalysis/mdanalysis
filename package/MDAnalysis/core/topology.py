"""The all singing all dancing new Topology system"""

class Topology(object):

    def __init__(self, *topologyattrs):

        # attach the TopologyAttrs
        for topologyattr in topologyattrs:
            self.__setattr__(topologyattr.attrname, topologyattr)
