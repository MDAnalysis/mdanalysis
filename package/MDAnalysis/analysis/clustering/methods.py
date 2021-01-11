from . import cluster_methods

class AffinityPropagation:

    def __init__(self, preference=-1.0, damping=0.9, max_iter=500,
                 conv_threshold=50, add_noise=True):
        self.preference = preference
        self.damping = damping
        self.max_iter = max_iter
        self.conv_threshold = conv_threshold
        self.add_noise = add_noise
    
    def fit_predict(self, data):
        cl = cluster_methods.affinity_propagation(data, self.preference,
                                                  self.damping, self.max_iter,
                                                  self.conv_threshold,
                                                  int(self.add_noise))
        return cl