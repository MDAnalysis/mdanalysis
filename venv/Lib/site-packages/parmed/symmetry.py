import numpy as np

class Symmetry:
    """Symmetry operations

    Parameters
    ----------
    tensor : np.array, shape=(3*n_symmetry, 4)
    """
    def __init__(self, tensor):
        self.data = np.asarray(tensor, dtype='f8')
        if self.data.ndim != 2 or self.data.shape[1] != 4:
            raise ValueError('Tensor shape must be (3*n_symmetry, 4)')
