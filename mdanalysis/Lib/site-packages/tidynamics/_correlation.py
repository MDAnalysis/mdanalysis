from __future__ import division
import numpy as np
from .core import autocorrelation_1d, correlation_1d
import itertools

def acf(data):
    """Autocorrelation of the input data using the Fast Correlation Algorithm.

    Computes the autocorrelation for all time lags in the input data. The numerical results for
    large lags contain fewer samples than for short lags and are not accurate. This is intrinsic to
    the computation and not a limitation of the algorithm.

    For D-dimensional time series, a sum is performed on the last dimension.

    Args:
        data (array-like): The input signal, of shape (N,) or (N,D).

    Returns:
        ndarray of shape (N,) with the autocorrelation for successive linearly
        spaced time delays

    """

    data = np.asarray(data)
    if data.ndim==1:
        return autocorrelation_1d(data)
    elif data.ndim>1:
        result = autocorrelation_1d(data[:,0])
        for j in range(1, data.shape[1]):
            result += autocorrelation_1d(data[:,j])
        return result

def msd(pos):
    """Mean-squared displacement (MSD) of the input trajectory using the Fast
    Correlation Algorithm.

    Computes the MSD for all possible time deltas in the trajectory. The numerical results for large
    time deltas contain fewer samples than for small time times and are less accurate. This is
    intrinsic to the computation and not a limitation of the algorithm.

    Args:
        pos (array-like): The input trajectory, of shape (N,) or (N,D).

    Returns:
        : ndarray of shape (N,) with the MSD for successive linearly spaced time
        delays.

    """

    pos = np.asarray(pos)
    if pos.ndim==1:
        pos = pos.reshape((-1,1))
    N = len(pos)
    rsq = np.sum(pos**2, axis=1)
    MSD = np.zeros(N, dtype=float)

    SAB = autocorrelation_1d(pos[:,0])
    for i in range(1, pos.shape[1]):
        SAB += autocorrelation_1d(pos[:,i])

    SUMSQ = 2*np.sum(rsq)

    m = 0
    MSD[m] = SUMSQ - 2*SAB[m]*N

    MSD[1:] = (SUMSQ - np.cumsum(rsq)[:-1] - np.cumsum(rsq[1:][::-1])) / (N-1-np.arange(N-1))
    MSD[1:] -= 2*SAB[1:]

    return MSD

def cross_displacement(pos):
    """Cross displacement of the components of the input trajectory.

    Args:
        pos (array-like): The input trajectory, of shape (N, D).

    Returns:
        : list of lists of times series, where the fist two indices [i][j]
        denote the coordinates for the cross displacement: "(Delta pos[:,i]) (Delta pos[:,j])".

    """

    pos = np.asarray(pos)
    if pos.ndim != 2:
        raise ValueError("Incorrect input data for cross_displacement")
    D = pos.shape[1]

    # Precompute the component-wise MSD
    split_msd = [msd(pos_i) for pos_i in pos.T]

    # Create list of lists for the output
    result = [[] for i in range(D)]
    for i, j in itertools.product(range(D), range(D)):
        result[i].append([])

    for i, j in itertools.product(range(D), range(D)):
        if i==j:
            result[i][j] = split_msd[i]
        else:
            sum_of_pos = msd(pos[:,i]+pos[:,j])
            result[i][j] = 0.5*(sum_of_pos - split_msd[i] - split_msd[j])

    return result

def correlation(data1, data2):
    """Correlation between the input data using the Fast Correlation Algorithm.

    For D-dimensional time series, a sum is performed on the last dimension.

    Args:
        data1 (array-like): The first input signal, of shape (N,) or (N,D).

        data2 (array-like): The first input signal, of equal shape as data1.

    Returns:
        : ndarray of shape (2*N-1,) with the correlation for
        "data1*data2[tau]" where tau is the lag in units of the timestep in the
        input data. The correlation is given from time -N to time N.

    """

    data1 = np.asarray(data1)
    data2 = np.asarray(data2)
    if data1.shape != data2.shape:
        raise ValueError('Incompatible shapes for data1 and data2')

    if data1.ndim==1:
        return correlation_1d(data1, data2)
    elif data1.ndim>1:
        result = correlation_1d(data1[:,0], data2[:,0])
        for j in range(1, data1.shape[1]):
            result += correlation_1d(data1[:,j], data2[:,j])
        return result
