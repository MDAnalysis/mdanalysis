from __future__ import division
import numpy as np

def select_power_of_two(n):
    """
    Select the closest i such that n<=2**i
    """
    current_exp = int(np.ceil(np.log2(n+1)))
    if n == 2**current_exp:
        n_fft = n
    if n < 2**current_exp:
        n_fft = 2**current_exp
    elif n > 2**current_exp:
        n_fft = 2**(current_exp+1)

    return n_fft


def autocorrelation_1d(data):
    """
    Compute the autocorrelation of a scalar time series.
    """

    N = len(data)
    n_fft = select_power_of_two(N)

    # Pad the signal with zeros to avoid the periodic images.

    R_data = np.zeros(2*n_fft)
    R_data[:N] = data

    F_data = np.fft.fft(R_data)

    result = np.fft.ifft(F_data*F_data.conj())[:N].real/(N-np.arange(N))

    return result[:N]


def correlation_1d(data1, data2):
    """
    Compute the correlation of two scalar time series.

    Args:
        data1 (array-like): Input time series of shape (N,)
        data2 (array-like): Input time series of shape (N,)

    Returns:
        : ndarray of shape (2*N-1,) with the correlation for
        "data1*data2[tau]" where tau is the lag in units of the timestep in the
        input data. The correlation is given from time -N to time N.
    """

    N = len(data1)
    assert N == len(data2)
    n_fft = select_power_of_two(N)

    # Pad the signal with zeros to avoid the periodic images.
    R_data1 = np.zeros(2*n_fft)
    R_data1[:N] = data1
    R_data2 = np.zeros(2*n_fft)
    R_data2[:N] = data2
    F_data1 = np.fft.fft(R_data1)
    F_data2 = np.fft.fft(R_data2)
    result = np.fft.ifft(F_data1.conj()*F_data2)
    positive_time = result[:N].real/(N-np.arange(N))
    negative_time = result[-N+1:][::-1].real/(N-1-np.arange(N-1))

    return np.concatenate((negative_time[::-1], positive_time))
