import numpy as np
import scipy

EPS = 1e-15
xlogy = np.vectorize(
    lambda x, y: 0.0 if (x <= EPS and y <= EPS) else x * np.log(y))


def discrete_kl_div(a, b):
    return np.sum(xlogy(a, a / b))
    # return scipy.special.rel_entr(a, b).sum()


def discrete_js_div(a, b):
    ab = (a+b) * 0.5
    mask = ab > 0
    a = a[mask]
    b = b[mask]
    ab = ab[mask]
    return 0.5 * (discrete_kl_div(a, ab) + discrete_kl_div(b, ab))


# def discrete_kl_matrix(matrix):
#     log = np.log(matrix)
#     cross_entropy = -(matrix @ log.T)
#     # print(cross_entropy)
#     entropy = -np.einsum('ij,ij->i', matrix, log)
#     # print(entropy)
#     kl = cross_entropy - entropy[..., None]
#     # kl[np.isnan(kl)] = 0
#     return kl


def discrete_kl_matrix(matrix):
    x, y = matrix.shape
    matrix = np.array([matrix] * x)
    div = matrix / matrix.transpose((1, 0, 2))
    # div[~np.isfinite(div)] = 0
    log = np.log(div)
    log[~np.isfinite(div)] = 0
    log[log == -np.inf] = 0
    div[~np.isfinite(div)] = 0
    div[div == -np.inf] = 0
    
    # logab = np.broadcast_to(log, (x, x, y))
    # logdiff = logab - logab.transpose((1, 0, 2))
    prod = matrix * log

    # prod[(div < EPS) & (matrix < EPS)] = 0
    # prod[prod < EPS] = 0
    return prod.sum(axis=-1)
    






def discrete_js_matrix(matrix):
    value = np.zeros((len(matrix), len(matrix)))
    for i, row1 in enumerate(matrix):
        for j, row2 in enumerate(matrix[i:], i):
            value[i][j] = value[j][i] = discrete_js_div(row1, row2)
    return value
    # x = matrix.shape[0]
    # matrix = np.array([matrix] * x)
    # matrixT = matrix.transpose((1, 0, 2))
    # half = (matrix + matrixT) * 0.5
    
    # pB = matrix / half
    # loghalf = np.log(pB)
    # kl = (matrix * loghalf)
    # mask = (half < EPS) | (matrix < EPS) | (pB < EPS)
    # kl[mask] = 0
    # # kl[half < EPS] = 0
    # klsum = kl.sum(axis=-1)
    # print(klsum.shape)
    # print(klsum)

    # # kl = discrete_kl_matrix(matrix)
    # return 0.5 * (klsum + klsum.T)


def max_likelihood_covariance(coordinates, center=True):
    if center:
        coordinates -= coordinates.mean(axis=0)

    cov = np.einsum('ij,ik->jk', coordinates, coordinates)
    cov /= coordinates.shape[0]
    return cov


def shrinkage_covariance(coordinates, shrinkage_parameter=None):
    offset = coordinates - coordinates.mean(axis=0)
    t, n = coordinates.shape
    xmkt = offset.mean(axis=1)[..., None]
    coords = np.hstack([offset, xmkt])
    sample = max_likelihood_covariance(coords, False)
    sample *= (t-1)/float(t)

    # Split covariance matrix into components
    covmkt = sample[:n, n]
    varmkt = sample[n, n]
    sample = sample[:n, :n]

    covouter = np.outer(covmkt, covmkt)
    prior = covouter / varmkt
    dix = np.diag_indices(n)
    prior[dix] = sample[dix]

    if shrinkage_parameter is None:
        c = np.linalg.norm(sample - prior, ord="fro") ** 2
        inv_t = 1/len(offset)
        y = offset ** 2
        p = inv_t * (y.T @ y).sum() - (sample ** 2).sum()
        rdiag = inv_t * (y ** 2).sum() - np.trace(sample ** 2)
        z = offset * np.repeat(xmkt, n, axis=1)
        covmkt_ = np.repeat(covmkt[..., None], n, axis=1)
        v1 = inv_t * (y.T @ z) - covmkt_ * sample
        roff1 = ((v1 * covmkt_.T).sum() - np.trace(v1*covmkt)) / varmkt
        v3 = inv_t * (z.T @ z) - (varmkt * sample)
        roff3 = ((v3 * covouter).sum() - np.trace(v3 * covmkt ** 2))
        roff3 /= varmkt ** 2
        roff = 2 * roff1 - roff3
        r = rdiag + roff
        k = (p - r) / c
        shrinkage_parameter = max(0, min(1, k * inv_t))
    
    return shrinkage_parameter * prior + (1 - shrinkage_parameter) * sample


def get_bootstrap_frames(frames, n_samples=100):
    n = len(frames)
    indices = [np.random.randit(0, high=n, size=n) for i in range(n_samples)]
    return [frames[x] for x in indices]