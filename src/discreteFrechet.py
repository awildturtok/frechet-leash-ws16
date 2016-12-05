import numpy as np
import numpy.linalg as linalg


# Euclidean distance.
def euc_dist(p, q):
    dt = (p[0] - p[0], p[1] - q[1])
    return np.sqrt(dt[0] ** 2 + dt[1] ** 2)


def _frechet_dist(ca, i, j, P, Q):
    if ca[i, j] > -1:
        return ca[i, j]
    elif i == 0 and j == 0:
        ca[i, j] = euc_dist(P[0], Q[0])
    elif i > 0 and j == 0:
        ca[i, j] = max(_frechet_dist(ca, i - 1, 0, P, Q), euc_dist(P[i], Q[0]))
    elif i == 0 and j > 0:
        ca[i, j] = max(_frechet_dist(ca, 0, j - 1, P, Q), euc_dist(P[0], Q[j]))
    elif i > 0 and j > 0:
        ca[i, j] = max(min(_frechet_dist(ca, i - 1, j, P, Q), _frechet_dist(ca, i - 1, j - 1, P, Q), _frechet_dist(ca, i, j - 1, P, Q)),
                       euc_dist(P[i], Q[j]))
    else:
        ca[i, j] = float("inf")
    return ca[i, j]


""" Computes the discrete frechet distance between two polygonal lines
Algorithm: http://www.kr.tuwien.ac.at/staff/eiter/et-archive/cdtr9464.pdf
P and Q are arrays of 2-element arrays (points)
"""


def frechet_dist(P, Q):
    ca = np.multiply(np.ones((len(P), len(Q))), -1)

    dist = _frechet_dist(ca, len(P) - 1, len(Q) - 1, P, Q)

    return dist, ca
