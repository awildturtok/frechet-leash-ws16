import numpy as np
import numpy.linalg as linalg
from src.discreteFrechet import *
from typing import Tuple, List, Dict
import math as math

Point = Tuple[float, float]
Curve = List[Point]
Matrix = Dict[Tuple[int, int], float]

'''
dists = frechetDistance Matrix
preds = predecessors of Points
'''


def euc_dist(p: Point, q: Point) -> float:
    dt = (p[0] - p[0], p[1] - q[1])
    return math.sqrt(dt[0] ** 2 + dt[1] ** 2)


def _lexicographic_path(ca: Matrix, i: int, j: int, P: Curve, Q: Curve) -> Tuple[float, Curve]:
    if i == 0 and j == 0:
        ca[i, j] = euc_dist(P[0], Q[0])
        return ca[i, j], [(0, 0)]

    # All neighbouring cells inside the grid on a monotone path.
    neighbours = [p
                  for p in [(i - 1, j), (i - 1, j - 1), (i, j - 1)]
                  if p[0] >= 0 and p[1] >= 0]

    val = float("+inf")
    preds = []

    for (pi, pj) in neighbours:
        dist, preds = _lexicographic_path(ca, pi, pj, P, Q)
        if dist < val:
            val = dist

    dst = euc_dist(P[i], Q[j])

    ca[i, j] = max([dst, val])

    return ca[i, j], preds + [(i, j)]


def lexicographic_path(P: Curve, Q: Curve) -> Tuple[float, Curve, np.matrix]:
    ca = np.multiply(np.ones((len(P), len(Q))), -1)
    dist, indices = _lexicographic_path(ca, len(P) - 1, len(Q) - 1, P, Q)

    return dist, indices, ca


ps = list(zip(np.arange(0, 3), np.zeros(6)))
ps += reversed(ps)
qs = list(zip(np.arange(0, 6), np.arange(0, 6)))

dist, indices, terrain = lexicographic_path(ps, qs)


print(dist)
print(indices)
print(terrain)

