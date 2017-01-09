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
    dt = (float(p[0]) - float(q[0]), float(p[1]) - float(q[1]))
    return math.sqrt(dt[0] ** 2 + dt[1] ** 2)


def _lexicographic_path(ca: Matrix, i: int, j: int, P: Curve, Q: Curve) -> Tuple[float, Curve]:

    ca[i, j] = euc_dist(P[i], Q[j])

    if i == len(P) - 1 and j == len(Q) - 1:
        return ca[i, j], [(len(P)- 1, len(Q) -1)]

    # All neighbouring cells inside the grid on a monotone path.
    neighbours = [p
                  for p in [(i + 1, j), (i + 1, j + 1), (i, j + 1)]
                  if p[0] < len(P) and p[1] < len(Q)]

    val = float("+inf")
    preds = []
    _val = 0

    for (pi, pj) in neighbours:
        _val, mypreds = _lexicographic_path(ca, pi, pj, P, Q)
        if _val < val:
            preds = mypreds
            val = _val

    #todo change to bfs: recursion after loop
    return max([val + ca[i, j]]), preds + [(i, j)]


def lexicographic_path(P: Curve, Q: Curve) -> Tuple[float, Curve, np.matrix]:
    ca = np.multiply(np.ones((len(P), len(Q))), -1)
    dist, indices = _lexicographic_path(ca, 0, 0, P, Q)

    return dist, list(reversed(indices)), ca


P = [(0,0), (1,2), (2, 3), (0,7)]
Q = list(zip(np.arange(0, 4), np.zeros(4)))

dist, indices, terrain = lexicographic_path(P, Q)


print(P)
print(Q)

print(dist)
print(indices)
print(terrain)

"""
	from matlotlib import pyplot as plt
	plt.imshow(ca,'gray'),plt.title('ORIGINAL')
	
	plt.imshow(ca, cmap='hot', interpolation='nearest')
	plt.show()
	
	"""