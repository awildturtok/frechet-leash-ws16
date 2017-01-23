from matplotlib import pyplot as plt
from src.discreteFrechet import *
from typing import Tuple, List, Dict, Any
import math as math


Point = Tuple[float, float]
Curve = List[Point]
Matrix = Dict[Point, float]

'''
dists = frechetDistance Matrix
preds = predecessors of Points
'''

def unzip(zipped : List[Tuple[Any, Any]]) -> Tuple[List[Any], Tuple[List[Any]]] :
    unzipped = zip(*zipped)

    second = next(unzipped)
    first = next(unzipped)

    return first, second


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

    for (pi, pj) in neighbours:
        _val, mypreds = _lexicographic_path(ca, pi, pj, P, Q)
        if _val < val:
            preds = mypreds
            val = _val

    #todo change to bfs: recursion after loop
    return max([val + ca[i, j]]), [(i, j)] + preds


def lexicographic_path(P: Curve, Q: Curve) -> Tuple[float, Curve, np.matrix]:
    ca = np.multiply(np.ones((len(P), len(Q))), -1)
    dist, indices = _lexicographic_path(ca, 0, 0, P, Q)

    return dist, list(indices), ca


P = [(0,0), (1,2), (2, 3), (0,7)]
Q = [(0,0), (1,1), (1, 2), (1,1), (0,7)]

dist, indices, terrain = lexicographic_path(P, Q)


print(P)
print(Q)

print(dist)
print(indices)
print(terrain)

normalized_terrain = terrain / terrain.max(axis=1)[:,np.newaxis]

x, y = unzip(indices)

# plt.imshow(ca,'gray'),plt.title('ORIGINAL')
plt.plot(x, y, color="red")
plt.imshow(normalized_terrain, cmap='inferno', interpolation='bicubic')

plt.xlim([0, len(Q) - 1])
plt.ylim([0, len(P) - 1])

plt.show()
