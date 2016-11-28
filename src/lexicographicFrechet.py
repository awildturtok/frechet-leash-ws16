import numpy as np
import numpy.linalg as linalg

'''
	dists = frechetDistance Matrix
	preds = predecessors of Points
'''
def _lexicographicPath(dists, i, j):
	
	if i == 0 and j == 0:
		return [(0, 0)]

	neighbours = [(p,dists[p]) for p in [(i-1,j), (i-1,j-1),(i,j-1)] if p[0] >= 0 and p[1] >= 1]
	(pi, pj), dist = neighbours[np.array([p[1] for p in neighbours]).argmin()]

	return _lexicographicPath(dists, pi, pj).append((i, j))
