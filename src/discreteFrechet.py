import numpy as np
import plotly as py
import numpy.linalg as linalg

# Euclidean distance.
def euc_dist(pt1,pt2):
	return linalg.norm(pt1 - pt2)

def _c(ca,i,j,P,Q):
	if ca[i,j] > -1:
		return ca[i,j]
	elif i == 0 and j == 0:
		ca[i,j] = euc_dist(P[0],Q[0])
	elif i > 0 and j == 0:
		ca[i,j] = max(_c(ca,i-1,0,P,Q),euc_dist(P[i],Q[0]))
	elif i == 0 and j > 0:
		ca[i,j] = max(_c(ca,0,j-1,P,Q),euc_dist(P[0],Q[j]))
	elif i > 0 and j > 0:
		ca[i,j] = max(min(_c(ca,i-1,j,P,Q),_c(ca,i-1,j-1,P,Q),_c(ca,i,j-1,P,Q)),euc_dist(P[i],Q[j]))
	else:
		ca[i,j] = float("inf")
	return ca[i,j]

""" Computes the discrete frechet distance between two polygonal lines
Algorithm: http://www.kr.tuwien.ac.at/staff/eiter/et-archive/cdtr9464.pdf
P and Q are arrays of 2-element arrays (points)
"""
def frechetDist(P,Q):

	ca = np.multiply(np.ones((len(P),len(Q))), -1)
	frechetDistance= _c(ca,len(P)-1,len(Q)-1,P,Q)
	
	return (frechetDistance , ca)


_p = np.asarray([(1,0), (2,2), (-1,3), (2,4)])
_q = np.asarray([(-3, 2), (-2,0), (0.5, 1)])

py.offline.iplot(data, filename='heatmap_frechet')

frechet, dists = frechetDist(_p, _q)
