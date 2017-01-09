import numpy as np
import numpy.linalg as linalg
import cv2 as cv
import numpy as np
from matplotlib import pyplot as plt

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

	"""P = [[1,2],[1,3],[1,4]]
	Q = [[2,1],[3,1],[4,1]]"""

	ca = np.multiply(np.ones((len(P),len(Q))), -1)
	frechetDistance= _c(ca,len(P)-1,len(Q)-1,P,Q)
	
	
	plt.imshow(ca,'gray')
	plt.title('ORIGINAL')
	
	"""plt.imshow(ca, cmap='hot', interpolation='nearest')"""
	plt.show()
	
	return (frechetDistance , ca)


"""
how to run the code
open ipython


	run C:/....
	
	_p = np.asarray([(1,1), (2,2), (3,3), (4,4)])
	_q = np.asarray([(1,0),(2,0),(3,0),(4,0)])
	frechet(_p,_q)
"""