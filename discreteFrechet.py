import numpy as np
import matplotlib.pyplot as plt
import plotly as py
import plotly.graph_objs as go

# Euclidean distance.
def euc_dist(pt1,pt2):
	return np.sqrt((pt2[0]-pt1[0])*(pt2[0]-pt1[0])+(pt2[1]-pt1[1])*(pt2[1]-pt1[1]))

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
	ca = np.ones((len(P),len(Q)))
	ca = np.multiply(ca,-1)
	leash = _c(ca,len(P)-1,len(Q)-1,P,Q)

	return ca
P = [(1,1),(9,3),(3,7)]
Q = [(2,3),(1,7),(0,5)]
	
data = [go.Heatmap(z=frechetDist(P,Q))]

py.offline.iplot(data, filename='heatmap_frechet')




