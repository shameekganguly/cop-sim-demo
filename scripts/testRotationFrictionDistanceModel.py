import numpy as np
from matplotlib import pyplot as plt
import sys

xlen = 2
ylen = 1

xvec = xlen*np.linspace(0, 1, 100)
yvec = ylen*np.linspace(0, 1, 100)
xv, yv = np.meshgrid(xvec, yvec, indexing='ij')

v1 = [0, 0]
v2 = [0, ylen]
v3 = [xlen, ylen]
v4 = [xlen, 0]
e1 = [v1, v2]
e2 = [v2, v3]
e3 = [v3, v4]
e4 = [v1, v4]

edges = [e1,e2,e3,e4]

def distanceToVertex(pt, v):
	return np.linalg.norm(pt - np.array(v))

def edgeLength(e):
	temp = np.array(e[1]) - np.array(e[0])
	return np.linalg.norm(temp)

def distanceToEdge(pt, e):
	temp = np.array(e[1]) - np.array(e[0])
	elen = edgeLength(e)
	temp2 = np.array(e[1]) - pt
	temp3 = temp2 - temp*(np.dot(temp2, temp)/(elen**2))
	return np.linalg.norm(temp3)

def edgeNormal(e):
	

def isPtInPolygon(pt, edge_list, interior_pt):


# print xvec[xvec.size/2], yvec[0]
# reg = 1e-4
# pt = np.array([1, 0])
# lengths = []
# devec = []
# for edge in edges:
# 	dv0 = distanceToVertex(pt, edge[0])
# 	dv1 = distanceToVertex(pt, edge[1])
# 	de = distanceToEdge(pt, edge)
# 	lengths.append(dv0*dv1/(de + edgeLength(edge)))
# 	devec.append(de)
# weighted_lengths = []
# for k in range(len(edges)):
# 	ew = 1
# 	for l in range(len(edges)):
# 		if l != k:
# 			temp = 1/(reg + devec[l]) - 1/(reg + devec[k])
# 			if temp < 1e3:
# 				ew += np.exp(1/(reg + devec[l]) - 1/(reg + devec[k]))
# 			else:
# 				ew = np.inf
# 				break
# 	weighted_lengths.append(1.0/ew*lengths[k])
# print np.sum(weighted_lengths)

# sys.exit(0)

dist_map = np.zeros((xvec.size, yvec.size))
reg = 1e-4
for i in range(xvec.size):
	for j in range(yvec.size):
		pt = np.array([xvec[i], yvec[j]])
		lengths = []
		devec = []
		for edge in edges:
			dv0 = distanceToVertex(pt, edge[0])
			dv1 = distanceToVertex(pt, edge[1])
			de = distanceToEdge(pt, edge)
			lengths.append(dv0*dv1/(de + edgeLength(edge)))
			devec.append(de)
		weighted_lengths = []
		for k in range(len(edges)):
			ew = 1
			for l in range(len(edges)):
				if l != k:
					temp = 1/(reg + devec[l]) - 1/(reg + devec[k])
					if temp < 1e3:
						ew += np.exp(1/(reg + devec[l]) - 1/(reg + devec[k]))
					else:
						ew = np.inf
						break
			weighted_lengths.append(1.0/ew*lengths[k])
		dist_map[i,j] = np.sum(weighted_lengths)
# print dist_map
print('Max distance: ', dist_map[xvec.size/2, yvec.size/2])
fig = plt.figure()
CS = plt.contourf(xv, yv, dist_map, levels=np.linspace(0,dist_map[xvec.size/2, yvec.size/2],20))
plt.gca().set_aspect('equal', 'box')
CB = fig.colorbar(CS, shrink=0.8)
plt.show()