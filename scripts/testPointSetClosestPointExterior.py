# testPointSetClosestPointExterior.py

import numpy as np
from matplotlib import pyplot as plt
from testPointSetBoundaryAlongRay import ptSetDetectBoundaryPtAlongRay, intersectLineSegRay

# test pt is assumed to be OUTSIDE the convex hull of the planar points
def ptSetDetectBoundaryClosestPt(test_pt, points):
	# first detect the closest point
	ind_closest = np.argmin(np.linalg.norm(points - test_pt, 2, 1))
	inside_pt = points[ind_closest]
	ray = test_pt - inside_pt # pointing from inside pt to test_pt
	ray_dist = np.linalg.norm(ray)
	if ray_dist < 1e-5:
		return (inside_pt, inside_pt, inside_pt)
	ray /= ray_dist
	(status, ret_pt, bdry_pt1, bdry_pt2, _) = ptSetDetectBoundaryPtAlongRay(inside_pt, ray, points)
	# PT_ON_BDRY indicates that the returned edge is just a boundary vertex
	# in this case, there are two possibilities. Either the returned point
	# is the closest to the test point, or the closest point lies on a boundary
	# edge which has the returned point as a vertex
	if np.linalg.norm(ret_pt - bdry_pt1) < 1e-5 or np.linalg.norm(ret_pt - bdry_pt2) < 1e-5:
		# check if least subtended angle is acute, then this is the other bdry point
		diffs = points - ret_pt
		# ignore points clustered around ret_pt
		norms = np.linalg.norm(diffs, 2, 1)
		inds = np.where(norms < 1e-5)[0]
		norms[inds] = 1e8
		diffs /= norms[:,None]
		projs = np.dot(diffs, ray)
		projs[inds] = -1
		ind_max = np.argmax(projs)
		if projs[ind_max] > 0:
			bdry_pt1 = ret_pt
			# there might be collinear points, so take the point furthest from ret_pt
			# to be the second boundary point
			collinear_inds = np.where(projs > (projs[ind_max] - 1e-6))[0]
			bdry_pt2 = points[collinear_inds[np.argmax(norms[collinear_inds])]]
		else:
			return (ret_pt, ret_pt, ret_pt)
	# find perpendicular point
	edge_dir = bdry_pt2 - bdry_pt1
	if np.dot(edge_dir, ray) < 0:
		temp = bdry_pt2
		bdry_pt2 = bdry_pt1
		bdry_pt1 = temp
		edge_dir *= -1
	edge_dir /= np.linalg.norm(edge_dir)
	perp_edge_dir = np.array([-edge_dir[1], edge_dir[0]])
	ret_pt = intersectLineSegRay(bdry_pt1, bdry_pt2, test_pt, -perp_edge_dir)
	return (ret_pt, bdry_pt1, bdry_pt2)

def testMain():
	SCALE = 10.0
	SHOULD_PLOT = False
	# we work with a predefined set of points here
	points = np.array([
		[0, 0],
		[1, 0],
		[2, 0],
		[0, 2],
		[1, 2],
		[2, 2],
		[-2,1],
		[4,1],
	], dtype=np.double)/SCALE

	test_pt = np.array([0.1,0.3])
	(ret_pt, bdry1, bdry2) = ptSetDetectBoundaryClosestPt(test_pt, points)
	plt.figure()
	plt.plot(points[:,0], points[:,1], '.m')
	plt.plot(test_pt[0], test_pt[1], 'or')
	plt.plot(ret_pt[0], ret_pt[1], 'og')
	edge_pts = np.zeros((2,2))
	edge_pts[0,:] = bdry1
	edge_pts[1,:] = bdry2
	plt.plot(edge_pts[:,0], edge_pts[:,1], '--k')
	plt.gca().set_aspect('equal')
	plt.show()

if __name__ == "__main__":
	testMain()