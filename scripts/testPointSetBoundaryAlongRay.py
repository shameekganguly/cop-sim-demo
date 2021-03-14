import numpy as np
import numpy.random as rand
from scipy.spatial import ConvexHull

import matplotlib.pyplot as plt
import sys

NUM_POINTS = 40
NUM_POINT_SETS = 2
NUM_RAYS = 2

PT_INSIDE = 0
PT_OUTSIDE = 1
PT_ON_BDRY = 2

FLAG_PLOT = True

def projectPtsOnRay(test_pt, ray, points):
	diffs = points - test_pt
	projs = np.dot(diffs, ray)
	return projs

def partitionPtsByRay(test_pt, ray, points):
	ray_rot = np.array([-ray[1], ray[0]])
	diffs = points - test_pt
	projs = np.dot(diffs, ray_rot)
	return projs > 0

PARTITION_SIGN_NAN = -1e100
def argMaxWithPartitionSign(data, data_sign, positives=True):
	mod_data = data.copy()
	if positives:
		mod_data[~data_sign] = PARTITION_SIGN_NAN
	else:
		mod_data[data_sign] = PARTITION_SIGN_NAN
	ind = np.argmax(mod_data)
	return (ind, mod_data[ind])

def intersectLineSegRay(pt1, pt2, test_pt, ray):
	pt1_d = pt1 - test_pt
	pt2_d = pt2 - test_pt
	theta = np.arctan2(ray[1], ray[0])
	rot = np.array([[np.cos(theta), np.sin(theta)],[-np.sin(theta), np.cos(theta)]])
	pt1_d = np.dot(rot,pt1_d)
	pt2_d = np.dot(rot,pt2_d)
	# print abs(pt1_d[1] - pt2_d[1])
	if abs(pt1_d[1] - pt2_d[1]) < 1e-5:
		print "New code path4"
		# indicates that the test pt lies on the boundary edge and the ray dir is
		# parallel to the edge.
		line_seg = pt2 - pt1
		line_seg_dot_ray = np.dot(line_seg, ray)
		if line_seg_dot_ray > 0:
			return pt2
		else:
			return pt1
	alp = pt2_d[0] - pt2_d[1]*(pt1_d[0] - pt2_d[0])/(pt1_d[1] - pt2_d[1])
	return test_pt + np.dot(rot.T,np.array([alp, 0]))

# def implReturnTestPtOnBdryRayOutwards(test_pt, ray_projs, ray_partition_signs, points, iters):
# 	# determine is pt is vertex or not
# 	ind_pos = argMaxWithPartitionSign(ray_projs, ray_partition_signs, True)
# 	ind_neg = argMaxWithPartitionSign(ray_projs, ray_partition_signs, False)
# 	pos_dist = np.linalg.norm(points[ind_pos] - test_pt)
# 	neg_dist = np.linalg.norm(points[ind_neg] - test_pt)
# 	if pos_dist < 1e-5 or neg_dist < 1e-5: # indicates vertex
# 		return (PT_ON_BDRY, test_pt, test_pt, test_pt, iters)
# 	else:
# 		return (PT_ON_BDRY, test_pt, points[ind_pos], points[ind_neg], iters)

def ptSetDetectBoundaryPtAlongRay(test_pt, ray, points):
	proj1 = projectPtsOnRay(test_pt, ray, points)
	if max(proj1) < -1e-8 or min(proj1) > 1e-8:
		return (PT_OUTSIDE, None, None, None, 0)
	partition_signs = partitionPtsByRay(test_pt, ray, points)
	# TODO: check if proj1 max point is on the ray. Then the boundary point
	# is this vertex itself
	(ind_pos, pos_max) = argMaxWithPartitionSign(proj1, partition_signs, True)
	(ind_neg, neg_max) = argMaxWithPartitionSign(proj1, partition_signs, False)
	if abs(pos_max - PARTITION_SIGN_NAN) < 1e-10 or abs(neg_max - PARTITION_SIGN_NAN) < 1e-10:
		print "New code path1"
		# indicates that either test_pt is outside the convex hull or test_pt lies on a
		# bdry edge and the ray is parallel to the edge or test_pt lies on a bdry vertex
		# and the ray is outwards
		# we assume that the point is not outside
		pos_only = abs(neg_max - PARTITION_SIGN_NAN) < 1e-10
		if pos_only:
			ray_perp_dir = np.array([ray[1], -ray[0]])
		else:
			ray_perp_dir = np.array([-ray[1], ray[0]])
		perp_proj = projectPtsOnRay(test_pt, ray_perp_dir, points)
		mod_projs = proj1.copy()
		mod_projs[abs(perp_proj) > 1e-5] = PARTITION_SIGN_NAN
		ret_vertex_ind = np.argmax(mod_projs)
		if abs(mod_projs[ret_vertex_ind] - PARTITION_SIGN_NAN) < 1e-10:
			raise RuntimeError("Should not be here!")
		ret_pt = points[ret_vertex_ind]
		return (PT_ON_BDRY, ret_pt, ret_pt, ret_pt, 0)

	iters = 0
	MAX_ITERS = 10
	while True:
		iters += 1
		bdry_dir = points[ind_pos] - points[ind_neg]
		bdry_len = np.linalg.norm(bdry_dir)
		if bdry_len < 1e-5:
			print "New code path2"
			# in this case, either the boundary edge is ill-formed (super short) or
			# the ray intersects a boundary vertex, and there is is an interior point
			# that we are marking mistakenly as a boundary point.
			# assuming the edges are well formed, return the furthest point
			pos_dist = np.linalg.norm(points[ind_pos] - test_pt)
			neg_dist = np.linalg.norm(points[ind_neg] - test_pt)
			ret_pt = points[ind_pos] if pos_dist > neg_dist else points[ind_neg]
			return (PT_INSIDE, ret_pt, ret_pt, ret_pt, iters)
		bdry_dir /= bdry_len
		perp_dir = np.array([bdry_dir[1], -bdry_dir[0]])
		perp_proj = projectPtsOnRay(points[ind_neg], perp_dir, points)
		ind = np.argmax(perp_proj)
		max_perp_proj = perp_proj[ind]
		spread = max_perp_proj - np.min(perp_proj)
		if spread < 0.01*bdry_len:
			raise RuntimeError("Not a surface contact. Likely a line contact")
		if iters > MAX_ITERS:
			raise RuntimeError("Max iters exceeded")
		if max_perp_proj < 1e-5:
			# we have the right edge, but there is still some ambiguity in case
			# collinear points are present on this edge. we want to return the
			# true maximal points in order to transition to a nice line contact
			collinear_edge_inds = np.where(perp_proj > -1e-5)[0]
			assert collinear_edge_inds.size >= 2
			if collinear_edge_inds.size > 2:
				print "New code path3"
				edge_dir_proj = projectPtsOnRay(points[ind_neg], bdry_dir, points[collinear_edge_inds])
				ind_pos = collinear_edge_inds[np.argmax(edge_dir_proj)]
				ind_neg = collinear_edge_inds[np.argmin(edge_dir_proj)]
			break
		if partition_signs[ind] > 0:
			ind_pos = ind
		else:
			ind_neg = ind
	# print "Iters: ", iters
	ret_point = intersectLineSegRay(points[ind_pos], points[ind_neg], test_pt, ray)
	return (PT_INSIDE, ret_point, points[ind_pos], points[ind_neg], iters)

def convexHullBoundaryPtAlongRay(ray, hull):
	# TODO: this only works if test_pt is at 0,0
    eq=hull.equations.T
    V,b=eq[:-1],eq[-1]
    alpha=-b/np.dot(V.T,ray)
    return np.min(alpha[alpha>0])*ray

def bothSameSide(pt1, pt2, line_pts):
	line_perp = line_pts[1,:] - line_pts[0,:]
	line_perp = np.array([line_perp[1], -line_perp[0]])
	return np.dot(line_perp, pt1-line_pts[0,:])*np.dot(line_perp, pt2-line_pts[0,:]) > 0

def inHull(pt, centroid, hull):
	for s in hull.simplices:
		line_pts = hull.points[s,:]
		if not bothSameSide(pt, centroid, line_pts):
			return False
	return True

def postProcessResults(points, ray, test_pt, truth_pt, result_tuple, fplot):
	(result, bdry_pt, ept1, ept2, algo_iters) = result_tuple
	if fplot:
			plt.figure()
			plt.plot(points[:,0], points[:,1], '.m')
			plt.plot(test_pt[0], test_pt[1], 'oc')
			ray_pts = np.zeros((2,2))
			ray_pts[0,:] = test_pt
			ray_pts[1,:] = test_pt + ray * 0.5
			plt.plot(ray_pts[:,0], ray_pts[:,1], 'k')
	if result == PT_OUTSIDE:
		print "test point ", test_pt, "is outside the convex hull of the point set"
	else:
		if result == PT_ON_BDRY:
			print "test point ", test_pt, "is on the boundary of the convex hull of the point set"
		if fplot:
			plt.plot(bdry_pt[0], bdry_pt[1], 'og')
			edge_pts = np.zeros((2,2))
			edge_pts[0,:] = ept1
			edge_pts[1,:] = ept2
			plt.plot(edge_pts[:,0], edge_pts[:,1], '--k')
	# compare with convex hull solution
	if truth_pt:
		if np.linalg.norm(truth_pt - bdry_pt) > 1e-4:
			print "MISMATCH"
			print points
			print "ray ", ray
			sys.exit(0)

		if fplot:
			plt.plot(truth_pt[0], truth_pt[1], 'xr')
	if fplot:
		plt.show()

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

	print "Test 1"
	test_pt = np.array([0, 0], dtype=np.double)/SCALE
	ray = np.array([1, 1], dtype=np.double)
	ray = ray/np.linalg.norm(ray)
	result_tuple = ptSetDetectBoundaryPtAlongRay(test_pt, ray, points)
	postProcessResults(points, ray, test_pt, None, result_tuple, SHOULD_PLOT)

	print  "Test 2"
	test_pt = np.array([0, 0], dtype=np.double)/SCALE
	ray = np.array([1, 0], dtype=np.double)
	ray = ray/np.linalg.norm(ray)
	result_tuple = ptSetDetectBoundaryPtAlongRay(test_pt, ray, points)
	postProcessResults(points, ray, test_pt, None, result_tuple, SHOULD_PLOT)

	print  "Test 3"
	test_pt = np.array([0, 0], dtype=np.double)/SCALE
	ray = np.array([0, 1], dtype=np.double)
	ray = ray/np.linalg.norm(ray)
	result_tuple = ptSetDetectBoundaryPtAlongRay(test_pt, ray, points)
	postProcessResults(points, ray, test_pt, None, result_tuple, SHOULD_PLOT)

	print  "Test 4"
	test_pt = np.array([0, 0], dtype=np.double)/SCALE
	ray = np.array([-2, 1], dtype=np.double)
	ray = ray/np.linalg.norm(ray)
	result_tuple = ptSetDetectBoundaryPtAlongRay(test_pt, ray, points)
	postProcessResults(points, ray, test_pt, None, result_tuple, SHOULD_PLOT)

	print  "Test 5"
	test_pt = np.array([0, 0], dtype=np.double)/SCALE
	ray = np.array([-1, 0], dtype=np.double)
	ray = ray/np.linalg.norm(ray)
	result_tuple = ptSetDetectBoundaryPtAlongRay(test_pt, ray, points)
	postProcessResults(points, ray, test_pt, None, result_tuple, SHOULD_PLOT)

	print  "Test 6"
	test_pt = np.array([-1, 0.5], dtype=np.double)/SCALE
	ray = np.array([-2, 1], dtype=np.double)
	ray = ray/np.linalg.norm(ray)
	result_tuple = ptSetDetectBoundaryPtAlongRay(test_pt, ray, points)
	postProcessResults(points, ray, test_pt, None, result_tuple, SHOULD_PLOT)

	print  "Test 7"
	test_pt = np.array([-1, 0.5], dtype=np.double)/SCALE
	ray = np.array([2, -1], dtype=np.double)
	ray = ray/np.linalg.norm(ray)
	result_tuple = ptSetDetectBoundaryPtAlongRay(test_pt, ray, points)
	postProcessResults(points, ray, test_pt, None, result_tuple, SHOULD_PLOT)

	print  "Test 8"
	test_pt = np.array([1.9999, 0], dtype=np.double)/SCALE
	ray = np.array([2.0002, 1], dtype=np.double)
	ray = ray/np.linalg.norm(ray)
	result_tuple = ptSetDetectBoundaryPtAlongRay(test_pt, ray, points)
	postProcessResults(points, ray, test_pt, None, result_tuple, SHOULD_PLOT)

def main():
	point_set_iter = 0
	algo_iter_set = np.zeros(NUM_POINT_SETS*NUM_RAYS)
	algo_iter_set_ind = 0
	num_sets_skipped = 0
	while point_set_iter < NUM_POINT_SETS:
		point_set_iter += 1
		points = rand.rand(NUM_POINTS*2)
		points = (points - 0.5).reshape((NUM_POINTS, 2))
		rand_translation = (rand.rand(2) - 0.5)*4.0/5.0
		points += rand_translation

		hull = ConvexHull(points)
		test_pt = np.array([0, 0]) # has to be zeros for the hull method to work currently
		if not inHull(test_pt, np.average(points, 0), hull):
			# print "Test point is not in convex hull of points, skipping"
			num_sets_skipped += 1
			continue

		for rayi in range(NUM_RAYS):
			theta = rand.rand(1)[0]*2.0*np.pi
			ray = np.array([np.cos(theta), np.sin(theta)])
			# print "ray: ", ray

			try:
				result_tuple = ptSetDetectBoundaryPtAlongRay(test_pt, ray, points)
				algo_iters = result_tuple[4]
				algo_iter_set[algo_iter_set_ind] = algo_iters
				algo_iter_set_ind += 1
				if algo_iter_set_ind % 100 == 0:
					print algo_iter_set_ind
			except:
				print "Detect failed for points: "
				# print points
				print "Test pt: ", test_pt
				raise
			truth_pt = convexHullBoundaryPtAlongRay(ray, hull)
			postProcessResults(points, ray, test_pt, truth_pt, result_tuple, FLAG_PLOT)

	print "Num points: ", len(points)
	print "Max iters: ", max(algo_iter_set)
	print "Num sets skipped: ", num_sets_skipped

if __name__ == "__main__":
	testMain()
	# main()