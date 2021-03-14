// ContactPatch.cpp

#include <iostream>
#include "ContactPatch.h"

using namespace Eigen;

namespace Sai2COPSim {

Vector3d ContactPatch::patchToLineSegmentCoordinates(uint line_segment_id) {
	//TODO: implement
	return Vector3d::Zero();
}

PointTestResult ContactPatch::testPoint(const Eigen::Vector2d& point) const {
	PointTestResult ret_result;
	double curve_min_dist = std::numeric_limits<double>::max();

	// check each curve
	for(const auto curve: _intersection_curves) {
		double dist = curve->distanceToPoint(point);
		if(dist < curve_min_dist) {
			curve_min_dist = dist;
		}
	}
	if(abs(curve_min_dist) < ContactPatchAlgorithmicConstants::POINT_ON_BOUNDARY_DISTANCE_THRESHOLD) {
		ret_result.f_is_in_patch = false;
		ret_result.f_is_on_vertex = true;
		ret_result.min_dist_to_boundary = curve_min_dist;
		return ret_result;
	}

	// check line segments
	double line_min_dist = curve_min_dist;
	if(_line_segments.size() == 1) {
		line_min_dist = _line_segments[0].distanceToPoint(point);
		if(abs(line_min_dist) < ContactPatchAlgorithmicConstants::POINT_ON_BOUNDARY_DISTANCE_THRESHOLD) {
			ret_result.f_is_in_patch = false;
			ret_result.f_is_on_line_seg = true;
			ret_result.line_seg_id = 0;
			ret_result.min_dist_to_boundary = curve_min_dist;
			return ret_result;
		}
	} else if (_line_segments.size() > 1) {
		Vector2d min_2dists_line_segs;
		min_2dists_line_segs << _line_segments[0].distanceToPoint(point),
								_line_segments[1].distanceToPoint(point);
		uint line_seg_id = 0;
		if(min_2dists_line_segs[0] > min_2dists_line_segs[1]) {
			line_seg_id = 1;
			double temp = min_2dists_line_segs[1];
			min_2dists_line_segs[1] = min_2dists_line_segs[0];
			min_2dists_line_segs[0] = temp;
		}
		for(uint i = 2; i < _line_segments.size(); i++) {
			double dist = _line_segments[i].distanceToPoint(point);
			if(dist < min_2dists_line_segs[0]) {
				line_seg_id = i;
				double temp = min_2dists_line_segs[0];
				min_2dists_line_segs << dist, temp;
			} else if (dist < min_2dists_line_segs[1]) {
				min_2dists_line_segs[1] = dist;
			}
		}
		line_min_dist = min_2dists_line_segs[0];
		if(abs(line_min_dist) < ContactPatchAlgorithmicConstants::POINT_ON_BOUNDARY_DISTANCE_THRESHOLD) {
			// check for vertex
			if(abs(min_2dists_line_segs[1]) < ContactPatchAlgorithmicConstants::POINT_ON_BOUNDARY_DISTANCE_THRESHOLD) {
				ret_result.f_is_in_patch = false;
				ret_result.f_is_on_vertex = true;
				ret_result.min_dist_to_boundary = line_min_dist;
				return ret_result;
			} else {
				ret_result.f_is_in_patch = false;
				ret_result.f_is_on_line_seg = true;
				ret_result.line_seg_id = line_seg_id;
				ret_result.min_dist_to_boundary = line_min_dist;
				return ret_result;
			}
		}
	}
	// return inside patch
	ret_result.min_dist_to_boundary = fmin(curve_min_dist, line_min_dist);
	ret_result.f_is_in_patch = (ret_result.min_dist_to_boundary > 0);
	return ret_result;
}

double ContactPatch::distanceFromBoundaryAlongRay(const Vector2d& point, const Vector2d& direction) const {
	double ret_dist = std::numeric_limits<double>::max();
	// check each curve
	for(const auto curve: _intersection_curves) {
		double dist = curve->distanceFromPointAlongRay(point, direction);
		if(dist > 1e-10 && dist < ret_dist) {
			ret_dist = dist;
		}
	}
	// NOTE: we use the 1e-10 THRESHOLD above to avoid numerical issues
	// The justification is that this call is only made in one of two cases:
	// 1. The point is in the interior of the patch. In this case, we can
	//		safely assume that the minimum distance to the boundary in the direction
	//		of the ray will be atleast THRESHOLD
	// 2. The point is on the contact patch, but the COP solver has forced a
	//		PatchCenter solution. As a result, due to numerical issues, it is
	//		possible that the distance to the line segment on which the test
	//		point is currently located is actually a small positive number < THRESHOLD
	//		e.g. 1e-17. Therefore, we should ignore this distance and look for a
	//		boundary curve/segment with a larger +ve distance along the ray
	// check each line segments
	for(const auto& lineseg: _line_segments) {
		double dist = lineseg.distanceFromPointAlongRay(point, direction);
		// std::cout << dist << " ";
		if(dist > 1e-10 && dist < ret_dist) {
			ret_dist = dist;
		}
	}
	// std::cout << std::endl;
	return ret_dist;
}


}