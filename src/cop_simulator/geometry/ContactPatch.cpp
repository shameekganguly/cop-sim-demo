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
	double min_dist = _intersection_curves[0]->distanceToPoint(point); // NOTE: we assume convex contact patch set
	// TODO: check each line segment

	// check each curve
	for(const auto curve: _intersection_curves) {
		double dist = curve->distanceToPoint(point);
		if(dist < min_dist) {
			min_dist = dist;
		}
	}
	ret_result.min_dist_to_boundary = min_dist;
	if(abs(min_dist) < ContactPatchAlgorithmicConstants::POINT_ON_BOUNDARY_DISTANCE_THRESHOLD) {
		ret_result.f_is_in_patch = false;
		ret_result.f_is_on_vertex = true;
	} else {
		ret_result.f_is_in_patch = (min_dist > 0);
	}
	return ret_result;
}

double ContactPatch::distanceFromBoundaryAlongRay(const Vector2d& point, const Vector2d& direction) const {
	// TODO: add line segment tests
	// check each curve
	for(const auto curve: _intersection_curves) {
		double dist = curve->distanceFromPointAlongRay(point, direction);
		if(dist > 0) {
			// since we assume that patch is convex, only one boundary curve can be
			// intersected
			return dist;
		}
	}
	return -1;
}


}