// GeometryUtils.cpp

#include "GeometryUtils.h"
#include <numeric>

namespace Sai2COPSim {

Eigen::Vector2d GeometryUtils::Point2DFromPoint3DPlaneTransform (
	const Eigen::Vector3d& point3d, // assumed to lie on plane
	const Eigen::Affine3d& plane_tf // normal assumed to be in +z
) {
	return (plane_tf.inverse()*point3d).head(2);
}

Eigen::Vector3d GeometryUtils::centroidOfPoints(const std::vector<Eigen::Vector3d>& points) {
	if (points.size() == 0) return Eigen::Vector3d::Zero();
	Eigen::Vector3d init = Eigen::Vector3d::Zero();
	return std::accumulate(points.begin(), points.end(), init)/points.size();
}

double GeometryUtils::distancePointToPlane(
	const Eigen::Vector3d& test_point, // assumed to be in same coord system as plane
	const Eigen::Vector3d& plane_point,
	const Eigen::Vector3d& plane_normal
) {
	return (test_point - plane_point).dot(plane_normal);
}

};