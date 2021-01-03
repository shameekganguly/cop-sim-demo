// GeometryUtils.h
#ifndef GEOMETRY_UTILS_H
#define GEOMETRY_UTILS_H

#include "ContactPatch.h"

namespace Sai2COPSim {

// class to encapsulate static utility functions
// TODO: move to individual primitives that are affected by this
class GeometryUtils {
public:
	static Eigen::Vector2d Point2DFromPoint3DPlaneTransform (
		const Eigen::Vector3d& point3d, // assumed to lie on plane
		const Eigen::Affine3d& plane_tf
	);

	static Eigen::Vector3d centroidOfPoints(const std::vector<Eigen::Vector3d>& points);

	static double distancePointToPlane(
		const Eigen::Vector3d& test_point, // assumed to be in same coord system as plane
		const Eigen::Vector3d& plane_point,
		const Eigen::Vector3d& plane_normal
	);
};
} // namespace Sai2COPSim

#endif //GEOMETRY_UTILS_H