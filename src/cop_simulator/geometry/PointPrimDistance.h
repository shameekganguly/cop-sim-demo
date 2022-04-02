// PointPrimDistance.h

#ifndef POINT_PRIM_DISTANCE_H
#define POINT_PRIM_DISTANCE_H

#include <Eigen/Dense>

#include "Primitive.h"
#include "Composite1PkN.h"

namespace Sai2COPSim {

struct PointPrimDistanceInfo {
	Eigen::Vector3d closest_prim_point;
	double distance;
};

// --------- distance computations. static functions ----------
class PointPrimDistance {
public:
	// returns contact info in the world frame
	static PointPrimDistanceInfo distancePointPrimitive(
		const Eigen::Vector3d& pointInWorld,
		const Primitive* prim, Eigen::Affine3d primInWorld
	);

	static PointPrimDistanceInfo distancePointPlane(
		const Eigen::Vector3d& pointInWorld,
		const PlanePrimitive& plane, Eigen::Affine3d primInWorld
	);

	static PointPrimDistanceInfo distancePointCapsule(
		const Eigen::Vector3d& pointInWorld,
		const CapsulePrimitive& capsule, Eigen::Affine3d primInWorld
	);

	static PointPrimDistanceInfo distancePointNegCapsule(
		const Eigen::Vector3d& pointInWorld,
		const NegCapsulePrimitive& neg_capsule, Eigen::Affine3d primInWorld
	);
};

}

#endif // POINT_PRIM_DISTANCE_H