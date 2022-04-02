// PrimNormal.h

#ifndef PRIM_NORMAL_H
#define PRIM_NORMAL_H

#include <Eigen/Dense>

#include "Primitive.h"
#include "Composite1PkN.h"

namespace Sai2COPSim {

// --------- normal computations. static functions ----------
class PrimNormal {
public:
	// Returns normal at the given point on the primitive in the world frame
	static bool primNormal(
		Eigen::Vector3d& ret_normal,
		const Eigen::Vector3d& pointOnPrimInWorld,
		const Primitive* prim, Eigen::Affine3d primInWorld
	);

	static bool planeNormal(
		Eigen::Vector3d& ret_normal,
		const Eigen::Vector3d& pointInWorld,
		const PlanePrimitive& plane, Eigen::Affine3d primInWorld
	);

	static bool capsuleNormal(
		Eigen::Vector3d& ret_normal,
		const Eigen::Vector3d& pointInWorld,
		const CapsulePrimitive& capsule, Eigen::Affine3d primInWorld
	);

	static bool negCapsuleNormal(
		Eigen::Vector3d& ret_normal,
		const Eigen::Vector3d& pointInWorld,
		const NegCapsulePrimitive& neg_capsule, Eigen::Affine3d primInWorld
	);
};

}

#endif // PRIM_NORMAL_H