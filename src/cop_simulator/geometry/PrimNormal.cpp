// PrimNormal.cpp

#include "PrimNormal.h"

namespace Sai2COPSim {

using namespace Eigen;

bool PrimNormal::primNormal(Eigen::Vector3d& ret_normal,
	const Eigen::Vector3d& pointOnPrimInWorld, const Primitive* prim,
	Eigen::Affine3d primInWorld
) {
	if(prim->_type == Primitive::GeometryType::Plane) {
		return planeNormal(
			ret_normal, pointOnPrimInWorld,
			*(dynamic_cast<const PlanePrimitive*>(prim)), primInWorld
		);
	} else if(prim->_type == Primitive::GeometryType::Capsule) {
		return capsuleNormal(
			ret_normal, pointOnPrimInWorld,
			*(dynamic_cast<const CapsulePrimitive*>(prim)), primInWorld
		);
	} else if(prim->_type == Primitive::GeometryType::NegCapsule) {
		return negCapsuleNormal(
			ret_normal, pointOnPrimInWorld,
			*(dynamic_cast<const NegCapsulePrimitive*>(prim)), primInWorld
		);
	} else {
		std::cerr << "Prim type: " << prim->_type << std::endl;
		throw(std::runtime_error("Normal computation at a point is not implemented for this primitive."));
	}
}

bool PrimNormal::planeNormal(
	Eigen::Vector3d& ret_normal,
	const Eigen::Vector3d& pointInWorld,
	const PlanePrimitive& plane, Eigen::Affine3d primInWorld
) {
	assert(plane._props != NULL);
	ret_normal = primInWorld.linear()*plane._props->normal;
	return true;
}

bool PrimNormal::capsuleNormal(
	Eigen::Vector3d& ret_normal,
	const Eigen::Vector3d& pointInWorld,
	const CapsulePrimitive& capsule, Eigen::Affine3d primInWorld
) {
	assert(capsule._props != NULL);
	Vector3d capsule_axis_world = primInWorld.linear().col(0);
	Vector3d capsule_center_world = primInWorld.translation();
	const double half_len = capsule._props->length*0.5;

	double zval = (pointInWorld - capsule_center_world).dot(capsule_axis_world);
	zval = fmin(half_len, fmax(-half_len, zval));
	Vector3d axis_pt = capsule_center_world + zval*capsule_axis_world;
	Vector3d dir = pointInWorld - axis_pt;
	if (dir.norm() < 1e-5) return false;
	ret_normal = dir / dir.norm();
	return true;
}

bool PrimNormal::negCapsuleNormal(
	Eigen::Vector3d& ret_normal,
	const Eigen::Vector3d& pointInWorld,
	const NegCapsulePrimitive& neg_capsule, Eigen::Affine3d primInWorld
) {
	assert(neg_capsule._props != NULL);
	Vector3d capsule_axis_world = primInWorld.linear().col(0);
	Vector3d capsule_center_world = primInWorld.translation();
	const double half_len = neg_capsule._props->length*0.5;

	double zval = (pointInWorld - capsule_center_world).dot(capsule_axis_world);
	zval = fmin(half_len, fmax(-half_len, zval));
	Vector3d axis_pt = capsule_center_world + zval*capsule_axis_world;
	Vector3d dir = axis_pt - pointInWorld;
	if (dir.norm() < 1e-5) return false;
	ret_normal = dir / dir.norm();
	return true;
}

}