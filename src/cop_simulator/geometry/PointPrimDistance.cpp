// PointPrimDistance.cpp

#include "PointPrimDistance.h"

namespace Sai2COPSim {

using namespace Eigen;

PointPrimDistanceInfo PointPrimDistance::distancePointPrimitive(
	const Vector3d& pointInWorld,
	const Primitive* prim, Eigen::Affine3d primInWorld
) {
	if(prim->_type == Primitive::GeometryType::Plane) {
		return distancePointPlane(
			pointInWorld,
			*(dynamic_cast<const PlanePrimitive*>(prim)), primInWorld
		);
	} else if(prim->_type == Primitive::GeometryType::Capsule) {
		return distancePointCapsule(
			pointInWorld,
			*(dynamic_cast<const CapsulePrimitive*>(prim)), primInWorld
		);
	} else if(prim->_type == Primitive::GeometryType::NegCapsule) {
		return distancePointNegCapsule(
			pointInWorld,
			*(dynamic_cast<const NegCapsulePrimitive*>(prim)), primInWorld
		);
	} else {
		std::cerr << "Prim type: " << prim->_type << std::endl;
		throw(std::runtime_error("Distance between point and this primitive is not implemented."));
	}
}

PointPrimDistanceInfo PointPrimDistance::distancePointPlane(
	const Vector3d& pointInWorld,
	const PlanePrimitive& plane, Eigen::Affine3d primInWorld
) {
	PointPrimDistanceInfo ret_info;
	assert(plane._props != NULL);

	Vector3d plane_point = plane._props->point;
	Vector3d plane_normal = plane._props->normal;

	Vector3d plane_point_world = primInWorld*plane_point;
	Vector3d plane_normal_world = primInWorld.linear()*plane_normal;

	double dist = (pointInWorld - plane_point).dot(plane_normal_world);
	ret_info.closest_prim_point = pointInWorld - dist*plane_normal_world;
	ret_info.distance = dist;
	return ret_info;
}

PointPrimDistanceInfo PointPrimDistance::distancePointCapsule(
	const Vector3d& pointInWorld,
	const CapsulePrimitive& capsule, Eigen::Affine3d primInWorld
) {
	PointPrimDistanceInfo ret_info;
	assert(capsule._props != NULL);

	Vector3d capsule_axis_world = primInWorld.linear().col(0);
	Vector3d capsule_center_world = primInWorld.translation();
	const double radius = capsule._props->radius;
	const double half_len = capsule._props->length*0.5;

	double zval = (pointInWorld - capsule_center_world).dot(capsule_axis_world);
	zval = fmin(half_len, fmax(-half_len, zval));
	Vector3d axis_pt = capsule_center_world + zval*capsule_axis_world;
	Vector3d dir = pointInWorld - axis_pt;
	double center_dist = dir.norm();
	ret_info.distance = center_dist - radius;
	if(center_dist > 1e-5) {
		ret_info.closest_prim_point = axis_pt + dir/center_dist*radius;
	} else {
		ret_info.closest_prim_point = axis_pt + primInWorld.linear().col(1)*radius;
	}
	return ret_info;
}

PointPrimDistanceInfo PointPrimDistance::distancePointNegCapsule(
	const Vector3d& pointInWorld,
	const NegCapsulePrimitive& neg_capsule, Eigen::Affine3d primInWorld
) {
	PointPrimDistanceInfo ret_info;
	assert(neg_capsule._props != NULL);

	Vector3d capsule_axis_world = primInWorld.linear().col(0);
	Vector3d capsule_center_world = primInWorld.translation();
	const double radius = neg_capsule._props->radius;
	const double half_len = neg_capsule._props->length*0.5;

	double zval = (pointInWorld - capsule_center_world).dot(capsule_axis_world);
	zval = fmin(half_len, fmax(-half_len, zval));
	Vector3d axis_pt = capsule_center_world + zval*capsule_axis_world;
	Vector3d dir = pointInWorld - axis_pt;
	double center_dist = dir.norm();
	ret_info.distance = radius - center_dist; // ONLY DIFFERENCE FROM POSITIVE CAPSULE
	if(center_dist > 1e-5) {
		ret_info.closest_prim_point = axis_pt + dir/center_dist*radius;
	} else {
		ret_info.closest_prim_point = axis_pt + primInWorld.linear().col(1)*radius;
	}
	return ret_info;
}

}