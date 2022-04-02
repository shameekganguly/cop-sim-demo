// IntersectionEdgeFactory.cpp

#include "IntersectionEdge.h"

namespace Sai2COPSim {

using namespace Eigen;

std::vector<IntersectionEdge*> IntersectionEdgeFactory::computeIntersectionEdges(
		const Primitive* positivePrim, const Eigen::Affine3d& positivePrimInWorld,
		const Primitive* negativePrim, const Eigen::Affine3d& negativePrimInWorld) {
	if(negativePrim->_type == Primitive::GeometryType::NegCapsule
		&& positivePrim->_type == Primitive::GeometryType::Plane) {
		return computeIntersectionEdgesPlaneNegCapsule(
			*(dynamic_cast<const PlanePrimitive*>(positivePrim)), positivePrimInWorld,
			*(dynamic_cast<const NegCapsulePrimitive*>(negativePrim)), negativePrimInWorld
		);
	} else {
		std::cerr << "Positive prim type: " << positivePrim->_type << std::endl;
		std::cerr << "Negative prim type: " << negativePrim->_type << std::endl;
		throw(std::runtime_error("Computing intersection edges between these primitives is not implemented."));
	}
}

std::vector<IntersectionEdge*> IntersectionEdgeFactory::computeIntersectionEdgesPlaneNegCapsule(
	const PlanePrimitive& plane, const Eigen::Affine3d& planeInWorld,
	const NegCapsulePrimitive& negCapsule, const Eigen::Affine3d& negCapsuleInWorld) {
	assert(plane._props != NULL);

	Vector3d plane_point = plane._props->point;
	Vector3d plane_normal = plane._props->normal;

	Vector3d plane_point_world = planeInWorld*plane_point;
	Vector3d plane_normal_world = planeInWorld.linear()*plane_normal;

	assert(negCapsule._props != NULL);

	Vector3d capsule_axis_world = negCapsuleInWorld.linear().col(0);
	Vector3d capsule_center_world = negCapsuleInWorld.translation();
	const double radius = negCapsule._props->radius;
	const double half_len = negCapsule._props->length*0.5;

	if(plane_normal_world.dot(capsule_axis_world) < 0.999) {
		throw(std::runtime_error("Intersection edge cannot be computed unless negative capsule is normal to the plane"));
	}

	std::vector<IntersectionEdge*> ret_edges;
	double zval = fabs((plane_point_world - capsule_center_world).dot(capsule_axis_world));
	if(zval >= half_len + radius) {
		// no intersection
		return ret_edges;
	}
	double circle_radius = radius;
	if(zval > half_len) {
		circle_radius = sqrt(radius*radius - (zval - half_len)*(zval - half_len));
	}
	Vector3d circle_pt_parent_frame = planeInWorld.inverse()*(capsule_center_world + zval*capsule_axis_world);
	Vector3d circle_normal_parent_frame = planeInWorld.linear().transpose()*capsule_axis_world;
	ret_edges.push_back(new Circle3D(circle_radius, circle_pt_parent_frame, circle_normal_parent_frame));
	return ret_edges;
}

}