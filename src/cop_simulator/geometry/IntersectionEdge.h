// IntersectionEdge.h

#ifndef INTERSECTION_EDGE_H
#define INTERSECTION_EDGE_H

#include "Primitive.h"
#include "PrimPrimContactInfo.h"

namespace Sai2COPSim {

// Parent is a Composite1PkN
struct IntersectionEdge {
	virtual PrimPrimContactInfo primDistance(Primitive* prim, const Eigen::Affine3d& primInParent) = 0;

	virtual ~IntersectionEdge() = default;
};

struct Circle3D: public IntersectionEdge {
	Circle3D(double in_radius, Eigen::Vector3d in_center, Eigen::Vector3d in_plane_normal)
	: radius(in_radius), center(in_center), plane_normal(in_plane_normal) {
		circle_frame_in_parent.translation() = center;
		Eigen::Vector3d circle_x (1, 0, 0);
		if (circle_x.dot(plane_normal) > 0.999) {
			circle_x << 0, 1, 0;
		}
		circle_x -= (circle_x.dot(plane_normal))*plane_normal;
		circle_x /= circle_x.norm();
		circle_frame_in_parent.linear().col(0) = circle_x;
		circle_frame_in_parent.linear().col(1) = plane_normal.cross(circle_x);
		circle_frame_in_parent.linear().col(2) = plane_normal;
	}

	PrimPrimContactInfo primDistance(Primitive* prim, const Eigen::Affine3d& primInParent) override;

	double radius;
	Eigen::Vector3d center;
	Eigen::Vector3d plane_normal;

	// center is origin, Z axis is circle normal
	Eigen::Affine3d circle_frame_in_parent;

protected:
	Circle3D() {} // disable default constructor
};

// Specialized functions for distance from different primitives
class Circle3DDistance {
public:
	static PrimPrimContactInfo capsuleDist(const CapsulePrimitive& cap, const Circle3D& circle, const Eigen::Affine3d& primInParent);
};

// TODO: add line segment chain intersection edge

}

#endif // INTERSECTION_EDGE_H