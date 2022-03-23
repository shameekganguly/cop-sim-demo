// Composite1PkN: composite shape with 1 positive primitive and k *non-intersecting*
// negative primitives

#ifndef COMPOSITE1PkN_H
#define COMPOSITE1PkN_H

#include "Primitive.h"
#include "PrimPrimContactInfo.h"

namespace Sai2COPSim {

// Parent is a Composite1PkN
struct IntersectionEdge {
	virtual PrimPrimContactInfo primDistance(Primitive* prim, const Eigen::Affine3d& primInParent) = 0;

	virtual ~IntersectionEdge() = default;
};

struct CircleIE: public IntersectionEdge {
	CircleIE(double in_radius, Eigen::Vector3d in_center, Eigen::Vector3d in_plane_normal)
		: radius(in_radius), center(in_center), plane_normal(in_plane_normal) {}

	PrimPrimContactInfo primDistance(Primitive* prim, const Eigen::Affine3d& primInParent) override;

	double radius;
	Eigen::Vector3d center;
	Eigen::Vector3d plane_normal;
protected:
	CircleIE() {} // disable default constructor
};

// TODO: add line segment chain intersection edge

// Parent is a Composite1PkN
struct NegativePrimitiveInfo {
	~NegativePrimitiveInfo() {
		delete prim;
	}

	Primitive* prim;
	std::vector<IntersectionEdge> intersection_edge;
	Eigen::Affine3d transform_in_parent;
};

// Shape coordinate system is same as the positive primitive.
class Composite1PkN: public Primitive {
public:
	// Takes ownership of positivePrimitive
	Composite1PkN(const std::string& name, Primitive* positivePrimitive);

	~Composite1PkN() {
		delete _positivePrimitive;
		_negativePrimitives.clear();
	}

	// Takes ownership of negativePrimitive
	bool addNegativePrimitive(Primitive* negativePrimitive, Eigen::Affine3d world) {
		// TODO: Handle transform and implement
		// Check if intersects with other negative prims
		// Compute intersections with surface of positive prim
		return true;
	}

public:
	Primitive* _positivePrimitive;
	std::vector<NegativePrimitiveInfo> _negativePrimitives;
};

}

#endif // COMPOSITE1PkN_H