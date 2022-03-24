// Composite1PkN: composite shape with 1 positive primitive and k *non-intersecting*
// negative primitives

#ifndef COMPOSITE1PkN_H
#define COMPOSITE1PkN_H

#include "Primitive.h"
#include "IntersectionEdge.h"

namespace Sai2COPSim {

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