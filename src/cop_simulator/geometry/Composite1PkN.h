// Composite1PkN: composite shape with 1 positive primitive and k *non-intersecting*
// negative primitives

#ifndef COMPOSITE1PkN_H
#define COMPOSITE1PkN_H

#include "Primitive.h"
#include "IntersectionEdge.h"

namespace Sai2COPSim {

// Parent is a Composite1PkN
struct NegativePrimitiveInfo {
	NegativePrimitiveInfo() = default;
	NegativePrimitiveInfo(NegativePrimitiveInfo&& other) {
		prim = other.prim;
		other.prim = NULL;
		transform_in_parent = std::move(other.transform_in_parent);
		intersection_edges.insert(intersection_edges.end(), other.intersection_edges.begin(), other.intersection_edges.end());
		other.intersection_edges.clear();
	}

	~NegativePrimitiveInfo() {
		if(prim != NULL) {
			delete prim;
		}
		for(uint i = 0; i < intersection_edges.size(); i++) {
			delete intersection_edges[i];
			intersection_edges[i] = NULL;
		}
	}

	Primitive* prim;
	std::vector<IntersectionEdge*> intersection_edges;
	Eigen::Affine3d transform_in_parent;

	// delete copy constructor to prevent bad memory release
	NegativePrimitiveInfo(const NegativePrimitiveInfo&) =delete;
	NegativePrimitiveInfo& operator=(const NegativePrimitiveInfo&) =delete;
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
	bool addNegativePrimitive(Primitive* negativePrimitive, Eigen::Affine3d primInPositivePrim);

public:
	Primitive* _positivePrimitive;
	std::vector<NegativePrimitiveInfo> _negativePrimitives;
};

}

#endif // COMPOSITE1PkN_H