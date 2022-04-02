// Composite1PkN.cpp

#include "Composite1PkN.h"

namespace Sai2COPSim {

Composite1PkN::Composite1PkN(const std::string& name, Primitive* positivePrimitive) {
	if(name.empty()) throw(std::runtime_error("Name cannot be empty"));
	_name = name;
	_type = GeometryType::Composite1PkN;
	_positivePrimitive = positivePrimitive;
}

bool Composite1PkN::addNegativePrimitive(Primitive* negativePrimitive, Eigen::Affine3d primInPositivePrim) {
	assert(negativePrimitive != NULL);
	// TODO: Check if intersects with other negative prims
	// Compute intersections with surface of positive prim
	std::vector<IntersectionEdge*> i_edges = IntersectionEdgeFactory::computeIntersectionEdges(
		_positivePrimitive, Eigen::Affine3d::Identity(), negativePrimitive, primInPositivePrim);
	NegativePrimitiveInfo info;
	info.prim=negativePrimitive;
	info.intersection_edges=i_edges;
	info.transform_in_parent=primInPositivePrim;
	_negativePrimitives.push_back(std::move(info));
	return true;
}


}