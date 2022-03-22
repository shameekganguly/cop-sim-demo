// Composite1PkN.cpp

#include "Composite1PkN.h"

namespace Sai2COPSim {

PrimPrimContactInfo CircleIE::primDistance(Primitive* prim, const Eigen::Affine3d& primInParent) {
	if (prim->_type != Primitive::GeometryType::Capsule) {
		std::cerr << "Prim type for circleIE distance: " << prim->_type << std::endl;
		throw(std::runtime_error("Distance between circleIE and this primitive is not implemented."));
	}

	// TODO: move to different static function
	PrimPrimContactInfo ret_info;
	// TODO: compute
	return ret_info;
}

Composite1PkN::Composite1PkN(const std::string& name, Primitive* positivePrimitive) {
	if(name.empty()) throw(std::runtime_error("Name cannot be empty"));
	_name = name;
	_type = GeometryType::Composite1PkN;
	_positivePrimitive = positivePrimitive;
}

}