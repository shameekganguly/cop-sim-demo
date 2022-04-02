// IntersectionEdge.cpp

#include "IntersectionEdge.h"

using namespace Eigen;

namespace Sai2COPSim {

PrimPrimContactInfo Circle3D::primDistance(const Primitive* prim, const Eigen::Affine3d& primInParent) const {
	if(prim->_type == Primitive::GeometryType::Capsule) {
		return Circle3DDistance::capsuleDist(*(dynamic_cast<const CapsulePrimitive*>(prim)), *this, primInParent);
	} else {
		std::cerr << "Prim type: " << prim->_type << std::endl;
		throw(std::runtime_error("Distance between Circle3D and this primitive is not implemented."));
	}
}

}