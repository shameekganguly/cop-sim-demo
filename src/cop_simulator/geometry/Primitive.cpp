// Primitive.cpp

#include <iostream>
#include "Primitive.h"

using namespace Eigen;

namespace Sai2COPSim {

PlanePrimitive::PlanePrimitive(const std::string& name, const Eigen::Vector3d& planeNormal, const Eigen::Vector3d& planePoint)
{
	if(name.empty()) throw(std::runtime_error("Name cannot be empty"));
	if(abs(planeNormal.norm() - 1.0) > 1e-10) {
		std::cerr << planeNormal << std::endl;
		throw(std::runtime_error("Plane normal not unit vector"));
	}
	// // TODO: remove the below once the contact solver can handle non-Z normals
	// if(planeNormal.dot(Vector3d(0.0, 0.0, 1.0)) < 0.999999999) {
	// 	std::cerr << planeNormal << std::endl;
	// 	throw(std::runtime_error("Plane normal not +Z"));
	// }
	_name = name;
	_type = GeometryType::Plane;
	_props = new PlaneProperties();
	_props->point = planePoint;
	_props->normal = planeNormal;
}

void PlanePrimitive::setPlaneProperties(PlaneProperties* props) {
	if(props == NULL) {
		throw(std::runtime_error("Props is null."));
	}
	if(abs(props->normal.norm() - 1.0) > 1e-10) {
		std::cerr << props->normal << std::endl;
		throw(std::runtime_error("Plane normal not unit vector"));
	}
	if(_props != NULL) {
		delete _props;
	}
	_props = props;
}

CapsulePrimitive::CapsulePrimitive(const std::string& name, double radius, double length)
{
	if(name.empty()) throw(std::runtime_error("Name cannot be empty"));
	if(radius < 0.001) {
		std::cerr << radius << std::endl;
		throw(std::runtime_error("Radius too small."));
	}
	if(length < radius*0.01 || length < 0.001) {
		std::cerr << length << std::endl;
		throw(std::runtime_error("Length too small."));
	}
	_name = name;
	_type = GeometryType::Capsule;
	_props = new CapsuleProperties();
	_props->radius = radius;
	_props->length = length;
}

void CapsulePrimitive::setCapsuleProperties(CapsuleProperties* props) {
	if(props == NULL) {
		throw(std::runtime_error("Props is null."));
	}
	if(props->radius < 0.001) {
		std::cerr << props->radius << std::endl;
		throw(std::runtime_error("Radius too small."));
	}
	if(props->length < props->radius*0.01 || props->length < 0.001) {
		std::cerr << props->length << std::endl;
		throw(std::runtime_error("Length too small."));
	}
	if(_props != NULL) {
		delete _props;
	}
	_props = props;
}

}