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

CylinderPrimitive::CylinderPrimitive(const std::string& name, double radius, double height, uint num_points)
{
	if(name.empty()) throw(std::runtime_error("Name cannot be empty"));
	if(radius < 0.001) {
		std::cerr << radius << std::endl;
		throw(std::runtime_error("Radius too small."));
	}
	if(height < radius*0.01 || height < 0.001) {
		std::cerr << height << std::endl;
		throw(std::runtime_error("Height too small."));
	}
	if(num_points < 3) {
		std::cerr << num_points << std::endl;
		throw(std::runtime_error("Num points too small."));
	}
	_name = name;
	_type = GeometryType::Cylinder;
	_props = new CylinderProperties();
	_props->radius = radius;
	_props->height = height;
	_props->num_points = num_points;

	computeFacePoints();
}

void CylinderPrimitive::setCylinderProperties(CylinderProperties* props) {
	if(props == NULL) {
		throw(std::runtime_error("Props is null."));
	}
	if(props->radius < 0.001) {
		std::cerr << props->radius << std::endl;
		throw(std::runtime_error("Radius too small."));
	}
	if(props->height < props->radius*0.01 || props->height < 0.001) {
		std::cerr << props->height << std::endl;
		throw(std::runtime_error("Height too small."));
	}
	if(props->num_points < 3) {
		std::cerr << props->num_points << std::endl;
		throw(std::runtime_error("Num points too small."));
	}
	if(_props != NULL) {
		delete _props;
	}
	_props = props;

	computeFacePoints();
}

void CylinderPrimitive::computeFacePoints() {
	_faceA_points.clear();
	_faceB_points.clear();
	// point one in both face is assumed to be on the X axis
	Vector3d point0FaceA(_props->radius, 0, 0);
	Vector3d faceAtoFaceB(0, 0, _props->height);
	double angle = (2*M_PI)/_props->num_points;
	Matrix3d rot;
	for(uint i = 0; i < _props->num_points; i++) {
		rot << cos(angle*i), -sin(angle*i), 0,
			   sin(angle*i), cos(angle*i),  0,
			   		0,				0,		1.0;
		_faceA_points.push_back(rot*point0FaceA);
		_faceB_points.push_back(rot*point0FaceA + faceAtoFaceB);
	}
}

}
