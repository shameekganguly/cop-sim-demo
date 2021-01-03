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

BoxPrimitive::BoxPrimitive(const std::string& name, double xlength, double ylength, double zlength)
{
	if(name.empty()) throw(std::runtime_error("Name cannot be empty"));
	if(xlength < 0.001) {
		std::cerr << xlength << std::endl;
		throw(std::runtime_error("xlength too small."));
	}
	if(ylength < 0.001) {
		std::cerr << ylength << std::endl;
		throw(std::runtime_error("ylength too small."));
	}
	if(zlength < 0.001) {
		std::cerr << zlength << std::endl;
		throw(std::runtime_error("zlength too small."));
	}
	double rat1 = xlength/ylength;
	double rat2 = xlength/zlength;
	double rat3 = ylength/zlength;
	if(rat1 < 0.01 || rat1 > 100) {
		std::cerr << rat1 << std::endl;
		throw(std::runtime_error("xlength/ylength out of bounds."));
	}
	if(rat2 < 0.01 || rat2 > 100) {
		std::cerr << rat2 << std::endl;
		throw(std::runtime_error("xlength/zlength out of bounds."));
	}
	if(rat3 < 0.01 || rat3 > 100) {
		std::cerr << rat3 << std::endl;
		throw(std::runtime_error("ylength/zlength out of bounds."));
	}
	_name = name;
	_type = GeometryType::Box;
	_props = new BoxProperties();
	_props->xlength = xlength;
	_props->ylength = ylength;
	_props->zlength = zlength;
}

void BoxPrimitive::setBoxProperties(BoxProperties* props) {
	if(props == NULL) {
		throw(std::runtime_error("Props is null."));
	}
	if(props->xlength < 0.001) {
		std::cerr << props->xlength << std::endl;
		throw(std::runtime_error("xlength too small."));
	}
	if(props->ylength < 0.001) {
		std::cerr << props->ylength << std::endl;
		throw(std::runtime_error("ylength too small."));
	}
	if(props->zlength < 0.001) {
		std::cerr << props->zlength << std::endl;
		throw(std::runtime_error("zlength too small."));
	}
	double rat1 = props->xlength/props->ylength;
	double rat2 = props->xlength/props->zlength;
	double rat3 = props->ylength/props->zlength;
	if(rat1 < 0.01 || rat1 > 100) {
		std::cerr << rat1 << std::endl;
		throw(std::runtime_error("xlength/ylength out of bounds."));
	}
	if(rat2 < 0.01 || rat2 > 100) {
		std::cerr << rat2 << std::endl;
		throw(std::runtime_error("xlength/zlength out of bounds."));
	}
	if(rat3 < 0.01 || rat3 > 100) {
		std::cerr << rat3 << std::endl;
		throw(std::runtime_error("ylength/zlength out of bounds."));
	}
	if(_props != NULL) {
		delete _props;
	}
	_props = props;
}

PyramidPrimitive::PyramidPrimitive(const std::string& name, uint num_sides_base, double length_base_side, double height)
{
	if(name.empty()) throw(std::runtime_error("Name cannot be empty"));
	if(num_sides_base < 3) {
		std::cerr << num_sides_base << std::endl;
		throw(std::runtime_error("num_sides_base must be >= 3."));
	}
	if(length_base_side < 0.001) {
		std::cerr << length_base_side << std::endl;
		throw(std::runtime_error("length_base_side too small."));
	}
	if(height < 0.001) {
		std::cerr << height << std::endl;
		throw(std::runtime_error("height too small."));
	}
	double rat = length_base_side/height;
	if(rat < 0.01 || rat > 100) {
		std::cerr << rat << std::endl;
		throw(std::runtime_error("length_base_side/height out of bounds."));
	}
	_name = name;
	_type = GeometryType::Pyramid;
	_props = new PyramidProperties();
	_props->num_sides_base = num_sides_base;
	_props->length_base_side = length_base_side;
	_props->height = height;

	computeInternalProperties();
	computeBasePoints();
}

void PyramidPrimitive::setPyramidProperties(PyramidProperties* props) {
	if(props == NULL) {
		throw(std::runtime_error("Props is null."));
	}
	if(props->num_sides_base < 3) {
		std::cerr << props->num_sides_base << std::endl;
		throw(std::runtime_error("num_sides_base must be >= 3."));
	}
	if(props->length_base_side < 0.001) {
		std::cerr << props->length_base_side << std::endl;
		throw(std::runtime_error("length_base_side too small."));
	}
	if(props->height < 0.001) {
		std::cerr << props->height << std::endl;
		throw(std::runtime_error("height too small."));
	}
	double rat = props->length_base_side/props->height;
	if(rat < 0.01 || rat > 100) {
		std::cerr << rat << std::endl;
		throw(std::runtime_error("length_base_side/height out of bounds."));
	}
	if(_props != NULL) {
		delete _props;
	}
	_props = props;

	computeInternalProperties();
	computeBasePoints();
}

void PyramidPrimitive::computeInternalProperties() {
	_included_angle = 2*M_PI/_props->num_sides_base;
	_incircle_radius = _props->length_base_side/2/tan(_included_angle/2);
	_side_face_angle = atan2(_incircle_radius, _props->height);
	_circum_radius = _props->length_base_side/2/sin(_included_angle/2);
	_side_edge_angle = atan2(_circum_radius, _props->height);
}

void PyramidPrimitive::computeBasePoints() {
	Vector3d vertex0(_circum_radius, 0, 0);
	for(uint i = 0; i < _props->num_sides_base; i++) {
		_base_points.push_back( Eigen::AngleAxisd(_included_angle*i, Vector3d::UnitZ()) * vertex0);
	}
}

}
