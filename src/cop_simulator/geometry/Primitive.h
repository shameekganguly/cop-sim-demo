// Primitive.h

#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include <string>
#include <vector>
#include <Eigen/Dense>

namespace Sai2COPSim {

// forward definitions of primitive properties
struct CapsuleProperties;  // used for positive as well as negative capsule volumes, type is different
struct PlaneProperties;
struct SphereProperties;
struct CylinderProperties;
struct BoxProperties;
struct PyramidProperties;

class Primitive {
public:
	enum GeometryType {
		Undefined,
		Plane,
		Sphere,
		Capsule,
		Cylinder,
		Box,
		Pyramid,
		NegCapsule,
		Composite1PkN
	};

public:
	// base constructor
	// TODO: add transform to body frame
	Primitive(const std::string& name) {
		if(name.empty()) throw(std::runtime_error("Name cannot be empty"));
		_type = GeometryType::Undefined;
		_name = name;
	}

	// base default constructor
	Primitive(): _articulated_body_name(""), _link_name(""), _is_static(true), _transform_in_link(Eigen::Affine3d::Identity())
	{
		// Nothing to do
	}

	// virtual capsule props
	virtual void setCapsuleProperties(CapsuleProperties* props) { }
	virtual CapsuleProperties* getCapsuleProperties() { return NULL; }

	// virtual plane props
	virtual void setPlaneProperties(PlaneProperties* props) { }
	virtual PlaneProperties* getPlaneProperties() { return NULL; }

	// virtual sphere props
	virtual void setSphereProperties(SphereProperties* props) { }
	virtual SphereProperties* getSphereProperties() { return NULL; }

	// virtual cylinder props
	virtual void setCylinderProperties(CylinderProperties* props) { }
	virtual CylinderProperties* getCylinderProperties() { return NULL; }

	// virtual box props
	virtual void setBoxProperties(BoxProperties* props) { }
	virtual BoxProperties* getBoxProperties() { return NULL; }

	// virtual pyramid props
	virtual void setPyramidProperties(PyramidProperties* props) { }
	virtual PyramidProperties* getPyramidProperties() { return NULL; }

public:
	// ownership info. TODO: consider something more robust
	std::string _articulated_body_name;
	std::string _link_name;
	bool _is_static; // belongs to world and not to a body

public:
	GeometryType _type;
	std::string _name;
	Eigen::Affine3d _transform_in_link;
};

// ------- specialized primitives --------
// TODO: move to individual files

// plane primitive
struct PlaneProperties {
	Eigen::Vector3d point;
	Eigen::Vector3d normal;
};

class PlanePrimitive: public Primitive {
public:
	PlanePrimitive(const std::string& name, const Eigen::Vector3d& planeNormal, const Eigen::Vector3d& planePoint);

	~PlanePrimitive() {
		delete _props;
	}

	virtual void setPlaneProperties(PlaneProperties* props);

	virtual PlaneProperties* getPlaneProperties() { return _props; }

public:
	PlaneProperties* _props; // owned by this primitive
};

// sphere primitive
// local frame at the center of the sphere.
struct SphereProperties {
	double radius;
};

class SpherePrimitive: public Primitive {
public:
	SpherePrimitive(const std::string& name, double radius);

	~SpherePrimitive() {
		delete _props;
	}

	virtual void setSphereProperties(SphereProperties* props);

	virtual SphereProperties* getSphereProperties() { return _props; }

public:
	SphereProperties* _props;
};


// capsule primitive
// capsule is assumed to be aligned with X axis. local frame at the center of the capsule.
struct CapsuleProperties {
	double radius;
	double length;
};

class CapsulePrimitive: public Primitive {
public:
	CapsulePrimitive(const std::string& name, double radius, double length);

	virtual ~CapsulePrimitive() {
		delete _props;
	}

	virtual void setCapsuleProperties(CapsuleProperties* props);

	virtual CapsuleProperties* getCapsuleProperties() { return _props; }

public:
	CapsuleProperties* _props;
};

class NegCapsulePrimitive: public CapsulePrimitive {
public:
	NegCapsulePrimitive(const std::string& name, double radius, double length);

	static NegCapsulePrimitive* CreateFromCapsule(const std::string& name, const CapsulePrimitive& capsule);
};

// cylinder primitive
// cylinder is assumed to be aligned with +Z axis to be consistent with chai.
// TODO: consider changing this in future.
// Also, local frame is at faceA of the cylinder.
struct CylinderProperties {
	double radius;
	double height;
	uint num_points; // number of points on a face for collision testing
};

class CylinderPrimitive: public Primitive {
public:
	CylinderPrimitive(const std::string& name, double radius, double height, uint num_points);

	~CylinderPrimitive() {
		delete _props;
	}

	virtual void setCylinderProperties(CylinderProperties* props);

	virtual CylinderProperties* getCylinderProperties() { return _props; }

public:
	CylinderProperties* _props;
	// FaceA and FaceB points in local frame
	std::vector<Eigen::Vector3d> _faceA_points;
	std::vector<Eigen::Vector3d> _faceB_points;

private:
	void computeFacePoints();
};

// box primitive
// box center is at 0. so x-extent is ± xlength/2
struct BoxProperties {
	double xlength;
	double ylength;
	double zlength;
};

class BoxPrimitive: public Primitive {
public:
	BoxPrimitive(const std::string& name, double xlength, double ylength, double zlength);

	~BoxPrimitive() {
		delete _props;
	}

	virtual void setBoxProperties(BoxProperties* props);

	virtual BoxProperties* getBoxProperties() { return _props; }

public:
	BoxProperties* _props;
};

// pyramid primitive
// pyramid base center is at 0. x-axis aligns with one base vertex
struct PyramidProperties {
	uint num_sides_base;
	double length_base_side;
	double height;
};

class PyramidPrimitive: public Primitive {
public:
	PyramidPrimitive(const std::string& name, uint num_sides_base, double length_base_side, double height);

	~PyramidPrimitive() {
		delete _props;
	}

	virtual void setPyramidProperties(PyramidProperties* props);

	virtual PyramidProperties* getPyramidProperties() { return _props; }

	double sideFaceAngle() const {
		return _side_face_angle;
	}

	double sideEdgeAngle() const {
		return _side_edge_angle;
	}

	double includedAngle() const {
		return _included_angle;
	}

	double circumRadius() const {
		return _circum_radius;
	}

	double incircleRadius() const {
		return _incircle_radius;
	}

public:
	PyramidProperties* _props;
	std::vector<Eigen::Vector3d> _base_points;

private:
	double _incircle_radius;
	double _circum_radius;
	double _side_face_angle;
	double _side_edge_angle;
	double _included_angle;
	void computeInternalProperties();
	void computeBasePoints();
};

}

#endif //PRIMITIVE_H