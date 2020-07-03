// Primitive.h

#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include <string>
#include <vector>
#include <Eigen/Dense>

namespace Sai2COPSim {

// Constants used
namespace PrimitiveAlgorithmicConstants {
	const double MULTI_POINT_HIGHER_PAIR_CONTACT_DISTANCE_DIFF_THRESHOLD = 6e-5; //m
	const double MIN_HIGHER_PAIR_CONTACT_EXTENT_ANY_DIR = 0.001; //m
};

enum ContactType {
	UNDEFINED,
	POINT,
	LINE,
	SURFACE
};

class PrimPrimContactInfo {
public:
	double min_distance;
	Eigen::Vector3d normal_dir; // in world frame. 
	// ^ Note that this is directed from prim A to prim B, in the order that the 
	// prim-prim distance computation function was called
	Eigen::Vector3d constraint_dir1; // in world frame. for line contacts and surface contacts
	Eigen::Vector3d constraint_dir2; // in world frame. for surface contacts only
	// TODO: for concave objects, the each contact point can have its own normal and 
	// constraint dirs. Consider if we want to support that, or assume that each primitive
	// is strictly convex
	std::vector<Eigen::Vector3d> contact_points; // closest points on either prim A or prim B in world frame
	//TODO: do we need to be consistent about which primitive the points lie on?
	// or do we need to return points on both bodies?

	// Note: contact_points might not be set for a surface-surface contact. e.g. cylinder
	// on plane
	// TODO: add contact patch info.
	ContactType type;

public:
	PrimPrimContactInfo(): min_distance(0), type(ContactType::UNDEFINED) { }
	PrimPrimContactInfo(double adist, ContactType atype): min_distance(adist), type(atype) { }

	// clear
	void clear() {
		contact_points.clear();
		type = ContactType::UNDEFINED;
	}

	// flip normal
	void flipNormal();
};

// forward definitions of primitive properties
struct CapsuleProperties;
struct PlaneProperties;
struct SphereProperties;
struct CylinderProperties;

class Primitive {
public:
	enum GeometryType { Undefined, Plane, Sphere, Capsule, Cylinder };

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

// capsule primitive
// capsule is assumed to be aligned with X axis. local frame at the center of the capsule.
struct CapsuleProperties {
	double radius;
	double length;
};

class CapsulePrimitive: public Primitive {
public:
	CapsulePrimitive(const std::string& name, double radius, double length);

	~CapsulePrimitive() {
		delete _props;
	}

	virtual void setCapsuleProperties(CapsuleProperties* props);

	virtual CapsuleProperties* getCapsuleProperties() { return _props; }

public:
	CapsuleProperties* _props;
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

// --------- distance computations. static functions ----------
class PrimPrimDistance {
public:
	// distance computation. TODO: Think about template specializations for speed up
	// returns contact info in the world frame
	// NOTE: the order of the primitives matter if both are associated with ARBs
	// the contact normal returned is directed from primA to primB
	// TODO: add flag for whether clues should be used from the passed info to compute
	// delta distance updates
	// NOTE: currently, each of these calls clears the existing prim_prim_info
	static void distancePrimitivePrimitive(
		PrimPrimContactInfo& prim_prim_info,
		const Primitive* primA, Eigen::Affine3d primAinWorld,
		const Primitive* primB, Eigen::Affine3d primBinWorld
	);

	static void distancePlaneCapsule(
		PrimPrimContactInfo& prim_prim_info,
		const PlanePrimitive& plane, Eigen::Affine3d planeInWorld,
		const CapsulePrimitive& capsule, Eigen::Affine3d capsuleInWorld
	);

	static void distanceCapsuleCapsule(
		PrimPrimContactInfo& prim_prim_info,
		const CapsulePrimitive& capsuleA, Eigen::Affine3d capsuleAInWorld,
		const CapsulePrimitive& capsuleB, Eigen::Affine3d capsuleBInWorld
	);

	static void distancePlaneCylinder(
		PrimPrimContactInfo& prim_prim_info,
		const PlanePrimitive& plane, Eigen::Affine3d planeInWorld,
		const CylinderPrimitive& cylinder, Eigen::Affine3d cylinderInWorld
	);

	// TODO: Cylinder capsule
	// TODO: Cylinder cylinder
};

}

#endif //PRIMITIVE_H