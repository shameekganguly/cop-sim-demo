// Primitive.h

#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include <string>
#include <vector>
#include <Eigen/Dense>

#include "ContactPatch.h"

namespace Sai2COPSim {

// Constants used
namespace PrimitiveAlgorithmicConstants {
	const double MULTI_POINT_HIGHER_PAIR_CONTACT_DISTANCE_DIFF_THRESHOLD = 6e-4; //m
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
	// TODO: for geometries with different min and max curvatures at the contact point, we need
	// to ensure that constraint_dir1 and constraint_dir2 are aligned with the max and min
	// curvature planes respectively of primA.
	// TODO: for concave objects, the each contact point can have its own normal and
	// constraint dirs. Consider if we want to support that, or assume that each primitive
	// is strictly convex
	std::vector<Eigen::Vector3d> contact_points; // closest points on either prim A or prim B in world frame
	//TODO: do we need to be consistent about which primitive the points lie on?
	// or do we need to return points on both bodies?

	// signed curvature for each body
	// sign is positive is center of curvature lies in the positive normal direction
	// and negative otherwise
	// TODO: this works only for point contacts. need to extend to line contacts and surface contacts
	double primA_max_radius;
	double primA_min_radius;
	double primB_max_radius;
	double primB_min_radius;

	// angle from max curvature plane of primA to max curvature plane of primB
	double inter_prim_max_curvature_plane_angle;

	// Note: contact_points might not be set for a surface-surface contact. e.g. cylinder
	// on plane
	// contact patch info
	ContactPatch contact_patch;

	ContactType type;

public:
	PrimPrimContactInfo():
		min_distance(0),
		primA_max_radius(0),
		primA_min_radius(0),
		primB_max_radius(0),
		primB_min_radius(0),
		inter_prim_max_curvature_plane_angle(0),
		type(ContactType::UNDEFINED) { }
	PrimPrimContactInfo(double adist, ContactType atype):
		min_distance(adist),
		primA_max_radius(0),
		primA_min_radius(0),
		primB_max_radius(0),
		primB_min_radius(0),
		inter_prim_max_curvature_plane_angle(0),
		type(atype) { }

	// clear
	void clear() {
		contact_points.clear();
		type = ContactType::UNDEFINED;
		contact_patch.clear();
		primA_max_radius = 0;
		primA_min_radius = 0;
		primB_max_radius = 0;
		primB_min_radius = 0;
		inter_prim_max_curvature_plane_angle = 0;
	}

	// flip normal
	void flipNormal();
};

// forward definitions of primitive properties
struct CapsuleProperties;
struct PlaneProperties;
struct SphereProperties;
struct CylinderProperties;
struct BoxProperties;
struct PyramidProperties;

class Primitive {
public:
	enum GeometryType { Undefined, Plane, Sphere, Capsule, Cylinder, Box, Pyramid };

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

	static void distancePlaneSphere(
		PrimPrimContactInfo& prim_prim_info,
		const PlanePrimitive& plane, Eigen::Affine3d planeInWorld,
		const SpherePrimitive& sphere, Eigen::Affine3d sphereInWorld
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

	static void distancePlaneBox(
		PrimPrimContactInfo& prim_prim_info,
		const PlanePrimitive& plane, Eigen::Affine3d planeInWorld,
		const BoxPrimitive& box, Eigen::Affine3d boxInWorld
	);

	// TODO: Box capsule
	// TODO: Box cylinder
	// TODO: Box box

	static void distancePlanePyramid(
		PrimPrimContactInfo& prim_prim_info,
		const PlanePrimitive& plane, Eigen::Affine3d planeInWorld,
		const PyramidPrimitive& pyramid, Eigen::Affine3d pyramidInWorld
	);

	// TODO: Pyramid capsule
	// TODO: Pyramid cylinder
	// TODO: Pyramid box
	// TODO: Pyramid pyramid
};

}

#endif //PRIMITIVE_H