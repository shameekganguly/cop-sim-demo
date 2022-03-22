// PrimPrimDistance.h

#ifndef PRIM_PRIM_DISTANCE_H
#define PRIM_PRIM_DISTANCE_H

#include "Primitive.h"
#include "Composite1PkN.h"

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
	SURFACE,
	// Currently, CONCAVE implies that the contact consists of a set of points, each with
	// a different normal, curvatures and tangent directions
	CONCAVE
};

// TODO: reorganize this to reflect the core assumption that two primitives can be in
// contact across a set of contact patches, where each contact patch is *planar*.
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

	static void distanceComposite1PkNCapsule(
		PrimPrimContactInfo& prim_prim_info,
		const Composite1PkN& composite, Eigen::Affine3d compositeInWorld,
		const CapsulePrimitive& capsule, Eigen::Affine3d capsuleInWorld
	);

	// TODO: Composite1PkN - other primitives
	// TODO: Composite1PkN - Composite1PkN

	static void distanceNegCapsuleCapsule(
		PrimPrimContactInfo& prim_prim_info,
		const NegCapsulePrimitive& negCapsule, Eigen::Affine3d negCapsuleInWorld,
		const CapsulePrimitive& capsule, Eigen::Affine3d capsuleInWorld
	);
};

}

#endif // PRIM_PRIM_DISTANCE_H