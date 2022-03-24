// PrimPrimDistance.h

#ifndef PRIM_PRIM_DISTANCE_H
#define PRIM_PRIM_DISTANCE_H

#include "Primitive.h"
#include "Composite1PkN.h"
#include "PrimPrimContactInfo.h"

namespace Sai2COPSim {

// Constants used
namespace PrimitiveAlgorithmicConstants {
	const double MULTI_POINT_HIGHER_PAIR_CONTACT_DISTANCE_DIFF_THRESHOLD = 6e-4; //m
	const double MIN_HIGHER_PAIR_CONTACT_EXTENT_ANY_DIR = 0.001; //m
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

	// NOTE: capsule is treated as a positive capsule here.
	static void distanceNegCapsuleCapsule(
		PrimPrimContactInfo& prim_prim_info,
		const NegCapsulePrimitive& negCapsule, Eigen::Affine3d negCapsuleInWorld,
		const CapsulePrimitive& capsule, Eigen::Affine3d capsuleInWorld
	);
};

}

#endif // PRIM_PRIM_DISTANCE_H