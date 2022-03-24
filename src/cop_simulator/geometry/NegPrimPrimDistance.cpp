// NegPrimPrimDistance.cpp

#include <algorithm>
#include <iostream>
#include "Primitive.h"
#include "PrimPrimDistance.h"
#include "GeometryUtils.h"

using namespace Eigen;

namespace Sai2COPSim {

void PrimPrimDistance::distanceNegCapsuleCapsule(
	PrimPrimContactInfo& prim_prim_info,
	const NegCapsulePrimitive& negCapsule, Eigen::Affine3d negCapsuleInWorld,
	const CapsulePrimitive& capsule, Eigen::Affine3d capsuleInWorld
) {
	assert(negCapsule._props != NULL);
	assert(capsule._props != NULL);

	prim_prim_info.clear();

	double neg_cap_rad = negCapsule._props->radius;
	double neg_cap_half_len = negCapsule._props->length * 0.5;

	PrimPrimContactInfo pos_dist_info;
	CapsulePrimitive pos_cap("tdncc", neg_cap_rad, 2*neg_cap_half_len);
	distanceCapsuleCapsule(pos_dist_info, pos_cap, negCapsuleInWorld, capsule, capsuleInWorld);
	// TODO: this check can be generalized to all negative volumes. The positive volume acts
	// as a cover for the negative volume
	if(pos_dist_info.min_distance > 0) {
		// no intersection with negative capsule, so just return
		prim_prim_info.type = ContactType::CONCAVE;
		prim_prim_info.min_distance = PrimitiveAlgorithmicConstants::DISTANCE_FULL_PENETRATION;
		return;
	}

	double pos_cap_rad = capsule._props->radius;
	double pos_cap_half_len = capsule._props->length * 0.5;

	// TODO: think a little bit more about this. There is a chance that duplicate
	// contact points are reported from here and from the distance check with the
	// intersection edge between the negative and the parent positive primitives.
	// TODO: make this ratio configurable
	// TODO: move this check to a canIntersect function instead, specific to
	// negative primitives.
	if(pos_cap_rad > neg_cap_rad*0.9) {
		prim_prim_info.type = ContactType::UNDEFINED;
		return;
	}

	// compute relative transform
	Affine3d pos_cap_local_tf = negCapsuleInWorld.inverse()*capsuleInWorld;
	Vector3d pos_cap_axis_local = pos_cap_local_tf.linear().col(0);
	Vector3d pos_cap_center_local = pos_cap_local_tf.translation();
	Vector3d left_end = pos_cap_center_local - pos_cap_axis_local*pos_cap_half_len;
	Vector3d right_end = pos_cap_center_local + pos_cap_axis_local*pos_cap_half_len;
	auto check_pos_end_lies_on_axis = [&neg_cap_half_len](const Vector3d& pos_end) -> bool {
		if(fabs(pos_end(0)) > neg_cap_half_len) {
			// end pos is closer to the end cap of the negative capsule than the side
			// walls
			return false;
		}
		return (pos_end.segment<2>(1).norm() < 1e-3);
	};

	auto closest_pt_neg_cap_axis = [&neg_cap_half_len](Vector3d& test_pt) -> Vector3d {
		return Vector3d(fmin(fmax(test_pt[0], -neg_cap_half_len), neg_cap_half_len), 0, 0);
	};

	// Call ONLY ONCE to add both tangent directions.
	auto add_constraint_dirs = [&prim_prim_info](Vector3d& normal) {
		if(abs(normal(0)) < 0.999) {
			// use neg capsule x axis
			Vector3d dir1(1.0, 0.0, 0.0);
			dir1 -= normal(0)*normal;
			prim_prim_info.constraint_dir1s.push_back(dir1);
		} else {
			// use neg capsule y axis
			Vector3d dir1(0.0, 1.0, 0.0);
			dir1 -= normal(1)*normal;
			prim_prim_info.constraint_dir1s.push_back(dir1);
		}
		prim_prim_info.constraint_dir2s.push_back(normal.cross(prim_prim_info.constraint_dir1s.back()));
	};

	auto get_normal_axisymmetric = [&neg_cap_rad](Vector3d& side_wall_pt) -> Vector3d {
		return -Vector3d(0, side_wall_pt(1), side_wall_pt(2))/neg_cap_rad;
	};

	for (uint end_id = 0; end_id < 2; end_id++) {
		Vector3d end_pt = (end_id == 0)? left_end: right_end;
		if(check_pos_end_lies_on_axis(end_pt)) {
			Vector3d ptA1(end_pt(0), neg_cap_rad, 0);
			Vector3d ptA2(end_pt(0), -neg_cap_rad*cos(M_PI/3), neg_cap_rad*sin(M_PI/3));
			Vector3d ptA3(end_pt(0), -neg_cap_rad*cos(M_PI/3), -neg_cap_rad*sin(M_PI/3));

			Vector3d normal1(get_normal_axisymmetric(ptA1));
			prim_prim_info.contact_points.push_back(ptA1);
			prim_prim_info.normal_dirs.push_back(normal1);
			prim_prim_info.distances.push_back(neg_cap_rad - pos_cap_rad);
			add_constraint_dirs(normal1);

			Vector3d normal2(get_normal_axisymmetric(ptA2));
			prim_prim_info.contact_points.push_back(ptA2);
			prim_prim_info.normal_dirs.push_back(normal2);
			prim_prim_info.distances.push_back(neg_cap_rad - pos_cap_rad);
			add_constraint_dirs(normal2);

			Vector3d normal3(get_normal_axisymmetric(ptA3));
			prim_prim_info.contact_points.push_back(ptA3);
			prim_prim_info.normal_dirs.push_back(normal3);
			prim_prim_info.distances.push_back(neg_cap_rad - pos_cap_rad);
			add_constraint_dirs(normal3);
		} else {
			Vector3d ptA = closest_pt_neg_cap_axis(end_pt);
			Vector3d normal = ptA - end_pt;
			double dist = normal.norm();
			// This is safe since we handled the axisymmetric case already
			normal = normal/dist;
			ptA += (-normal)*neg_cap_rad; // project from axis to surface
			prim_prim_info.contact_points.push_back(ptA);
			prim_prim_info.normal_dirs.push_back(normal);
			prim_prim_info.distances.push_back(neg_cap_rad - dist - pos_cap_rad);
			add_constraint_dirs(normal);
		}
	}

	prim_prim_info.type = ContactType::CONCAVE;
	prim_prim_info.min_distance = *(std::min_element(prim_prim_info.distances.begin(), prim_prim_info.distances.end()));

	// transfer from negCapsule coordinates to world coordinates
	for(auto& normal_dir: prim_prim_info.normal_dirs) {
		normal_dir = negCapsuleInWorld.linear()*normal_dir;
	}
	for(auto& constraint_dir: prim_prim_info.constraint_dir1s) {
		constraint_dir = negCapsuleInWorld.linear()*constraint_dir;
	}
	for(auto& constraint_dir: prim_prim_info.constraint_dir2s) {
		constraint_dir = negCapsuleInWorld.linear()*constraint_dir;
	}
	for(auto& pt: prim_prim_info.contact_points) {
		pt = negCapsuleInWorld*pt;
	}
}

}