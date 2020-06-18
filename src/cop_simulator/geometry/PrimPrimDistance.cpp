// PrimPrimDistance.cpp

#include <iostream>
#include "Primitive.h"

using namespace Eigen;

namespace Sai2COPSim {
	void PrimPrimContactInfo::flipNormal() {
		normal_dir *= -1.0;
		constraint_dir2 *= -1.0;
	}

	PrimPrimContactInfo PrimPrimDistance::distancePrimitivePrimitive(
		const Primitive* primA, Eigen::Affine3d primAinWorld,
		const Primitive* primB, Eigen::Affine3d primBinWorld
	) {
		if(primA->_type == Primitive::GeometryType::Capsule && primB->_type == Primitive::GeometryType::Plane) {
			auto ret_info = distancePlaneCapsule(
				*(dynamic_cast<const PlanePrimitive*>(primB)), primBinWorld,
				*(dynamic_cast<const CapsulePrimitive*>(primA)), primAinWorld
			);
			ret_info.flipNormal();
			return ret_info;
		} else if(primB->_type == Primitive::GeometryType::Capsule && primA->_type == Primitive::GeometryType::Plane) {
			return distancePlaneCapsule(
				*(dynamic_cast<const PlanePrimitive*>(primA)), primAinWorld,
				*(dynamic_cast<const CapsulePrimitive*>(primB)), primBinWorld
			);
		} else if(primA->_type == Primitive::GeometryType::Capsule && primB->_type == Primitive::GeometryType::Capsule) {
			return distanceCapsuleCapsule(
				*(dynamic_cast<const CapsulePrimitive*>(primA)), primAinWorld,
				*(dynamic_cast<const CapsulePrimitive*>(primB)), primBinWorld
			);
		} else {
			std::cerr << "PrimA type: " << primA->_type << std::endl;
			std::cerr << "PrimB type: " << primB->_type << std::endl;
			throw(std::runtime_error("Distance between these primitives is not implemented."));
		}
	}



	// returns surface normal pointing outward from plane towards capsule
	PrimPrimContactInfo PrimPrimDistance::distancePlaneCapsule(
		const PlanePrimitive& plane, Eigen::Affine3d planeInWorld,
		const CapsulePrimitive& capsule, Eigen::Affine3d capsuleInWorld
	) {
		assert(capsule._props != NULL);
		assert(plane._props != NULL);

		double capsule_len = capsule._props->length;
		double capsule_radius = capsule._props->radius;
		Vector3d plane_point = plane._props->point;
		Vector3d plane_normal = plane._props->normal;

		Vector3d end1_world, end2_world;
		Vector3d end1_local(-capsule_len/2, 0.0, 0.0);
		Vector3d end2_local(-capsule_len/2, 0.0, 0.0);
		double end1_distance, end2_distance;
		end1_world = capsuleInWorld*end1_local;
		end2_world = capsuleInWorld*end2_local;

		Vector3d plane_point_world = planeInWorld*plane_point;
		Vector3d plane_normal_world = planeInWorld.linear()*plane_normal;
		
		end1_distance = (end1_world - plane_point_world).dot(plane_normal_world) - capsule_radius;
		end2_distance = (end2_world - plane_point_world).dot(plane_normal_world) - capsule_radius;

		PrimPrimContactInfo ret_info;
		ret_info.normal_dir = plane_normal_world;
		// determine tangent directions
		// first try the capsule axis
		Vector3d inter_end_axis = capsuleInWorld.linear().col(0); // local x axis
		double testaxis1_normal_proj = inter_end_axis.dot(plane_normal_world);
		Vector3d worldx_axis(1.0,0.0,0.0);
		double testaxis2_normal_proj = worldx_axis.dot(plane_normal_world);
		Vector3d worldy_axis(0.0,1.0,0.0);
		double testaxis3_normal_proj = worldy_axis.dot(plane_normal_world);
		if(testaxis1_normal_proj < 0.999) {
			// we can use the inter_end axis as contact_direction1, i.e. the x-axis in 
			// the tangent plane
			ret_info.constraint_dir1 = inter_end_axis - testaxis1_normal_proj*plane_normal_world;
			ret_info.constraint_dir1 /= ret_info.constraint_dir1.norm();
		} else if(testaxis2_normal_proj < 0.999) { // test world x axis next
			ret_info.constraint_dir1 = worldx_axis - testaxis2_normal_proj*plane_normal_world;
			ret_info.constraint_dir1 /= ret_info.constraint_dir1.norm();
		} else { // use world y axis
			ret_info.constraint_dir1 = worldy_axis - testaxis3_normal_proj*plane_normal_world;
			ret_info.constraint_dir1 /= ret_info.constraint_dir1.norm();
		}
		ret_info.constraint_dir2 = ret_info.normal_dir.cross(ret_info.constraint_dir1);
		if(abs(end1_distance - end2_distance) < 
							PrimitiveAlgorithmicConstants::MULTI_POINT_HIGHER_PAIR_CONTACT_DISTANCE_DIFF_THRESHOLD
		) {
			ret_info.type = ContactType::LINE;
		// 	// if(LOG_DEBUG) cout << "Line contact " << endl;
			ret_info.contact_points.push_back(end1_world - capsule_radius * plane_normal_world);
			ret_info.contact_points.push_back(end2_world - capsule_radius * plane_normal_world);
		} else if(end1_distance < end2_distance) {
			ret_info.type = ContactType::POINT;
		// 	// if(LOG_DEBUG) cout << "End 1 contact " << endl;
			ret_info.contact_points.push_back(end1_world - capsule_radius * plane_normal_world);
		} else {
			ret_info.type = ContactType::POINT;
		// 	// if(LOG_DEBUG) cout << "End 2 contact " << endl;
			ret_info.contact_points.push_back(end2_world - capsule_radius * plane_normal_world);
		}
		ret_info.min_distance = fmin(end1_distance, end2_distance);
		return ret_info;
	}

	PrimPrimContactInfo PrimPrimDistance::distanceCapsuleCapsule(
		const CapsulePrimitive& capsuleA, Eigen::Affine3d capsuleAInWorld,
		const CapsulePrimitive& capsuleB, Eigen::Affine3d capsuleBInWorld
	) {
		throw(std::runtime_error("Not implemented"));
		return PrimPrimContactInfo();
	}
}