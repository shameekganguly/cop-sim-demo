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
		double testaxis1_normal_proj = abs(inter_end_axis.dot(plane_normal_world));
		Vector3d worldx_axis(1.0,0.0,0.0);
		double testaxis2_normal_proj = abs(worldx_axis.dot(plane_normal_world));
		Vector3d worldy_axis(0.0,1.0,0.0);
		double testaxis3_normal_proj = abs(worldy_axis.dot(plane_normal_world));
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
		assert(capsuleA._props != NULL);
		assert(capsuleB._props != NULL);

		double capAradius = capsuleA._props->radius;
		double capBradius = capsuleB._props->radius;
		double capAlength = capsuleA._props->length;
		double capBlength = capsuleB._props->length;

		// compute relative transform
		Affine3d capsuleBToCapsuleA = capsuleAInWorld.inverse()*capsuleBInWorld;
		Vector3d capBaxisToA = capsuleBToCapsuleA.linear().col(0);
		Vector3d capBcenterToA = capsuleBToCapsuleA.translation();

		PrimPrimContactInfo ret_info;
		Vector3d ptA, ptB; // for point contact

		// case 1: capsules are parallel to each other and close enough for line contact
		if(abs(capBaxisToA(0)) > 0.99) { // we allow some almost parallel cases
			if(abs(capBcenterToA(0)) < 0.5*(capAlength + capBlength)) {
				// case 1a: line contact
				assert(capBcenterToA.norm() > 0.5*(capAradius + capBradius)); // this should not happen in practice
				ret_info.type = ContactType::LINE;
				ret_info.min_distance = capBcenterToA.tail(2).norm() - capAradius - capBradius;
				assert(capBcenterToA.tail(2).norm() > 0.5*(capAradius + capBradius)); // this should not happen in practice
				// compute normal and constraint directions
				ret_info.normal_dir << 0.0, capBcenterToA(1), capBcenterToA(2);
				ret_info.normal_dir /= ret_info.normal_dir.norm();
				ret_info.constraint_dir1 << 1.0, 0.0, 0.0;
				ret_info.constraint_dir2 = ret_info.normal_dir.cross(ret_info.constraint_dir1);
				// compute overlap points for line contact
				if(capBcenterToA(0) - capBlength/2.0 > -capAlength/2.0) {
					Vector3d pt = capBcenterToA - Vector3d(capBlength/2.0, 0.0, 0.0);
					ret_info.contact_points.push_back(pt - ret_info.normal_dir*capBradius);
				} else {
					Vector3d pt = Vector3d(-capAlength/2.0, 0.0, 0.0);
					ret_info.contact_points.push_back(pt + ret_info.normal_dir*capAradius);
				}
				if (capBcenterToA(0) + capBlength/2.0 < capAlength/2.0) {
					Vector3d pt = capBcenterToA + Vector3d(capBlength/2.0, 0.0, 0.0);
					ret_info.contact_points.push_back(pt - ret_info.normal_dir*capBradius);
				} else {
					Vector3d pt = Vector3d(capAlength/2.0, 0.0, 0.0);
					ret_info.contact_points.push_back(pt + ret_info.normal_dir*capAradius);
				}
			} else {
				// case 1b: end contact
				ret_info.type = ContactType::POINT;
				if(capBcenterToA(0) >= 0) {
					ptA << capAlength/2.0, 0.0, 0.0;
					ptB << -capBlength/2.0, 0.0, 0.0;
					ptB += capBcenterToA;
				} else {
					ptA << -capAlength/2.0, 0.0, 0.0;
					ptB << capBlength/2.0, 0.0, 0.0;
					ptB += capBcenterToA;
				}
			}
		} else {	// case 2: point contact
			ret_info.type = ContactType::POINT;
			// compute minimum distance between the inner line segments of the two capsules
			double temp = capBaxisToA(0);
			double Delta = temp*temp - 1.0;
			double p2p1d1 = capBcenterToA(0);
			double p2p1d2 = capBcenterToA.dot(capBaxisToA);
			double z1 = -1.0/Delta*(p2p1d1 - temp*p2p1d2);
			double z2 = -1.0/Delta*(temp*p2p1d1 - p2p1d2);

			if(abs(z1) > 0.5*capAlength) {
				// clamp z1
				z1 = z1/abs(z1) * 0.5*capAlength;
			}
			if(abs(z2) > 0.5*capBlength) {
				// clamp z2
				z2 = z2/abs(z2) * 0.5*capBlength;
			}
			ptA << z1, 0.0, 0.0;
			ptB = capBcenterToA + z2*capBaxisToA;
		}

		// set normal and constraint directions for point contact
		if(ret_info.type == ContactType::POINT) {
			ret_info.normal_dir = ptB - ptA;
			ret_info.min_distance = ret_info.normal_dir.norm() - capAradius - capBradius;
			assert(ret_info.normal_dir.norm() > 0.25*(capAradius + capBradius));
			ret_info.normal_dir /= ret_info.normal_dir.norm();
			// calculate constraint_dir
			if(abs(ret_info.normal_dir(0)) < 0.999) {
				// use capsule A x axis
				ret_info.constraint_dir1 << 1.0, 0.0, 0.0;
				ret_info.constraint_dir1 -= ret_info.normal_dir(0)*ret_info.normal_dir;
			} else {
				// use capsule A y axis
				ret_info.constraint_dir1 << 0.0, 1.0, 0.0;
				ret_info.constraint_dir1 -= ret_info.normal_dir(1)*ret_info.normal_dir;
			}
			ret_info.constraint_dir1 /= ret_info.constraint_dir1.norm();
			ret_info.constraint_dir2 = ret_info.normal_dir.cross(ret_info.constraint_dir1);
			ret_info.contact_points.push_back(ptA + ret_info.normal_dir*capAradius);
		} else {
			// check line length for line contact
			if((ret_info.contact_points[0] - ret_info.contact_points[1]).norm() <
					PrimitiveAlgorithmicConstants::MIN_HIGHER_PAIR_CONTACT_EXTENT_ANY_DIR) {
				ret_info.type = ContactType::POINT;
				ret_info.contact_points.resize(1);
			}
		}

		// transfer from capsuleA coordinates to world coordinates
		ret_info.normal_dir = capsuleAInWorld.linear()*ret_info.normal_dir;
		ret_info.constraint_dir1 = capsuleAInWorld.linear()*ret_info.constraint_dir1;
		ret_info.constraint_dir2 = capsuleAInWorld.linear()*ret_info.constraint_dir2;
		for(auto& pt: ret_info.contact_points) {
			pt = capsuleAInWorld*pt;
		}
		return ret_info;
	}
}