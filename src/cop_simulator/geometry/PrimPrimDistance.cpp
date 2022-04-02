// PrimPrimDistance.cpp

#include <algorithm>
#include <iostream>
#include "PointPrimDistance.h"
#include "Primitive.h"
#include "PrimNormal.h"
#include "PrimPrimDistance.h"
#include "GeometryUtils.h"

using namespace Eigen;

namespace Sai2COPSim {
	void PrimPrimDistance::distancePrimitivePrimitive(
		PrimPrimContactInfo& prim_prim_info,
		const Primitive* primA, Eigen::Affine3d primAinWorld,
		const Primitive* primB, Eigen::Affine3d primBinWorld
	) {
		if(primA->_type == Primitive::GeometryType::Capsule && primB->_type == Primitive::GeometryType::Plane) {
			distancePlaneCapsule(
				prim_prim_info,
				*(dynamic_cast<const PlanePrimitive*>(primB)), primBinWorld,
				*(dynamic_cast<const CapsulePrimitive*>(primA)), primAinWorld
			);
			prim_prim_info.flipNormal();
			return;
		} else if(primA->_type == Primitive::GeometryType::Sphere && primB->_type == Primitive::GeometryType::Plane) {
			distancePlaneSphere(
				prim_prim_info,
				*(dynamic_cast<const PlanePrimitive*>(primB)), primBinWorld,
				*(dynamic_cast<const SpherePrimitive*>(primA)), primAinWorld
			);
			prim_prim_info.flipNormal();
			return;
		} else if(primB->_type == Primitive::GeometryType::Sphere && primA->_type == Primitive::GeometryType::Plane) {
			distancePlaneSphere(
				prim_prim_info,
				*(dynamic_cast<const PlanePrimitive*>(primA)), primAinWorld,
				*(dynamic_cast<const SpherePrimitive*>(primB)), primBinWorld
			);
			return;
		} else if(primB->_type == Primitive::GeometryType::Capsule && primA->_type == Primitive::GeometryType::Plane) {
			distancePlaneCapsule(
				prim_prim_info,
				*(dynamic_cast<const PlanePrimitive*>(primA)), primAinWorld,
				*(dynamic_cast<const CapsulePrimitive*>(primB)), primBinWorld
			);
			return;
		} else if(primA->_type == Primitive::GeometryType::Capsule && primB->_type == Primitive::GeometryType::Capsule) {
			distanceCapsuleCapsule(
				prim_prim_info,
				*(dynamic_cast<const CapsulePrimitive*>(primA)), primAinWorld,
				*(dynamic_cast<const CapsulePrimitive*>(primB)), primBinWorld
			);
			return;
		} else if(primA->_type == Primitive::GeometryType::Cylinder && primB->_type == Primitive::GeometryType::Plane) {
			distancePlaneCylinder(
				prim_prim_info,
				*(dynamic_cast<const PlanePrimitive*>(primB)), primBinWorld,
				*(dynamic_cast<const CylinderPrimitive*>(primA)), primAinWorld
			);
			prim_prim_info.flipNormal();
			return;
		} else if(primB->_type == Primitive::GeometryType::Cylinder && primA->_type == Primitive::GeometryType::Plane) {
			distancePlaneCylinder(
				prim_prim_info,
				*(dynamic_cast<const PlanePrimitive*>(primA)), primAinWorld,
				*(dynamic_cast<const CylinderPrimitive*>(primB)), primBinWorld
			);
			return;
		} else if(primA->_type == Primitive::GeometryType::Box && primB->_type == Primitive::GeometryType::Plane) {
			distancePlaneBox(
				prim_prim_info,
				*(dynamic_cast<const PlanePrimitive*>(primB)), primBinWorld,
				*(dynamic_cast<const BoxPrimitive*>(primA)), primAinWorld
			);
			prim_prim_info.flipNormal();
			return;
		} else if(primB->_type == Primitive::GeometryType::Box && primA->_type == Primitive::GeometryType::Plane) {
			distancePlaneBox(
				prim_prim_info,
				*(dynamic_cast<const PlanePrimitive*>(primA)), primAinWorld,
				*(dynamic_cast<const BoxPrimitive*>(primB)), primBinWorld
			);
			return;
		} else if(primA->_type == Primitive::GeometryType::Pyramid && primB->_type == Primitive::GeometryType::Plane) {
			distancePlanePyramid(
				prim_prim_info,
				*(dynamic_cast<const PlanePrimitive*>(primB)), primBinWorld,
				*(dynamic_cast<const PyramidPrimitive*>(primA)), primAinWorld
			);
			prim_prim_info.flipNormal();
			return;
		} else if(primB->_type == Primitive::GeometryType::Pyramid && primA->_type == Primitive::GeometryType::Plane) {
			distancePlanePyramid(
				prim_prim_info,
				*(dynamic_cast<const PlanePrimitive*>(primA)), primAinWorld,
				*(dynamic_cast<const PyramidPrimitive*>(primB)), primBinWorld
			);
			return;
		} else if(primA->_type == Primitive::GeometryType::Capsule && primB->_type == Primitive::GeometryType::Composite1PkN) {
			distanceComposite1PkNCapsule(
				prim_prim_info,
				*(dynamic_cast<const Composite1PkN*>(primB)), primBinWorld,
				*(dynamic_cast<const CapsulePrimitive*>(primA)), primAinWorld
			);
			prim_prim_info.flipNormal();
			return;
		} else if(primB->_type == Primitive::GeometryType::Capsule && primA->_type == Primitive::GeometryType::Composite1PkN) {
			distanceComposite1PkNCapsule(
				prim_prim_info,
				*(dynamic_cast<const Composite1PkN*>(primA)), primAinWorld,
				*(dynamic_cast<const CapsulePrimitive*>(primB)), primBinWorld
			);
			return;
		} else if(primA->_type == Primitive::GeometryType::NegCapsule && primB->_type == Primitive::GeometryType::Capsule) {
			// TODO: what if primA is pos capsule and primB is neg capsule?
			distanceNegCapsuleCapsule(
				prim_prim_info,
				*(dynamic_cast<const NegCapsulePrimitive*>(primA)), primAinWorld,
				*(dynamic_cast<const CapsulePrimitive*>(primB)), primBinWorld
			);
			return;
		} else {
			std::cerr << "PrimA type: " << primA->_type << std::endl;
			std::cerr << "PrimB type: " << primB->_type << std::endl;
			throw(std::runtime_error("Distance between these primitives is not implemented."));
		}
	}

	void PrimPrimDistance::distancePlaneSphere(
		PrimPrimContactInfo& prim_prim_info,
		const PlanePrimitive& plane, Eigen::Affine3d planeInWorld,
		const SpherePrimitive& sphere, Eigen::Affine3d sphereInWorld
	) {
		assert(sphere._props != NULL);
		assert(plane._props != NULL);

		prim_prim_info.clear();

		Vector3d plane_point = plane._props->point;
		Vector3d plane_normal = plane._props->normal;

		Vector3d plane_point_world = planeInWorld*plane_point;
		Vector3d plane_normal_world = planeInWorld.linear()*plane_normal;
		prim_prim_info.normal_dir = plane_normal_world;

		Vector3d worldx_axis(1.0,0.0,0.0);
		double testaxis2_normal_proj = worldx_axis.dot(plane_normal_world);
		Vector3d worldy_axis(0.0,1.0,0.0);
		double testaxis3_normal_proj = worldy_axis.dot(plane_normal_world);
		if(abs(testaxis2_normal_proj) < 0.999) { // use world x axis
			prim_prim_info.constraint_dir1 = worldx_axis - testaxis2_normal_proj*plane_normal_world;
			prim_prim_info.constraint_dir1 /= prim_prim_info.constraint_dir1.norm();
		} else { // use world y axis
			prim_prim_info.constraint_dir1 = worldy_axis - testaxis3_normal_proj*plane_normal_world;
			prim_prim_info.constraint_dir1 /= prim_prim_info.constraint_dir1.norm();
		}
		prim_prim_info.constraint_dir2 = prim_prim_info.normal_dir.cross(prim_prim_info.constraint_dir1);

		double radius = sphere._props->radius;
		Vector3d center = sphereInWorld.translation();
		prim_prim_info.type = ContactType::POINT;
		prim_prim_info.contact_points.push_back(center - radius * plane_normal_world);
		prim_prim_info.primB_max_radius = radius;
		prim_prim_info.primB_min_radius = radius;

		double distance = (center - plane_point_world).dot(plane_normal_world) - radius;
		prim_prim_info.min_distance = distance;
	}

	// returns surface normal pointing outward from plane towards capsule
	void PrimPrimDistance::distancePlaneCapsule(
		PrimPrimContactInfo& prim_prim_info,
		const PlanePrimitive& plane, Eigen::Affine3d planeInWorld,
		const CapsulePrimitive& capsule, Eigen::Affine3d capsuleInWorld
	) {
		assert(capsule._props != NULL);
		assert(plane._props != NULL);

		prim_prim_info.clear();

		Vector3d plane_point = plane._props->point;
		Vector3d plane_normal = plane._props->normal;

		Vector3d plane_point_world = planeInWorld*plane_point;
		Vector3d plane_normal_world = planeInWorld.linear()*plane_normal;

		double capsule_len = capsule._props->length;
		double capsule_radius = capsule._props->radius;

		// TODO: this should move out to the sphere hierarchy broad phase
		double cover_dist = (capsuleInWorld.translation() - plane_point_world).dot(plane_normal_world)
				- 1.5*0.5*(capsule_len+2*capsule_radius);
		if(cover_dist > 0) {
			// std::cout << "Plane capsule cover: " << cover_dist << std::endl;
			prim_prim_info.min_distance = (capsuleInWorld.translation() - plane_point_world).dot(plane_normal_world)
				- 0.5*(capsule_len+2*capsule_radius);
			return;
		}

		Vector3d end1_world, end2_world;
		Vector3d end1_local(-capsule_len/2, 0.0, 0.0);
		Vector3d end2_local(capsule_len/2, 0.0, 0.0);
		double end1_distance, end2_distance;
		end1_world = capsuleInWorld*end1_local;
		end2_world = capsuleInWorld*end2_local;

		end1_distance = (end1_world - plane_point_world).dot(plane_normal_world) - capsule_radius;
		end2_distance = (end2_world - plane_point_world).dot(plane_normal_world) - capsule_radius;

		prim_prim_info.normal_dir = plane_normal_world;
		// determine tangent directions
		// first try the capsule axis
		// std::cout <<"Cap rotation " << capsuleInWorld.linear() << std::endl;
		Vector3d inter_end_axis = capsuleInWorld.linear().col(0); // local x axis
		double testaxis1_normal_proj = inter_end_axis.dot(plane_normal_world);
		Vector3d worldx_axis(1.0,0.0,0.0);
		double testaxis2_normal_proj = worldx_axis.dot(plane_normal_world);
		Vector3d worldy_axis(0.0,1.0,0.0);
		double testaxis3_normal_proj = worldy_axis.dot(plane_normal_world);
		if(abs(testaxis1_normal_proj) < 0.999) {
			// we can use the inter_end axis as contact_direction1, i.e. the x-axis in
			// the tangent plane
			prim_prim_info.constraint_dir1 = inter_end_axis - testaxis1_normal_proj*plane_normal_world;
			prim_prim_info.constraint_dir1 /= prim_prim_info.constraint_dir1.norm();
		} else if(abs(testaxis2_normal_proj) < 0.999) { // test world x axis next
			prim_prim_info.constraint_dir1 = worldx_axis - testaxis2_normal_proj*plane_normal_world;
			prim_prim_info.constraint_dir1 /= prim_prim_info.constraint_dir1.norm();
		} else { // use world y axis
			prim_prim_info.constraint_dir1 = worldy_axis - testaxis3_normal_proj*plane_normal_world;
			prim_prim_info.constraint_dir1 /= prim_prim_info.constraint_dir1.norm();
		}
		prim_prim_info.constraint_dir2 = prim_prim_info.normal_dir.cross(prim_prim_info.constraint_dir1);
		if(abs(end1_distance - end2_distance) <
							PrimitiveAlgorithmicConstants::MULTI_POINT_HIGHER_PAIR_CONTACT_DISTANCE_DIFF_THRESHOLD
		) {
			prim_prim_info.type = ContactType::LINE;
		// 	// if(LOG_DEBUG) cout << "Line contact " << endl;
			prim_prim_info.contact_points.push_back(end1_world - capsule_radius * plane_normal_world);
			prim_prim_info.contact_points.push_back(end2_world - capsule_radius * plane_normal_world);
		} else if(end1_distance < end2_distance) {
			prim_prim_info.type = ContactType::POINT;
		// 	// if(LOG_DEBUG) cout << "End 1 contact " << endl;
			prim_prim_info.contact_points.push_back(end1_world - capsule_radius * plane_normal_world);
		} else {
			prim_prim_info.type = ContactType::POINT;
		// 	// if(LOG_DEBUG) cout << "End 2 contact " << endl;
			prim_prim_info.contact_points.push_back(end2_world - capsule_radius * plane_normal_world);
		}
		prim_prim_info.min_distance = fmin(end1_distance, end2_distance);
	}

	void PrimPrimDistance::distanceCapsuleCapsule(
		PrimPrimContactInfo& prim_prim_info,
		const CapsulePrimitive& capsuleA, Eigen::Affine3d capsuleAInWorld,
		const CapsulePrimitive& capsuleB, Eigen::Affine3d capsuleBInWorld
	) {
		assert(capsuleA._props != NULL);
		assert(capsuleB._props != NULL);

		prim_prim_info.clear();

		double capAradius = capsuleA._props->radius;
		double capBradius = capsuleB._props->radius;
		double capAlength = capsuleA._props->length;
		double capBlength = capsuleB._props->length;

		// compute relative transform
		Affine3d capsuleBToCapsuleA = capsuleAInWorld.inverse()*capsuleBInWorld;
		Vector3d capBaxisToA = capsuleBToCapsuleA.linear().col(0);
		Vector3d capBcenterToA = capsuleBToCapsuleA.translation();

		// TODO: this should move out to the sphere hierarchy broad phase
		double cover_dist = capBcenterToA.norm()
				- 1.5*0.5*(capAlength+2*capAradius)
				- 1.5*0.5*(capBlength+2*capBradius);
		if(cover_dist > 0) {
			// std::cout << "Capsule capsule cover: " << cover_dist << std::endl;
			prim_prim_info.min_distance = capBcenterToA.norm() - 0.5*(capAlength+2*capAradius) - 0.5*(capBlength+2*capBradius);
			return;
		}

		Vector3d ptA, ptB; // for point contact

		// case 1: capsules are parallel to each other and close enough for line contact
		if(abs(capBaxisToA(0)) > 0.99) { // we allow some almost parallel cases
			if(abs(capBcenterToA(0)) < 0.5*(capAlength + capBlength)) {
				// case 1a: line contact
				assert(capBcenterToA.norm() > 0.5*(capAradius + capBradius)); // this should not happen in practice
				prim_prim_info.type = ContactType::LINE;
				prim_prim_info.min_distance = capBcenterToA.tail(2).norm() - capAradius - capBradius;
				assert(capBcenterToA.tail(2).norm() > 0.5*(capAradius + capBradius)); // this should not happen in practice
				// compute normal and constraint directions
				prim_prim_info.normal_dir << 0.0, capBcenterToA(1), capBcenterToA(2);
				prim_prim_info.normal_dir /= prim_prim_info.normal_dir.norm();
				prim_prim_info.constraint_dir1 << 1.0, 0.0, 0.0;
				prim_prim_info.constraint_dir2 = prim_prim_info.normal_dir.cross(prim_prim_info.constraint_dir1);
				// compute overlap points for line contact
				if(capBcenterToA(0) - capBlength/2.0 > -capAlength/2.0) {
					Vector3d pt = capBcenterToA - Vector3d(capBlength/2.0, 0.0, 0.0);
					prim_prim_info.contact_points.push_back(pt - prim_prim_info.normal_dir*capBradius);
				} else {
					Vector3d pt = Vector3d(-capAlength/2.0, 0.0, 0.0);
					prim_prim_info.contact_points.push_back(pt + prim_prim_info.normal_dir*capAradius);
				}
				if (capBcenterToA(0) + capBlength/2.0 < capAlength/2.0) {
					Vector3d pt = capBcenterToA + Vector3d(capBlength/2.0, 0.0, 0.0);
					prim_prim_info.contact_points.push_back(pt - prim_prim_info.normal_dir*capBradius);
				} else {
					Vector3d pt = Vector3d(capAlength/2.0, 0.0, 0.0);
					prim_prim_info.contact_points.push_back(pt + prim_prim_info.normal_dir*capAradius);
				}
			} else {
				// case 1b: end contact
				prim_prim_info.type = ContactType::POINT;
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
			prim_prim_info.type = ContactType::POINT;
			// compute minimum distance between the inner line segments of the two capsules
			double temp = capBaxisToA(0);
			double Delta = temp*temp - 1.0;
			double p2p1d1 = capBcenterToA(0);
			double p2p1d2 = capBcenterToA.dot(capBaxisToA);
			double z1 = -1.0/Delta*(p2p1d1 - temp*p2p1d2);
			double z2 = -1.0/Delta*(temp*p2p1d1 - p2p1d2);
			if(abs(z1) > 0.5*capAlength || abs(z2) > 0.5*capBlength) {
				double z1v = abs(z1) - 0.5*capAlength;
				double z2v = abs(z2) - 0.5*capBlength;
				if(z1v >= z2v) {
					// capsule A should be end point, capsule B can be anywhere
					z1 = z1/abs(z1) * 0.5*capAlength;
					z2 = -(capBcenterToA - Vector3d(z1, 0.0, 0.0)).dot(capBaxisToA);
				} else {
					// capsule B should be end point, capsule A can be anywhere
					z2 = z2/abs(z2) * 0.5*capBlength;
					z1 = capBcenterToA(0) + z2*capBaxisToA(0);
				}
			}
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
			// std::cout << ptA.transpose() << " - " << ptB.transpose() << std::endl;
		}

		// set normal and constraint directions for point contact
		if(prim_prim_info.type == ContactType::POINT) {
			prim_prim_info.normal_dir = ptB - ptA;
			prim_prim_info.min_distance = prim_prim_info.normal_dir.norm() - capAradius - capBradius;
			assert(prim_prim_info.normal_dir.norm() > 0.25*(capAradius + capBradius));
			prim_prim_info.normal_dir /= prim_prim_info.normal_dir.norm();
			// calculate constraint_dir
			if(abs(prim_prim_info.normal_dir(0)) < 0.999) {
				// use capsule A x axis
				prim_prim_info.constraint_dir1 << 1.0, 0.0, 0.0;
				prim_prim_info.constraint_dir1 -= prim_prim_info.normal_dir(0)*prim_prim_info.normal_dir;
			} else {
				// use capsule A y axis
				prim_prim_info.constraint_dir1 << 0.0, 1.0, 0.0;
				prim_prim_info.constraint_dir1 -= prim_prim_info.normal_dir(1)*prim_prim_info.normal_dir;
			}
			prim_prim_info.constraint_dir1 /= prim_prim_info.constraint_dir1.norm();
			prim_prim_info.constraint_dir2 = prim_prim_info.normal_dir.cross(prim_prim_info.constraint_dir1);
			prim_prim_info.contact_points.push_back(ptA + prim_prim_info.normal_dir*capAradius);
		} else {
			// check line length for line contact
			if((prim_prim_info.contact_points[0] - prim_prim_info.contact_points[1]).norm() <
					PrimitiveAlgorithmicConstants::MIN_HIGHER_PAIR_CONTACT_EXTENT_ANY_DIR) {
				prim_prim_info.type = ContactType::POINT;
				prim_prim_info.contact_points.resize(1);
			}
		}

		// transfer from capsuleA coordinates to world coordinates
		prim_prim_info.normal_dir = capsuleAInWorld.linear()*prim_prim_info.normal_dir;
		prim_prim_info.constraint_dir1 = capsuleAInWorld.linear()*prim_prim_info.constraint_dir1;
		prim_prim_info.constraint_dir2 = capsuleAInWorld.linear()*prim_prim_info.constraint_dir2;
		for(auto& pt: prim_prim_info.contact_points) {
			pt = capsuleAInWorld*pt;
		}
	}

	void PrimPrimDistance::distancePlaneCylinder(
		PrimPrimContactInfo& prim_prim_info,
		const PlanePrimitive& plane, Eigen::Affine3d planeInWorld,
		const CylinderPrimitive& cylinder, Eigen::Affine3d cylinderInWorld
	) {
		assert(cylinder._props != NULL);
		assert(plane._props != NULL);

		prim_prim_info.clear();

		double cylinder_len = cylinder._props->height;
		// std::cout << "Height: " << cylinder_len << std::endl;
		double cylinder_radius = cylinder._props->radius;
		Vector3d plane_point = plane._props->point;
		Vector3d plane_normal = plane._props->normal;

		Vector3d endA_world, endB_world;
		Vector3d endA_local(0.0, 0.0, 0.0);
		Vector3d endB_local(0.0, 0.0, cylinder_len);
		endA_world = cylinderInWorld*endA_local;
		endB_world = cylinderInWorld*endB_local;
		Vector3d cylinder_axis = cylinderInWorld.linear().col(2);

		Vector3d plane_point_world = planeInWorld*plane_point;
		Vector3d plane_normal_world = planeInWorld.linear()*plane_normal;

		prim_prim_info.normal_dir = plane_normal_world;

		// determine tangent directions
		// first try the cylinder axis
		// std::cout <<"Cylinder rotation " << cylinderInWorld.linear() << std::endl;
		Vector3d inter_end_axis = cylinder_axis; // local z axis
		double testaxis1_normal_proj = inter_end_axis.dot(plane_normal_world);
		Vector3d worldx_axis(1.0,0.0,0.0);
		double testaxis2_normal_proj = worldx_axis.dot(plane_normal_world);
		Vector3d worldy_axis(0.0,1.0,0.0);
		double testaxis3_normal_proj = worldy_axis.dot(plane_normal_world);
		if(abs(testaxis1_normal_proj) < 0.999) {
			// we can use the inter_end axis as contact_direction1, i.e. the x-axis in
			// the tangent plane
			prim_prim_info.constraint_dir1 = inter_end_axis - testaxis1_normal_proj*plane_normal_world;
			prim_prim_info.constraint_dir1 /= prim_prim_info.constraint_dir1.norm();
		} else if(abs(testaxis2_normal_proj) < 0.999) { // test world x axis next
			prim_prim_info.constraint_dir1 = worldx_axis - testaxis2_normal_proj*plane_normal_world;
			prim_prim_info.constraint_dir1 /= prim_prim_info.constraint_dir1.norm();
		} else { // use world y axis
			prim_prim_info.constraint_dir1 = worldy_axis - testaxis3_normal_proj*plane_normal_world;
			prim_prim_info.constraint_dir1 /= prim_prim_info.constraint_dir1.norm();
		}
		prim_prim_info.constraint_dir2 = prim_prim_info.normal_dir.cross(prim_prim_info.constraint_dir1);

		// check for line contact
		if(abs(cylinder_axis.dot(plane_normal_world)) <
			PrimitiveAlgorithmicConstants::MULTI_POINT_HIGHER_PAIR_CONTACT_DISTANCE_DIFF_THRESHOLD
		) {
			prim_prim_info.type = ContactType::LINE;
			// if(LOG_DEBUG) cout << "Line contact " << endl;
			prim_prim_info.contact_points.push_back(endA_world - cylinder_radius * plane_normal_world);
			prim_prim_info.contact_points.push_back(endB_world - cylinder_radius * plane_normal_world);
			double endA_distance = (endA_world - plane_point_world).dot(plane_normal_world) - cylinder_radius;
			double endB_distance = (endB_world - plane_point_world).dot(plane_normal_world) - cylinder_radius;
			prim_prim_info.min_distance = fmin(endA_distance, endB_distance);
			return;
		}
		// check for end surface contact
		if(abs(cylinder_axis.dot(plane_normal_world)) > 0.999) {
			prim_prim_info.type = ContactType::SURFACE;
			double endA_distance = (endA_world - plane_point_world).dot(plane_normal_world);
			double endB_distance = (endB_world - plane_point_world).dot(plane_normal_world);
			//TODO: add contact patch
			// if(LOG_DEBUG) cout << "Line contact " << endl;
			Vector3d interior_pt;
			if(endA_distance < endB_distance) {
				// add FaceA points to contact points
				for(auto& pt: cylinder._faceA_points) {
					prim_prim_info.contact_points.push_back(cylinderInWorld*pt);
				}
				interior_pt = endA_world;
			} else {
				// add FaceB points to contact points
				for(auto& pt: cylinder._faceB_points) {
					prim_prim_info.contact_points.push_back(cylinderInWorld*pt);
				}
				interior_pt = endB_world;
			}
			prim_prim_info.contact_patch._interior_point = interior_pt;
			prim_prim_info.contact_patch.max_extent = 2*cylinder_radius;
			auto * c = new Circle2D();
			c->center = Vector2d::Zero();
			c->radius = cylinder_radius;
			prim_prim_info.contact_patch._intersection_curves.push_back(c);

			prim_prim_info.min_distance = fmin(endA_distance, endB_distance);
			return;
		}
		// last case: point contact with circular edge
		{
			prim_prim_info.type = ContactType::POINT;
			double endA_distance = (endA_world - plane_point_world).dot(plane_normal_world);
			double endB_distance = (endB_world - plane_point_world).dot(plane_normal_world);
			Vector3d pt_dir = (plane_normal_world.cross(cylinder_axis)).cross(cylinder_axis);
			pt_dir /= pt_dir.norm();
			Vector3d close_pt;
			if(endA_distance < endB_distance) {
				// pt contact on face A
				close_pt = endA_world + cylinder_radius * pt_dir;
			} else {
				// pt contact on face B
				close_pt = endB_world + cylinder_radius * pt_dir;
			}
			prim_prim_info.contact_points.push_back(close_pt);
			prim_prim_info.min_distance = (close_pt - plane_point_world).dot(plane_normal_world);
			return;
		}
	}

	// utility function for planeBox distance computation
	static void addBoxPlaneSurfaceContactPatch(
		PrimPrimContactInfo& prim_prim_info,
		Eigen::Affine3d& contact_patch_tf,
		const Eigen::Vector3d& vertex1, // edges are assumed to be 1-2, 1-3, 2-4, 3-4
		const Eigen::Vector3d& vertex2,
		const Eigen::Vector3d& vertex3,
		const Eigen::Vector3d& vertex4
	) {
		prim_prim_info.contact_points.push_back(vertex1);
		prim_prim_info.contact_points.push_back(vertex2);
		prim_prim_info.contact_points.push_back(vertex3);
		prim_prim_info.contact_points.push_back(vertex4);
		prim_prim_info.contact_patch._interior_point =
					GeometryUtils::centroidOfPoints(prim_prim_info.contact_points);
		contact_patch_tf.translation() = prim_prim_info.contact_patch._interior_point;
		prim_prim_info.contact_patch.max_extent = (vertex1 - vertex4).norm();
		prim_prim_info.contact_patch._line_segments.push_back(LineSegment2D(
				GeometryUtils::Point2DFromPoint3DPlaneTransform(vertex1, contact_patch_tf),
				GeometryUtils::Point2DFromPoint3DPlaneTransform(vertex2, contact_patch_tf)
		));
		prim_prim_info.contact_patch._line_segments.push_back(LineSegment2D(
				GeometryUtils::Point2DFromPoint3DPlaneTransform(vertex1, contact_patch_tf),
				GeometryUtils::Point2DFromPoint3DPlaneTransform(vertex3, contact_patch_tf)
		));
		prim_prim_info.contact_patch._line_segments.push_back(LineSegment2D(
				GeometryUtils::Point2DFromPoint3DPlaneTransform(vertex2, contact_patch_tf),
				GeometryUtils::Point2DFromPoint3DPlaneTransform(vertex4, contact_patch_tf)
		));
		prim_prim_info.contact_patch._line_segments.push_back(LineSegment2D(
				GeometryUtils::Point2DFromPoint3DPlaneTransform(vertex3, contact_patch_tf),
				GeometryUtils::Point2DFromPoint3DPlaneTransform(vertex4, contact_patch_tf)
		));
	}

	void PrimPrimDistance::distancePlaneBox(
		PrimPrimContactInfo& prim_prim_info,
		const PlanePrimitive& plane, Eigen::Affine3d planeInWorld,
		const BoxPrimitive& box, Eigen::Affine3d boxInWorld
	) {
		assert(box._props != NULL);
		assert(plane._props != NULL);

		prim_prim_info.clear();

		// convert plane point and normal to box frame
		Vector3d plane_point = plane._props->point;
		Vector3d plane_normal = plane._props->normal;
		Vector3d plane_point_world = planeInWorld*plane_point;
		Vector3d plane_normal_world = planeInWorld.linear()*plane_normal;
		Vector3d plane_point_box = boxInWorld.inverse()*plane_point_world;
		Vector3d plane_normal_box = boxInWorld.inverse().linear()*plane_normal_world;

		prim_prim_info.normal_dir = plane_normal_world;

		// determine tangent directions
		// first try the cylinder axis
		// std::cout <<"Cylinder rotation " << cylinderInWorld.linear() << std::endl;
		Vector3d worldx_axis(1.0,0.0,0.0);
		double testaxis2_normal_proj = worldx_axis.dot(plane_normal_world);
		Vector3d worldy_axis(0.0,1.0,0.0);
		double testaxis3_normal_proj = worldy_axis.dot(plane_normal_world);
		if(abs(testaxis2_normal_proj) < 0.999) { // test world x axis next
			prim_prim_info.constraint_dir1 = worldx_axis - testaxis2_normal_proj*plane_normal_world;
			prim_prim_info.constraint_dir1 /= prim_prim_info.constraint_dir1.norm();
		} else { // use world y axis
			prim_prim_info.constraint_dir1 = worldy_axis - testaxis3_normal_proj*plane_normal_world;
			prim_prim_info.constraint_dir1 /= prim_prim_info.constraint_dir1.norm();
		}
		prim_prim_info.constraint_dir2 = prim_prim_info.normal_dir.cross(prim_prim_info.constraint_dir1);

		// local points list
		Vector3d vppp(box._props->xlength/2, box._props->ylength/2, box._props->zlength/2);
		Vector3d vppn(box._props->xlength/2, box._props->ylength/2, -box._props->zlength/2);
		Vector3d vpnp(box._props->xlength/2, -box._props->ylength/2, box._props->zlength/2);
		Vector3d vpnn(box._props->xlength/2, -box._props->ylength/2, -box._props->zlength/2);
		Vector3d vnpp(-box._props->xlength/2, box._props->ylength/2, box._props->zlength/2);
		Vector3d vnpn(-box._props->xlength/2, box._props->ylength/2, -box._props->zlength/2);
		Vector3d vnnp(-box._props->xlength/2, -box._props->ylength/2, box._props->zlength/2);
		Vector3d vnnn(-box._props->xlength/2, -box._props->ylength/2, -box._props->zlength/2);
		// get plane alignment measures
		double pnx = plane_normal_box.dot(Vector3d(1.0, 0, 0));
		double pny = plane_normal_box.dot(Vector3d(0.0, 1.0, 0));
		double pnz = plane_normal_box.dot(Vector3d(0.0, 0.0, 1.0));

		// check for surface contact:
		{
			Affine3d contact_patch_tf;
			contact_patch_tf.linear().col(0) = prim_prim_info.constraint_dir1;
			contact_patch_tf.linear().col(1) = prim_prim_info.constraint_dir2;
			contact_patch_tf.linear().col(2) = plane_normal_world;

			// TODO: add contact patch
			// static_assert(false, "TODO: add contact patch");
			if(pnx < -0.999) {// contact with +x surface
				prim_prim_info.type = ContactType::SURFACE;
				prim_prim_info.min_distance = plane_point_box[0] - box._props->xlength/2;
				addBoxPlaneSurfaceContactPatch(
					prim_prim_info,
					contact_patch_tf,
					boxInWorld*vppp,
					boxInWorld*vppn,
					boxInWorld*vpnp,
					boxInWorld*vpnn
				);
				return;
			} else if(pnx > 0.999) {// contact with -x surface
				prim_prim_info.type = ContactType::SURFACE;
				prim_prim_info.min_distance = - box._props->xlength/2 - plane_point_box[0];
				addBoxPlaneSurfaceContactPatch(
					prim_prim_info,
					contact_patch_tf,
					boxInWorld*vnpp,
					boxInWorld*vnpn,
					boxInWorld*vnnp,
					boxInWorld*vnnn
				);
				return;
			} else if(pny < -0.999) {// contact with +y surface
				prim_prim_info.type = ContactType::SURFACE;
				prim_prim_info.min_distance = plane_point_box[1] - box._props->ylength/2;
				addBoxPlaneSurfaceContactPatch(
					prim_prim_info,
					contact_patch_tf,
					boxInWorld*vppp,
					boxInWorld*vppn,
					boxInWorld*vnpp,
					boxInWorld*vnpn
				);
				return;
			} else if(pny > 0.999) {// contact with -y surface
				prim_prim_info.type = ContactType::SURFACE;
				prim_prim_info.min_distance = - box._props->ylength/2 - plane_point_box[1];
				addBoxPlaneSurfaceContactPatch(
					prim_prim_info,
					contact_patch_tf,
					boxInWorld*vpnp,
					boxInWorld*vpnn,
					boxInWorld*vnnp,
					boxInWorld*vnnn
				);
				return;
			} else if(pnz < -0.999) {// contact with +z surface
				prim_prim_info.type = ContactType::SURFACE;
				prim_prim_info.min_distance = plane_point_box[2] - box._props->zlength/2;
				addBoxPlaneSurfaceContactPatch(
					prim_prim_info,
					contact_patch_tf,
					boxInWorld*vppp,
					boxInWorld*vpnp,
					boxInWorld*vnpp,
					boxInWorld*vnnp
				);
				return;
			} else if(pnz > 0.999) {// contact with -z surface
				prim_prim_info.type = ContactType::SURFACE;
				prim_prim_info.min_distance = - box._props->zlength/2 - plane_point_box[2];
				addBoxPlaneSurfaceContactPatch(
					prim_prim_info,
					contact_patch_tf,
					boxInWorld*vppn,
					boxInWorld*vpnn,
					boxInWorld*vnpn,
					boxInWorld*vnnn
				);
				return;
			}
		}
		// check for edge contact
		{
			if(abs(pny) < 0.001) { //TODO: should be relative to the edge size
				// edge contact with plane normal in X-Z plane
				if(pnx < 0 && pnz < 0) {
					prim_prim_info.type = ContactType::LINE;
					prim_prim_info.min_distance = -plane_normal_box.dot(
							Vector3d(plane_point_box[0], 0, plane_point_box[2]) -
							Vector3d(box._props->xlength/2, 0, box._props->zlength/2)
						);
					prim_prim_info.contact_points.push_back(boxInWorld*vppp);
					prim_prim_info.contact_points.push_back(boxInWorld*vpnp);
					return;
				} else if (pnx < 0 && pnz > 0) {
					prim_prim_info.type = ContactType::LINE;
					prim_prim_info.min_distance = -plane_normal_box.dot(
							Vector3d(plane_point_box[0], 0, plane_point_box[2]) -
							Vector3d(box._props->xlength/2, 0, -box._props->zlength/2)
						);
					prim_prim_info.contact_points.push_back(boxInWorld*vppn);
					prim_prim_info.contact_points.push_back(boxInWorld*vpnn);
					return;
				} else if (pnx > 0 && pnz > 0) {
					prim_prim_info.type = ContactType::LINE;
					prim_prim_info.min_distance = -plane_normal_box.dot(
							Vector3d(plane_point_box[0], 0, plane_point_box[2]) -
							Vector3d(-box._props->xlength/2, 0, -box._props->zlength/2)
						);
					prim_prim_info.contact_points.push_back(boxInWorld*vnpn);
					prim_prim_info.contact_points.push_back(boxInWorld*vnnn);
					return;
				} else if (pnx > 0 && pnz < 0) {
					prim_prim_info.type = ContactType::LINE;
					prim_prim_info.min_distance = -plane_normal_box.dot(
							Vector3d(plane_point_box[0], 0, plane_point_box[2]) -
							Vector3d(-box._props->xlength/2, 0, box._props->zlength/2)
						);
					prim_prim_info.contact_points.push_back(boxInWorld*vnpp);
					prim_prim_info.contact_points.push_back(boxInWorld*vnnp);
					return;
				} else {
					std::cerr << "pny = 0; Should not be here" << std::endl;
				}
			} else if(abs(pnx) < 0.001) {
				// edge contact with plane normal in Y-Z plane
				if(pny < 0 && pnz < 0) {
					prim_prim_info.type = ContactType::LINE;
					prim_prim_info.min_distance = -plane_normal_box.dot(
							Vector3d(0, plane_point_box[1], plane_point_box[2]) -
							Vector3d(0, box._props->ylength/2, box._props->zlength/2)
						);
					prim_prim_info.contact_points.push_back(boxInWorld*vppp);
					prim_prim_info.contact_points.push_back(boxInWorld*vnpp);
					return;
				} else if (pny < 0 && pnz > 0) {
					prim_prim_info.type = ContactType::LINE;
					prim_prim_info.min_distance = -plane_normal_box.dot(
							Vector3d(0, plane_point_box[1], plane_point_box[2]) -
							Vector3d(0, box._props->ylength/2, -box._props->zlength/2)
						);
					prim_prim_info.contact_points.push_back(boxInWorld*vppn);
					prim_prim_info.contact_points.push_back(boxInWorld*vnpn);
					return;
				} else if (pny > 0 && pnz > 0) {
					prim_prim_info.type = ContactType::LINE;
					prim_prim_info.min_distance = -plane_normal_box.dot(
							Vector3d(0, plane_point_box[1], plane_point_box[2]) -
							Vector3d(0, -box._props->ylength/2, -box._props->zlength/2)
						);
					prim_prim_info.contact_points.push_back(boxInWorld*vpnn);
					prim_prim_info.contact_points.push_back(boxInWorld*vnnn);
					return;
				} else if (pny > 0 && pnz < 0) {
					prim_prim_info.type = ContactType::LINE;
					prim_prim_info.min_distance = -plane_normal_box.dot(
							Vector3d(0, plane_point_box[1], plane_point_box[2]) -
							Vector3d(0, -box._props->ylength/2, box._props->zlength/2)
						);
					prim_prim_info.contact_points.push_back(boxInWorld*vpnp);
					prim_prim_info.contact_points.push_back(boxInWorld*vnnp);
					return;
				} else {
					std::cerr << "pnx = 0; Should not be here" << std::endl;
				}
			} else if(abs(pnz) < 0.001) {
				// edge contact with plane normal in Y-Z plane
				if(pnx < 0 && pny < 0) {
					prim_prim_info.type = ContactType::LINE;
					prim_prim_info.min_distance = -plane_normal_box.dot(
							Vector3d(plane_point_box[0], plane_point_box[1], 0) -
							Vector3d(box._props->xlength/2, box._props->ylength/2, 0)
						);
					prim_prim_info.contact_points.push_back(boxInWorld*vppp);
					prim_prim_info.contact_points.push_back(boxInWorld*vppn);
					return;
				} else if (pnx < 0 && pny > 0) {
					prim_prim_info.type = ContactType::LINE;
					prim_prim_info.min_distance = -plane_normal_box.dot(
							Vector3d(plane_point_box[0], plane_point_box[1], 0) -
							Vector3d(box._props->xlength/2, -box._props->ylength/2, 0)
						);
					prim_prim_info.contact_points.push_back(boxInWorld*vpnp);
					prim_prim_info.contact_points.push_back(boxInWorld*vpnn);
					return;
				} else if (pnx > 0 && pny > 0) {
					prim_prim_info.type = ContactType::LINE;
					prim_prim_info.min_distance = -plane_normal_box.dot(
							Vector3d(plane_point_box[0], plane_point_box[1], 0) -
							Vector3d(-box._props->xlength/2, -box._props->ylength/2, 0)
						);
					prim_prim_info.contact_points.push_back(boxInWorld*vnnp);
					prim_prim_info.contact_points.push_back(boxInWorld*vnnn);
					return;
				} else if (pnx > 0 && pny < 0) {
					prim_prim_info.type = ContactType::LINE;
					prim_prim_info.min_distance = -plane_normal_box.dot(
							Vector3d(plane_point_box[0], plane_point_box[1], 0) -
							Vector3d(-box._props->xlength/2, box._props->ylength/2, 0)
						);
					prim_prim_info.contact_points.push_back(boxInWorld*vnpp);
					prim_prim_info.contact_points.push_back(boxInWorld*vnpn);
					return;
				} else {
					std::cerr << "pnz = 0; Should not be here" << std::endl;
				}
			}
		}
		// check for point contact
		{
			if(pnx < 0 && pny < 0 && pnz < 0) {
				prim_prim_info.type = ContactType::POINT;
				prim_prim_info.min_distance = -plane_normal_box.dot(plane_point_box - vppp);
				prim_prim_info.contact_points.push_back(boxInWorld*vppp);
				return;
			} else if(pnx < 0 && pny > 0 && pnz < 0) {
				prim_prim_info.type = ContactType::POINT;
				prim_prim_info.min_distance = -plane_normal_box.dot(plane_point_box - vpnp);
				prim_prim_info.contact_points.push_back(boxInWorld*vpnp);
				return;
			} else if(pnx > 0 && pny > 0 && pnz < 0) {
				prim_prim_info.type = ContactType::POINT;
				prim_prim_info.min_distance = -plane_normal_box.dot(plane_point_box - vnnp);
				prim_prim_info.contact_points.push_back(boxInWorld*vnnp);
				return;
			} else if(pnx > 0 && pny < 0 && pnz < 0) {
				prim_prim_info.type = ContactType::POINT;
				prim_prim_info.min_distance = -plane_normal_box.dot(plane_point_box - vnpp);
				prim_prim_info.contact_points.push_back(boxInWorld*vnpp);
				return;
			} else if(pnx < 0 && pny < 0 && pnz > 0) {
				prim_prim_info.type = ContactType::POINT;
				prim_prim_info.min_distance = -plane_normal_box.dot(plane_point_box - vppn);
				prim_prim_info.contact_points.push_back(boxInWorld*vppn);
				return;
			} else if(pnx < 0 && pny > 0 && pnz > 0) {
				prim_prim_info.type = ContactType::POINT;
				prim_prim_info.min_distance = -plane_normal_box.dot(plane_point_box - vpnn);
				prim_prim_info.contact_points.push_back(boxInWorld*vpnn);
				return;
			} else if(pnx > 0 && pny > 0 && pnz > 0) {
				prim_prim_info.type = ContactType::POINT;
				prim_prim_info.min_distance = -plane_normal_box.dot(plane_point_box - vnnn);
				prim_prim_info.contact_points.push_back(boxInWorld*vnnn);
				return;
			} else if(pnx > 0 && pny < 0 && pnz > 0) {
				prim_prim_info.type = ContactType::POINT;
				prim_prim_info.min_distance = -plane_normal_box.dot(plane_point_box - vnpn);
				prim_prim_info.contact_points.push_back(boxInWorld*vnpn);
				return;
			}
		}
		std::cerr << "Unknown contact type for plane:box" << std::endl;
	}


	// utility function for adding contact patch between plane and polygon
	// DOES NOT SET THE MAXIMUM EXTENT
	static void addPolygonPlaneSurfaceContactPatch(
		PrimPrimContactInfo& prim_prim_info,
		Eigen::Affine3d& contact_patch_tf,
		const std::vector<Eigen::Vector3d> vertices // assumed to be in order, clockwise or anti-clockwise
	) {
		for(auto vertex: vertices) {
			prim_prim_info.contact_points.push_back(vertex);
		}
		prim_prim_info.contact_patch._interior_point =
					GeometryUtils::centroidOfPoints(prim_prim_info.contact_points);
		contact_patch_tf.translation() = prim_prim_info.contact_patch._interior_point;
		for(uint i = 0; i < vertices.size()-1; i++) {
			prim_prim_info.contact_patch._line_segments.push_back(LineSegment2D(
				GeometryUtils::Point2DFromPoint3DPlaneTransform(vertices[i], contact_patch_tf),
				GeometryUtils::Point2DFromPoint3DPlaneTransform(vertices[i+1], contact_patch_tf)
			));
		}
		// add last edge
		prim_prim_info.contact_patch._line_segments.push_back(LineSegment2D(
			GeometryUtils::Point2DFromPoint3DPlaneTransform(vertices[vertices.size()-1], contact_patch_tf),
			GeometryUtils::Point2DFromPoint3DPlaneTransform(vertices[0], contact_patch_tf)
		));
	}

	void PrimPrimDistance::distancePlanePyramid(
		PrimPrimContactInfo& prim_prim_info,
		const PlanePrimitive& plane, Eigen::Affine3d planeInWorld,
		const PyramidPrimitive& pyramid, Eigen::Affine3d pyramidInWorld
	) {
		assert(pyramid._props != NULL);
		assert(plane._props != NULL);

		prim_prim_info.clear();

		// convert plane point and normal to box frame
		Vector3d plane_point = plane._props->point;
		Vector3d plane_normal = plane._props->normal;
		Vector3d plane_point_world = planeInWorld*plane_point;
		Vector3d plane_normal_world = planeInWorld.linear()*plane_normal;
		Vector3d plane_point_pyramid = pyramidInWorld.inverse()*plane_point_world;
		Vector3d plane_normal_pyramid = pyramidInWorld.inverse().linear()*plane_normal_world;

		prim_prim_info.normal_dir = plane_normal_world;

		// determine tangent directions
		Vector3d worldx_axis(1.0,0.0,0.0);
		double testaxis2_normal_proj = worldx_axis.dot(plane_normal_world);
		Vector3d worldy_axis(0.0,1.0,0.0);
		double testaxis3_normal_proj = worldy_axis.dot(plane_normal_world);
		if(abs(testaxis2_normal_proj) < 0.999) { // test world x axis next
			prim_prim_info.constraint_dir1 = worldx_axis - testaxis2_normal_proj*plane_normal_world;
			prim_prim_info.constraint_dir1 /= prim_prim_info.constraint_dir1.norm();
		} else { // use world y axis
			prim_prim_info.constraint_dir1 = worldy_axis - testaxis3_normal_proj*plane_normal_world;
			prim_prim_info.constraint_dir1 /= prim_prim_info.constraint_dir1.norm();
		}
		prim_prim_info.constraint_dir2 = prim_prim_info.normal_dir.cross(prim_prim_info.constraint_dir1);

		Affine3d surface_contact_patch_tf;
		surface_contact_patch_tf.linear().col(0) = prim_prim_info.constraint_dir1;
		surface_contact_patch_tf.linear().col(1) = prim_prim_info.constraint_dir2;
		surface_contact_patch_tf.linear().col(2) = plane_normal_world;

		// compute plane normal product = dot product between base normal and plane normal
		Vector3d apex_in_pyramid(0, 0, pyramid._props->height);
		Vector3d base_normal(0, 0, -1);
		double plane_norm_dot = plane_normal_pyramid.dot(base_normal);
		double plane_norm_angle = acos(plane_norm_dot);

		// compute apex distance to plane
		double apex_plane_dist = GeometryUtils::distancePointToPlane(
			apex_in_pyramid,
			plane_point_pyramid,
			plane_normal_pyramid
		);

		// - if acos(plane normal product) < side angle, return point contact with apex
		// std::cout << plane_norm_angle << " " << pyramid.sideEdgeAngle() << std::endl;
		if(plane_norm_angle < M_PI/2 - pyramid.sideEdgeAngle()) {
			//return point contact with apex
			prim_prim_info.type = ContactType::POINT;
			prim_prim_info.min_distance = apex_plane_dist;
			prim_prim_info.contact_points.push_back(pyramidInWorld*apex_in_pyramid);
			return;
		}
		// - if plane normal == -base normal, return plane contact with base
		if(1+plane_norm_dot < 0.001) {
			//return plane contact with base
			prim_prim_info.type = ContactType::SURFACE;
			prim_prim_info.min_distance = apex_plane_dist - pyramid._props->height;
			prim_prim_info.contact_patch.max_extent =
				(pyramid._props->num_sides_base % 2 == 0)? 2*pyramid.circumRadius():
					(pyramid.circumRadius() + pyramid.incircleRadius());
			std::vector<Eigen::Vector3d> verticesInWorld;
			for(const auto& vertex: pyramid._base_points) {
				verticesInWorld.push_back(pyramidInWorld*vertex);
			}
			addPolygonPlaneSurfaceContactPatch(
				prim_prim_info,
				surface_contact_patch_tf,
				verticesInWorld
			);
			return;
		}

		// project plane normal on pyramid X-Y plane
		Vector3d plane_normal_base_proj = plane_normal_pyramid - base_normal*(base_normal.dot(plane_normal_pyramid));
		plane_normal_base_proj /= plane_normal_base_proj.norm();

		// set plane angle = acos(- x component of projected plane normal)
		double plane_angle = atan2(-plane_normal_base_proj[1], -plane_normal_base_proj[0]);
		if(plane_angle < 0) plane_angle += 2*M_PI;
		// std::cout << "plane_angle " << plane_angle << std::endl;

		// get side_number = (plane angle - 0.5*included angle) % included angle
		int side_number = round((plane_angle - 0.5*pyramid.includedAngle()) / pyramid.includedAngle());
		int vertex_number = round((plane_angle) / pyramid.includedAngle());
		double side_residual_angle = (plane_angle - 0.5*pyramid.includedAngle()) - pyramid.includedAngle()*side_number;
		double vertex_residual_angle = (plane_angle) - pyramid.includedAngle()*vertex_number;

		double vertex_plane_dist = GeometryUtils::distancePointToPlane(
			pyramid._base_points[vertex_number],
			plane_point_pyramid,
			plane_normal_pyramid
		);
		// std::cout << "Side residual angle: " << side_residual_angle << std::endl;
		// std::cout << "Side number: " << side_number << std::endl;
		// std::cout << "Vertex residual angle: " << vertex_residual_angle << std::endl;
		// std::cout << "Vertex number: " << vertex_number << std::endl;
		// if side_number < thresh,
		if(abs(side_residual_angle) < 0.001) {
			if(abs(plane_norm_angle - (M_PI/2 - pyramid.sideFaceAngle())) < 0.001) {
				//return plane contact with corresponding side face
				prim_prim_info.type = ContactType::SURFACE;
				prim_prim_info.min_distance = apex_plane_dist;
				prim_prim_info.contact_patch.max_extent =
					sqrt(pyramid.incircleRadius()*pyramid.incircleRadius() +
						pyramid._props->height*pyramid._props->height);
				std::vector<Eigen::Vector3d> verticesInWorld;
				verticesInWorld.push_back(pyramidInWorld*apex_in_pyramid);
				verticesInWorld.push_back(pyramidInWorld*pyramid._base_points[side_number]);
				if(side_number < pyramid._props->num_sides_base - 1) {
					verticesInWorld.push_back(pyramidInWorld*pyramid._base_points[side_number+1]);
				} else {
					verticesInWorld.push_back(pyramidInWorld*pyramid._base_points[0]);
				}
				addPolygonPlaneSurfaceContactPatch(
					prim_prim_info,
					surface_contact_patch_tf,
					verticesInWorld
				);
				return;
			} else {
				// return edge contact with corresponding base edge
				prim_prim_info.type = ContactType::LINE;
				prim_prim_info.min_distance = vertex_plane_dist;
				prim_prim_info.contact_points.push_back(pyramidInWorld*pyramid._base_points[side_number]);
				if(side_number < pyramid._props->num_sides_base - 1) {
					prim_prim_info.contact_points.push_back(pyramidInWorld*pyramid._base_points[side_number+1]);
				} else {
					prim_prim_info.contact_points.push_back(pyramidInWorld*pyramid._base_points[0]);
				}
			}
		} else if(abs(vertex_residual_angle) < 0.001) {
			if(abs(plane_norm_angle - (M_PI/2 - pyramid.sideEdgeAngle())) < 0.001) {
				// return edge contact with corresponding vertical edge
				prim_prim_info.type = ContactType::LINE;
				prim_prim_info.min_distance = apex_plane_dist;
				prim_prim_info.contact_points.push_back(pyramidInWorld*apex_in_pyramid);
				prim_prim_info.contact_points.push_back(pyramidInWorld*pyramid._base_points[vertex_number]);
			} else {
				// return point contact with corresponding base vertex
				prim_prim_info.type = ContactType::POINT;
				prim_prim_info.min_distance = vertex_plane_dist;
				prim_prim_info.contact_points.push_back(pyramidInWorld*pyramid._base_points[vertex_number]);
				return;
			}
		} else {
			if(abs(vertex_plane_dist) < abs(apex_plane_dist)) {
				//return point contact with base vertex
				prim_prim_info.type = ContactType::POINT;
				prim_prim_info.min_distance = vertex_plane_dist;
				prim_prim_info.contact_points.push_back(pyramidInWorld*pyramid._base_points[vertex_number]);
				return;
			} else {
				//return point contact with apex
				prim_prim_info.type = ContactType::POINT;
				prim_prim_info.min_distance = apex_plane_dist;
				prim_prim_info.contact_points.push_back(pyramidInWorld*apex_in_pyramid);
				return;
			}
		}
	}

	void PrimPrimDistance::distanceComposite1PkNCapsule(
		PrimPrimContactInfo& prim_prim_info,
		const Composite1PkN& composite, Eigen::Affine3d compositeInWorld,
		const CapsulePrimitive& capsule, Eigen::Affine3d capsuleInWorld
	) {
		distancePrimitivePrimitive(prim_prim_info,
									composite._positivePrimitive,
									compositeInWorld,
									&capsule,
									capsuleInWorld);
		const uint pos_prim_contact_pts_size = prim_prim_info.contact_points.size();

		auto set_keep_inds_from_exclude_inds = [](std::vector<uint>& keep_inds,
			const std::vector<uint>& exclude_inds, uint size) {
			uint eind = 0;
			for(uint i = 0; i < size; i++) {
				if(eind < exclude_inds.size() && i == exclude_inds[eind]) {
					eind++;
				} else {
					keep_inds.push_back(i);
				}
			}
		};
		// Handle negative primitives
		std::vector<uint> pos_inds_to_remove;
		for(const auto& neg_prim_info: composite._negativePrimitives) {
			PrimPrimContactInfo neg_dist_info;
			Affine3d neg_prim_in_world = compositeInWorld*neg_prim_info.transform_in_parent;
			distancePrimitivePrimitive(neg_dist_info, neg_prim_info.prim,
									neg_prim_in_world, &capsule, capsuleInWorld);
			// If all negative primitives report full penetration, return the distance info
			// from the positive primitive check.
			if(fabs(neg_dist_info.min_distance - PrimitiveAlgorithmicConstants::DISTANCE_FULL_PENETRATION) < 1e-3) {
				continue;
			}

			// Filter out the contact points from positive primitive which lie inside the
			// negative primitive
			for(uint i = 0; i < pos_prim_contact_pts_size; i++) {
				auto pos_pt = prim_prim_info.contact_points[i];
				PointPrimDistanceInfo pos_pt_in_neg_prim_info =
					PointPrimDistance::distancePointPrimitive(
						pos_pt, neg_prim_info.prim, neg_prim_in_world);
				if(pos_pt_in_neg_prim_info.distance > 0) {
					pos_inds_to_remove.push_back(i);
				}
			}

			// If negative primitive reports undefined contact type, ignore the distance
			// info result. This can happen if it is physically impossible for the positive
			// primitive to fit in the negative primitive.
			if(neg_dist_info.type != ContactType::UNDEFINED) {
				// If a negative primitive reports partial penetration, filter contact
				// points which lie inside the positive primitive
				std::vector<uint> keep_neg_pts_indices;
				for(uint j = 0; j < neg_dist_info.contact_points.size(); j++) {
					PointPrimDistanceInfo neg_pt_in_pos_prim_info =
						PointPrimDistance::distancePointPrimitive(
							neg_dist_info.contact_points[j],
							composite._positivePrimitive, compositeInWorld);
					if(neg_pt_in_pos_prim_info.distance <= -1e-6) {
						keep_neg_pts_indices.push_back(j);
					} else {
						// std::cout << "Remove point " << neg_dist_info.contact_points[j].transpose() << std::endl;
					}
				}
				neg_dist_info.filterContactPoints(keep_neg_pts_indices);

				// set contact type to CONCAVE and insert points from neg prim contact
				prim_prim_info.addOtherContactInfo(neg_dist_info);
			}

			// check for contact with intersection edge
			// TODO: if the test primitive is completely contained inside the positive,
			// primitive, the intersection check can be omitted
			for(const auto* i_edge: neg_prim_info.intersection_edges) {
				PrimPrimContactInfo edge_dist_info =
					i_edge->primDistance(&capsule, compositeInWorld.inverse()*capsuleInWorld);
				// filter points from contact with intersection edge where the normal points
				// outwards from the negative primitive or inwards from the positive primitive
				// std::cout << "Edge points " << edge_dist_info.contact_points.size() << std::endl;
				std::vector<uint> edge_pts_to_keep;
				for(uint pt_ind = 0; pt_ind < edge_dist_info.contact_points.size(); pt_ind++) {
					Vector3d normal;
					// check with positive primitive
					if(!PrimNormal::primNormal(normal,
							edge_dist_info.contact_points[pt_ind],
							composite._positivePrimitive, compositeInWorld)) {
						throw(std::runtime_error("Pt not on positive prim"));
					}
					if(normal.dot(edge_dist_info.normal_dirs[pt_ind]) < -1e-3) {
						// std::cout << "Remove edge point " << edge_dist_info.contact_points[pt_ind].transpose()
						//           << " normal(1): " << edge_dist_info.normal_dirs[pt_ind].transpose() << std::endl;
						continue;
					}
					// check with negative primitive
					if(!PrimNormal::primNormal(normal,
							edge_dist_info.contact_points[pt_ind],
							neg_prim_info.prim, neg_prim_in_world)) {
						throw(std::runtime_error("Pt not on negative prim"));
					}
					if(normal.dot(edge_dist_info.normal_dirs[pt_ind]) < -1e-3) {
						// std::cout << "Remove edge point " << edge_dist_info.contact_points[pt_ind].transpose()
						//           << " normal(2): " << edge_dist_info.normal_dirs[pt_ind].transpose() << std::endl;
						continue;
					}
					edge_pts_to_keep.push_back(pt_ind);
				}
				edge_dist_info.filterContactPoints(edge_pts_to_keep);

				// add remaining points to prim_prim_info
				prim_prim_info.addOtherContactInfo(edge_dist_info);
			}
		}

		// remove all positive contact points in pos_inds_to_remove
		// have to be careful, the minimum distance should be adjusted accordingly
		// e.g. capsule hovering right above the hole
		// hopefully, this is handled by the distance from the intersection edge?
		if(!pos_inds_to_remove.empty()) {
			std::vector<uint> pos_inds_to_keep;
			set_keep_inds_from_exclude_inds(pos_inds_to_keep, pos_inds_to_remove, prim_prim_info.contact_points.size());
			prim_prim_info.filterContactPoints(pos_inds_to_keep);
		}
	}
}
