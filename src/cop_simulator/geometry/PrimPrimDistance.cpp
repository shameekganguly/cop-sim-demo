// PrimPrimDistance.cpp

#include <iostream>
#include "Primitive.h"

using namespace Eigen;

namespace Sai2COPSim {
	void PrimPrimContactInfo::flipNormal() {
		normal_dir *= -1.0;
		constraint_dir2 *= -1.0;
	}

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
		} else {
			std::cerr << "PrimA type: " << primA->_type << std::endl;
			std::cerr << "PrimB type: " << primB->_type << std::endl;
			throw(std::runtime_error("Distance between these primitives is not implemented."));
		}
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
			if(abs(z1) > 0.5*capAlength && abs(z2) > 0.5*capBlength) {
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
			Circle* c = new Circle();
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
}
