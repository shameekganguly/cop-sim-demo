// IntersectionEdgeDistance.cpp

#include <gte/Mathematics/UIntegerAP32.h>
#include <gte/Mathematics/BSRational.h>
#include <gte/Mathematics/RootsPolynomial.h>
#include "IntersectionEdge.h"

using namespace Eigen;

namespace Sai2COPSim {

Vector3d Circle3DDistance::pointDist(const Vector3d& pointInParent, const Circle3D& circle) {
	Vector3d pointInCircle = circle.circle_frame_in_parent.inverse()*pointInParent;
	if(pointInCircle.segment<2>(0).norm() < 1e-3) {
		// point lies on the circle's axis
		return circle.circle_frame_in_parent*Vector3d(circle.radius, 0, 0);
	}
	double gamma = atan2(pointInCircle(1), pointInCircle(0));
	return circle.circle_frame_in_parent*(circle.radius*Vector3d(cos(gamma), sin(gamma), 0));
}

PrimPrimContactInfo Circle3DDistance::capsuleDist(const CapsulePrimitive& cap,
												const Circle3D& circle,
												const Eigen::Affine3d& primInParent) {
	assert(cap._props != NULL);

	PrimPrimContactInfo ret_info;
	Affine3d prim_in_circle = circle.circle_frame_in_parent.inverse()*primInParent;
	Vector3d cap_axis_local = prim_in_circle.linear().col(0);
	Vector3d cap_center_local = prim_in_circle.translation();
	Vector3d local_z(0,0,1);
	double tn = cap_axis_local.dot(local_z);

	double circle_rad = circle.radius;
	double cap_rad = cap._props->radius;
	double cap_half_len = cap._props->length*0.5;

	// early return if distance is very large
	if (cap_center_local.norm() > 1.5*(circle_rad + cap_half_len + cap_rad)) {
		return ret_info;
	}

	// Call ONLY ONCE to add both tangent directions.
	auto add_constraint_dirs = [](Vector3d& normal, PrimPrimContactInfo& add_info) {
		if(fabs(normal(2)) < 0.999) {
			// use local z axis
			Vector3d dir1(0.0, 0.0, 1.0);
			dir1 -= normal(2)*normal;
			dir1 /= dir1.norm();
			add_info.constraint_dir1s.push_back(dir1);
		} else {
			// use local x axis
			Vector3d dir1(1.0, 0.0, 0.0);
			dir1 -= normal(0)*normal;
			dir1 /= dir1.norm();
			add_info.constraint_dir1s.push_back(dir1);
		}
		add_info.constraint_dir2s.push_back(normal.cross(add_info.constraint_dir1s.back()));
	};

	auto add_prim_prim_info_for_line_circle_pt =
		[&add_constraint_dirs, &cap_rad](Vector3d& circle_pt, Vector3d& line_pt, PrimPrimContactInfo& add_info) {
			// Add line point
			Vector3d normal = line_pt - circle_pt;
			double dist = normal.norm();
			assert(fabs(dist) > 1e-3);
			normal /= dist;
			dist -= cap_rad;
			add_info.contact_points.push_back(circle_pt);
			add_info.normal_dirs.push_back(normal);
			add_constraint_dirs(normal, add_info);
			add_info.distances.push_back(dist);
	};

	auto append_vectorB_to_vectorA_double = [](std::vector<double>& vectorA, std::vector<double>& vectorB) {
		vectorA.insert(vectorA.end(), vectorB.begin(), vectorB.end());
	};

	auto append_vectorB_to_vectorA_Vector3d = [](std::vector<Vector3d>& vectorA, std::vector<Vector3d>& vectorB) {
		vectorA.insert(vectorA.end(), vectorB.begin(), vectorB.end());
	};

	auto append_prim_info_into_ret_info =
		[&ret_info, &append_vectorB_to_vectorA_double, &append_vectorB_to_vectorA_Vector3d](PrimPrimContactInfo& add_info) {
			append_vectorB_to_vectorA_Vector3d(ret_info.contact_points, add_info.contact_points);
			append_vectorB_to_vectorA_Vector3d(ret_info.normal_dirs, add_info.normal_dirs);
			append_vectorB_to_vectorA_Vector3d(ret_info.constraint_dir1s, add_info.constraint_dir1s);
			append_vectorB_to_vectorA_Vector3d(ret_info.constraint_dir2s, add_info.constraint_dir2s);
			append_vectorB_to_vectorA_double(ret_info.distances, add_info.distances);
	};

	auto clamped_cap_z = [&cap_half_len] (double z) -> double {
		return fmax(-cap_half_len, fmin(cap_half_len, z));
	};

	// Get distance from left and right ends
	Vector3d left_end_local = cap_center_local - cap_half_len*cap_axis_local;
	Vector3d right_end_local = cap_center_local + cap_half_len*cap_axis_local;
	std::vector<PrimPrimContactInfo> end_contact_info; // 0 is left end, 1 is right end
	Vector3d circle_axissym_pt1(circle_rad, 0, 0);
	Vector3d circle_axissym_pt2(-circle_rad*cos(M_PI/3), circle_rad*sin(M_PI/3), 0);
	Vector3d circle_axissym_pt3(-circle_rad*cos(M_PI/3), -circle_rad*sin(M_PI/3), 0);
	auto add_prim_prim_info_for_line_circle_pt_cap_end =
		[&add_prim_prim_info_for_line_circle_pt, &cap_axis_local](Vector3d& circle_pt,
																Vector3d& end_pt,
																PrimPrimContactInfo& add_info,
																bool right_end) {
			// only add the contact point if the normal from the circle towards the end
			// point is directed towards the capsule axis. i.e. the ray from the circle
			// pt towards the end point should not intersect the side wall of the capsule
			if((end_pt - circle_pt).dot(cap_axis_local)*(right_end? -1: 1) > -2e-3) {
				add_prim_prim_info_for_line_circle_pt(circle_pt, end_pt, add_info);
			}
	};
	for(uint i = 0; i < 2; i++) {
		PrimPrimContactInfo end_info;
		Vector3d test_end_local = (i == 0)? left_end_local: right_end_local;
		if(test_end_local.segment<2>(0).norm() < 1e-3) {
			// axisymmetric case, add 3 points
			add_prim_prim_info_for_line_circle_pt_cap_end(circle_axissym_pt1, test_end_local, end_info, i);
			add_prim_prim_info_for_line_circle_pt_cap_end(circle_axissym_pt2, test_end_local, end_info, i);
			add_prim_prim_info_for_line_circle_pt_cap_end(circle_axissym_pt3, test_end_local, end_info, i);
		} else {
			Vector3d test_end_parent = circle.circle_frame_in_parent*test_end_local;
			Vector3d closest_circle_pt =
				circle.circle_frame_in_parent.inverse()*pointDist(test_end_parent, circle);
			// only add the contact point if the normal from the circle towards the end
			// point is directed towards the capsule axis. i.e. the ray from the circle
			// pt towards the end point should not intersect the side wall of the capsule
			add_prim_prim_info_for_line_circle_pt_cap_end(closest_circle_pt, test_end_local, end_info, i);
		}
		end_contact_info.push_back(end_info);
	}
	auto left_end_dist_info = end_contact_info.begin();
	auto right_end_dist_info = end_contact_info.rbegin();


	// -- Find closest point from line passing through capsule axis to circle --

	if(fabs(tn) < 1e-3) {
		// Capsule is parallel to circle plane
		double a_val = cap_center_local(2);
		double b_val = -cap_center_local.dot(cap_axis_local);
		Vector3d closest_line_pt_zaxis = cap_center_local + b_val*cap_axis_local;
		double line_zaxis_dist = closest_line_pt_zaxis.segment<2>(0).norm();
		std::vector<double> z_vals;
		std::vector<double> gamma_vals;

		if(line_zaxis_dist > circle_rad || fabs(line_zaxis_dist - circle_rad) < 1e-3) {
			// only one closest pt on circle
			z_vals.push_back(b_val);
		} else {
			double circle_chord_half_len = sqrt(circle_rad*circle_rad - line_zaxis_dist*line_zaxis_dist);
			double z1 = b_val + circle_chord_half_len;
			double clamped_z1 = clamped_cap_z(z1);
			double z2 = b_val - circle_chord_half_len;
			double clamped_z2 = clamped_cap_z(z2);
			// label: redundant_contact_pt_handling
			if((fabs(z1) > cap_half_len && fabs(z2) <= cap_half_len && fabs(clamped_z1 - z2) < circle_chord_half_len)) {
				// only add z2 to avoid reporting duplicate points
				z_vals.push_back(z2);
			} else if ((fabs(z1) <= cap_half_len && fabs(z2) > cap_half_len && fabs(clamped_z2 - z1) < circle_chord_half_len)) {
				// only add z2 to avoid reporting duplicate points
				z_vals.push_back(z1);
			} else {
				z_vals.push_back(b_val + circle_chord_half_len);
				z_vals.push_back(b_val - circle_chord_half_len);
			}
		}

		bool added_right_pt = false;
		bool added_left_pt = false;
		for(auto& z_val: z_vals) {
			if(z_val > cap_half_len || fabs(z_val - cap_half_len) < 1e-3) {
				if(!added_right_pt) {
					added_right_pt = true;
					append_prim_info_into_ret_info(*right_end_dist_info);
				}
			} else if(z_val < -cap_half_len || fabs(z_val + cap_half_len) < 1e-3) {
				if(!added_left_pt) {
					added_left_pt = true;
					append_prim_info_into_ret_info(*left_end_dist_info);
				}
			} else {
				Vector3d line_pt = cap_center_local + z_val*cap_axis_local;
				double gamma = atan2(line_pt(1), line_pt(0));
				Vector3d circle_pt(circle_rad*cos(gamma), circle_rad*sin(gamma), 0);
				add_prim_prim_info_for_line_circle_pt(circle_pt, line_pt, ret_info);
			}
		}

	} else {
		// Capsule is NOT parallel to circle plane
		Vector3d line_plane_pt = cap_center_local - (cap_center_local.dot(local_z)/tn)*cap_axis_local;

		double tq = line_plane_pt.norm();
		double tm = line_plane_pt.dot(cap_axis_local);

		bool axisymmetric = false;

		// Each 'a' value is the Z-intercept for a common normal from the line to the circle
		std::vector<double> sol_a_vec;

		double tnsqr_comp = (1 - tn*tn);
		double tnsqr = tn*tn;
		double tmsqr = tm*tm;
		double tqsqr = tq*tq;
		double Rsqr = circle_rad*circle_rad;

		// special case 1: line is parallel to Z axis
		if(fabs(tn) > 0.999) {
			sol_a_vec.push_back(0);
			axisymmetric = (fabs(tq) < 1e-3);

		// special case 2: line_plane_pt is very close to origin
		} else if(fabs(tq) < 1e-3) {
			// TODO: handle edge case where 2 identical contact points might be reported
			// one along the capsule axis, and one at the capsule end, similar to
			// \ref{redundant_contact_pt_handling} above
			sol_a_vec.push_back(circle_rad*sqrt(1 - tn*tn)/tn);
			sol_a_vec.push_back(-circle_rad*sqrt(1 - tn*tn)/tn);

		// special case 3: tm = 0
		} else if(fabs(tm) < 1e-3) {
			double temp = tqsqr - Rsqr*tnsqr_comp*tnsqr_comp;
			if(temp > 0 || fabs(temp) < 1e-3) {
				sol_a_vec.push_back(0);
			} else {
				// TODO: handle edge case where 2 identical contact points might be
				// reported, one along the capsule axis, and one at the capsule end,
				// similar to \ref{redundant_contact_pt_handling} above
				double tempa = -temp/(tnsqr*tnsqr_comp);
				sol_a_vec.push_back(sqrt(tempa));
				sol_a_vec.push_back(-sqrt(tempa));
				// Technically, a_val = 0 can also be a valid solution in this case.
				sol_a_vec.push_back(0);
			}
		} else {
			// std::cout << "Quartic\n";
			// set up quartic equation
			std::map<double, int32_t> roots;
			// NOTE: Using rationals is more accurate, but leads to 40x slowdown!
			// using gte::BSRational;
			// using gte::UIntegerAP32;
			// BSRational<UIntegerAP32> p4(tnsqr*tnsqr_comp);
			// BSRational<UIntegerAP32> p3(2*tm * tn*tn*tn);
			// BSRational<UIntegerAP32> p2(tqsqr - tmsqr*(1 + tnsqr) - Rsqr*tnsqr_comp*tnsqr_comp);
			// BSRational<UIntegerAP32> p1(-2*Rsqr*tm*tn*tnsqr_comp);
			// BSRational<UIntegerAP32> p0(-Rsqr*tmsqr*tnsqr);

			double p4(tnsqr*tnsqr_comp);
			double p3(2*tm * tn*tn*tn);
			double p2(tqsqr - tmsqr*(1 + tnsqr) - Rsqr*tnsqr_comp*tnsqr_comp);
			double p1(-2*Rsqr*tm*tn*tnsqr_comp);
			double p0(-Rsqr*tmsqr*tnsqr);

			// only returns real roots
			gte::RootsPolynomial<double>::SolveQuartic(p0, p1, p2, p3, p4, roots);
			if (roots.empty()) throw(std::runtime_error("Could not solve quartic."));

			// All real roots correspond to either minimum distance or maximum distance
			// from the line as the circle is traversed.
			// Ignore all solutions which are close to zero since they are explicitly
			// handled above.
			// Technically, repeated roots are not possible for the circle-line distance
			// problem unless a = 0. So we disregard the multiplicity of the solution.
			// TODO: handle edge case where 2 identical contact points might be reported
			// one along the capsule axis, and one at the capsule end, similar to
			// \ref{redundant_contact_pt_handling} above
			for (auto& it: roots) {
				double root = it.first;
				if(fabs(root) > 1e-3) {
					// std::cout << "Root a: " << root << std::endl;
					sol_a_vec.push_back(root);
				}
			}
		}

		bool added_right_pt = false;
		bool added_left_pt = false;
		if(axisymmetric) {
			double z_val_origin = -cap_center_local(2)*tn;
			if(z_val_origin > cap_half_len || fabs(z_val_origin - cap_half_len) < 1e-3) {
				if(!added_right_pt) {
					added_right_pt = true;
					append_prim_info_into_ret_info(*right_end_dist_info);
				}
			} else if(z_val_origin < -cap_half_len || fabs(z_val_origin + cap_half_len) < 1e-3) {
				if(!added_left_pt) {
					added_left_pt = true;
					append_prim_info_into_ret_info(*left_end_dist_info);
				}
			} else {
				Vector3d line_pt = cap_center_local + z_val_origin*cap_axis_local;
				add_prim_prim_info_for_line_circle_pt(circle_axissym_pt1, line_pt, ret_info);
				add_prim_prim_info_for_line_circle_pt(circle_axissym_pt2, line_pt, ret_info);
				add_prim_prim_info_for_line_circle_pt(circle_axissym_pt3, line_pt, ret_info);
			}
		} else {
			for(auto& a_val: sol_a_vec) {
				double b_val = a_val*tn - tm;
				Vector3d closest_line_pt = line_plane_pt + b_val*cap_axis_local;

				// Ignore a_val's where the closest line point lies on the Z-axis
				// This case should have been handled in the axisymmetric case above.
				// Usually these are additional extremum points returned from the quartic
				// solver. But they don't actually correspond to minimum or maximum
				// distance. They correspond to inflection points in the distance
				// function as we traverse the circle.
				if(closest_line_pt.segment<2>(0).norm() < 1e-3) {
					continue;
				}

				// std::cout << "closest_line_pt " << closest_line_pt.transpose() << std::endl;
				double z_val = (closest_line_pt - cap_center_local).dot(cap_axis_local);
				if(z_val > cap_half_len || fabs(z_val - cap_half_len) < 1e-3) {
					if(!added_right_pt) {
						added_right_pt = true;
						append_prim_info_into_ret_info(*right_end_dist_info);
					}
				} else if(z_val < -cap_half_len || fabs(z_val + cap_half_len) < 1e-3) {
					if(!added_left_pt) {
						added_left_pt = true;
						append_prim_info_into_ret_info(*left_end_dist_info);
					}
				} else {
					Vector3d circle_pt;
					if(fabs(a_val) < 1e-3) {
						// when a_val is zero, it is guaranteed that tm is also zero. So
						// c_val cannot be calculated as below
						circle_pt = line_plane_pt/line_plane_pt.norm()*circle_rad;
					} else {
						double c_val = a_val/(a_val*tnsqr_comp + tm*tn);
						// std::cout << "a: " << a_val << " b: " << b_val << " c: " << c_val << std::endl;
						circle_pt = (1 - c_val)*Vector3d(0, 0, a_val) + c_val*closest_line_pt;
					}
					// std::cout << "Circle pt " << circle_pt.transpose() << std::endl;
					assert(fabs(circle_pt.norm() - circle_rad) < 1e-3);
					add_prim_prim_info_for_line_circle_pt(circle_pt, closest_line_pt, ret_info);
				}
			}
			assert(ret_info.contact_points.size() < 3);
		}
	}
	if(ret_info.contact_points.size() == 0) {
		// throw(std::runtime_error("Did not generate any contact points"));
		return ret_info;
	}

	ret_info.type = ContactType::CONCAVE;
	ret_info.setMinDistanceFromDistances();

	// Transform points, normals, constraint dirs to parent frame
	for(auto& normal_dir: ret_info.normal_dirs) {
		normal_dir = circle.circle_frame_in_parent.linear()*normal_dir;
	}
	for(auto& constraint_dir: ret_info.constraint_dir1s) {
		constraint_dir = circle.circle_frame_in_parent.linear()*constraint_dir;
	}
	for(auto& constraint_dir: ret_info.constraint_dir2s) {
		constraint_dir = circle.circle_frame_in_parent.linear()*constraint_dir;
	}
	for(auto& pt: ret_info.contact_points) {
		pt = circle.circle_frame_in_parent*pt;
	}

	return ret_info;
}

}