// COPSolverExtended.cpp

#include <iostream>
#include "COPSolverExtended.h"
#include "lcp_solvers/LCPSolver.h"

using namespace Eigen;

namespace Sai2COPSim {

namespace {
bool isSurfaceContactSeparating(const Vector3d& rhs_lin_acc,
							    const Vector3d& rhs_ang_acc,
							    const Vector3d& omegaA,
							    const Vector3d& omegaB,
							    const ContactPatch& patch) {
	if(rhs_lin_acc[2] > -1e-8) {
		if(rhs_ang_acc.head<2>(0).norm() < 1e-8) {
			return true;
		}

		// pick a worse case point to test for penetration due to angular acceleration
		Vector3d faraway_test_pt = rhs_ang_acc.cross(Vector3d(0, 0, 1));
		faraway_test_pt *= patch.max_extent / faraway_test_pt.norm();
		if(!COPSolver::isAccelerationPenetrating(
			faraway_test_pt,
			rhs_lin_acc,
			rhs_ang_acc,
			omegaA,
			omegaB
		)) {
			return true;
		}
	}
	return false;
}

void getDisplacedMatricesForMultipleSurfaceContacts(
	Eigen::MatrixXd& A_disp, 			// 6n x 6n
	Eigen::VectorXd& rhs_disp, 			// 6n
	std::vector<Eigen::Vector3d>& lin_vels_disp,
	const std::vector<Eigen::Vector3d>& displacements,
	const Eigen::MatrixXd& A,			// 6n x 6n
	const Eigen::VectorXd& rhs,			// 6n
	const std::vector<Eigen::Vector3d>& omegaAs,
	const std::vector<Eigen::Vector3d>& omegaBs,
	const std::vector<Eigen::Vector3d>& lin_vels
) {
	const uint num_surfaces = displacements.size();
	const uint A_size = A.rows();
	A_disp = A;
	rhs_disp = rhs;
	for (uint i = 0; i < num_surfaces; i++) {
		if (displacements[i].norm() < 1e-6) {
			continue;
		}
		Eigen::Matrix3d cM = -crossMat(displacements[i]);
		A_disp.block(0, 6*i, A_size, 3) +=
			A_disp.block(0, 6*i + 3, A_size, 3)*cM.transpose();

		A_disp.block(6*i, 0, 3, A_size) +=
			cM*A_disp.block(6*i + 3, 0, 3, A_size);

		rhs_disp.segment<3>(6*i) += cM*rhs_disp.segment<3>(6*i + 3);
		rhs_disp.segment<3>(6*i) += omegaBs[i].cross(omegaBs[i].cross(displacements[i]))
				- omegaAs[i].cross(omegaAs[i].cross(displacements[i]));

		lin_vels_disp[i] = lin_vels[i] +
						   		(omegaBs[i] - omegaAs[i]).cross(displacements[i]);
	}
}

}  // namespace

ContactCOPSolution COPSolver::solveStartWithLastCOP(
	double friction_coeff,
	const Eigen::MatrixXd& A_constraint,
	const Eigen::VectorXd& rhs_constraint,
	const std::vector<std::vector<Eigen::Vector3d>>& boundary_points,
	const std::vector<uint>& patch_indices,
	const std::vector<ContactType>& contact_types,
	const std::vector<Eigen::Vector3d>& omega_bodyA,
	const std::vector<Eigen::Vector3d>& omega_bodyB,
	const std::vector<Eigen::Vector3d>& linear_contact_velocity,
	const std::vector<Eigen::Vector3d>& last_cop_point,
	ContactCOPSolution last_COP_sol
) {
	// TODO: implement, with careful consideration of how to transform the last COP solution to the
	// current active contact geometry
	// For now, we simply ignore the last COP sol
	return solveStartWithPatchCentroid(
		friction_coeff,
		A_constraint,
		rhs_constraint,
		boundary_points,
		patch_indices,
		contact_types,
		omega_bodyA,
		omega_bodyB,
		linear_contact_velocity
	);
}

void COPSolver::getCOPLineContactDisplacedMatricesExtended(
	Eigen::MatrixXd& A_disp,
	Eigen::VectorXd& rhs_disp,
	Eigen::Vector3d& lin_vel_disp,
	double signed_dist,
	const Eigen::MatrixXd& A,
	const Eigen::VectorXd& rhs,
	const Eigen::Vector3d& omegaA,
	const Eigen::Vector3d& omegaB,
	const Eigen::Vector3d& lin_vel
) {
	// if(COP_LOG_DEBUG) std::cout << "COP solver: getCOPLineContactDisplacedMatricesExtended: signed_distance_last_point: " << signed_distance_last_point << std::endl;
	MatrixXd cross_mat(5,5);
	cross_mat.setZero();
	cross_mat.block<3,3>(0,0) = Eigen::Matrix3d::Identity();
	cross_mat.block<3,2>(0,3) << 0,	0,
				0, signed_dist,
				-signed_dist, 0;
	cross_mat.block<2,2>(3,3) = Eigen::Matrix2d::Identity();
	A_disp = cross_mat*A*cross_mat.transpose();
	rhs_disp = cross_mat*rhs;
	// add the omega x (omega x r) component
	Matrix3d rcross;
	rcross << 0, 0, 0,
			  0, 0, -signed_dist,
			  0, signed_dist, 0;
	rhs_disp.segment<3>(0) += omegaB.cross(-rcross * omegaB) - omegaA.cross(-rcross * omegaA);
	lin_vel_disp = lin_vel + rcross*(omegaB - omegaA);
}

ContactCOPSolution COPSolver::solveStartWithPatchCentroid(
	double friction_coeff,
	const Eigen::MatrixXd& A_constraint,
	const Eigen::VectorXd& rhs_constraint, // 6 defined at point 0. Includes Jdot_qdot
	const std::vector<std::vector<Eigen::Vector3d>>& boundary_points, // points are in the COP frame with point0 being origin
	const std::vector<uint>& patch_indices, // ith entry gives starting row of ith contact patch in A_constraint, rhs_constraint
	const std::vector<ContactType>& contact_types, // ith entry gives contact type of ith contact patch
	const std::vector<Eigen::Vector3d>& omega_bodyA, // vector of body A angular velocities, in COP frame, one entry per contact patch
	const std::vector<Eigen::Vector3d>& omega_bodyB, // vector of body B angular velocities, in COP frame, one entry per contact patch
	const std::vector<Eigen::Vector3d>& linear_contact_velocity // 3 dof relative translation velocity at point 0, in COP frame, one entry per contact patch
	//^ expected to be zero in the z direction, but we don't explicitly check
) {
	// if(COP_LOG_DEBUG) {
	// 	std::cout << "COPSolver: resolveCOPLineContactWithLastCOPSol: inputs: " << std::endl;
	// 	std::cout << "A_full: " << A_full << std::endl;
	// 	std::cout << "rhs_nonlin_full: " << rhs_nonlin_full.transpose() << std::endl;
	// 	std::cout << "omega: " << omega.transpose() << std::endl;
	// 	std::cout << "r_last_cop: " << r_last_cop.transpose() << std::endl;
	// 	std::cout << "linear velocity: " << linear_contact_velocity.transpose() << std::endl;
	// }

	ContactCOPSolution ret_sol;
	if(boundary_points.size() != 1 && boundary_points[0].size() != 2) {
		ret_sol.result = COPSolResult::UnimplementedCase;
		return ret_sol;
	}

	assert(boundary_points[0][0].norm() < 1e-15); // should be (0,0,0)

	return solveOneLineStartWithPatchCentroidWithLCP(
		friction_coeff,
		A_constraint,
		rhs_constraint,
		boundary_points[0],
		omega_bodyA[0],
		omega_bodyB[0],
		linear_contact_velocity[0]);

	// ----------------------------------------------------------------------------------

	 // we assume that contact line segment is along x-axis, but do not know the direction
	Eigen::Vector3d point0_to_point1_dir = boundary_points[0][1];
	point0_to_point1_dir /= point0_to_point1_dir.norm(); // is either (1,0,0) or (-1,0,0)

	// compute the patch centroid
	Vector3d patch_centroid_point = boundary_points[0][1]/2.0;
	double signed_distance_from_point0 = patch_centroid_point.dot(Vector3d(1, 0, 0));
	// ^+ve in +x dir, -ve in -x dir

	// compute all matrices in the global COP frame, displaced to the patch centroid point
	double max_signed_dist = 0.0;
	double min_signed_dist = 0.0;
	if(signed_distance_from_point0 > 0) {
		max_signed_dist = 2*signed_distance_from_point0;
	} else {
		min_signed_dist = 2*signed_distance_from_point0;
	}

	// compute sign of rotational slip
	const double rotation_slip = omega_bodyB[0][2] - omega_bodyA[0][2];
	const double rotation_sign = (rotation_slip > 0)? 1: -1;

	// internal variables
	const int max_iters = 20;
	const double mu = friction_coeff;
	int iter_ind = 0;
	Eigen::MatrixXd A_disp(5, 5);
	Eigen::VectorXd rhs_disp(5);
	Vector3d lin_vel_disp;
	VectorXd constraint_vel_for_test = VectorXd::Zero(5);
	constraint_vel_for_test(4) = rotation_slip;
	Eigen::VectorXd f_sol(5), f_sol_full_roll(5);
	Eigen::Vector3d local_sol_point;
	Eigen::VectorXd a_sol(5);
	double mu_rot = 0.0;
	FeasibleCOPLineResult it_test_res;
	bool fForceLineCenter = false;
	double last_signed_distance = 0.0;
	double last_signed_distance_weight = 0.0;
	bool fForceIgnoreSlip = false;
	bool fBisectionSearch = false;
	double bisection_dist_bound_upper = 0.0;
	double bisection_dist_bound_lower = 0.0;
	double last_zero_moment_error = 0.0;

	// search for solution starting from centroid
	while (iter_ind < max_iters) {
		if(COP_LOG_DEBUG) std::cout << "COPSolver: iters: " << iter_ind << std::endl;
		iter_ind++;
		local_sol_point << signed_distance_from_point0, 0, 0;
		double dist1 = (local_sol_point).norm();
		double dist2 = (local_sol_point - boundary_points[0][1]).norm();
		// check whether we are on the boundary or inside the contact segment
		if(fForceLineCenter) {
			ret_sol.cop_type = COPContactType::LineCenter;
		} else {
			if(dist1 < 1e-15 || dist2 < 1e-15) {
				ret_sol.cop_type = COPContactType::LineEnd;
			} else {
				ret_sol.cop_type = COPContactType::LineCenter;
			}
		}
		// if(COP_LOG_DEBUG) std::cout << "COPSolver: contact type: " << ((ret_sol.cop_type == COPContactType::LineCenter)? "Center": "End") << std::endl;
		// if(COP_LOG_DEBUG) std::cout << "COPSolver: distance from last pt: " << signed_distance_from_point0 << std::endl;
		// compute contact matrices
		getCOPLineContactDisplacedMatricesExtended(
			A_disp, rhs_disp, lin_vel_disp,
			signed_distance_from_point0,
			A_constraint, rhs_constraint,
			omega_bodyA[0], omega_bodyB[0],
			linear_contact_velocity[0]
		);
		constraint_vel_for_test.segment<3>(0) = lin_vel_disp;
		// if(COP_LOG_DEBUG) {
		// 	std::cout << "COPSolver: displaced matrices: " << std::endl;
		// 	std::cout << "A: " << A_disp << std::endl;
		// 	std::cout << "rhs_nonlin: " << rhs_disp.transpose() << std::endl;
		// 	std::cout << "lin_vel_disp: " << lin_vel_disp.transpose() << std::endl;
		// }

		// compute mu_rotation
		mu_rot = getMuRotation(
			mu,
			(local_sol_point).norm(),
			(local_sol_point - boundary_points[0][1]).norm()
		);
		// if(COP_LOG_DEBUG) std::cout << "COPSolver: computed mu rot: " << mu_rot << std::endl;
		// ------------- test for rolling -----------------
		f_sol_full_roll = A_disp.ldlt().solve(-rhs_disp);
		// if(COP_LOG_DEBUG) std::cout << "COPSolver: computed full roll sol: " << f_sol_full_roll.transpose() << std::endl;
		a_sol.setZero();
		Eigen::MatrixXd tA;
		Eigen::VectorXd tsol, trhs;
		Eigen::Vector2d slip_vel = lin_vel_disp.segment<2>(0);
		// if(COP_LOG_DEBUG) std::cout << "COPSolver: translation slip: " << slip_vel.transpose() << " norm: " << slip_vel.norm() << std::endl;
		// if(COP_LOG_DEBUG) std::cout << "COPSolver: rotation slip: " << rotation_slip << std::endl;
		if(ret_sol.cop_type == COPContactType::LineCenter) {
			if(!fForceIgnoreSlip && slip_vel.norm() > 1e-6 && abs(rotation_slip) > 1e-8) {
				// slip velocity on both translation and rotation
				trhs.setZero(2);
				trhs[0] = -rhs_disp[2];
				trhs[1] = -rhs_disp[3];
				tA.setZero(2,2);
				tA(0,0) = A_disp(2,2)
							- mu*(A_disp.block<1,2>(2, 0).dot(slip_vel))/slip_vel.norm()
							- mu_rot*rotation_sign*A_disp(2,4);
				tA(0,1) = A_disp(2,3);
				tA(1,0) = A_disp(3,2)
							- mu*(A_disp.block<1,2>(3, 0).dot(slip_vel))/slip_vel.norm()
							- mu_rot*rotation_sign*A_disp(3,4);
				tA(1,1) = A_disp(3,3);
				tsol = tA.partialPivLu().solve(trhs);
				f_sol << - mu*tsol[0]*slip_vel/slip_vel.norm(),
							tsol,
							- mu_rot*rotation_sign*tsol[0];
				if(COP_LOG_DEBUG) std::cout << "COPSolver: trans & rot sliding force sol : " << f_sol.transpose() << std::endl;

			} else if (!fForceIgnoreSlip && slip_vel.norm() > 1e-6 && abs(rotation_slip) <= 1e-8) {
				// slip velocity on translation only. rotation velocity is zero
				trhs.setZero(3);
				trhs = -rhs_disp.segment<3>(2);
				tA.setZero(3,3);
				tA.block<3,1>(0,0) = A_disp.block<3,1>(2,2) - mu*A_disp.block<3,2>(2, 0)*slip_vel/slip_vel.norm();
				tA.block<3,2>(0,1) = A_disp.block<3,2>(2,3);
				tsol = tA.partialPivLu().solve(trhs);
				f_sol << - mu*tsol[0]*slip_vel/slip_vel.norm(),
							tsol;
				if(COP_LOG_DEBUG) std::cout << "COPSolver: trans slide & rot rolling force sol : " << f_sol.transpose() << std::endl;
			} else if ((fForceIgnoreSlip || slip_vel.norm() <= 1e-6) && abs(rotation_slip) > 1e-8) {
				// slip velocity on rotation only. translation slip velocity is zero
				fForceIgnoreSlip = true; // TODO: think more about this
				// TODO: Instead of this approach, we can consider the slip direction at the next time step.
				// i.e. slip_vel = current_slip_vel + tangential_acc*dt
				// Have to think about whether this can introduce an increase in energy, like in
				// Kane's example for frictional collision
				trhs.setZero(4);
				trhs = -rhs_disp.segment<4>(0);
				tA.setZero(4,4);
				tA.block<4,2>(0,0) = A_disp.block<4,2>(0,0);
				tA.block<4,1>(0,2) = A_disp.block<4,1>(0,2)
							- mu_rot*rotation_sign*A_disp.block<4,1>(0,4);
				tA.block<4,1>(0,3) = A_disp.block<4,1>(0,3);
				tsol = tA.partialPivLu().solve(trhs);
				f_sol << tsol.segment<4>(0), -mu_rot*rotation_sign*tsol[2];
				if(COP_LOG_DEBUG) std::cout << "COPSolver: trans rolling & rot slide force sol : " << f_sol.transpose() << std::endl;
			} else {
				// perfect rolling in both translation and rotation
				fForceIgnoreSlip = true; // TODO: think more about this
				f_sol = f_sol_full_roll;
			}
		} else { // basically a single point contact
			// TODO: call the one pt LCP solver?
			if(slip_vel.norm() > 1e-6) {
				// slip velocity on both translation and rotation
				double trhs_scalar = -rhs_disp[2];
				double tA_scalar = A_disp(2,2)
							- mu*(A_disp.block<1,2>(2, 0).dot(slip_vel))/slip_vel.norm();
				double tsol_scalar = trhs_scalar/tA_scalar;
				f_sol << - mu*tsol_scalar*slip_vel/slip_vel.norm(),
							tsol_scalar,
							0,
							0;
			} else {
				// slip velocity on rotation only. translation slip velocity is zero
				trhs.setZero(3);
				trhs = -rhs_disp.segment<3>(0);
				tA.setZero(3,3);
				tA = A_disp.block<3,3>(0,0);
				tsol = tA.ldlt().solve(trhs);
				f_sol << tsol, 0, 0;
			}
		}
		a_sol = A_disp*f_sol + rhs_disp;
		double dist_to_other_end = 0;
		if(ret_sol.cop_type == COPContactType::LineEnd) {
			if(max_signed_dist - signed_distance_from_point0 > 1e-15) {
				dist_to_other_end = max_signed_dist - signed_distance_from_point0;
			} else {
				dist_to_other_end = min_signed_dist - signed_distance_from_point0;
			}
		}
		it_test_res = testFeasibleCOPLine(f_sol, a_sol, constraint_vel_for_test, ret_sol.cop_type, mu, mu_rot, dist_to_other_end);
		if(COP_LOG_DEBUG) std::cout << "COPSolver: rolling solution result: " << static_cast<int>(it_test_res) << std::endl;
		// std::cout << "COPSolver: rolling solution result: " << static_cast<int>(it_test_res) << std::endl;
		if(it_test_res == FeasibleCOPLineResult::FRSuccess) {
			ret_sol.local_cop_pos = local_sol_point;
			ret_sol.force_sol = f_sol;
			ret_sol.result = COPSolResult::Success;
			return ret_sol;
		}
		if(fForceIgnoreSlip && it_test_res == FeasibleCOPLineResult::FRTranslationFrictionNonDissipative) {
			if(COP_LOG_DEBUG) std::cout << "Force slip, therefore ignore FRTranslationFrictionNonDissipative" << std::endl;
			ret_sol.local_cop_pos = local_sol_point;
			ret_sol.force_sol = f_sol;
			ret_sol.result = COPSolResult::Success;
			return ret_sol;
		}
		if(it_test_res == FeasibleCOPLineResult::FRNegativeNormalForce) {
			if(rhs_disp[2] > -1e-15) {
				ret_sol.force_sol.setZero(5);
				ret_sol.local_cop_pos = local_sol_point;
				ret_sol.result = COPSolResult::Success;
				return ret_sol;
			}
			// else: maybe there is a sliding solution? so continue checking
		}
		if(it_test_res == FeasibleCOPLineResult::FRZeroMomentViolation) {
			double moment_error = f_sol[3];
			if(fBisectionSearch && !fForceIgnoreSlip) {
				last_signed_distance = signed_distance_from_point0;
				if(moment_error > 0) {
					bisection_dist_bound_upper = signed_distance_from_point0;
				} else {
					bisection_dist_bound_lower = signed_distance_from_point0;
				}
				signed_distance_from_point0 = 0.5*bisection_dist_bound_lower + 0.5*bisection_dist_bound_upper;
				if(COP_LOG_DEBUG) std::cout << "Update signed distance last point: " << signed_distance_from_point0 << std::endl;
			}
			else if(moment_error*last_zero_moment_error < -1e-10 && !fForceIgnoreSlip) {
				// sign change in moment detected. Search using bisection
				// if force ignore slip is on, then there is a change in the moment curve
				// so we cannot search using bisection
				fBisectionSearch = true;
				// if(COP_LOG_DEBUG) std::cout << "Start bisection search for COP" << std::endl;
				// std::cout << "Start bisection search for COP" << std::endl;
				if(moment_error > 0) {
					bisection_dist_bound_lower = last_signed_distance;
					bisection_dist_bound_upper = signed_distance_from_point0;
				} else {
					bisection_dist_bound_lower = signed_distance_from_point0;
					bisection_dist_bound_upper = last_signed_distance;
				}
				last_signed_distance = signed_distance_from_point0;
				signed_distance_from_point0 = 0.5*bisection_dist_bound_lower + 0.5*bisection_dist_bound_upper;
				if(COP_LOG_DEBUG) std::cout << "Last signed distance: " << last_signed_distance << std::endl;
				if(COP_LOG_DEBUG) std::cout << "Set new signed distance last point: " << signed_distance_from_point0 << std::endl;
			} else {
				double zero_moment_dist = -f_sol[3]/fmax(f_sol[2], 1e-5);
				last_signed_distance = signed_distance_from_point0;
				signed_distance_from_point0 += zero_moment_dist;
				fBisectionSearch = false;
			}
			last_zero_moment_error = moment_error;
			// signed_distance_from_point0 = signed_distance_from_point0*(1.0 - last_signed_distance_weight) + last_signed_distance*last_signed_distance_weight;
			// if(last_signed_distance_weight < 1e-15) {
			// 	last_signed_distance_weight = 0.5;
			// }
			if(signed_distance_from_point0 > 0) {
				signed_distance_from_point0 = fmin(signed_distance_from_point0, max_signed_dist);
			} else {
				signed_distance_from_point0 = fmax(signed_distance_from_point0, min_signed_dist);
			}
			fForceLineCenter = false;
			continue;
		} else {
			fBisectionSearch = false;
		}
		if(it_test_res == FeasibleCOPLineResult::FROtherEndPenetration) {
			fForceLineCenter = true;
			continue;
		}
		if(it_test_res == FeasibleCOPLineResult::FRTranslationFrictionConeViolation) {
			// this case can only be possible when the translational slip velocity is zero
			// and the rolling constraint was violated.
			Eigen::Vector2d roll_slip_direction = -f_sol_full_roll.segment<2>(0);
			if(roll_slip_direction.norm() < 1e-8) {
				// no slip direction
				ret_sol.result = COPSolResult::NoSolution;
				return ret_sol;
			}
			if(ret_sol.cop_type == COPContactType::LineCenter) {
				// we need to test two conditions here, one with rotational sticking,
				// and another with rotational sliding
				if(abs(rotation_slip) > 1e-8) {
					// slip velocity on both translation and rotation
					trhs.setZero(2);
					trhs[0] = -rhs_disp[2];
					trhs[1] = -rhs_disp[3];
					tA.setZero(2,2);
					tA(0,0) = A_disp(2,2)
								- mu*(A_disp.block<1,2>(2, 0).dot(roll_slip_direction))
								- mu_rot*rotation_sign*A_disp(2,4);
					tA(0,1) = A_disp(2,3);
					tA(1,0) = A_disp(3,2)
								- mu*(A_disp.block<1,2>(3, 0).dot(roll_slip_direction))
								- mu_rot*rotation_sign*A_disp(3,4);
					tA(1,1) = A_disp(3,3);
					tsol = tA.partialPivLu().solve(trhs);
					f_sol << - mu*tsol[0]*roll_slip_direction,
								tsol,
								- mu_rot*rotation_sign*tsol[0];
				} else {
					// test rotational rolling. In case it is not feasible, the rotation friction cone constraint
					// will be violated. we will test this in the next step.
					trhs.setZero(3);
					trhs = -rhs_disp.segment<3>(2);
					tA.setZero(3,3);
					tA.block<3,1>(0,0) = A_disp.block<3,1>(2,2) - mu*A_disp.block<3,2>(2, 0)*roll_slip_direction;
					tA.block<3,2>(0,1) = A_disp.block<3,2>(2,3);
					tsol = tA.partialPivLu().solve(trhs);
					f_sol << - mu*tsol[0]*roll_slip_direction,
								tsol;
				}
			} else {
				// slip velocity on both translation and rotation
				double trhs_scalar = -rhs_disp[2];
				double tA_scalar = A_disp(2,2)
							- mu*(A_disp.block<1,2>(2, 0).dot(roll_slip_direction));
				double tsol_scalar = trhs_scalar/tA_scalar;
				f_sol << - mu*tsol_scalar*roll_slip_direction,
							tsol_scalar,
							0,
							0;
			}
			a_sol = A_disp*f_sol + rhs_disp;
			it_test_res = testFeasibleCOPLine(f_sol, a_sol, constraint_vel_for_test, ret_sol.cop_type, mu, mu_rot, dist_to_other_end);
			if(it_test_res == FeasibleCOPLineResult::FRSuccess) {
				ret_sol.local_cop_pos = local_sol_point;
				ret_sol.force_sol = f_sol;
				ret_sol.result = COPSolResult::Success;
				return ret_sol;
			} else {
				// if(COP_LOG_DEBUG) std::cout << "Translation sliding solution did not work " << static_cast<int>(it_test_res) << std::endl;
				std::cout << "Translation sliding solution did not work " << static_cast<int>(it_test_res) << std::endl;
			}
		}
		if(it_test_res == FeasibleCOPLineResult::FRRotationFrictionViolation) {
			// we only hit this if rotational slip is zero, but slip is unavoidable
			double roll_slip_direction = (f_sol_full_roll[4] > 0)? -1: 1;
			if(abs(f_sol_full_roll[4]) < 1e-8) {
				// no slip direction
				ret_sol.result = COPSolResult::NoSolution;
				return ret_sol;
			}
			if(ret_sol.cop_type == COPContactType::LineCenter) {
				if (slip_vel.norm() > 1e-8) {
					trhs.setZero(2);
					trhs[0] = -rhs_disp[2];
					trhs[1] = -rhs_disp[3];
					tA.setZero(2,2);
					tA(0,0) = A_disp(2,2)
								- mu*(A_disp.block<1,2>(2, 0).dot(slip_vel))/slip_vel.norm()
								- mu_rot*roll_slip_direction*A_disp(2,4);
					tA(0,1) = A_disp(2,3);
					tA(1,0) = A_disp(3,2)
								- mu*(A_disp.block<1,2>(3, 0).dot(slip_vel))/slip_vel.norm()
								- mu_rot*roll_slip_direction*A_disp(3,4);
					tA(1,1) = A_disp(3,3);
					tsol = tA.partialPivLu().solve(trhs);
					f_sol << - mu*tsol[0]*slip_vel/slip_vel.norm(),
								tsol,
								- mu_rot*roll_slip_direction*tsol[0];
				} else {
					trhs.setZero(4);
					trhs = -rhs_disp.segment<4>(0);
					tA.setZero(4,4);
					tA.block<4,2>(0,0) = A_disp.block<4,2>(0,0);
					tA.block<4,1>(0,2) = A_disp.block<4,1>(0,2)
								- mu_rot*roll_slip_direction*A_disp.block<4,1>(0,4);
					tA.block<4,1>(0,3) = A_disp.block<4,1>(0,3);
					tsol = tA.partialPivLu().solve(trhs);
					f_sol << tsol.segment<4>(0), -mu_rot*roll_slip_direction*tsol[2];
				}
			}
			// else, not possible since LineEnd case does not allow any moments at all
			a_sol = A_disp*f_sol + rhs_disp;
			it_test_res = testFeasibleCOPLine(f_sol, a_sol, constraint_vel_for_test, ret_sol.cop_type, mu, mu_rot, dist_to_other_end);
			if(it_test_res == FeasibleCOPLineResult::FRSuccess) {
				ret_sol.local_cop_pos = local_sol_point;
				ret_sol.force_sol = f_sol;
				ret_sol.result = COPSolResult::Success;
				return ret_sol;
			} else {
				// if(COP_LOG_DEBUG) std::cout << "Rotation sliding solution did not work " << static_cast<int>(it_test_res) << std::endl;
				std::cout << "Rotation sliding solution did not work " << static_cast<int>(it_test_res) << std::endl;
			}
		}
		// lastly check if we can find a sliding solution for both translation and rotation
		if(it_test_res == FeasibleCOPLineResult::FRTranslationFrictionConeViolation) {
			if(ret_sol.cop_type == COPContactType::LineCenter) {
				double rotation_slip_direction = (f_sol_full_roll[4] > 0)? -1: 1;
				Eigen::Vector2d roll_slip_direction = -f_sol_full_roll.segment<2>(0);
				if(abs(f_sol_full_roll[4]) >= 1e-8 && roll_slip_direction.norm() >= 1e-8) {
					trhs.setZero(2);
					trhs[0] = -rhs_disp[2];
					trhs[1] = -rhs_disp[3];
					tA.setZero(2,2);
					tA(0,0) = A_disp(2,2)
								- mu*(A_disp.block<1,2>(2, 0).dot(roll_slip_direction))
								- mu_rot*rotation_slip_direction*A_disp(2,4);
					tA(0,1) = A_disp(2,3);
					tA(1,0) = A_disp(3,2)
								- mu*(A_disp.block<1,2>(3, 0).dot(roll_slip_direction))
								- mu_rot*rotation_slip_direction*A_disp(3,4);
					tA(1,1) = A_disp(3,3);
					tsol = tA.partialPivLu().solve(trhs);
					f_sol << - mu*tsol[0]*roll_slip_direction,
								tsol,
								- mu_rot*rotation_slip_direction*tsol[0];
					a_sol = A_disp*f_sol + rhs_disp;
					it_test_res = testFeasibleCOPLine(f_sol, a_sol, constraint_vel_for_test, ret_sol.cop_type, mu, mu_rot, dist_to_other_end);
					if(it_test_res == FeasibleCOPLineResult::FRSuccess) {
						ret_sol.local_cop_pos = local_sol_point;
						ret_sol.force_sol = f_sol;
						ret_sol.result = COPSolResult::Success;
						return ret_sol;
					}  else {
						// if(COP_LOG_DEBUG) std::cout << "Rotation+translation sliding solution did not work " << static_cast<int>(it_test_res) << std::endl;
						std::cout << "Rotation+translation sliding solution did not work " << static_cast<int>(it_test_res) << std::endl;
					}
				}
			} // else condition is handled above already
		}
		// TODO: should we loop again?
		if(it_test_res != FeasibleCOPLineResult::FRSuccess) {
			break;
		}
		// TODO: handle Painleve's problem?
	}
	ret_sol.result = COPSolResult::NoSolution;
	if(COP_LOG_DEBUG) std::cout << "No cop sol!" << std::endl;
	return ret_sol;
}

ContactCOPSolution COPSolver::solveOneLineStartWithPatchCentroidWithLCP(
		double friction_coeff,
		const Eigen::MatrixXd& A_constraint, // defined at point0
		const Eigen::VectorXd& rhs_constraint, // defined at point 0. Includes Jdot_qdot
		const std::vector<Eigen::Vector3d>& boundary_points, // points are in the COP frame with point0 being origin
		const Eigen::Vector3d& omega_bodyA, // vector of body A angular velocities, in COP frame, one entry per contact patch
		const Eigen::Vector3d& omega_bodyB, // vector of body B angular velocities, in COP frame, one entry per contact patch
		const Eigen::Vector3d& linear_contact_velocity // 3 dof relative translation velocity at point 0, in COP frame, one entry per contact patch
		//^ expected to be zero in the z direction, but we don't explicitly check
) {
	ContactCOPSolution ret_sol;

	assert(boundary_points[0].norm() < 1e-15); // should be (0,0,0)

	 // we assume that contact line segment is along x-axis, but do not know the direction
	Eigen::Vector3d point0_to_point1_dir = boundary_points[1];
	point0_to_point1_dir /= point0_to_point1_dir.norm(); // is either (1,0,0) or (-1,0,0)

	// compute the patch centroid
	Vector3d patch_centroid_point = boundary_points[1]/2.0;
	double signed_distance_from_point0 = patch_centroid_point.dot(Vector3d(1, 0, 0));
	// ^+ve in +x dir, -ve in -x dir

	// compute all matrices in the global COP frame, displaced to the patch centroid point
	double max_signed_dist = 0.0;
	double min_signed_dist = 0.0;
	if(signed_distance_from_point0 > 0) {
		max_signed_dist = 2*signed_distance_from_point0;
	} else {
		min_signed_dist = 2*signed_distance_from_point0;
	}

	// // compute sign of rotational slip
	const double rotation_slip = omega_bodyB[2] - omega_bodyA[2];
	// const double rotation_sign = (rotation_slip > 0)? 1: -1;

	// internal variables
	const int max_iters = 20;
	const double mu = friction_coeff;
	int iter_ind = 0;
	Eigen::MatrixXd A_disp(5, 5);
	Eigen::VectorXd rhs_disp(5);
	Vector3d lin_vel_disp;
	VectorXd constraint_vel_for_test = VectorXd::Zero(5);
	constraint_vel_for_test(4) = rotation_slip;
	Eigen::VectorXd f_sol(5);
	Eigen::VectorXd a_sol(5);
	Eigen::Vector3d local_sol_point;
	double mu_rot = 0.0;
	FeasibleCOPLineResult it_test_res;
	bool fForceLineCenter = false;
	double last_signed_distance = 0.0;
	bool fForceIgnoreSlip = false;
	bool fBisectionSearch = false;
	double bisection_dist_bound_upper = 0.0;
	double bisection_dist_bound_lower = 0.0;
	double last_zero_moment_error = 0.0;
	auto lcp_solver = Sai2LCPSolver::LCPSolver();

	// search for solution starting from centroid
	while (iter_ind < max_iters) {
		if(COP_LOG_DEBUG) std::cout << "COPSolver: iters: " << iter_ind << std::endl;
		iter_ind++;
		local_sol_point << signed_distance_from_point0, 0, 0;
		double dist1 = (local_sol_point).norm();
		double dist2 = (local_sol_point - boundary_points[1]).norm();
		// check whether we are on the boundary or inside the contact segment
		if(fForceLineCenter) {
			ret_sol.cop_type = COPContactType::LineCenter;
		} else {
			if(dist1 < 1e-15 || dist2 < 1e-15) {
				ret_sol.cop_type = COPContactType::LineEnd;
			} else {
				ret_sol.cop_type = COPContactType::LineCenter;
			}
		}
		// if(COP_LOG_DEBUG) std::cout << "COPSolver: contact type: " << ((ret_sol.cop_type == COPContactType::LineCenter)? "Center": "End") << std::endl;
		// if(COP_LOG_DEBUG) std::cout << "COPSolver: distance from last pt: " << signed_distance_from_point0 << std::endl;
		// compute contact matrices
		getCOPLineContactDisplacedMatricesExtended(
			A_disp, rhs_disp, lin_vel_disp,
			signed_distance_from_point0,
			A_constraint, rhs_constraint,
			omega_bodyA, omega_bodyB,
			linear_contact_velocity
		);
		constraint_vel_for_test.segment<3>(0) = lin_vel_disp;
		if(COP_LOG_DEBUG) {
			std::cout << "COPSolver: displaced matrices: " << std::endl;
			std::cout << "A: " << A_disp << std::endl;
			std::cout << "rhs_nonlin: " << rhs_disp.transpose() << std::endl;
			std::cout << "lin_vel_disp: " << lin_vel_disp.transpose() << std::endl;
		}

		// Compute viscous rotational friction
		mu_rot = getMuRotation(
			mu,
			(local_sol_point).norm(),
			(local_sol_point - boundary_points[1]).norm()
		);
		const double visc_rot_torque = -mu_rot*rotation_slip;
		rhs_disp += A_disp.col(4)*visc_rot_torque;

		if(COP_LOG_DEBUG) {
			std::cout << "COPSolver: computed mu rot: " << mu_rot << std::endl;
			std::cout << "COPSolver: computed viscous rot torque: " << visc_rot_torque << std::endl;
		}

		Eigen::Vector2d slip_vel = lin_vel_disp.segment<2>(0);
		if (slip_vel.norm() <= 1e-6 && abs(rotation_slip) > 1e-8) {
			// Ignore slip in all future iterations since there is a pt on the line
			// where friction is purely rotational. So we apply just rotational friction
			// everywhere.
			fForceIgnoreSlip = true;
		}

		// Call LCP solver
		auto lcp_moment_constraint = Sai2LCPSolver::MomentConstraints::NoMoment;
		if (ret_sol.cop_type == COPContactType::LineCenter) {
			// Technically moment is about Y axis, but the way we constructed
			// the constraint matrix puts the Y moment at index 3 of the solution
			// force vector. From the LCP's pov, this is the X moment.
			lcp_moment_constraint = Sai2LCPSolver::MomentConstraints::XMoment;
		}

		// LCP solver expects pre_v to be the same size as vector b. So we add 0's to pad.
		VectorXd lcp_lin_vel_disp(5);
		lcp_lin_vel_disp << lin_vel_disp, 0, 0;
		if (fForceIgnoreSlip) {
			lcp_lin_vel_disp[0] = 0;
			lcp_lin_vel_disp[1] = 0;
		}
		auto lcp_result = lcp_solver.solveWithMoments(
			A_disp,
			rhs_disp,
			lcp_lin_vel_disp,
			{lcp_moment_constraint},
			mu);
		if (lcp_result.result != Sai2LCPSolver::LCPSolResult::Success) {
			break;
		}
		f_sol = lcp_result.p_sol;
		// std::cout << f_sol.transpose() << std::endl;

		a_sol = A_disp*f_sol + rhs_disp;
		double dist_to_other_end = 0;
		if(ret_sol.cop_type == COPContactType::LineEnd) {
			if(max_signed_dist - signed_distance_from_point0 > 1e-15) {
				dist_to_other_end = max_signed_dist - signed_distance_from_point0;
			} else {
				dist_to_other_end = min_signed_dist - signed_distance_from_point0;
			}
		}
		it_test_res = testFeasibleCOPLine(f_sol, a_sol, constraint_vel_for_test, ret_sol.cop_type, mu, mu_rot, dist_to_other_end);
		if(COP_LOG_DEBUG) std::cout << "COPSolver: rolling solution result: " << static_cast<int>(it_test_res) << std::endl;
		if(it_test_res == FeasibleCOPLineResult::FRSuccess) {
			ret_sol.local_cop_pos = local_sol_point;
			ret_sol.force_sol = f_sol;
			ret_sol.result = COPSolResult::Success;
			return ret_sol;
		}
		if(it_test_res == FeasibleCOPLineResult::FRZeroMomentViolation) {
			double moment_error = f_sol[3];
			if(fBisectionSearch && !fForceIgnoreSlip) {
				last_signed_distance = signed_distance_from_point0;
				if(moment_error > 0) {
					bisection_dist_bound_upper = signed_distance_from_point0;
				} else {
					bisection_dist_bound_lower = signed_distance_from_point0;
				}
				signed_distance_from_point0 = 0.5*bisection_dist_bound_lower + 0.5*bisection_dist_bound_upper;
				if(COP_LOG_DEBUG) std::cout << "Update signed distance last point: " << signed_distance_from_point0 << std::endl;
			}
			else if(moment_error*last_zero_moment_error < -1e-10 && !fForceIgnoreSlip) {
				// sign change in moment detected. Search using bisection
				// if force ignore slip is on, then there is a change in the moment curve
				// so we cannot search using bisection
				fBisectionSearch = true;
				// if(COP_LOG_DEBUG) std::cout << "Start bisection search for COP" << std::endl;
				// std::cout << "Start bisection search for COP" << std::endl;
				if(moment_error > 0) {
					bisection_dist_bound_lower = last_signed_distance;
					bisection_dist_bound_upper = signed_distance_from_point0;
				} else {
					bisection_dist_bound_lower = signed_distance_from_point0;
					bisection_dist_bound_upper = last_signed_distance;
				}
				last_signed_distance = signed_distance_from_point0;
				signed_distance_from_point0 = 0.5*bisection_dist_bound_lower + 0.5*bisection_dist_bound_upper;
				if(COP_LOG_DEBUG) std::cout << "Last signed distance: " << last_signed_distance << std::endl;
				if(COP_LOG_DEBUG) std::cout << "Set new signed distance last point: " << signed_distance_from_point0 << std::endl;
			} else {
				double zero_moment_dist = -f_sol[3]/fmax(f_sol[2], 1e-5);
				last_signed_distance = signed_distance_from_point0;
				signed_distance_from_point0 += zero_moment_dist;
				fBisectionSearch = false;
			}
			last_zero_moment_error = moment_error;
			if(signed_distance_from_point0 > 0) {
				signed_distance_from_point0 = fmin(signed_distance_from_point0, max_signed_dist);
			} else {
				signed_distance_from_point0 = fmax(signed_distance_from_point0, min_signed_dist);
			}
			fForceLineCenter = false;
			continue;
		} else {
			fBisectionSearch = false;
		}
		if(it_test_res == FeasibleCOPLineResult::FROtherEndPenetration) {
			fForceLineCenter = true;
			continue;
		}
		// TODO: should we loop again?
		if(it_test_res != FeasibleCOPLineResult::FRSuccess) {
			break;
		}
	}
	ret_sol.result = COPSolResult::NoSolution;
	if(COP_LOG_DEBUG) std::cout << "No cop sol!" << std::endl;
	return ret_sol;
}

ContactCOPSolution COPSolver::solvePtOnly(
		double friction_coeff,
		const Eigen::MatrixXd& A_constraint, // 3x3 constraint inertia matrix
		const Eigen::Vector3d& rhs_constraint,
		const Eigen::Vector3d& linear_contact_velocity
) {
	auto lcp_solver = Sai2LCPSolver::LCPSolver();

	// TODO: optimize the LCP solver for 1 pt contact with fixed size matrices.
	// This generic call is almost 30% slower than a specialized solver for
	// 1 pt.
	auto lcp_result = lcp_solver.solve(
		A_constraint,
		rhs_constraint,
		linear_contact_velocity,
		/*epsilon=*/ 0,
		friction_coeff,
		/*force_sliding_if_pre_slip=*/ true);

	ContactCOPSolution ret_sol;
	if (lcp_result.result == Sai2LCPSolver::LCPSolResult::Success) {
		ret_sol.force_sol = lcp_result.p_sol;
		ret_sol.local_cop_pos.setZero();
		ret_sol.result = COPSolResult::Success;
		return ret_sol;
	} else {
		ret_sol.result = COPSolResult::NoSolution;
		return ret_sol;
	}
}


static void getCOPLineAndPtContactDisplacedMatricesExtended (
	Eigen::MatrixXd& A_disp,
	Eigen::VectorXd& rhs_disp,
	Eigen::Vector3d& lin_vel_disp,
	double signed_dist,
	const Eigen::MatrixXd& A,
	const Eigen::VectorXd& rhs,
	const Eigen::Vector3d& line_omega_bodyA,
	const Eigen::Vector3d& line_omega_bodyB,
	const Eigen::Vector3d& line_lin_vel,
	uint line_start_row_id,
	uint pt_start_row_id
) {
	// if(COP_LOG_DEBUG) std::cout << "COP solver: getCOPLineContactDisplacedMatricesExtended: signed_distance_last_point: " << signed_distance_last_point << std::endl;
	MatrixXd cross_mat(5,5);
	cross_mat.setZero();
	cross_mat.block<3,3>(0,0) = Eigen::Matrix3d::Identity();
	cross_mat.block<3,2>(0,3) << 0,	0,
				0, signed_dist,
				-signed_dist, 0;
	cross_mat.block<2,2>(3,3) = Eigen::Matrix2d::Identity();

	// compute A matrix at new position
	A_disp.setZero(A.rows(), A.cols());
	A_disp.block<5,5>(line_start_row_id, line_start_row_id) =
		cross_mat*A.block<5,5>(line_start_row_id, line_start_row_id)*cross_mat.transpose();
	A_disp.block<3,3>(pt_start_row_id, pt_start_row_id) =
		A.block<3,3>(pt_start_row_id, pt_start_row_id);
	A_disp.block<5,3>(line_start_row_id, pt_start_row_id) =
		cross_mat*A.block<5,3>(line_start_row_id, pt_start_row_id);
	A_disp.block<3,5>(pt_start_row_id, line_start_row_id) =
		A.block<3,5>(pt_start_row_id, line_start_row_id)*cross_mat.transpose();

	// compute RHS at new position
	rhs_disp.setZero(rhs.size());
	rhs_disp.segment<5>(line_start_row_id) = cross_mat*rhs.segment<5>(line_start_row_id);
	rhs_disp.segment<3>(pt_start_row_id) = rhs.segment<3>(pt_start_row_id);

	// add the omega x (omega x r) component
	Matrix3d rcross;
	rcross << 0, 0, 0,
			  0, 0, -signed_dist,
			  0, signed_dist, 0;
	rhs_disp.segment<3>(line_start_row_id) +=
		line_omega_bodyB.cross(-rcross * line_omega_bodyB)
		- line_omega_bodyA.cross(-rcross * line_omega_bodyA);
	lin_vel_disp = line_lin_vel + rcross*(line_omega_bodyB - line_omega_bodyA);
}

// Solve for COPs with one point contact and one line contact:
ContactCOPSolution COPSolver::solveStartWithPatchCentroidOnePointOneLine(
	double friction_coeff,
	const Eigen::MatrixXd& A_constraint,
	const Eigen::VectorXd& rhs_constraint, // 5+3/ 3+5 defined at point 0. Includes Jdot_qdot
	const std::vector<std::vector<Eigen::Vector3d>>& boundary_points, // points are in the COP frame with point0 being origin
	const std::vector<uint>& patch_indices, // ith entry gives starting row of ith contact patch in A_constraint, rhs_constraint
	const std::vector<ContactType>& contact_types, // ith entry gives contact type of ith contact patch
	const std::vector<Eigen::Vector3d>& omega_bodyA, // vector of body A angular velocities, in COP frame, one entry per contact patch
	const std::vector<Eigen::Vector3d>& omega_bodyB, // vector of body B angular velocities, in COP frame, one entry per contact patch
	const std::vector<Eigen::Vector3d>& linear_contact_velocity // 3 dof relative translation velocity at point 0, in COP frame, one entry per contact patch
	//^ expected to be zero in the z direction, but we don't explicitly check
) {
	assert(contact_types.size() == 2);

	pt_start_row_id = 0;
	line_start_row_id = 0;
	uint vector_line_ind = 0;
	uint vector_pt_ind = 0;

	// std::cout << "New solver call " << rhs_constraint.transpose() << std::endl;

	if(contact_types[0] == ContactType::POINT) {
		assert(contact_types[1] == ContactType::LINE);
		pt_start_row_id = patch_indices[0];
		line_start_row_id = patch_indices[1];
		vector_line_ind = 1;
		vector_pt_ind = 0;
	} else if(contact_types[0] == ContactType::LINE) {
		assert(contact_types[1] == ContactType::POINT);
		pt_start_row_id = patch_indices[1];
		line_start_row_id = patch_indices[0];
		vector_line_ind = 0;
		vector_pt_ind = 1;
	} else {
		throw(std::runtime_error("Unsupported types"));
	}

	// check if all normal accelerations are > 0. if so, we are done
	double pt1_end_dist = boundary_points[vector_line_ind][1](0);
	if(rhs_constraint(pt_start_row_id+2) > -1e-8 &&
		rhs_constraint(line_start_row_id+2) > -1e-8 &&
		(rhs_constraint(line_start_row_id+2) - rhs_constraint(line_start_row_id+3)*pt1_end_dist) > -1e-8
	) {
		if(COP_LOG_DEBUG) {
			std::cout << "OnePointOneLine: Both contacts separating.\n";
		}
		ContactCOPSolution ret_sol;
		ret_sol.result = COPSolResult::Success;
		ret_sol.force_sol = VectorXd::Zero(8);
		ret_sol.local_cop_pos.setZero();
		return ret_sol;
	}

	// if b[pt contact normal] > 0, test with pt contact inactive.
	// - simply call the 1 pt line contact solver with reduced matrices
	// - check if pt contact is accelerating as a result. if not, we are done
	if(rhs_constraint(pt_start_row_id+2) > -1e-8) {
		std::vector<std::vector<Eigen::Vector3d>> lb_pts;
		lb_pts.push_back(boundary_points[vector_line_ind]);
		std::vector<uint> ptch_inds;
		ptch_inds.push_back(patch_indices[vector_line_ind]);
		std::vector<ContactType> ctype_inds;
		ctype_inds.push_back(contact_types[vector_line_ind]);
		std::vector<Eigen::Vector3d> omAs, omBs, lvs;
		omAs.push_back(omega_bodyA[vector_line_ind]);
		omBs.push_back(omega_bodyB[vector_line_ind]);
		lvs.push_back(linear_contact_velocity[vector_line_ind]);
		ContactCOPSolution line_only_sol = solveStartWithPatchCentroid(
												friction_coeff,
												A_constraint.block<5,5>(line_start_row_id, line_start_row_id),
												rhs_constraint.segment<5>(line_start_row_id),
												lb_pts,
												ptch_inds,
												ctype_inds,
												omAs,
												omBs,
												lvs
											);
		if(line_only_sol.result == COPSolResult::Success) {
			// check if line sol causes penetration at pt
			VectorXd fsol_pt0 = VectorXd::Zero(8);
			fsol_pt0.segment<5>(line_start_row_id) = line_only_sol.force_sol;
			fsol_pt0(line_start_row_id+3) -= line_only_sol.local_cop_pos(0)*line_only_sol.force_sol(2);
			fsol_pt0(line_start_row_id+4) += line_only_sol.local_cop_pos(0)*line_only_sol.force_sol(1);
			double pta = A_constraint.row(pt_start_row_id+2).dot(fsol_pt0) + rhs_constraint(pt_start_row_id+2);
			if(pta >= -1e-8) { // pt contact can be inactive
				if(COP_LOG_DEBUG) {
					std::cout << "OnePointOneLine: Line sol only succeeded.\n";
				}
				ContactCOPSolution ret_sol;
				ret_sol.force_sol.setZero(8);
				ret_sol.force_sol.segment<5>(line_start_row_id) = line_only_sol.force_sol;
				ret_sol.local_cop_pos = line_only_sol.local_cop_pos;
				ret_sol.result = COPSolResult::Success;
				ret_sol.cop_type = line_only_sol.cop_type;
				return ret_sol;
			}
		} else {
			// how can this fail? TODO: maybe the line+pt sol can still work?
			if(COP_LOG_DEBUG) {
				std::cout << "OnePointOneLine: Line sol only failed but still "
					      << "using it.\n";
			}
			return line_only_sol;
		}
	}

	// if b[line contact normal] > 0, test with line contact inactive.
	// - check if penetration acceleration occurs at any of the line contact test points
	// - if not, we are done
	if(rhs_constraint(line_start_row_id+2) > -1e-8 &&
			(rhs_constraint(line_start_row_id+2) - rhs_constraint(line_start_row_id+3)*pt1_end_dist) > -1e-8
			// TODO: this does not take into account the w x (w x r) term shift to the other point
	) {
		ContactCOPSolution pt_only_sol = solvePtOnly(
												friction_coeff,
												A_constraint.block<3,3>(pt_start_row_id, pt_start_row_id),
												rhs_constraint.segment<3>(pt_start_row_id),
												linear_contact_velocity[vector_pt_ind]
											);
		if(pt_only_sol.result == COPSolResult::Success) {
			Vector2d line_a = A_constraint.block<2,3>(line_start_row_id+2,pt_start_row_id)*(pt_only_sol.force_sol)
								+ rhs_constraint.segment<2>(line_start_row_id+2);
			if(line_a(0) > -1e-8 &&
				(line_a(0) - line_a(1)*pt1_end_dist) > -1e-8
			) {
				if(COP_LOG_DEBUG) {
					std::cout << "OnePointOneLine: Point sol only succeeded.\n";
				}
				ContactCOPSolution ret_sol;
				ret_sol.force_sol.setZero(8);
				ret_sol.force_sol.segment<3>(pt_start_row_id) = pt_only_sol.force_sol;
				ret_sol.local_cop_pos = Vector3d::Zero();
				ret_sol.result = COPSolResult::Success;
				return ret_sol;
			}
		} else {
			// how can this fail? TODO: maybe the line+pt sol can still work?
			if(COP_LOG_DEBUG) {
				std::cout << "OnePointOneLine: Point sol only failed but still "
					      << "using it.\n";
			}
			return pt_only_sol;
		}
	}

	// select starting point for line contact
	// set initial line-COP state to LineCenter
	 // we assume that contact line segment is along x-axis, but do not know the direction
	Eigen::Vector3d point0_to_point1_dir = boundary_points[vector_line_ind][1];
	point0_to_point1_dir /= point0_to_point1_dir.norm(); // is either (1,0,0) or (-1,0,0)

	// compute the patch centroid
	Vector3d patch_centroid_point = boundary_points[vector_line_ind][1]/2.0;
	double signed_distance_from_point0 = patch_centroid_point.dot(Vector3d(1, 0, 0));
	// ^+ve in +x dir, -ve in -x dir

	// compute all matrices in the global COP frame, displaced to the patch centroid point
	double max_signed_dist = 0.0;
	double min_signed_dist = 0.0;
	if(signed_distance_from_point0 > 0) {
		max_signed_dist = 2*signed_distance_from_point0;
	} else {
		min_signed_dist = 2*signed_distance_from_point0;
	}

	// compute line rotational slip
	line_rotation_slip_dir = omega_bodyB[vector_line_ind][2] - omega_bodyA[vector_line_ind][2];
	if(abs(line_rotation_slip_dir) > 1e-8) {
		line_rotational_state = FrictionState::Sliding;
		line_rotation_slip_dir /= abs(line_rotation_slip_dir);
	} else {
		// initialize to rolling
		line_rotational_state = FrictionState::Rolling;
	}

	// compute point translation slip
	// save rolling slip and translation slip behavior in point contact.
	point_slip_dir = linear_contact_velocity[vector_pt_ind].segment<2>(0);
	if (point_slip_dir.norm() > 1e-6) {
		point_state = FrictionState::Sliding;
		point_slip_dir /= point_slip_dir.norm();
	} else {
		// initialize to rolling
		point_state = FrictionState::Rolling;
	}

	const int max_iters = 40;
	int iter_ind = 0;
	Eigen::MatrixXd A_disp(8, 8);
	Eigen::VectorXd rhs_disp(8);
	Vector3d lin_vel_disp;
	bool fForceLineCenter = false;
	bool fForceIgnoreSlip = false;
	Eigen::VectorXd full_f_sol(8), full_a_sol(8);
	ContactCOPSolution ret_sol;
	double mu_rot;
	const double mu = friction_coeff;
	Vector3d local_line_cop_sol_point;

	// initial internal matrices
	TA = A_constraint;
	Trhs = rhs_constraint;
	TA_size = 0;

	bool did_update_cop_pos = true;
	while (iter_ind < max_iters) {
		iter_ind++;

		// if line-cop position was updated, update displaced terms
		if(did_update_cop_pos) {
			did_update_cop_pos = false;

			local_line_cop_sol_point << signed_distance_from_point0, 0, 0;
			double dist1 = (local_line_cop_sol_point).norm();
			double dist2 = (local_line_cop_sol_point - boundary_points[vector_line_ind][1]).norm();
			// check whether we are on the boundary or inside the contact segment
			if(fForceLineCenter) {
				line_cop_type = COPContactType::LineCenter;
			} else {
				if(dist1 < 1e-15 || dist2 < 1e-15) {
					line_cop_type = COPContactType::LineEnd;
				} else {
					line_cop_type = COPContactType::LineCenter;
				}
			}

			// get displaced matrices
			getCOPLineAndPtContactDisplacedMatricesExtended(
				A_disp, rhs_disp, lin_vel_disp,
				signed_distance_from_point0,
				A_constraint, rhs_constraint,
				omega_bodyA[vector_line_ind], omega_bodyB[vector_line_ind],
				linear_contact_velocity[vector_line_ind],
				line_start_row_id, pt_start_row_id
			);

			// compute mu_rotation
			mu_rot = getMuRotation(mu, dist1, dist2);

			// get line translational slip velocities at current COP
			if(!fForceIgnoreSlip) {
				line_translation_slip_dir = lin_vel_disp.segment<2>(0);
				if (line_translation_slip_dir.norm() > 1e-6) {
					line_translation_state = FrictionState::Sliding;
					line_translation_slip_dir /= line_translation_slip_dir.norm();
				} else {
					// initialize to rolling
					line_translation_state = FrictionState::Rolling;
					// if rolling at line-COP, set fForceIgnoreSlip true
					fForceIgnoreSlip = true;
				}
			}
		}

		// solve with current state
		computeMatrices(A_disp, rhs_disp, mu, mu_rot);
		solve(mu, mu_rot, full_f_sol);
		full_a_sol = A_disp * full_f_sol + rhs_disp;

		// check for normal force violation
		if(full_f_sol(line_start_row_id+2) < -1e-8 || full_f_sol(pt_start_row_id+2) < 1e-8) {
			// TODO: think about disabling a contact
			ret_sol.result = COPSolResult::UnimplementedCase;
			return ret_sol;
		}

		// point: check for friction cone violation
		if(point_state == FrictionState::Rolling) {
			Vector2d rolling_force = full_f_sol.segment<2>(pt_start_row_id);
			if(rolling_force.norm() > mu*full_f_sol(pt_start_row_id+2) + 1e-15) {
				// save impending direction
				// TODO: proper impending slip direction
				point_rolling_dir = -rolling_force/rolling_force.norm();
				// switch point state to impending dir
				point_state = FrictionState::Impending;
				continue;
			}
		}
		// point: check for impending slip wrong rolling direction
		if(point_state == FrictionState::Impending) {
			Vector2d friction_force = full_f_sol.segment<2>(pt_start_row_id);
			Vector2d slip_accel = full_a_sol.segment<2>(pt_start_row_id);
			if(slip_accel.dot(friction_force) > 1e-8) {
				// switch back to rolling so that we recompute the rolling force
				point_state = FrictionState::Rolling;
				continue;
			}
		}
		// line: check for translation friction cone violation
		if(line_translation_state == FrictionState::Rolling) {
			Vector2d rolling_force = full_f_sol.segment<2>(line_start_row_id);
			if(rolling_force.norm() > mu*full_f_sol(line_start_row_id+2) + 1e-15) {
				// save impending direction
				// TODO: proper impending slip direction
				line_translation_rolling_dir = -rolling_force/rolling_force.norm();
				// switch point state to impending dir
				line_translation_state = FrictionState::Impending;
				continue;
			}
		}
		// line: check for impending translation slip wrong rolling direction
		if(line_translation_state == FrictionState::Impending) {
			Vector2d friction_force = full_f_sol.segment<2>(line_start_row_id);
			Vector2d slip_accel = full_a_sol.segment<2>(line_start_row_id);
			if(slip_accel.dot(friction_force) > 1e-8) {
				// switch back to rolling so that we recompute the rolling force
				line_translation_state = FrictionState::Rolling;
				continue;
			}
		}
		// line: check for rotation friction limit violation
		if(line_rotational_state == FrictionState::Rolling) {
			double rolling_force = full_f_sol(line_start_row_id+4);
			if(abs(rolling_force) > mu_rot*full_f_sol(line_start_row_id+2) + 1e-15) {
				// save impending direction
				line_rotation_rolling_dir = -rolling_force/abs(rolling_force);
				// switch point state to impending dir
				line_rotational_state = FrictionState::Impending;
				continue;
			}
		}
		// line: check for impending rotational slip wrong rolling direction
		if(line_rotational_state == FrictionState::Impending) {
			double rolling_force = full_f_sol(line_start_row_id+4);
			double rolling_acc = full_a_sol(line_start_row_id+4);
			if(rolling_acc*rolling_force > 1e-8) {
				// switch back to rolling so that we recompute the rolling force
				line_rotational_state = FrictionState::Rolling;
				continue;
			}
		}

		// - check for other end penetration. if so, force line center solution. print to cerr.
		// - we dont expect this case since we start from the patch center point currently.
		if(line_cop_type == COPContactType::LineEnd) {
			// TODO: handle w x (w x r) term correctly
			double acc_end1 = full_a_sol[line_start_row_id+2] +
								full_a_sol[line_start_row_id+3]*signed_distance_from_point0;
			double acc_end2 = full_a_sol[line_start_row_id+2] -
								full_a_sol[line_start_row_id+3]*(boundary_points[vector_line_ind][1](0) - signed_distance_from_point0);
			if(acc_end1 < -1e-8 || acc_end2 < -1e-8) {
				std::cerr << "Other end penetration occurred." << std::endl;
				fForceLineCenter = true;
				continue;
			}
		}

		// - check if zero moment condition is violated for line contact,
		if(line_cop_type == COPContactType::LineCenter) {
			if(abs(full_f_sol[line_start_row_id+3]) > 1e-3) {
				// std::cout << "Full F sol line: " << full_f_sol.segment<5>(line_start_row_id).transpose() << std::endl;
				// - - if so, compute distance to new line-COP.
				double dist_move_cop = -full_f_sol[line_start_row_id+3]/fmax(full_f_sol[line_start_row_id+2], 1e-5);
				// - - if distance to new line-COP is very small, continue to check other violations
				if(abs(dist_move_cop) > 0.02) {// TODO: integrate bisection search and lower this threshold
					signed_distance_from_point0 += dist_move_cop;
					signed_distance_from_point0 = fmin(signed_distance_from_point0, max_signed_dist);
					signed_distance_from_point0 = fmax(signed_distance_from_point0, min_signed_dist);
					did_update_cop_pos = true;
					// std::cout << "New cop pos " << signed_distance_from_point0 << std::endl;
					// std::cout << "Update dist " << dist_move_cop << std::endl;
					// std::cout << "Max dist " << max_signed_dist << std::endl;
					// std::cout << "Min dist " << min_signed_dist << std::endl;
					//TODO: if we have fForceLineCenter set to true, consider setting it false
					continue;
				}
				// - - if distance is small, check first if we are currently doing a binary search
				// - - if so, continue binary search with point state fixed.
				// - - if not, check last line-COP moment. if we switch moment directions, start binary
				// - - search with point state fixed.
			}
		}

		// std::cout << "Success: COP pos " << local_line_cop_sol_point.transpose() << std::endl;
		ret_sol.force_sol = full_f_sol;
		ret_sol.result = COPSolResult::Success;
		ret_sol.local_cop_pos = local_line_cop_sol_point;
		ret_sol.cop_type = line_cop_type;
		return ret_sol;
	}
	std::cerr << "solveStartWithPatchCentroidOnePointOneLine failed to converge" << std::endl;
	ret_sol.result = COPSolResult::NoSolution;
	return ret_sol;
}

void COPSolver::computeMatrices(const Eigen::MatrixXd& A_disp, const Eigen::VectorXd& rhs_disp, double mu, double mu_rot) {
	TA_size = 0;
	uint lrid = line_start_row_id;
	uint prid = pt_start_row_id;

	//TODO: code simplifications (expected ~50% code reduction):
	// - compute slip directions for point, line-translation and line-rotation only once
	// at the beginning
	// - compute the independent blocks of the matrices and the rhs vectors first
	// - - might need to compute the block sizes first for this
	// - compute the TA coupling blocks

	if(point_state == FrictionState::Sliding || point_state == FrictionState::Impending) {
		Vector2d pt_slip_dir;
		if(point_state == FrictionState::Sliding) {
			pt_slip_dir = point_slip_dir;
		} else {
			pt_slip_dir = point_rolling_dir;
		}

		if(line_cop_type == COPContactType::LineCenter) {
			// line: both translation and rolling friction forces are limited
			if((line_translation_state == FrictionState::Sliding ||
				line_translation_state == FrictionState::Impending) &&
				(line_rotational_state == FrictionState::Sliding ||
				line_rotational_state == FrictionState::Impending)
			) {
				Vector2d ltslip_dir;
				double lrslip_dir;

				if(line_translation_state == FrictionState::Sliding) {
					ltslip_dir = line_translation_slip_dir;
				} else {
					ltslip_dir = line_translation_rolling_dir;
				}
				if(line_rotational_state == FrictionState::Sliding) {
					lrslip_dir = line_rotation_slip_dir;
				} else {
					lrslip_dir = line_rotation_rolling_dir;
				}

				Trhs.segment<2>(0) = -rhs_disp.segment<2>(lrid+2);
				TA.block<2,1>(0,0) = A_disp.block<2,1>(lrid+2,lrid+2)
							- mu*A_disp.block<2,2>(lrid+2, lrid+0)*(ltslip_dir)
							- mu_rot*lrslip_dir*A_disp.block<2,1>(lrid+2,lrid+4);
				TA.block<2,1>(0,1) = A_disp.block<2,1>(lrid+2,lrid+3);

				// pt stuff
				Trhs(2) = -rhs_disp(prid+2);
				TA(2,2) = A_disp(prid+2, prid+2)
							- mu*A_disp.block<1,2>(prid+2,prid+0)*(pt_slip_dir);

				// coupling stuff
				TA.block<2,1>(0,2) = A_disp.block<2,1>(lrid+2, prid+2)
							- mu*A_disp.block<2,2>(lrid+2,prid+0)*(pt_slip_dir);
				TA(2,0) = A_disp(prid+2, lrid+2)
							- mu*A_disp.block<1,2>(prid+2,lrid+0)*(ltslip_dir)
							- mu_rot*lrslip_dir*A_disp(prid+2,lrid+4);
				TA(2,1) = A_disp(prid+2, lrid+3);
				TA_size = 3;
			}
			// line: translation forces are limited, rolling forces are not
			else if((line_translation_state == FrictionState::Sliding ||
				line_translation_state == FrictionState::Impending) &&
				line_rotational_state == FrictionState::Rolling
			) {
				Vector2d ltslip_dir;

				if(line_translation_state == FrictionState::Sliding) {
					ltslip_dir = line_translation_slip_dir;
				} else {
					ltslip_dir = line_translation_rolling_dir;
				}

				// line stuff
				Trhs.segment<3>(0) = -rhs_disp.segment<3>(lrid+2);
				TA.block<3,1>(0,0) = A_disp.block<3,1>(lrid+2,lrid+2)
							- mu*A_disp.block<3,2>(lrid+2, lrid+0)*(ltslip_dir);
				TA.block<3,2>(0,1) = A_disp.block<3,2>(lrid+2,lrid+3);

				// pt stuff
				Trhs(3) = -rhs_disp(prid+2);
				TA(3,3) = A_disp(prid+2, prid+2)
							- mu*A_disp.block<1,2>(prid+2,prid+0)*(pt_slip_dir);

				// coupling stuff
				TA.block<3,1>(0,3) = A_disp.block<3,1>(lrid+2, prid+2)
							- mu*A_disp.block<3,2>(lrid+2,prid+0)*(pt_slip_dir);
				TA(3,0) = A_disp(prid+2, lrid+2)
							- mu*A_disp.block<1,2>(prid+2,lrid+0)*(ltslip_dir);
				TA.block<1,2>(3,1) = A_disp.block<1,2>(prid+2,lrid+3);
				TA_size = 4;
			}
			// line: rotational forces are limited, translation forces are not
			else if((line_rotational_state == FrictionState::Sliding ||
				line_rotational_state == FrictionState::Impending) &&
				line_translation_state == FrictionState::Rolling
			) {
				double lrslip_dir;

				if(line_rotational_state == FrictionState::Sliding) {
					lrslip_dir = line_rotation_slip_dir;
				} else {
					lrslip_dir = line_rotation_rolling_dir;
				}

				// line stuff
				Trhs.segment<4>(0) = -rhs_disp.segment<4>(lrid);
				TA.block<4,2>(0,0) = A_disp.block<4,2>(lrid+0,lrid+0);
				TA.block<4,1>(0,2) = A_disp.block<4,1>(lrid+0,lrid+2)
									- mu_rot*lrslip_dir*A_disp.block<4,1>(lrid+0,lrid+4);
				TA.block<4,1>(0,3) = A_disp.block<4,1>(lrid+0,lrid+3);

				// pt stuff
				Trhs(4) = -rhs_disp(prid+2);
				TA(4,4) = A_disp(prid+2, prid+2)
							- mu*A_disp.block<1,2>(prid+2,prid+0)*(pt_slip_dir);

				// coupling stuff
				TA.block<4,1>(0,4) = A_disp.block<4,1>(lrid+0, prid+2)
							- mu*A_disp.block<4,2>(lrid+0,prid+0)*(pt_slip_dir);
				TA.block<1,2>(4,0) = A_disp.block<1,2>(prid+2, lrid+0);
				TA(4,2) = A_disp(prid+2, lrid+2)
							- mu_rot*lrslip_dir*A_disp(prid+2,lrid+4);
				TA(4,3) = A_disp(prid+2, lrid+3);
				TA_size = 5;
			}
			// line: translational and rotational forces are not limited
			else if(line_rotational_state == FrictionState::Rolling &&
				line_translation_state == FrictionState::Rolling
			) {
				// line stuff
				Trhs.segment<5>(0) = -rhs_disp.segment<5>(lrid);
				TA.block<5,5>(0,0) = A_disp.block<5,5>(lrid+0,lrid+0);

				// pt stuff
				Trhs(5) = -rhs_disp(prid+2);
				TA(5,5) = A_disp(prid+2, prid+2)
							- mu*A_disp.block<1,2>(prid+2,prid+0)*(pt_slip_dir);

				// coupling stuff
				TA.block<5,1>(0,5) = A_disp.block<5,1>(lrid+0, prid+2)
							- mu*A_disp.block<5,2>(lrid+0,prid+0)*(pt_slip_dir);
				TA.block<1,5>(5,0) = A_disp.block<1,5>(prid+2, lrid+0);
				TA_size = 6;
			} else {
				throw(std::runtime_error("Unimplemented case"));
			}
		} else if(line_cop_type == COPContactType::LineEnd) {
			// TODO: combine with simple 2-pt solution
			if((line_translation_state == FrictionState::Sliding ||
				line_translation_state == FrictionState::Impending)
			) {
				Vector2d ltslip_dir;

				if(line_translation_state == FrictionState::Sliding) {
					ltslip_dir = line_translation_slip_dir;
				} else {
					ltslip_dir = line_translation_rolling_dir;
				}

				Trhs(0) = -rhs_disp(lrid+2);
				TA(0,0) = A_disp(lrid+2,lrid+2)
							- mu*A_disp.block<1,2>(lrid+2, lrid+0)*(ltslip_dir);

				// pt stuff
				Trhs(1) = -rhs_disp(prid+2);
				TA(1,1) = A_disp(prid+2, prid+2)
							- mu*A_disp.block<1,2>(prid+2,prid+0)*(pt_slip_dir);

				// coupling stuff
				TA(0,1) = A_disp(lrid+2, prid+2)
							- mu*A_disp.block<1,2>(lrid+2,prid+0)*(pt_slip_dir);
				TA(1,0) = A_disp(prid+2, lrid+2)
							- mu*A_disp.block<1,2>(prid+2,lrid+0)*(ltslip_dir);
				TA_size = 2;
			} else { // rolling at the line end point
				// line stuff
				Trhs.segment<3>(0) = -rhs_disp.segment<3>(lrid);
				TA.block<3,3>(0,0) = A_disp.block<3,3>(lrid+0,lrid+0);

				// pt stuff
				Trhs(3) = -rhs_disp(prid+2);
				TA(3,3) = A_disp(prid+2, prid+2)
							- mu*A_disp.block<1,2>(prid+2,prid+0)*(pt_slip_dir);

				// coupling stuff
				TA.block<3,1>(0,3) = A_disp.block<3,1>(lrid+0, prid+2)
							- mu*A_disp.block<3,2>(lrid+0,prid+0)*(pt_slip_dir);
				TA.block<1,3>(3,0) = A_disp.block<1,3>(prid+2, lrid+0);
				TA_size = 4;
			}
		}
	} else if(point_state == FrictionState::Rolling) {
		if(line_cop_type == COPContactType::LineCenter) {
			// line: both translation and rolling friction forces are limited
			if((line_translation_state == FrictionState::Sliding ||
				line_translation_state == FrictionState::Impending) &&
				(line_rotational_state == FrictionState::Sliding ||
				line_rotational_state == FrictionState::Impending)
			) {
				Vector2d ltslip_dir;
				double lrslip_dir;

				if(line_translation_state == FrictionState::Sliding) {
					ltslip_dir = line_translation_slip_dir;
				} else {
					ltslip_dir = line_translation_rolling_dir;
				}
				if(line_rotational_state == FrictionState::Sliding) {
					lrslip_dir = line_rotation_slip_dir;
				} else {
					lrslip_dir = line_rotation_rolling_dir;
				}

				// line stuff
				Trhs.segment<2>(0) = -rhs_disp.segment<2>(lrid+2);
				TA.block<2,1>(0,0) = A_disp.block<2,1>(lrid+2,lrid+2)
							- mu*A_disp.block<2,2>(lrid+2, lrid+0)*(ltslip_dir)
							- mu_rot*lrslip_dir*A_disp.block<2,1>(lrid+2,lrid+4);
				TA.block<2,1>(0,1) = A_disp.block<2,1>(lrid+2,lrid+3);

				// pt stuff
				Trhs.segment<3>(2) = -rhs_disp.segment<3>(prid+0);
				TA.block<3,3>(2,2) = A_disp.block<3,3>(prid+0, prid+0);

				// coupling stuff
				TA.block<2,3>(0,2) = A_disp.block<2,3>(lrid+2, prid+0);
				TA.block<3,1>(2,0) = A_disp.block<3,1>(prid+0, lrid+2)
							- mu*A_disp.block<3,2>(prid+0,lrid+0)*(ltslip_dir)
							- mu_rot*lrslip_dir*A_disp.block<3,1>(prid+0,lrid+4);
				TA.block<3,1>(2,1) = A_disp.block<3,1>(prid+0, lrid+3);
				TA_size = 5;
			}
			// line: translation forces are limited, rolling forces are not
			else if((line_translation_state == FrictionState::Sliding ||
				line_translation_state == FrictionState::Impending) &&
				line_rotational_state == FrictionState::Rolling
			) {
				Vector2d ltslip_dir;

				if(line_translation_state == FrictionState::Sliding) {
					ltslip_dir = line_translation_slip_dir;
				} else {
					ltslip_dir = line_translation_rolling_dir;
				}

				// line stuff
				Trhs.segment<3>(0) = -rhs_disp.segment<3>(lrid+2);
				TA.block<3,1>(0,0) = A_disp.block<3,1>(lrid+2,lrid+2)
							- mu*A_disp.block<3,2>(lrid+2, lrid+0)*(ltslip_dir);
				TA.block<3,2>(0,1) = A_disp.block<3,2>(lrid+2,lrid+3);

				// pt stuff
				Trhs.segment<3>(3) = -rhs_disp.segment<3>(prid+0);
				TA.block<3,3>(3,3) = A_disp.block<3,3>(prid+0, prid+0);

				// coupling stuff
				TA.block<3,3>(0,3) = A_disp.block<3,3>(lrid+2, prid+0);
				TA.block<3,1>(3,0) = A_disp.block<3,1>(prid+0, lrid+2)
							- mu*A_disp.block<3,2>(prid+0,lrid+0)*(ltslip_dir);
				TA.block<3,2>(3,1) = A_disp.block<3,2>(prid+0,lrid+3);
				TA_size = 6;
			}
			// line: rotational forces are limited, translation forces are not
			else if((line_rotational_state == FrictionState::Sliding ||
				line_rotational_state == FrictionState::Impending) &&
				line_translation_state == FrictionState::Rolling
			) {
				double lrslip_dir;

				if(line_rotational_state == FrictionState::Sliding) {
					lrslip_dir = line_rotation_slip_dir;
				} else {
					lrslip_dir = line_rotation_rolling_dir;
				}

				// line stuff
				Trhs.segment<4>(0) = -rhs_disp.segment<4>(lrid);
				TA.block<4,2>(0,0) = A_disp.block<4,2>(lrid+0,lrid+0);
				TA.block<4,1>(0,2) = A_disp.block<4,1>(lrid+0,lrid+2)
									- mu_rot*lrslip_dir*A_disp.block<4,1>(lrid+0,lrid+4);
				TA.block<4,1>(0,3) = A_disp.block<4,1>(lrid+0,lrid+3);

				// pt stuff
				Trhs.segment<3>(4) = -rhs_disp.segment<3>(prid+0);
				TA.block<3,3>(4,4) = A_disp.block<3,3>(prid+0, prid+0);

				// coupling stuff
				TA.block<4,3>(0,4) = A_disp.block<4,3>(lrid+0, prid+0);
				TA.block<3,2>(4,0) = A_disp.block<3,2>(prid+0, lrid+0);
				TA.block<3,1>(4,2) = A_disp.block<3,1>(prid+0, lrid+2)
							- mu_rot*lrslip_dir*A_disp.block<3,1>(prid+0,lrid+4);
				TA.block<3,1>(4,3) = A_disp.block<3,1>(prid+0, lrid+3);
				TA_size = 7;
			}
			// line: translational and rotational forces are not limited
			else if(line_rotational_state == FrictionState::Rolling &&
				line_translation_state == FrictionState::Rolling
			) {
				// line stuff
				Trhs.segment<5>(0) = -rhs_disp.segment<5>(lrid);
				TA.block<5,5>(0,0) = A_disp.block<5,5>(lrid+0,lrid+0);

				// pt stuff
				Trhs.segment<3>(5) = -rhs_disp.segment<3>(prid+0);
				TA.block<3,3>(5,5) = A_disp.block<3,3>(prid+0, prid+0);

				// coupling stuff
				TA.block<5,3>(0,5) = A_disp.block<5,3>(lrid+0, prid+0);
				TA.block<3,5>(5,0) = A_disp.block<3,5>(prid+0, lrid+0);
				TA_size = 8;
			} else {
				throw(std::runtime_error("Unimplemented case"));
			}
		} else if(line_cop_type == COPContactType::LineEnd) {
			// TODO: combine with simple 2-pt solution
			if((line_translation_state == FrictionState::Sliding ||
				line_translation_state == FrictionState::Impending)
			) {
				Vector2d ltslip_dir;

				if(line_translation_state == FrictionState::Sliding) {
					ltslip_dir = line_translation_slip_dir;
				} else {
					ltslip_dir = line_translation_rolling_dir;
				}

				Trhs(0) = -rhs_disp(lrid+2);
				TA(0,0) = A_disp(lrid+2,lrid+2)
							- mu*A_disp.block<1,2>(lrid+2, lrid+0)*(ltslip_dir);

				// pt stuff
				Trhs.segment<3>(1) = -rhs_disp.segment<3>(prid+0);
				TA.block<3,3>(1,1) = A_disp.block<3,3>(prid+0, prid+0);

				// coupling stuff
				TA.block<1,3>(0,1) = A_disp.block<1,3>(lrid+2, prid+0);
				TA.block<3,1>(1,0) = A_disp.block<3,1>(prid+0, lrid+2)
							- mu*A_disp.block<3,2>(prid+0,lrid+0)*(ltslip_dir);
				TA_size = 4;
			} else { // rolling at the line end point
				// line stuff
				Trhs.segment<3>(0) = -rhs_disp.segment<3>(lrid);
				TA.block<3,3>(0,0) = A_disp.block<3,3>(lrid+0,lrid+0);

				// pt stuff
				Trhs.segment<3>(3) = -rhs_disp.segment<3>(prid+0);
				TA.block<3,3>(3,3) = A_disp.block<3,3>(prid+0, prid+0);

				// coupling stuff
				TA.block<3,3>(0,3) = A_disp.block<3,3>(lrid+0, prid+0);
				TA.block<3,3>(3,0) = A_disp.block<3,3>(prid+0, lrid+0);
				TA_size = 6;
			}
		}
	} else {
		throw(std::runtime_error("Unimplemented case"));
	}
}

void COPSolver::solve(double mu, double mu_rot, Eigen::VectorXd& full_f_sol) {
	Tf_sol = TA.block(0,0,TA_size,TA_size).partialPivLu().solve(Trhs.segment(0,TA_size));

	uint lrid = line_start_row_id;
	uint prid = pt_start_row_id;


	// get line solution
	Vector2d ltslip_dir;
	double lrslip_dir;

	if(line_translation_state == FrictionState::Sliding) {
		ltslip_dir = line_translation_slip_dir;
	} else {
		ltslip_dir = line_translation_rolling_dir;
	}
	if(line_rotational_state == FrictionState::Sliding) {
		lrslip_dir = line_rotation_slip_dir;
	} else {
		lrslip_dir = line_rotation_rolling_dir;
	}

	uint tp_ind = 0;

	if(line_cop_type == COPContactType::LineCenter) {
		if(line_translation_state == FrictionState::Sliding ||
			line_translation_state == FrictionState::Impending
		) {
			full_f_sol.segment<2>(lrid+0) = -mu*Tf_sol(0)*ltslip_dir;
			full_f_sol.segment<2>(lrid+2) = Tf_sol.segment<2>(0);
			if(line_rotational_state == FrictionState::Sliding ||
				line_rotational_state == FrictionState::Impending
			) {
				full_f_sol(lrid+4) = -mu_rot*Tf_sol(0)*lrslip_dir;
				tp_ind = 2;
			} else { // rotation rolling
				full_f_sol(lrid+4) = Tf_sol(2);
				tp_ind = 3;
			}
		} else { // translation rolling
			full_f_sol.segment<4>(lrid+0) = Tf_sol.segment<4>(0);
			if(line_rotational_state == FrictionState::Sliding ||
				line_rotational_state == FrictionState::Impending
			) {
				full_f_sol(lrid+4) = -mu_rot*Tf_sol(2)*lrslip_dir;
				tp_ind = 4;
			} else { // rotation rolling
				full_f_sol(lrid+4) = Tf_sol(4);
				tp_ind = 5;
			}
		}
	} else if(line_cop_type == COPContactType::LineEnd) {
		full_f_sol.segment<2>(lrid+3) = Vector2d::Zero();
		if(line_translation_state == FrictionState::Sliding ||
			line_translation_state == FrictionState::Impending
		) {
			full_f_sol.segment<2>(lrid+0) = -mu*Tf_sol(0)*ltslip_dir;
			full_f_sol(lrid+2) = Tf_sol(0);
			tp_ind = 1;
		} else { // translation rolling
			full_f_sol.segment<3>(lrid+0) = Tf_sol.segment<3>(0);
			tp_ind = 3;
		}
	}

	// get point solution
	Vector2d pt_slip_dir;
	if(point_state == FrictionState::Sliding) {
		pt_slip_dir = point_slip_dir;
	} else {
		pt_slip_dir = point_rolling_dir;
	}
	if(point_state == FrictionState::Sliding ||
		point_state == FrictionState::Impending
	) {
		full_f_sol.segment<2>(prid+0) = -mu*Tf_sol(tp_ind+0)*pt_slip_dir;
		full_f_sol(prid+2) = Tf_sol(tp_ind+0);
	} else {// pt rolling
		full_f_sol.segment<3>(prid+0) = Tf_sol.segment<3>(tp_ind+0);
	}
}

bool COPSolver::isAccelerationPenetrating(
	const Eigen::Vector3d& cop_to_test_pos,
	const Eigen::Vector3d& cop_lin_acc,
	const Eigen::Vector3d& cop_ang_acc,
	const Eigen::Vector3d& omegaA,
	const Eigen::Vector3d& omegaB
) {
	Vector3d shift_nonlin_acc = omegaB.cross(omegaB.cross(cop_to_test_pos))
									- omegaA.cross(omegaA.cross(cop_to_test_pos));
	double testpt_z_acc = cop_lin_acc(2)
							- cop_ang_acc(1)*(cop_to_test_pos)(0)
							+ cop_ang_acc(0)*(cop_to_test_pos)(1)
							+ shift_nonlin_acc(2);
	return (testpt_z_acc < -1e-8);
}

void COPSolver::getCOPSurfaceContactDisplacedMatrices(
	Eigen::MatrixXd& A_disp,
	Eigen::VectorXd& rhs_disp,
	Eigen::Vector3d& lin_vel_disp,
	const Eigen::Vector3d& displacement,
	const Eigen::MatrixXd& A,
	const Eigen::VectorXd& rhs,
	const Eigen::Vector3d& omegaA,
	const Eigen::Vector3d& omegaB,
	const Eigen::Vector3d& lin_vel
) {
	MatrixXd transform_mat = MatrixXd::Identity(6,6);
	transform_mat.block<3,3>(0,3) = -crossMat(displacement);
	A_disp = transform_mat*A*transform_mat.transpose();
	rhs_disp = transform_mat*rhs;
	rhs_disp.segment<3>(0) += omegaB.cross(omegaB.cross(displacement))
				- omegaA.cross(omegaA.cross(displacement));
	lin_vel_disp = lin_vel + (omegaB - omegaA).cross(displacement);
}

static double getMuRotationSurface(double mu, double min_boundary_dist, double max_patch_extent) {
	return 2*mu*min_boundary_dist*(max_patch_extent - min_boundary_dist)/max_patch_extent;
}

ContactCOPSolution COPSolver::solveOneSurfaceContactWithLCP(
		double friction_coeff,
		const Eigen::MatrixXd& A_constraint, // defined at Contact Patch Interior point.
		const Eigen::VectorXd& rhs_constraint, // defined at Contact Patch Interior point. Includes Jdot_qdot
		const ContactPatch& contact_patch,
		const Eigen::Vector3d& omegaA, // vector of body A angular velocities, in COP frame, one entry per contact patch
		const Eigen::Vector3d& omegaB, // vector of body B angular velocities, in COP frame, one entry per contact patch
		const Eigen::Vector3d& linear_contact_velocity // 3 dof relative translation velocity at point 0, in COP frame, one entry per contact patch
		//^ expected to be zero in the z direction, but we don't explicitly check
) {
	// No need to check for constraints being satisfied, this optimization will be called
	// before this method is called.
	ContactCOPSolution ret_sol;

	// Check if the center of rotation of the object in the contact plane lies
	// within the contact patch.
	// This is required to have a proper estimate of the Coulomb friction.
	double rot_slip = omegaB(2) - omegaA(2);
	bool fForceIgnoreSlip = false;
	if (abs(rot_slip) > 1e-8) {
		Vector2d rot_point;
		rot_point << -linear_contact_velocity(1), linear_contact_velocity(0);
		rot_point /= rot_slip;
		if(rot_point.norm() < contact_patch.max_extent/2.0) {
			fForceIgnoreSlip = true;
		}
	}

	// Iteration variables
	const int max_iters = 40;
	int iter_ind = 0;
	Eigen::MatrixXd A_disp(6, 6);
	Eigen::VectorXd rhs_disp(6);
	Vector3d lin_vel_disp;
	const Vector3d lin_vel = linear_contact_velocity;
	Eigen::VectorXd full_f_sol(6), full_a_sol(6);
	double mu_rot;
	const double mu = friction_coeff;

	// Set cop_point initially to patch interior point.
	Vector3d cop_point = Vector3d::Zero();
	COPContactType patch_cop_type = COPContactType::PatchCenter;

	auto lcp_solver = Sai2LCPSolver::LCPSolver();
	while (iter_ind < max_iters) {
		iter_ind++;

		// -- Assemble LCP variables --
		getCOPSurfaceContactDisplacedMatrices(
			A_disp, rhs_disp, lin_vel_disp,
			cop_point,
			A_constraint, rhs_constraint,
			omegaA, omegaB,
			lin_vel
		);

		// Ignore slip in all future iterations if we are not slipping at the current
		// COP point.
		if(!fForceIgnoreSlip &&
		   lin_vel_disp.segment<2>(0).norm() < 1e-6) {
		   	fForceIgnoreSlip = true;
		}
		// LCP solver expects pre_v to be the same size as vector b. So we add 0's to pad.
		VectorXd lin_vel_disp_for_lcp(5);
		lin_vel_disp_for_lcp << lin_vel_disp, 0, 0;
		if (fForceIgnoreSlip) {
			lin_vel_disp_for_lcp.segment<2>(0).setZero();
		}

		auto lcp_moment_constraint = Sai2LCPSolver::MomentConstraints::XandYMoments;
		if (patch_cop_type != COPContactType::PatchCenter) {
			lcp_moment_constraint = Sai2LCPSolver::MomentConstraints::NoMoment;
		}

		PointTestResult test_result = contact_patch.testPoint(cop_point.head<2>());
		mu_rot = getMuRotationSurface(mu, test_result.min_dist_to_boundary,
									  contact_patch.max_extent);
		const double visc_rot_torque = -mu_rot*rot_slip;
		rhs_disp += A_disp.col(5)*visc_rot_torque;
		if(COP_LOG_DEBUG) {
			std::cout << "COPSolver Patch LCP: computed mu rot: " << mu_rot << std::endl;
			std::cout << "COPSolver Patch LCP: computed viscous rot torque: "
					  << visc_rot_torque << std::endl;
		}

		// -- Solve LCP --
		auto lcp_result = lcp_solver.solveWithMoments(
			A_disp.block<5, 5>(0, 0),
			rhs_disp.head<5>(),
			lin_vel_disp_for_lcp,
			{lcp_moment_constraint},
			mu);
		if (lcp_result.result != Sai2LCPSolver::LCPSolResult::Success) {
			break;
		}
		full_f_sol << lcp_result.p_sol, visc_rot_torque;
		// std::cout << f_sol.transpose() << std::endl;

		// Don't multiply A_disp with full_f_sol directly since we added
		// A_disp.col(5)*visc_rot_torque to rhs_disp earlier.
		full_a_sol = A_disp.block<6, 5>(0, 0)*full_f_sol.head<5>() +
						rhs_disp;

		// -- Check for COP and non-penetration constraints --
		if (patch_cop_type == COPContactType::PatchCurvePoint) {
			// Check for penetration at the interior point.
			Vector3d test_disp = -cop_point;
			if(isAccelerationPenetrating(test_disp,
				full_a_sol.segment<3>(0),
				full_a_sol.segment<3>(3),
				omegaA,
				omegaB
			)) {
				// We dont expect this case since we start from the patch center point currently.
				std::cerr << "Other end penetration occurred." << std::endl;
				patch_cop_type = COPContactType::PatchCenter;
				// No need to check for COP constraint since moment is zero.
				continue;
			}
		} else {
			// patch_cop_type = COPContactType::PatchCenter, check for COP constraint.
			Vector2d violation_moment = full_f_sol.segment<2>(3);
			if(violation_moment.norm() > 1e-3) {
				// Cap max COP displacement distance to violation_moment.norm()/1e-5 to
				// prevent divide by zero in case the normal force is really small.
				// Note: cop_disp_dist is non-negative.
				const double cop_disp_dist = violation_moment.norm()/fmax(full_f_sol[2], 1e-5);

				// - - if distance to new COP is very small, we are done
				if(cop_disp_dist > 0.01) {
					Vector2d cop_disp_dir;
					cop_disp_dir << -violation_moment(1), violation_moment(0);
					cop_disp_dir /= cop_disp_dir.norm();
					const double max_dist_to_bdry =
						contact_patch.distanceFromBoundaryAlongRay(cop_point.head<2>(),
																   cop_disp_dir);
					if(cop_disp_dist < max_dist_to_bdry) {
						cop_point.head<2>() += cop_disp_dir*cop_disp_dist;
					} else {
						// we are on the boundary
						cop_point.head<2>() += cop_disp_dir*max_dist_to_bdry;
						patch_cop_type = COPContactType::PatchCurvePoint;
					}
					continue;
				}
			}
		}
		// We are done, return success
		ret_sol.force_sol = full_f_sol;
		ret_sol.result = COPSolResult::Success;
		ret_sol.local_cop_pos = cop_point;
		ret_sol.cop_type = patch_cop_type;
		return ret_sol;
	}
	std::cerr << "solveSurfaceContact failed to converge" << std::endl;
	ret_sol.result = COPSolResult::NoSolution;
	return ret_sol;
}

std::vector<ContactCOPSolution> COPSolver::solveTwoSurfaceContactsWithLCP(
		double friction_coeff,
		const Eigen::MatrixXd& A_constraint,
		const Eigen::VectorXd& rhs_constraint,
		const std::vector<uint>& patch_indices, // ith entry gives starting row of ith contact patch in A_constraint, rhs_constraint
		const std::vector<ContactType>& contact_types, // ith entry gives contact type of ith contact patch
		const std::vector<ContactPatch>& contact_patches,
		const std::vector<Eigen::Vector3d>& omega_bodyA, // vector of body A angular velocities, in COP frame, one entry per contact patch
		const std::vector<Eigen::Vector3d>& omega_bodyB, // vector of body B angular velocities, in COP frame, one entry per contact patch
		const std::vector<Eigen::Vector3d>& linear_contact_velocity // 3 dof relative translation velocity at point 0, in COP frame, one entry per contact patch
		//^ expected to be zero in the z direction, but we don't explicitly check
	) {
	constexpr bool print_time_analytics = false;
	assert(contact_types.size() == 2);
	assert(patch_indices[0] == 0);
	assert(patch_indices[1] == 6);

	auto pre_time_pt1 = std::chrono::high_resolution_clock::now();
	if(contact_types[0] != ContactType::SURFACE ||
	   contact_types[1] != ContactType::SURFACE) {
		throw(std::runtime_error(
			"Unsupported types " + std::to_string(contact_types[0]) + " " +
			std::to_string(contact_types[1])));
	}

	// Check special cases
	auto is_ith_surface_separating = [&](uint i) -> bool {
		return isSurfaceContactSeparating(
					rhs_constraint.segment<3>(patch_indices[i]),
				    rhs_constraint.segment<3>(patch_indices[i] + 3),
				    omega_bodyA[i], omega_bodyB[i],
				    contact_patches[i]);
	};
	auto is_ith_surface_still_separating = [&](uint i,
											   const VectorXd& new_acc) -> bool {
		return isSurfaceContactSeparating(
					new_acc.segment<3>(0),
				    new_acc.segment<3>(3),
				    omega_bodyA[i], omega_bodyB[i],
				    contact_patches[i]);
	};
	auto solve_one_surface_contact = [&](uint i) -> ContactCOPSolution {
		return solveOneSurfaceContactWithLCP(
					friction_coeff,
					A_constraint.block<6, 6>(patch_indices[i], patch_indices[i]),
					rhs_constraint.segment<6>(patch_indices[i]),
					contact_patches[i],
					omega_bodyA[i], omega_bodyB[i],
					linear_contact_velocity[i]);
	};
	bool patch1_separating = is_ith_surface_separating(0);
	bool patch2_separating = is_ith_surface_separating(1);

	ContactCOPSolution no_force_sol;
	no_force_sol.force_sol = VectorXd::Zero(6);
	no_force_sol.local_cop_pos.setZero();
	no_force_sol.result = COPSolResult::Success;

	ContactCOPSolution failure_sol;
	no_force_sol.result = COPSolResult::NoSolution;

	if (patch1_separating && patch2_separating) {
		// - no contact on either surface
		return {no_force_sol, no_force_sol};
	} else if (patch1_separating) {
		// - check for contact on only patch2
		ContactCOPSolution patch2_sol = solve_one_surface_contact(1);
		VectorXd patch1_acc =
			A_constraint.block<6, 6>(patch_indices[0],
									 patch_indices[1])*patch2_sol.force_sol +
			rhs_constraint.segment<6>(patch_indices[0]);
		if (patch2_sol.result != COPSolResult::Success) {
			// unlikely that we'll find a solution
			return {failure_sol, failure_sol};
		} else if (is_ith_surface_still_separating(0, patch1_acc)) {
			return {no_force_sol, patch2_sol};
		}
	} else if (patch2_separating) {
		// - check for contact on only patch1
		ContactCOPSolution patch1_sol = solve_one_surface_contact(0);
		VectorXd patch2_acc =
			A_constraint.block<6, 6>(patch_indices[1],
									 patch_indices[0])*patch1_sol.force_sol +
			rhs_constraint.segment<6>(patch_indices[1]);
		if (patch1_sol.result != COPSolResult::Success) {
			// unlikely that we'll find a solution
			return {failure_sol, failure_sol};
		} else if (is_ith_surface_still_separating(1, patch2_acc)) {
			return {patch1_sol, no_force_sol};
		}
	}

	const double mu = friction_coeff;

	// Check if the center of rotation of the object in the contact plane lies
	// within the contact patch.
	// This is required to have a proper estimate of the Coulomb friction.
	std::vector<double> rot_slips = {omega_bodyB[0](2) - omega_bodyA[0](2),
									 omega_bodyB[1](2) - omega_bodyA[1](2)};
	auto ignore_slip = [&](uint i) -> bool {
		if (abs(rot_slips[i]) > 1e-8) {
			Vector2d rot_point;
			rot_point << -linear_contact_velocity[i](1), linear_contact_velocity[i](0);
			rot_point /= rot_slips[i];
			if(rot_point.norm() < contact_patches[i].max_extent/2.0) {
				return true;
			}
		}
		return false;
	};

	std::vector<bool> fForceIgnoreSlip = {ignore_slip(0), ignore_slip(1)};

	// Iteration variables
	const int max_iters = 40*patch_indices.size();
	int iter_ind = 0;
	Eigen::MatrixXd A_disp(12, 12);
	Eigen::VectorXd rhs_disp(12);
	std::vector<Vector3d> lin_vel_disp = linear_contact_velocity;
	const Vector3d lin_vel1 = linear_contact_velocity[0];
	const Vector3d lin_vel2 = linear_contact_velocity[1];
	Eigen::VectorXd full_f_sol(12), full_a_sol(12);

	// Set cop_point initially to patch interior point.
	std::vector<Vector3d> cop_points(2, Vector3d::Zero());
	Vector3d& cop_point1 = cop_points[0];
	Vector3d& cop_point2 = cop_points[1];
	std::vector<COPContactType> patches_cop_type(2, COPContactType::PatchCenter);
	COPContactType& patch1_cop_type = patches_cop_type[0];
	COPContactType& patch2_cop_type = patches_cop_type[1];

	// Helper lambdas.
	auto update_ignore_slip = [&](uint i) {
		// Ignore slip in all future iterations if we are not slipping at the current
		// COP point.
		if(!fForceIgnoreSlip[i] &&
		   lin_vel_disp[i].segment<2>(0).norm() < 1e-6) {
		   	fForceIgnoreSlip[i] = true;
		}
	};

	std::vector<double> visc_rot_torque(2, 0);
	auto update_viscous_torque = [&](uint i) {
		PointTestResult test_result =
			contact_patches[i].testPoint(cop_points[i].head<2>());
		double mu_rot = getMuRotationSurface(mu, test_result.min_dist_to_boundary,
									  contact_patches[i].max_extent);
		visc_rot_torque[i] = -mu_rot*rot_slips[i];
		rhs_disp += A_disp.col(6*i + 5)*visc_rot_torque[i];
		if(COP_LOG_DEBUG) {
			std::cout << "COPSolver 2-Patch LCP: patch " << i
					  << " computed mu rot: " << mu_rot << std::endl;
			std::cout << "COPSolver 2-Patch LCP: patch " << i
					  << " computed viscous rot torque: " << visc_rot_torque[i]
					  << std::endl;
		}
	};

	auto pre_time_pt2 = std::chrono::high_resolution_clock::now();

	auto lcp_solver = Sai2LCPSolver::LCPSolver();
	MatrixXd lcp_A_disp(10, 10);
	VectorXd lcp_rhs_disp(10);
	VectorXd lcp_lin_vel_disp(10);

	auto pre_time_pt3 = std::chrono::high_resolution_clock::now();
	double iter_time_total = 0;
	while (iter_ind < max_iters) {
		iter_ind++;

		// -- Assemble LCP variables --
		auto time_pt1 = std::chrono::high_resolution_clock::now();
		getDisplacedMatricesForMultipleSurfaceContacts(
			A_disp, rhs_disp, lin_vel_disp,
			cop_points,
			A_constraint, rhs_constraint,
			omega_bodyA, omega_bodyB,
			linear_contact_velocity);

		update_ignore_slip(0);
		update_ignore_slip(1);

		// LCP solver expects pre_v to be the same size as vector b. So we add 0's to pad
		// between linear velocity components for each point.
		lcp_lin_vel_disp << lin_vel_disp[0], 0, 0,
						    lin_vel_disp[1], 0, 0;
		// std::cout << "lcp_lin_vel_disp " << lcp_lin_vel_disp.transpose() << std::endl;
		if (fForceIgnoreSlip[0]) {
			lcp_lin_vel_disp.segment<2>(0).setZero();
		}
		if (fForceIgnoreSlip[1]) {
			lcp_lin_vel_disp.segment<2>(5).setZero();
		}

		auto lcp_moment_constraints =
			std::vector<Sai2LCPSolver::MomentConstraints>(
				2, Sai2LCPSolver::MomentConstraints::XandYMoments);
		if (patch1_cop_type != COPContactType::PatchCenter) {
			lcp_moment_constraints[0] = Sai2LCPSolver::MomentConstraints::NoMoment;
		}
		if (patch2_cop_type != COPContactType::PatchCenter) {
			lcp_moment_constraints[1] = Sai2LCPSolver::MomentConstraints::NoMoment;
		}

		auto time_pt2 = std::chrono::high_resolution_clock::now();
		update_viscous_torque(0);
		update_viscous_torque(1);

		auto time_pt3 = std::chrono::high_resolution_clock::now();
		// -- Solve LCP --
		// The LCP A matrix does not include row for rotational acc and col for
		// rotational torque for any patch. So we copy A_disp to lcp_A_disp without these
		// rows and columns.
		lcp_A_disp.block<5, 5>(0, 0) = A_disp.block<5, 5>(0, 0);
		lcp_A_disp.block<5, 5>(5, 5) = A_disp.block<5, 5>(6, 6);
		lcp_A_disp.block<5, 5>(0, 5) = A_disp.block<5, 5>(0, 6);
		lcp_A_disp.block<5, 5>(5, 0) = A_disp.block<5, 5>(6, 0);
		lcp_rhs_disp.segment<5>(0) = rhs_disp.segment<5>(0);
		lcp_rhs_disp.segment<5>(5) = rhs_disp.segment<5>(6);

		// std::cout << "lcp A " << std::endl << lcp_A_disp << std::endl;
		// std::cout << "lcp rhs " << lcp_rhs_disp.transpose() << std::endl;
		auto time_pt4 = std::chrono::high_resolution_clock::now();
		auto lcp_result = lcp_solver.solveWithMoments(
			lcp_A_disp,
			lcp_rhs_disp,
			lcp_lin_vel_disp,
			lcp_moment_constraints,
			mu);
		if (lcp_result.result != Sai2LCPSolver::LCPSolResult::Success) {
			break;
		}
		auto time_pt5 = std::chrono::high_resolution_clock::now();

		full_f_sol << lcp_result.p_sol.segment<5>(0), visc_rot_torque[0],
					  lcp_result.p_sol.segment<5>(5), visc_rot_torque[1];
		// std::cout << full_f_sol.transpose() << std::endl;

		// Don't multiply A_disp with full_f_sol directly since we added
		// the acc contributions from the viscous torques to rhs_disp earlier.
		full_a_sol = A_disp.block<12, 5>(0, 0)*full_f_sol.head<5>() +
		             	A_disp.block<12, 5>(0, 6)*full_f_sol.segment<5>(6) +
							rhs_disp;

		auto time_pt6 = std::chrono::high_resolution_clock::now();
		bool did_succeed = true;

		for (uint i = 0; i < patch_indices.size(); i++) {
			// -- Check for COP and non-penetration constraints --
			if (patches_cop_type[i] == COPContactType::PatchCurvePoint) {
				// Check for penetration at the interior point.
				Vector3d test_disp = -cop_points[i];
				if(isAccelerationPenetrating(test_disp,
					full_a_sol.segment<3>(6*i + 0),
					full_a_sol.segment<3>(6*i + 3),
					omega_bodyA[i],
					omega_bodyB[i]
				)) {
					// We dont expect this case since we start from the patch center point currently.
					std::cerr << "Other end penetration occurred for patch "
							  << i << std::endl;
					patches_cop_type[i] = COPContactType::PatchCenter;
					// No need to check for COP constraint since moment is zero.
					did_succeed = false;
					break;
				}
			} else {
				// patch_cop_type = COPContactType::PatchCenter, check for COP constraint.
				Vector2d violation_moment = full_f_sol.segment<2>(6*i + 3);
				if(violation_moment.norm() > 1e-3) {
					// Cap max COP displacement distance to violation_moment.norm()/1e-5 to
					// prevent divide by zero in case the normal force is really small.
					// Note: cop_disp_dist is non-negative.
					const double cop_disp_dist =
						violation_moment.norm()/fmax(full_f_sol[6*i + 2], 1e-5);

					// - - if distance to new COP is very small, we are done
					if(cop_disp_dist > 0.01) {
						Vector2d cop_disp_dir;
						cop_disp_dir << -violation_moment(1), violation_moment(0);
						cop_disp_dir /= cop_disp_dir.norm();
						const double max_dist_to_bdry =
							contact_patches[i].distanceFromBoundaryAlongRay(
								cop_points[i].head<2>(),cop_disp_dir);
						if(cop_disp_dist < max_dist_to_bdry) {
							cop_points[i].head<2>() += cop_disp_dir*cop_disp_dist;
						} else {
							// we are on the boundary
							cop_points[i].head<2>() += cop_disp_dir*max_dist_to_bdry;
							patches_cop_type[i] = COPContactType::PatchCurvePoint;
						}
						did_succeed = false;
						break;
					}
				}
			}
		}

		if (print_time_analytics) {
			auto time_pt7 = std::chrono::high_resolution_clock::now();
			iter_time_total += 1e3*std::chrono::duration<double>(time_pt7 - time_pt1).count();
			std::cout << "Iter took "
					  << 1e3*std::chrono::duration<double>(time_pt7 - time_pt1).count()
					  << "ms\n";
			std::cout << "Step 1: " << 1e3*std::chrono::duration<double>(time_pt2 - time_pt1).count()
			          << "Step 2: " << 1e3*std::chrono::duration<double>(time_pt3 - time_pt2).count()
			          << "Step 3: " << 1e3*std::chrono::duration<double>(time_pt4 - time_pt3).count()
			          << "Step 4: " << 1e3*std::chrono::duration<double>(time_pt5 - time_pt4).count()
			          << "Step 5: " << 1e3*std::chrono::duration<double>(time_pt6 - time_pt5).count()
			          << "Step 6: " << 1e3*std::chrono::duration<double>(time_pt7 - time_pt6).count()
			          << std::endl;
		}

		if (did_succeed) {
			// We are done, return success
		    auto copy_time_pt1 = std::chrono::high_resolution_clock::now();
			std::vector<ContactCOPSolution> ret_sols;
			for (uint i = 0; i < patch_indices.size(); i++) {
				ContactCOPSolution sol;
				sol.force_sol = full_f_sol.segment<6>(6*i);
				sol.result = COPSolResult::Success;
				sol.local_cop_pos = cop_points[i];
				sol.cop_type = patches_cop_type[i];
				ret_sols.push_back(std::move(sol));
			}
			auto copy_time_pt2 = std::chrono::high_resolution_clock::now();

			if (print_time_analytics) {
				double copy_time = 1e3*std::chrono::duration<double>(copy_time_pt2 - copy_time_pt1).count();
				double setup_time = 1e3*std::chrono::duration<double>(pre_time_pt3 - pre_time_pt1).count();
				std::cout << "Setup took " << setup_time << "ms\n";
				std::cout << "Result copy took " << copy_time << "ms\n";
				std::cout << "Two sol took "
						  << setup_time + copy_time + iter_time_total
						  << "ms\n";
			}
			return ret_sols;
		}
	}

	std::cerr << "solveSurfaceContact failed to converge" << std::endl;
	return {failure_sol, failure_sol};
}

ContactCOPSolution COPSolver::solveSurfaceContact(
	double friction_coeff,
	const Eigen::MatrixXd& A_constraint, // defined at Contact Patch Interior point, which is also origin for COP frame.
	const Eigen::VectorXd& rhs_constraint, // defined at Contact Patch Interior point. Includes Jdot_qdot
	const ContactPatch& contact_patch,
	const std::vector<uint>& patch_indices, // ith entry gives starting row of ith contact patch in A_constraint, rhs_constraint
	const std::vector<ContactType>& contact_types, // ith entry gives contact type of ith contact patch
	const std::vector<Eigen::Vector3d>& omega_bodyA, // vector of body A angular velocities, in COP frame, one entry per contact patch
	const std::vector<Eigen::Vector3d>& omega_bodyB, // vector of body B angular velocities, in COP frame, one entry per contact patch
	const std::vector<Eigen::Vector3d>& linear_contact_velocity // 3 dof relative translation velocity at point 0, in COP frame, one entry per contact patch
	//^ expected to be zero in the z direction, but we don't explicitly check
) {
	// NOTE: only single surface contact solver for now.
	assert(patch_indices.size() == 1);
	Vector3d omegaA = omega_bodyA[0];
	Vector3d omegaB = omega_bodyB[0];
	Vector3d lin_vel = linear_contact_velocity[0];

	ContactCOPSolution ret_sol;
	ret_sol.local_cop_pos.setZero();

	// check for separating contact
	if(rhs_constraint(2) > -1e-8) {
		if(rhs_constraint.segment<2>(3).norm() < 1e-8) {
			ret_sol.force_sol = VectorXd::Zero(6);
			ret_sol.result = COPSolResult::Success;
			return ret_sol;
		} else {
			// pick a worse case point to test for penetration due to angular acceleration
			Vector3d faraway_test_pt = rhs_constraint.segment<3>(3).cross(Vector3d(0, 0, 1));
			faraway_test_pt *= contact_patch.max_extent / faraway_test_pt.norm();
			if(!isAccelerationPenetrating(
				faraway_test_pt,
				rhs_constraint.segment<3>(0),
				rhs_constraint.segment<3>(3),
				omegaA,
				omegaB
			)) {
				ret_sol.force_sol = VectorXd::Zero(6);
				ret_sol.result = COPSolResult::Success;
				return ret_sol;
			}
		}
	}

	return solveOneSurfaceContactWithLCP(
		friction_coeff,
		A_constraint, rhs_constraint,
		contact_patch,
		omegaA, omegaB, lin_vel);

	// ----------------------------------------------------------------------------------

	// compute rotational slip
	bool fForceIgnoreSlip = false;

	patch_rotation_slip_dir = omegaB(2) - omegaA(2);
	if(abs(patch_rotation_slip_dir) > 1e-8) {
		patch_rotational_state = FrictionState::Sliding;
		double rot_slip = patch_rotation_slip_dir;
		patch_rotation_slip_dir /= abs(patch_rotation_slip_dir);

		// check if the center of rotation of the object in the contact plane lies
		// within the contact patch
		// This is required to have a proper estimate of the Coulomb friction
		Vector2d rot_point;
		rot_point << -lin_vel(1), lin_vel(0);
		rot_point /= rot_slip;
		if(rot_point.norm() < contact_patch.max_extent/2.0) { //TODO: switch to a smooth transition from 0 to 1 on mu
			// std::cout << "Rot point: " << rot_point.transpose() << std::endl;
			fForceIgnoreSlip = true;
			patch_translation_state = FrictionState::Rolling;
		}
	} else {
		// initialize to rolling
		patch_rotational_state = FrictionState::Rolling;
	}

	const int max_iters = 40;
	int iter_ind = 0;
	Eigen::MatrixXd A_disp(6, 6);
	Eigen::VectorXd rhs_disp(6);
	Vector3d lin_vel_disp;
	bool fForcePatchCenter = false;
	Eigen::VectorXd full_f_sol(6), full_a_sol(6);
	double mu_rot;
	const double mu = friction_coeff;
	Vector3d cop_point = Vector3d::Zero();
	// Vector3d cop_point;
	// cop_point << contact_patch._line_segments[0].point1, 0;
	// std::cout << cop_point.transpose() << std::endl;
	patch_cop_type = COPContactType::PatchCenter;

	// initial internal matrices
	TA = A_constraint;
	Trhs = rhs_constraint;
	TA_size = 0;

	bool did_update_cop_pos = true; // we set initially to true to compute mu_rot
	while (iter_ind < max_iters) {
		iter_ind++;

		// if patch-cop position was updated, update displaced terms
		if(did_update_cop_pos) {
			did_update_cop_pos = false;

			// get displaced matrices
			getCOPSurfaceContactDisplacedMatrices(
				A_disp, rhs_disp, lin_vel_disp,
				cop_point,
				A_constraint, rhs_constraint,
				omegaA, omegaB,
				lin_vel
			);

			// get min distance to patch boundary
			auto test_result = contact_patch.testPoint(cop_point.segment<2>(0));

			// debug stuff
			if(patch_cop_type == COPContactType::PatchCenter && !fForcePatchCenter) {
				assert(test_result.f_is_in_patch);
				assert(test_result.min_dist_to_boundary > -1e-4);
			} else {
				assert(test_result.f_is_on_vertex);
				assert(abs(test_result.min_dist_to_boundary) < 1e-4);
			}

			// compute mu_rotation //TODO: generalize getMuRotation function for patch and line
			mu_rot = getMuRotationSurface(mu, test_result.min_dist_to_boundary, contact_patch.max_extent);

			// get translational slip velocities at current COP
			if(!fForceIgnoreSlip) {
				patch_translation_slip_dir = lin_vel_disp.segment<2>(0);
				if (patch_translation_slip_dir.norm() > 1e-6) {
					patch_translation_state = FrictionState::Sliding;
					patch_translation_slip_dir /= patch_translation_slip_dir.norm();
				} else {
					// initialize to rolling
					patch_translation_state = FrictionState::Rolling;
					// if rolling at COP, set fForceIgnoreSlip true
					fForceIgnoreSlip = true;
				}
			}
		}

		// solve with current state
		computeMatricesForSurface(A_disp, rhs_disp, mu, mu_rot);
		solveForSurface(mu, mu_rot, full_f_sol);
		full_a_sol = A_disp * full_f_sol + rhs_disp;

		// check for normal force violation
		if(full_f_sol(2) < -1e-8) {
			std::cout << "Lin vel: " << lin_vel.transpose() << std::endl;
			std::cout << "Lin vel disp: " << lin_vel_disp.transpose() << std::endl;
			std::cout << "Displaced A mat: " << std::endl;
			std::cout << A_disp << std::endl;
			std::cout << "RHS: " << rhs_constraint.transpose() << std::endl;
			std::cout << "Trans state: " << patch_translation_state << std::endl;
			std::cout << "Rot state: " << patch_rotational_state << std::endl;
			std::cout << full_f_sol.transpose() << std::endl;
			std::cout << cop_point.transpose() << std::endl;
			ret_sol.result = COPSolResult::UnimplementedCase;
			return ret_sol;
		}

		// check for translation friction cone violation
		if(patch_translation_state == FrictionState::Rolling) {
			Vector2d rolling_force = full_f_sol.segment<2>(0);
			if(rolling_force.norm() > mu*full_f_sol(2) + 1e-15) {
				// save impending direction
				// TODO: proper impending slip direction
				patch_translation_rolling_dir = -rolling_force/rolling_force.norm();
				// switch point state to impending dir
				patch_translation_state = FrictionState::Impending;
				continue;
			}
		}
		// check for impending translation slip wrong rolling direction
		if(patch_translation_state == FrictionState::Impending) {
			Vector2d friction_force = full_f_sol.segment<2>(0);
			Vector2d slip_accel = full_a_sol.segment<2>(0);
			if(slip_accel.dot(friction_force) > 1e-8) {
				// switch back to rolling so that we recompute the rolling force
				patch_translation_state = FrictionState::Rolling;
				continue;
			}
		}
		// check for rotation friction limit violation
		if(patch_rotational_state == FrictionState::Rolling) {
			double rolling_force = full_f_sol(5);
			if(abs(rolling_force) > mu_rot*full_f_sol(2) + 1e-15) {
				// save impending direction
				patch_rotation_rolling_dir = -rolling_force/abs(rolling_force);
				// switch point state to impending dir
				patch_rotational_state = FrictionState::Impending;
				continue;
			}
		}
		// check for impending rotational slip wrong rolling direction
		if(patch_rotational_state == FrictionState::Impending) {
			double rolling_force = full_f_sol(5);
			double rolling_acc = full_a_sol(5);
			if(rolling_acc*rolling_force > 1e-8) {
				// switch back to rolling so that we recompute the rolling force
				patch_rotational_state = FrictionState::Rolling;
				continue;
			}
		}

		// - check for other end penetration. if so, force patch center solution. print to cerr.
		// - we dont expect this case since we start from the patch center point currently.
		if(patch_cop_type == COPContactType::PatchCurvePoint ||
			patch_cop_type == COPContactType::PatchVertex
		) {
			// TODO: perhaps we need to test at one more point that is not collinear
			// with the cop_point and the interior_point?
			Vector3d test_disp = -cop_point;// we are testing at the interior point
			if(isAccelerationPenetrating(test_disp,
				full_a_sol.segment<3>(0),
				full_a_sol.segment<3>(3),
				omegaA,
				omegaB
			)) {
				// std::cout << "Lin vel: " << lin_vel.transpose() << std::endl;
				// std::cout << "Omega: " << omegaB.transpose() << std::endl;
				// std::cout << "Lin vel disp: " << lin_vel_disp.transpose() << std::endl;
				// std::cout << "Displaced A mat: " << std::endl;
				// std::cout << A_disp << std::endl;
				// std::cout << "RHS: " << rhs_constraint.transpose() << std::endl;
				// std::cout << "Trans state: " << patch_translation_state << std::endl;
				// std::cout << "Rot state: " << patch_rotational_state << std::endl;
				// std::cout << full_f_sol.transpose() << std::endl;
				// std::cout << cop_point.transpose() << std::endl;
				std::cerr << "Other end penetration occurred." << std::endl;
				fForcePatchCenter = true;
				did_update_cop_pos = true; // set this to recompute the solver matrices
				patch_cop_type = COPContactType::PatchCenter;
				continue;
			}
		}

		// - check if zero moment condition is violated for patch contact,
		if(patch_cop_type == COPContactType::PatchCenter) {
			Vector2d violation_moment = full_f_sol.segment<2>(3);
			if(violation_moment.norm() > 1e-3) {
				fForcePatchCenter = false;
				// ^reset this because we want to be able to go to a different cop
				// position on the boundary if needed

				// std::cout << "Full F sol line: " << full_f_sol.segment<5>(line_start_row_id).transpose() << std::endl;
				// - - if so, compute distance to new line-COP.
				double cop_disp_dist = violation_moment.norm()/fmax(full_f_sol[2], 1e-5);

				// - - if distance to new line-COP is very small, we are done
				if(abs(cop_disp_dist) > 0.01) {// TODO: integrate bisection search and lower this threshold
					Vector2d cop_disp_dir;
					cop_disp_dir << -violation_moment(1), violation_moment(0);
					cop_disp_dir /= cop_disp_dir.norm();
					double max_dist_to_bdry = contact_patch.distanceFromBoundaryAlongRay(cop_point.segment<2>(0), cop_disp_dir);
					// std::cout << "Here " << cop_disp_dist << " max: " << max_dist_to_bdry << std::endl;
					if(abs(cop_disp_dist) < max_dist_to_bdry) {
						// we stay within the patch
						cop_point.segment<2>(0) += cop_disp_dir*cop_disp_dist;
					} else {
						// we are on the boundary
						// std::cout << "Here2" << std::endl;
						cop_point.segment<2>(0) += cop_disp_dir*max_dist_to_bdry;
						patch_cop_type = COPContactType::PatchCurvePoint;
					}
					did_update_cop_pos = true;
					// std::cout << "Full f sol " << full_f_sol.transpose() << std::endl;
					// std::cout << "Patch cop type " << (uint)(patch_cop_type) << std::endl;
					// std::cout << "New cop pos " << cop_point.transpose() << std::endl;
					continue;
				}
				// - - if distance is small, check first if we are currently doing a binary search
				// - - if so, continue binary search with point state fixed.
				// - - if not, check last COP moment. if we switch moment directions, start binary
				// - - search with point state fixed.
			}
		}

		ret_sol.force_sol = full_f_sol;
		ret_sol.result = COPSolResult::Success;
		ret_sol.local_cop_pos = cop_point;
		ret_sol.cop_type = patch_cop_type;
		return ret_sol;
	}
	std::cout << "Lin vel: " << lin_vel.transpose() << std::endl;
	std::cout << "Omega: " << omegaB.transpose() << std::endl;
	std::cout << "A mat: " << std::endl;
	std::cout << A_constraint << std::endl;
	std::cout << "RHS: " << rhs_constraint.transpose() << std::endl;
	std::cout << "Trans state: " << patch_translation_state << std::endl;
	std::cout << "Rot state: " << patch_rotational_state << std::endl;
	std::cerr << "solveSurfaceContact failed to converge" << std::endl;
	ret_sol.result = COPSolResult::NoSolution;
	return ret_sol;
}

void COPSolver::computeMatricesForSurface(const Eigen::MatrixXd& A_disp, const Eigen::VectorXd& rhs_disp, double mu, double mu_rot) {
	TA_size = 0;

	Vector2d ptslip_dir;
	double prslip_dir;
	if(patch_translation_state == FrictionState::Sliding) {
		ptslip_dir = patch_translation_slip_dir;
	} else {
		ptslip_dir = patch_translation_rolling_dir;
	}
	if(patch_rotational_state == FrictionState::Sliding) {
		prslip_dir = patch_rotation_slip_dir;
	} else {
		prslip_dir = patch_rotation_rolling_dir;
	}

	if(patch_cop_type == COPContactType::PatchCenter) {
		// both translation and rotation friction forces are limited
		if((patch_translation_state == FrictionState::Sliding ||
			patch_translation_state == FrictionState::Impending) &&
			(patch_rotational_state == FrictionState::Sliding ||
			patch_rotational_state == FrictionState::Impending)
		) {
			Trhs.segment<3>(0) = -rhs_disp.segment<3>(2);
			TA.block<3,1>(0,0) = A_disp.block<3,1>(2,2)
						- mu*A_disp.block<3,2>(2, 0)*ptslip_dir
						- mu_rot*prslip_dir*A_disp.block<3,1>(2,5);
			TA.block<3,2>(0,1) = A_disp.block<3,2>(2,3);

			TA_size = 3;
		}
		// translation forces are limited, rotation forces are not
		else if((patch_translation_state == FrictionState::Sliding ||
			patch_translation_state == FrictionState::Impending) &&
			patch_rotational_state == FrictionState::Rolling
		) {
			Trhs.segment<4>(0) = -rhs_disp.segment<4>(2);
			TA.block<3,1>(0,0) = A_disp.block<3,1>(2,2)
						- mu*A_disp.block<3,2>(2, 0)*ptslip_dir;
			TA.block<3,3>(0,1) = A_disp.block<3,3>(2,3);

			TA_size = 4;
		}
		// rotational forces are limited, translation forces are not
		else if((patch_rotational_state == FrictionState::Sliding ||
			patch_rotational_state == FrictionState::Impending) &&
			patch_translation_state == FrictionState::Rolling
		) {
			Trhs.segment<5>(0) = -rhs_disp.segment<5>(0);
			TA.block<5,2>(0,0) = A_disp.block<5,2>(0,0);
			TA.block<5,1>(0,2) = A_disp.block<5,1>(0,2)
								- mu_rot*prslip_dir*A_disp.block<5,1>(0,5);
			TA.block<5,2>(0,3) = A_disp.block<5,2>(0,3);

			TA_size = 5;
		}
		// translational and rotational forces are not limited
		else if(patch_rotational_state == FrictionState::Rolling &&
			patch_translation_state == FrictionState::Rolling
		) {
			Trhs.segment<6>(0) = -rhs_disp;
			TA.block<6,6>(0,0) = A_disp;

			TA_size = 6;
		}
	} else {
		// TODO: combine with simple 2-pt solution
		if((patch_translation_state == FrictionState::Sliding ||
			patch_translation_state == FrictionState::Impending)
		) {
			Trhs(0) = -rhs_disp(2);
			TA(0,0) = A_disp(2,2)
						- mu*A_disp.block<1,2>(2, 0)*(ptslip_dir);

			TA_size = 1;
		} else { // rolling at the line end point
			// line stuff
			Trhs.segment<3>(0) = -rhs_disp.segment<3>(0);
			TA.block<3,3>(0,0) = A_disp.block<3,3>(0,0);

			TA_size = 3;
		}
	}
}

void COPSolver::solveForSurface(double mu, double mu_rot, Eigen::VectorXd& full_f_sol) {
	// std::cout << "TA " << std::endl;
	// std::cout << TA.block(0,0,TA_size,TA_size) << std::endl;
	// std::cout << "Trhs: " << Trhs.segment(0,TA_size).transpose() << std::endl;
	if(TA_size > 1) {
		Tf_sol = TA.block(0,0,TA_size,TA_size).partialPivLu().solve(Trhs.segment(0,TA_size));
	} else {
		Tf_sol = VectorXd::Zero(1);
		Tf_sol(0) = Trhs(0)/TA(0,0);
	}

	Vector2d ptslip_dir;
	double prslip_dir;

	if(patch_translation_state == FrictionState::Sliding) {
		ptslip_dir = patch_translation_slip_dir;
	} else {
		ptslip_dir = patch_translation_rolling_dir;
	}
	if(patch_rotational_state == FrictionState::Sliding) {
		prslip_dir = patch_rotation_slip_dir;
	} else {
		prslip_dir = patch_rotation_rolling_dir;
	}

	uint tp_ind = 0;

	if(patch_cop_type == COPContactType::PatchCenter) {
		if(patch_translation_state == FrictionState::Sliding ||
			patch_translation_state == FrictionState::Impending
		) {
			full_f_sol.segment<2>(0) = -mu*Tf_sol(0)*ptslip_dir; // fx, fy
			full_f_sol.segment<3>(2) = Tf_sol.segment<3>(0); // fz, mx, my
			if(patch_rotational_state == FrictionState::Sliding ||
				patch_rotational_state == FrictionState::Impending
			) {
				full_f_sol(5) = -mu_rot*Tf_sol(0)*prslip_dir; // mz
			} else { // rotation rolling
				full_f_sol(5) = Tf_sol(3); // mz
			}
		} else { // translation rolling
			full_f_sol.segment<5>(0) = Tf_sol.segment<5>(0); // fx, fy, fz, mx, my
			if(patch_rotational_state == FrictionState::Sliding ||
				patch_rotational_state == FrictionState::Impending
			) {
				full_f_sol(5) = -mu_rot*Tf_sol(2)*prslip_dir; // mz
			} else { // rotation rolling
				full_f_sol(5) = Tf_sol(5); // mz
			}
		}
	} else {
		full_f_sol.segment<3>(3) = Vector3d::Zero(); // mx, my, mz
		if(patch_translation_state == FrictionState::Sliding ||
			patch_translation_state == FrictionState::Impending
		) {
			full_f_sol.segment<2>(0) = -mu*Tf_sol(0)*ptslip_dir; // fx, fy
			full_f_sol(2) = Tf_sol(0); // fz
		} else { // translation rolling
			full_f_sol.segment<3>(0) = Tf_sol.segment<3>(0); // fx, fy, fz
		}
	}
}

}