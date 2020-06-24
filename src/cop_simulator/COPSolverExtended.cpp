// COPSolverExtended.cpp

#include <iostream>
#include "COPSolverExtended.h"

using namespace Eigen;

namespace Sai2COPSim {

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

static void getCOPLineContactDisplacedMatricesExtended(
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
		// if(COP_LOG_DEBUG) std::cout << "COPSolver: iters: " << iter_ind << std::endl;
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
				// if(COP_LOG_DEBUG) std::cout << "COPSolver: trans & rot sliding force sol : " << f_sol.transpose() << std::endl;

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

}