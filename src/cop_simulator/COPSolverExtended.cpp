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

ContactCOPSolution COPSolver::solvePtOnly(
		double friction_coeff,
		const Eigen::MatrixXd& A_constraint, // 3x3 constraint inertia matrix
		const Eigen::Vector3d& rhs_constraint,
		const Eigen::Vector3d& linear_contact_velocity
) {
	double mu = friction_coeff;
	// check for separating contact
	ContactCOPSolution ret_sol;
	if(rhs_constraint[2] > -1e-8) {
		ret_sol.force_sol.setZero(3);
		ret_sol.local_cop_pos.setZero();
		ret_sol.result = COPSolResult::Success;
		return ret_sol;
	}
	Vector2d slip_dir = linear_contact_velocity.segment<2>(0);
	if(slip_dir.norm() > 1e-6) {
		point_state = FrictionState::Sliding;
		point_slip_dir = slip_dir/slip_dir.norm();
	} else {
		point_state = FrictionState::Rolling;
	}

	if(point_state == FrictionState::Sliding) {
		double trhs_scalar = -rhs_constraint[2];
		double tA_scalar = A_constraint(2,2) 
					- mu*A_constraint.block<1,2>(2, 0).dot(point_slip_dir);
		double tsol_scalar = trhs_scalar/tA_scalar;
		Tf_sol.setZero(3);
		Tf_sol << - mu*tsol_scalar*point_slip_dir,
					tsol_scalar;
	} else {
		// slip velocity on rotation only. translation slip velocity is zero
		Tf_sol = A_constraint.ldlt().solve(-rhs_constraint);
	}
	if(Tf_sol[2] < -1e-8) { // Painleves condition
		ret_sol.result = COPSolResult::NoSolution;
		return ret_sol;
	}
	if(Tf_sol.segment<2>(0).norm() > mu*abs(Tf_sol[2]) + 1e-15) {
		point_rolling_dir = -Tf_sol.segment<2>(0);
		point_rolling_dir /= point_rolling_dir.norm();
		// finally try an impending slip sol
		// TODO: use correct impending slip direction estimate
		double trhs_scalar = -rhs_constraint[2];
		double tA_scalar = A_constraint(2,2) 
					- mu*A_constraint.block<1,2>(2, 0).dot(point_rolling_dir);
		double tsol_scalar = trhs_scalar/tA_scalar;
		Tf_sol.setZero(3);
		Tf_sol << - mu*tsol_scalar*point_rolling_dir,
					tsol_scalar;
	}
	ret_sol.force_sol = Tf_sol;
	ret_sol.local_cop_pos.setZero();
	ret_sol.result = COPSolResult::Success;
	return ret_sol;
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
		// std::cout << "Both contacts separating, optimized sol" << std::endl;
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
		// std::cout << "Pt contacts separating, optimized sol" << std::endl;
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
		// std::cout << "Line contacts separating, optimized sol" << std::endl;
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
				ContactCOPSolution ret_sol;
				ret_sol.force_sol.setZero(8);
				ret_sol.force_sol.segment<3>(pt_start_row_id) = pt_only_sol.force_sol;
				ret_sol.local_cop_pos = Vector3d::Zero();
				ret_sol.result = COPSolResult::Success;
				return ret_sol;
			}
		} else {
			// how can this fail? TODO: maybe the line+pt sol can still work?
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
				if(abs(dist_move_cop) > 0.01) {// TODO: integrate bisection search and lower this threshold
					signed_distance_from_point0 += dist_move_cop;
					signed_distance_from_point0 = fmax(signed_distance_from_point0, max_signed_dist);
					signed_distance_from_point0 = fmin(signed_distance_from_point0, min_signed_dist);
					did_update_cop_pos = true;
					//TODO: if we have fForceLineCenter set to true, consider setting it false
					continue;
				}
				// - - if distance is small, check first if we are currently doing a binary search
				// - - if so, continue binary search with point state fixed.
				// - - if not, check last line-COP moment. if we switch moment directions, start binary 
				// - - search with point state fixed.
			}
		}

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

static void getCOPSurfaceContactDisplacedMatrices(
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
	Vector3d omegaA = omega_bodyA[0];
	Vector3d omegaB = omega_bodyB[0];
	Vector3d lin_vel = linear_contact_velocity[0];

	ContactCOPSolution ret_sol;

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

	// compute rotational slip
	patch_rotation_slip_dir = omegaB(2) - omegaA(2);
	if(abs(patch_rotation_slip_dir) > 1e-8) {
		patch_rotational_state = FrictionState::Sliding;
		patch_rotation_slip_dir /= abs(patch_rotation_slip_dir);
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
	bool fForceIgnoreSlip = false;
	Eigen::VectorXd full_f_sol(6), full_a_sol(6);
	double mu_rot;
	const double mu = friction_coeff;
	Vector3d cop_point = Vector3d::Zero();
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
					if(abs(cop_disp_dist) < max_dist_to_bdry) {
						// we stay within the patch
						cop_point.segment<2>(0) += cop_disp_dir*cop_disp_dist;
					} else {
						// we are on the boundary
						cop_point.segment<2>(0) += cop_disp_dir*max_dist_to_bdry;
						patch_cop_type = COPContactType::PatchCurvePoint;
					}
					did_update_cop_pos = true;
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