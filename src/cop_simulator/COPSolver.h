// COP solvers

#ifndef COPSOLVER_H
#define COPSOLVER_H

#include <vector>
#include <Eigen/Dense>
#include <iostream>
//Module TODO: Move all Vector/MatrixXd to static templated Vector5d, Matrix5d

#define COP_LOG_DEBUG false

enum struct COPContactType {
	Unknown,
	LineEnd,
	LineCenter
};

enum struct COPSolResult {
	None,
	Success,
	NoSolution,
	UnimplementedCase,
	FromCollisionLCP
};

struct ContactCOPSolution {
	Eigen::Vector3d local_cop_pos;
	COPContactType cop_type;
	Eigen::VectorXd force_sol;
	COPSolResult result; 

	//ctor
	ContactCOPSolution()
	:cop_type(COPContactType::Unknown), result(COPSolResult::None) {
		// nothing to do
	}
};

// obtain the resolved local_cop that can be used to instantiate the lcp solver
// in the stead contact force resolution step
inline ContactCOPSolution getLCPForInelasticCollResult(
	const Eigen::VectorXd& local_psol,
	const std::vector<Eigen::Vector3d>& local_point_list
	// const Eigen::Matrix3d& R_body_COP_frame
) {
	ContactCOPSolution ret_sol;
	ret_sol.force_sol.setZero(5);
	// ret_sol.force_sol.segment<3>(0) = local_psol.segment<3>(0) + local_psol.segment<3>(3);
	double end1_imp = fmax(0, local_psol[2]);
	double end2_imp = fmax(0, local_psol[5]);
	double cop_dist_pos0 = end2_imp/(end2_imp+end1_imp);
	if(COP_LOG_DEBUG) std::cout << "COPSolver: getLCPForInelasticCollResult: deduced cop distance from point 0: " << cop_dist_pos0 << std::endl;
	assert(cop_dist_pos0 < (1+1e15) && cop_dist_pos0 >= -1e-15);
	ret_sol.local_cop_pos = local_point_list[0]*(1.0 - cop_dist_pos0) + local_point_list[1]*(cop_dist_pos0);
	ret_sol.cop_type = (abs(abs(cop_dist_pos0 - 0.5) - 0.5) < 1e-15)? COPContactType::LineEnd : COPContactType::LineCenter;
	// if(ret_sol.cop_type == COPContactType::LineCenter) {
		// Eigen::Vector3d temp1 = R_body_COP_frame*(local_point_list[0] - ret_sol.local_cop_pos);
		// Eigen::Vector3d temp2 = R_body_COP_frame*(local_point_list[1] - ret_sol.local_cop_pos);
		// ret_sol.force_sol[4] = local_psol[1]*temp1[0] + local_psol[4]*temp2[0];
	// }
	ret_sol.result = COPSolResult::FromCollisionLCP;
	return ret_sol;
}

inline ContactCOPSolution getLCPForOnePtContactResult(
	const Eigen::VectorXd& local_psol,
	const std::vector<Eigen::Vector3d>& local_point_list,
	bool isActivePtZero
) {
	ContactCOPSolution ret_sol;
	ret_sol.force_sol.setZero(5);
	ret_sol.force_sol.segment<3>(0) = local_psol.segment<3>(0);
	// ret_sol.force_sol.segment<3>(0) = local_psol.segment<3>(0) + local_psol.segment<3>(3);
	double cop_dist_pos0 = (isActivePtZero)? 0: 1;
	if(COP_LOG_DEBUG) std::cout << "COPSolver: getLCPForOnePtContactResult: deduced cop distance from point 0: " << cop_dist_pos0 << std::endl;
	ret_sol.local_cop_pos = local_point_list[0]*(1.0 - cop_dist_pos0) + local_point_list[1]*(cop_dist_pos0);
	ret_sol.cop_type = COPContactType::LineEnd;
	ret_sol.result = COPSolResult::Success;
	return ret_sol;
}

inline Eigen::Matrix3d crossMat(const Eigen::Vector3d& r) {
	Eigen::Matrix3d ret_mat;
	ret_mat << 0, -r(2), r(1),
    r(2), 0, -r(0),
    -r(1), r(0), 0;
    return ret_mat;
}

// TODO: merge the below functions into one
inline void getCOPJ6FullFromPoint0J6Full(
	const Eigen::MatrixXd& point0_J_full, // in global cop frame
	const Eigen::Vector3d& r_point0_to_cop, // in global cop frame
	Eigen::MatrixXd& last_cop_J6
) {
	last_cop_J6.setZero();
	Eigen::Matrix3d part_cross_mat = crossMat(r_point0_to_cop);
	last_cop_J6.block(0, 0, 3, point0_J_full.cols()) = 
		point0_J_full.block(0, 0, 3, point0_J_full.cols()) - part_cross_mat*point0_J_full.block(3, 0, 3, point0_J_full.cols());
	last_cop_J6.block(3, 0, 3, point0_J_full.cols()) = point0_J_full.block(3, 0, 3, point0_J_full.cols());
}

inline void getCOPLambdaInv6FullFromPoint0LambdaInv6Full(
	const Eigen::MatrixXd& point0_lambdaInv_full,
	const Eigen::Vector3d& r_point0_to_cop,
	Eigen::MatrixXd& last_cop_lambdainv_full
) {
	if(COP_LOG_DEBUG) std::cout << "COPsolver: getCOPLambdaInv6FullFromPoint0LambdaInv6Full: " << std::endl;
	if(COP_LOG_DEBUG) std::cout << "r_point0_to_cop " << r_point0_to_cop.transpose() << std::endl;
	last_cop_lambdainv_full.setZero();
	Eigen::Matrix3d part_cross_mat = crossMat(r_point0_to_cop);
	Eigen::MatrixXd cross_mat = Eigen::MatrixXd::Identity(6,6);
	cross_mat.block<3,3>(0, 3) = -part_cross_mat;
	last_cop_lambdainv_full = cross_mat * point0_lambdaInv_full * cross_mat.transpose();
	if(COP_LOG_DEBUG) std::cout << "point0_lambdaInv_full " << point0_lambdaInv_full << std::endl;
	if(COP_LOG_DEBUG) std::cout << "last_cop_lambdainv_full " << last_cop_lambdainv_full << std::endl;
}

inline void getCOPRhsNonlinFullFromPoint0RhsNonlinFull(
	const Eigen::VectorXd& point0_rhs_nonlin_6,
	const Eigen::Vector3d& r_point0_to_cop,
	Eigen::VectorXd& last_cop_rhsnonlin_full
) {
	last_cop_rhsnonlin_full.setZero();
	Eigen::Matrix3d part_cross_mat = crossMat(r_point0_to_cop);
	Eigen::MatrixXd cross_mat = Eigen::MatrixXd::Identity(6,6);
	cross_mat.block<3,3>(0, 3) = -part_cross_mat;
	last_cop_rhsnonlin_full = cross_mat*point0_rhs_nonlin_6;
}


// // get contact Jacobian at given local point
// void getCOPLineContactJacobian(
// 	Eigen::MatrixXd& cop_jacobian,
// 	const Eigen::MatrixXd& ref_jacobian,
// 	double dist_ref_to_cop_pos // measured along the line of contact
// ) {
// 	Eigen::MatrixXd cross_mat(3,2);
// 	cross_mat << 0, 	0,
// 				0, -dist_ref_to_cop_pos,
// 				dist_ref_to_cop_pos, 0;
// 	 = crossMat(cop_point - ref_point);
// 	int dof = ref_jacobian.cols();
// 	cop_jacobian.block(0, 0, 3, dof) = ref_jacobian.block(0, 0, 3, dof) 
// 										- cross_mat * ref_jacobian.block(3, 0, 2, dof);
// 	cop_jacobian.block(3, 0, 2, dof) = ref_jacobian.block(3, 0, 2, dof);
// }

inline void getCOPLineContactDisplacedMatricesFromFull(
	Eigen::MatrixXd& A,
	Eigen::VectorXd& rhs_nonlin,
	Eigen::VectorXd& Jdot_qdot,
	Eigen::VectorXd& pre_vel,
	const Eigen::Vector3d& last_cop_to_cop_proj,
	const Eigen::MatrixXd& A_full,
	const Eigen::VectorXd& rhs_nonlin_full,
	const Eigen::Vector3d& omega,
	const Eigen::Vector3d& r_last_cop_proj,
	const Eigen::VectorXd& pre_vel_full)
{
	Eigen::MatrixXd cross_mat(5,6);
	Eigen::Matrix3d part_cross_mat = crossMat(last_cop_to_cop_proj);
	if(COP_LOG_DEBUG) std::cout << "COPsolver: getCOPLineContactDisplacedMatricesFromFull: last_cop_to_cop_proj: " << last_cop_to_cop_proj.transpose() << std::endl;
	cross_mat.setZero();
	cross_mat.block<3,3>(0,0) = Eigen::Matrix3d::Identity();
	cross_mat.block<3,3>(0,3) = -part_cross_mat;
	cross_mat.block<2,2>(3,4) = Eigen::Matrix2d::Identity();
	A = cross_mat*A_full*cross_mat.transpose();
	rhs_nonlin = cross_mat*rhs_nonlin_full;
	Jdot_qdot.setZero(5);
	Jdot_qdot.segment<3>(0) = omega.cross(omega.cross(r_last_cop_proj));
	pre_vel = cross_mat*pre_vel_full;
}

inline void getCOPLineContactDisplacedMatrices(
	Eigen::MatrixXd& A_disp,
	Eigen::VectorXd& rhs_nonlin_disp,
	Eigen::VectorXd& Jdot_qdot_disp,
	Eigen::VectorXd& pre_vel_disp,
	double signed_distance_last_point,
	const Eigen::MatrixXd& A,
	const Eigen::VectorXd& rhs_nonlin,
	const Eigen::Vector3d& omega,
	const Eigen::Vector3d& r_last_cop_proj,
	const Eigen::VectorXd& pre_vel)
{
	if(COP_LOG_DEBUG) std::cout << "COP solver: getCOPLineContactDisplacedMatrices: signed_distance_last_point: " << signed_distance_last_point << std::endl;
	Eigen::MatrixXd cross_mat(5,5);
	cross_mat.setZero();
	cross_mat.block<3,3>(0,0) = Eigen::Matrix3d::Identity();
	cross_mat.block<3,2>(0,3) << 0,	0,
				0, signed_distance_last_point,
				-signed_distance_last_point, 0;
	cross_mat.block<2,2>(3,3) = Eigen::Matrix2d::Identity();
	A_disp = cross_mat*A*cross_mat.transpose();
	rhs_nonlin_disp = cross_mat*rhs_nonlin;
	Jdot_qdot_disp.setZero(5);
	Jdot_qdot_disp.segment<3>(0) = omega.cross(omega.cross(r_last_cop_proj + signed_distance_last_point*Eigen::Vector3d(1.0,0.0,0.0)));
	pre_vel_disp = pre_vel;
	pre_vel_disp[1] += signed_distance_last_point*pre_vel[4];
	pre_vel_disp[2] += -signed_distance_last_point*pre_vel[3];
}

// compute mu_rot from mu and location of cop point
inline double getMuRotation(double mu, double dist_end1, double dist_end2) {
	assert(dist_end1 > -1e-15 && dist_end2 > -1e-15 && (dist_end1 + dist_end2) > 1e-15 );
	return mu * 2 * dist_end1*dist_end2/(dist_end1 + dist_end2);
}

enum struct FeasibleCOPLineResult {
	FRSuccess = 0,
	FRNegativeNormalForce,
	FROtherEndPenetration,
	FRTranslationFrictionConeViolation,
	FRRotationFrictionViolation,
	FRLineEndNonZeroRotationFriction,
	FRZeroMomentViolation,
	FRPenetrationAcceleration,
	FRPenetrationRotationAcceleration,
	FRTranslationFrictionNonDissipative,
	FRRotationFrictionNonDissipative
};

inline FeasibleCOPLineResult testFeasibleCOPLine(
	const Eigen::VectorXd& force_sol,
	const Eigen::VectorXd& acc_sol,
	const Eigen::VectorXd& pre_vel,
	COPContactType contact_type,
	double mu,
	double mu_rot,
	double signed_dist_opp_end
) {
	if(force_sol[2] < -1e-8) {
		return FeasibleCOPLineResult::FRNegativeNormalForce;
	}
	if(contact_type == COPContactType::LineCenter && abs(force_sol[3]) > 1e-3) { //TODO: this is too low!!
		return FeasibleCOPLineResult::FRZeroMomentViolation;
	}
	if(contact_type == COPContactType::LineEnd) {
		if(acc_sol[2] - signed_dist_opp_end*acc_sol[3] < -1e-15) {
			return FeasibleCOPLineResult::FROtherEndPenetration;
		}
	}
	if(force_sol.segment<2>(0).norm() > mu*force_sol[2] + 1e-15) {
		return FeasibleCOPLineResult::FRTranslationFrictionConeViolation;
	}
	if(contact_type == COPContactType::LineCenter && abs(force_sol[4]) > mu_rot*force_sol[2] + 1e-15) {
		return FeasibleCOPLineResult::FRRotationFrictionViolation;
	}
	if(contact_type == COPContactType::LineEnd && abs(force_sol[4]) > 1e-15) {
		return FeasibleCOPLineResult::FRLineEndNonZeroRotationFriction;
	}
	if(acc_sol[2] < -1e-15) {
		return FeasibleCOPLineResult::FRPenetrationAcceleration;
	}
	// if(contact_type == COPContactType::LineCenter && acc_sol[2] < 1e-15 && abs(acc_sol[3]) > 1e-15) {
	// 	return FeasibleCOPLineResult::FRPenetrationRotationAcceleration;
	// }
	if(pre_vel.segment<2>(0).dot(force_sol.segment<2>(0)) > 1e-8) {
		if(COP_LOG_DEBUG) std::cout << "Translational speed: " << pre_vel.segment<2>(0).transpose() << std::endl;
		if(COP_LOG_DEBUG) std::cout << "Tangent force: " << force_sol.segment<2>(0).transpose() << std::endl;
		return FeasibleCOPLineResult::FRTranslationFrictionNonDissipative;
	}
	if(contact_type == COPContactType::LineCenter && (pre_vel[4]*force_sol[4]) > 1e-8) {
		if(COP_LOG_DEBUG) std::cout << "Rotational speed: " << pre_vel[4] << std::endl;
		if(COP_LOG_DEBUG) std::cout << "Moment: " << force_sol[4] << std::endl;
		return FeasibleCOPLineResult::FRRotationFrictionNonDissipative;
	}
	return FeasibleCOPLineResult::FRSuccess;
}

// All quantities are assumed to be defined in the GLOBAL COP FRAME
// which is x = line direction
//          y = z cross x
//          z = normal to surface
// This INCLUDES last_COP_sol.force_sol
//
// EXCEPTION: local_point_list, and last_COP_sol.local_cop_pos, which 
// are defined in the BODY FRAME
inline ContactCOPSolution resolveCOPLineContactWithLastCOPSol(
	const Eigen::MatrixXd& A_full, // 6 x 6 defined at last_COP_sol local point
	const Eigen::VectorXd& rhs_nonlin_full, // 6 defined at last_COP_sol local point
	const Eigen::Vector3d& omega, // defined in global COP frame
	const Eigen::Vector3d& r_last_cop, // position vector from center of object to last_COP_sol local point, defined in global COP frame
	const Eigen::VectorXd& pre_vel_full, // full 6 dof at last_COP_sol local point
	const Eigen::Matrix3d& R_body_COP_frame, // rotation from body frame to new global COP frame
	std::vector<Eigen::Vector3d>& local_point_list, // points are in the local body frame, NOT in the global COP frame
	double mu,
	ContactCOPSolution& last_COP_sol
) {
	if(COP_LOG_DEBUG) {
		std::cout << "COPSolver: resolveCOPLineContactWithLastCOPSol: inputs: " << std::endl;
		std::cout << "A_full: " << A_full << std::endl;
		std::cout << "rhs_nonlin_full: " << rhs_nonlin_full.transpose() << std::endl;
		std::cout << "omega: " << omega.transpose() << std::endl;
		std::cout << "r_last_cop: " << r_last_cop.transpose() << std::endl;
		std::cout << "pre_vel_full: " << pre_vel_full.transpose() << std::endl;
		std::cout << "R_body_COP_frame: " << R_body_COP_frame << std::endl;		
	}
	
	ContactCOPSolution ret_sol;
	if(local_point_list.size() != 2) {
		ret_sol.result = COPSolResult::UnimplementedCase;
		return ret_sol;
	}

	Eigen::Vector3d local_sol_point, projected_last_cop_local_point;
	Eigen::Vector3d last_local_sol_point = last_COP_sol.local_cop_pos;
	
	Eigen::VectorXd a_sol(5);
	double mu_rot = 0.0;
	// project last_COP_sol.local_cop_pos to new contact line
	Eigen::Vector3d point0_to_point1_dir = local_point_list[1] - local_point_list[0];
	point0_to_point1_dir /= point0_to_point1_dir.norm();
	double projected_component = (last_local_sol_point - local_point_list[0]).dot(point0_to_point1_dir);
	if(projected_component < 0) {
		ret_sol.result = COPSolResult::UnimplementedCase;
		return ret_sol;
	}
	projected_last_cop_local_point = local_point_list[0] + point0_to_point1_dir*projected_component;
	if(COP_LOG_DEBUG) std::cout << "COP Solver: distance projected last cop to contact line: " << sqrt(pow((last_local_sol_point - local_point_list[0]).norm(), 2) - projected_component*projected_component) << std::endl;

	// compute all matrices in the global COP frame, displaced to the
	// projected COP point
	Eigen::MatrixXd A(5,5);
	Eigen::VectorXd rhs_nonlin(5);
	Eigen::VectorXd Jdot_qdot(5);
	Eigen::VectorXd pre_vel(5);
	Eigen::Vector3d r_last_cop_proj = r_last_cop + R_body_COP_frame*(projected_last_cop_local_point - last_local_sol_point);
	getCOPLineContactDisplacedMatricesFromFull(A, rhs_nonlin, Jdot_qdot, pre_vel,
			R_body_COP_frame*(projected_last_cop_local_point - last_local_sol_point), A_full, rhs_nonlin_full, omega, r_last_cop_proj, pre_vel_full);
	double max_signed_distance_last_point = 0.0; // distance along global COP frame X axis from projected last cop point
	double min_signed_distance_last_point = 0.0;
	if((R_body_COP_frame*point0_to_point1_dir).dot(Eigen::Vector3d(1,0,0)) > 0) { //this should be +1 or -1
		max_signed_distance_last_point = (local_point_list[1] - projected_last_cop_local_point).dot(point0_to_point1_dir);
		min_signed_distance_last_point = (local_point_list[0] - projected_last_cop_local_point).dot(point0_to_point1_dir);
	} else {
		max_signed_distance_last_point = (local_point_list[0] - projected_last_cop_local_point).dot(-point0_to_point1_dir);
		min_signed_distance_last_point = (local_point_list[1] - projected_last_cop_local_point).dot(-point0_to_point1_dir);
	}
	if(last_COP_sol.result == COPSolResult::Success) {
		// test if last_COP_sol is still valid 
		a_sol = A*last_COP_sol.force_sol + rhs_nonlin + Jdot_qdot;
		// compute max mu_rotation given position on contact line
		mu_rot = getMuRotation(mu,
				(projected_last_cop_local_point - local_point_list[0]).norm(),
				(projected_last_cop_local_point - local_point_list[1]).norm());
		double dist_to_other_end = 0;
		if(last_COP_sol.cop_type == COPContactType::LineEnd) {
			if(max_signed_distance_last_point > 1e-15) {
				dist_to_other_end = max_signed_distance_last_point;
			} else {
				dist_to_other_end = min_signed_distance_last_point;
			}
		}
		FeasibleCOPLineResult test_res = testFeasibleCOPLine(last_COP_sol.force_sol, a_sol, pre_vel, last_COP_sol.cop_type, mu, mu_rot, dist_to_other_end);
		if(test_res == FeasibleCOPLineResult::FRSuccess) {
			ret_sol.local_cop_pos = projected_last_cop_local_point;
			ret_sol.force_sol = last_COP_sol.force_sol;
			ret_sol.cop_type = last_COP_sol.cop_type;
			ret_sol.result = COPSolResult::Success;
			return ret_sol;
		}
		if(COP_LOG_DEBUG) {
			std::cout << "Mu rotation: " << mu << std::endl;
			std::cout << "Last cop sol force: " << last_COP_sol.force_sol.transpose() << std::endl;
			std::cout << "cop accel: " << a_sol.transpose() << std::endl;
			std::cout << "check fail result: " << static_cast<int>(test_res) <<std::endl;
		}
	}

	// if projected point is no longer a valid solution, try to search for one
	const int max_iters = 20;
	int iter_ind = 0;
	Eigen::MatrixXd A_disp(5, 5);
	Eigen::VectorXd rhs_nonlin_disp(5);
	Eigen::VectorXd Jdot_qdot_disp(5), pre_vel_disp(5);
	Eigen::VectorXd f_sol(5), f_sol_full_roll(5);
	double signed_distance_last_point = 0.0;
	FeasibleCOPLineResult it_test_res;
	bool fForceLineCenter = false;
	double last_signed_distance = 0.0;
	double last_signed_distance_weight = 0.0;
	bool fForceIgnoreSlip = false;
	bool fBisectionSearch = false;
	double bisection_dist_bound_upper = 0.0;
	double bisection_dist_bound_lower = 0.0;
	double last_zero_moment_error = 0.0;
	while (iter_ind < max_iters) {
		if(COP_LOG_DEBUG) std::cout << "COPSolver: iters: " << iter_ind << std::endl;
		iter_ind++;
		local_sol_point = projected_last_cop_local_point + R_body_COP_frame.transpose()*Eigen::Vector3d(signed_distance_last_point, 0, 0);
		double dist1 = (local_sol_point - local_point_list[0]).norm();
		double dist2 = (local_sol_point - local_point_list[1]).norm();
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
		if(COP_LOG_DEBUG) std::cout << "COPSolver: contact type: " << ((ret_sol.cop_type == COPContactType::LineCenter)? "Center": "End") << std::endl;
		if(COP_LOG_DEBUG) std::cout << "COPSolver: distance from last pt: " << signed_distance_last_point << std::endl;
		// compute contact matrices
		getCOPLineContactDisplacedMatrices(A_disp, rhs_nonlin_disp, Jdot_qdot_disp, pre_vel_disp,
			signed_distance_last_point, A, rhs_nonlin, omega, r_last_cop_proj, pre_vel);
		if(COP_LOG_DEBUG) {
			std::cout << "COPSolver: displaced matrices: " << std::endl;
			std::cout << "A: " << A_disp << std::endl;
			std::cout << "rhs_nonlin: " << rhs_nonlin_disp.transpose() << std::endl;
			std::cout << "Jdot_qdot: " << Jdot_qdot_disp.transpose() << std::endl;
			std::cout << "pre_vel: " << pre_vel_disp.transpose() << std::endl;
		}
		// compute mu_rotation
		mu_rot = getMuRotation(mu,
			(local_sol_point - local_point_list[0]).norm(),
			(local_sol_point - local_point_list[1]).norm());
		double rotation_sign = (signbit(pre_vel_disp[4]))? -1: 1;
		if(COP_LOG_DEBUG) std::cout << "COPSolver: computed mu rot: " << mu_rot << std::endl;
		// ------------- test for rolling -----------------
		f_sol_full_roll = A_disp.ldlt().solve(-rhs_nonlin_disp - Jdot_qdot_disp);
		if(COP_LOG_DEBUG) std::cout << "COPSolver: computed full roll sol: " << f_sol_full_roll.transpose() << std::endl;
		a_sol.setZero();
		Eigen::MatrixXd tA;
		Eigen::VectorXd tsol, trhs;
		Eigen::Vector2d slip_vel = pre_vel_disp.segment<2>(0);
		if(COP_LOG_DEBUG) std::cout << "COPSolver: translation slip: " << slip_vel.transpose() << " norm: " << slip_vel.norm() << std::endl;
		if(COP_LOG_DEBUG) std::cout << "COPSolver: rotation slip: " << pre_vel_disp[4] << std::endl;
		if(ret_sol.cop_type == COPContactType::LineCenter) {
			if(!fForceIgnoreSlip && slip_vel.norm() > 1e-6 && abs(pre_vel_disp[4]) > 1e-8) {
				// slip velocity on both translation and rotation
				trhs.setZero(2);
				trhs[0] = -rhs_nonlin_disp[2] - Jdot_qdot_disp[2];
				trhs[1] = -rhs_nonlin_disp[3] - Jdot_qdot_disp[3];
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

			} else if (!fForceIgnoreSlip && slip_vel.norm() > 1e-6 && abs(pre_vel_disp[4]) <= 1e-8) {
				// slip velocity on translation only. rotation velocity is zero
				trhs.setZero(3);
				trhs = -rhs_nonlin_disp.segment<3>(2) - Jdot_qdot_disp.segment<3>(2);
				tA.setZero(3,3);
				tA.block<3,1>(0,0) = A_disp.block<3,1>(2,2) - mu*A_disp.block<3,2>(2, 0)*slip_vel/slip_vel.norm();
				tA.block<3,2>(0,1) = A_disp.block<3,2>(2,3);
				tsol = tA.partialPivLu().solve(trhs);
				f_sol << - mu*tsol[0]*slip_vel/slip_vel.norm(),
							tsol;
				if(COP_LOG_DEBUG) std::cout << "COPSolver: trans slide & rot rolling force sol : " << f_sol.transpose() << std::endl;
			} else if ((fForceIgnoreSlip || slip_vel.norm() <= 1e-6) && abs(pre_vel_disp[4]) > 1e-8) {
				// slip velocity on rotation only. translation slip velocity is zero
				fForceIgnoreSlip = true; // TODO: think more about this
				// TODO: Instead of this approach, we can consider the slip direction at the next time step.
				// i.e. slip_vel = current_slip_vel + tangential_acc*dt
				// Have to think about whether this can introduce an increase in energy, like in 
				// Kane's example for frictional collision
				trhs.setZero(4);
				trhs = -rhs_nonlin_disp.segment<4>(0) - Jdot_qdot_disp.segment<4>(0);
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
				double trhs_scalar = -rhs_nonlin_disp[2] - Jdot_qdot_disp[2];
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
				trhs = -rhs_nonlin_disp.segment<3>(0) - Jdot_qdot_disp.segment<3>(0);
				tA.setZero(3,3);
				tA = A_disp.block<3,3>(0,0);
				tsol = tA.ldlt().solve(trhs);
				f_sol << tsol, 0, 0;
			}
		}
		a_sol = A_disp*f_sol + rhs_nonlin_disp + Jdot_qdot_disp;
		double dist_to_other_end = 0;
		if(ret_sol.cop_type == COPContactType::LineEnd) {
			if(max_signed_distance_last_point - signed_distance_last_point > 1e-15) {
				dist_to_other_end = max_signed_distance_last_point - signed_distance_last_point;
			} else {
				dist_to_other_end = min_signed_distance_last_point - signed_distance_last_point;
			}
		}
		it_test_res = testFeasibleCOPLine(f_sol, a_sol, pre_vel_disp, ret_sol.cop_type, mu, mu_rot, dist_to_other_end);
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
			if(rhs_nonlin_disp[2] + Jdot_qdot_disp[2] > -1e-15) {
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
				last_signed_distance = signed_distance_last_point;
				if(moment_error > 0) {
					bisection_dist_bound_upper = signed_distance_last_point;
				} else {
					bisection_dist_bound_lower = signed_distance_last_point;
				}
				signed_distance_last_point = 0.5*bisection_dist_bound_lower + 0.5*bisection_dist_bound_upper;
				if(COP_LOG_DEBUG) std::cout << "Update signed distance last point: " << signed_distance_last_point << std::endl;
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
					bisection_dist_bound_upper = signed_distance_last_point;
				} else {
					bisection_dist_bound_lower = signed_distance_last_point;
					bisection_dist_bound_upper = last_signed_distance;
				}
				last_signed_distance = signed_distance_last_point;
				signed_distance_last_point = 0.5*bisection_dist_bound_lower + 0.5*bisection_dist_bound_upper;
				if(COP_LOG_DEBUG) std::cout << "Last signed distance: " << last_signed_distance << std::endl;
				if(COP_LOG_DEBUG) std::cout << "Set new signed distance last point: " << signed_distance_last_point << std::endl;
			} else {
				double zero_moment_dist = -f_sol[3]/fmax(f_sol[2], 1e-5);
				last_signed_distance = signed_distance_last_point;
				signed_distance_last_point += zero_moment_dist;
				fBisectionSearch = false;
			}
			last_zero_moment_error = moment_error;
			// signed_distance_last_point = signed_distance_last_point*(1.0 - last_signed_distance_weight) + last_signed_distance*last_signed_distance_weight;
			// if(last_signed_distance_weight < 1e-15) {
			// 	last_signed_distance_weight = 0.5;
			// }
			if(signed_distance_last_point > 0) {
				signed_distance_last_point = fmin(signed_distance_last_point, max_signed_distance_last_point);
			} else {
				signed_distance_last_point = fmax(signed_distance_last_point, min_signed_distance_last_point);
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
				if(abs(pre_vel_disp[4]) > 1e-8) {
					// slip velocity on both translation and rotation
					trhs.setZero(2);
					trhs[0] = -rhs_nonlin_disp[2] - Jdot_qdot_disp[2];
					trhs[1] = -rhs_nonlin_disp[3] - Jdot_qdot_disp[3];
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
					trhs = -rhs_nonlin_disp.segment<3>(2) - Jdot_qdot_disp.segment<3>(2);
					tA.setZero(3,3);
					tA.block<3,1>(0,0) = A_disp.block<3,1>(2,2) - mu*A_disp.block<3,2>(2, 0)*roll_slip_direction;
					tA.block<3,2>(0,1) = A_disp.block<3,2>(2,3);
					tsol = tA.partialPivLu().solve(trhs);
					f_sol << - mu*tsol[0]*roll_slip_direction,
								tsol;
				}
			} else {
				// slip velocity on both translation and rotation
				double trhs_scalar = -rhs_nonlin_disp[2] - Jdot_qdot_disp[2];
				double tA_scalar = A_disp(2,2) 
							- mu*(A_disp.block<1,2>(2, 0).dot(roll_slip_direction));
				double tsol_scalar = trhs_scalar/tA_scalar;
				f_sol << - mu*tsol_scalar*roll_slip_direction,
							tsol_scalar,
							0,
							0;
			}
			a_sol = A_disp*f_sol + rhs_nonlin_disp + Jdot_qdot_disp;
			it_test_res = testFeasibleCOPLine(f_sol, a_sol, pre_vel_disp, ret_sol.cop_type, mu, mu_rot, dist_to_other_end);
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
					trhs[0] = -rhs_nonlin_disp[2] - Jdot_qdot_disp[2];
					trhs[1] = -rhs_nonlin_disp[3] - Jdot_qdot_disp[3];
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
					trhs = -rhs_nonlin_disp.segment<4>(0) - Jdot_qdot_disp.segment<4>(0);
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
			a_sol = A_disp*f_sol + rhs_nonlin_disp + Jdot_qdot_disp;
			it_test_res = testFeasibleCOPLine(f_sol, a_sol, pre_vel_disp, ret_sol.cop_type, mu, mu_rot, dist_to_other_end);
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
					trhs[0] = -rhs_nonlin_disp[2] - Jdot_qdot_disp[2];
					trhs[1] = -rhs_nonlin_disp[3] - Jdot_qdot_disp[3];
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
					a_sol = A_disp*f_sol + rhs_nonlin_disp + Jdot_qdot_disp;
					it_test_res = testFeasibleCOPLine(f_sol, a_sol, pre_vel_disp, ret_sol.cop_type, mu, mu_rot, dist_to_other_end);
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

#endif //COPSOLVER_H