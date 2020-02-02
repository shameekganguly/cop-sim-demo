// LCP solver
// Note that only spatial solvers are included
#ifndef LCPSOLVER_H
#define LCPSOLVER_H

#include <vector>
#include <Eigen/Dense>
#include <iostream>

// Collision LCP points solver solves a set of equations of the form
// v = Ap + b 					where A is assumed positive definite
// v_n >= -eps * pre_v_n  \perp  p_n >= 0 	where n is the normal direction. Assumed to
// 								be the [2] component of each v_i 3-velocity vector and 
//								p_i 3-force vector. pre_v is the velocity before coll.
//								eps is the coefficient of restitution.
// |v_t| >= 0  \perp  p_n * \mu >= |p_t| 	where \mu is the coefficient of friction
//											p_t is the tangential force component
// Note that the above does not yield a valid solution in the case that v_t ≠ 0
// In this case, a further resolution of the direction of p_t is required.
// We assume that in this case, p_t direction is opposite to initial slip direction
// TODO: the above assumption is only correct when the slip direction does not change
// Otherwise, the classic counter example is the Kane's pendulum problem where it is
// clear that the above assumption leads to an increase in energy for elastic collision.
// NOTE 2: pre_v_n must be ≤ 0 for all contacts. Otherwise an energetically correct
// solution is not guaranteed. If there are contacts where pre_v_n > 0, these contacts
// must not be included in the collision solution. They must be handled through
// sequential collision handling.

enum LCPSolResult {Success = 0, NoSolution, UnimplementedCase};

struct CollLCPPointSolution {
	LCPSolResult result;
	Eigen::VectorXd p_sol;
	// std::vector<uint> no_contact_list; //internal, for inspection only

	// ctor
	CollLCPPointSolution(LCPSolResult a_result): result(a_result) {
		// nothing to do
	}
};

enum FeasibleCollLCPPointResult {
	FRSuccess = 0,
	FRNegativeNormalImpulse,
	FRFrictionConeViolation,
	FRNonActivePenetration
};

FeasibleCollLCPPointResult testFeasibleCollLCPPoint(uint Nactive,
	uint Nnonactive,
	const Eigen::VectorXd& pre_v_nonactive,
	const Eigen::VectorXd& psol_active,
	const Eigen::VectorXd& vsol_nonactive,
	double epsilon,
	double mu
) {
	// check for non-negative contact forces and friction cone constraint.
	for(uint i = 0; i < Nactive; i++) {
		if(psol_active(i*3 + 2) < 0) {
			return FeasibleCollLCPPointResult::FRNegativeNormalImpulse;
		}
		if(psol_active.segment(i*3, 2).norm() > mu * psol_active(i*3 + 2)) {
			return FeasibleCollLCPPointResult::FRFrictionConeViolation;
		}
	} //TODO: mem optimization by putting all normal indices first, then all tangent indices
	if(Nnonactive) {
		// check for non-negative contact velocity.
		for(uint i = 0; i < Nnonactive; i++) {
			if(vsol_nonactive(i*3 + 2) < 0) {
				return FeasibleCollLCPPointResult::FRNonActivePenetration;
			}
		}
	}
	return FeasibleCollLCPPointResult::FRSuccess;
}

// NOTE: the method below assumes that the 1st Cartesian coordinate is aligned with the
// line joining the two contact points. For two contact points, there is a redundancy
// in this direction because equal and opposite forces will cancel each other out.

CollLCPPointSolution solveCollLCPPoint (uint Npoints, 
	const Eigen::MatrixXd& A,
	const Eigen::VectorXd& b,
	const Eigen::VectorXd& pre_v,
	double epsilon,
	double mu
) {
	// TODO: extend to beyond 2 points
	if(Npoints != 2) {
		return CollLCPPointSolution(LCPSolResult::UnimplementedCase);
	}
	Eigen::VectorXd psol(Npoints*3);
	Eigen::VectorXd vsol(Npoints*3);
	Eigen::Vector2d slip1_direction, slip2_direction;
	Eigen::Matrix2d A_slide_slide;
	Eigen::Vector2d psol_slide_slide, lhs_slide_slide;
	Eigen::Matrix4d A_slide_roll;
	Eigen::Vector4d psol_slide_roll, lhs_slide_roll;
	bool sliding_with_no_preslip = false;

	// --------- case1: contact 1 only ---------
	uint cind = 0;
	uint nind = 1;
	// --------- case 1a: rolling contact ---------
	psol.setZero();
	vsol.setZero();
	vsol(cind*3 + 2) = -epsilon * pre_v(cind*3 + 2);
	//TODO: ^ use composed A_contact for active contacts
	psol.segment(cind*3, 3) = A.block(cind*3, cind*3, 3, 3).ldlt().solve(
		vsol.segment(cind*3, 3) - b.segment(cind*3, 3)
	);
	// compute non active velocity
	vsol.segment(nind*3, 3) = A.block(nind*3, cind*3, 3, 3) * psol.segment(cind*3, 3) + b.segment(nind*3, 3);
	FeasibleCollLCPPointResult ret_case1a = testFeasibleCollLCPPoint(1,
		1,
		pre_v.segment(nind*3,3),
		psol.segment(cind*3, 3),
		vsol.segment(nind*3, 3),
		epsilon,
		mu
	);
	if(ret_case1a == FeasibleCollLCPPointResult::FRSuccess) {
		std::cout << "LCP Success: Contact 1 roll, No contact 2" << std::endl;
		CollLCPPointSolution ret_lcp_sol(LCPSolResult::Success);
		ret_lcp_sol.p_sol = psol;
		return ret_lcp_sol;
	}
	// --------- case 1b: sliding contact ---------
	vsol.setZero();
	slip1_direction = pre_v.segment(cind*3, 2);
	if(slip1_direction.norm() < 1e-5) {
		sliding_with_no_preslip = true;
	} else {
		slip1_direction /= slip1_direction.norm();
		psol.segment(cind*3, 3) << -mu*slip1_direction, 1.0;
		double normal_impulse = 1/(A.block(cind*3, cind*3, 3, 3).row(2).dot(psol.segment(cind*3, 3)))
							* (-epsilon * pre_v(cind*3 + 2) - b(cind*3 + 2));
		psol.segment(cind*3, 3) *= normal_impulse;
		// compute non active velocity
		vsol.segment(nind*3, 3) = A.block(nind*3, cind*3, 3, 3) * psol.segment(cind*3, 3) + b.segment(nind*3, 3);
		FeasibleCollLCPPointResult ret_case1b = testFeasibleCollLCPPoint(1,
			1,
			pre_v.segment(nind*3,3),
			psol.segment(cind*3, 3),
			vsol.segment(nind*3, 3),
			epsilon,
			mu
		);
		if(ret_case1b == FeasibleCollLCPPointResult::FRSuccess) {
			std::cout << "LCP Success: Contact 1 slide, No contact 2" << std::endl;
			CollLCPPointSolution ret_lcp_sol(LCPSolResult::Success);
			ret_lcp_sol.p_sol = psol;
			return ret_lcp_sol;
		}
	}
	// --------- case2: both contacts ---------
	// --------- case 2a: both contacts rolling
	psol.setZero();
	vsol.setZero();
	vsol(0*3 + 2) = -epsilon * pre_v(0*3 + 2);
	vsol(1*3 + 2) = -epsilon * pre_v(1*3 + 2);
	psol = A.ldlt().solve(vsol - b);
	// The solution above is one of many possible because A is positive semi-definite.
	// There is a line of force redundancy in the [0] coordinate in the tangent plane.
	FeasibleCollLCPPointResult ret_case2a = testFeasibleCollLCPPoint(2,
		0,
		pre_v, //unused
		psol,
		vsol, //unused
		epsilon,
		mu
	);
	if(ret_case2a == FeasibleCollLCPPointResult::FRSuccess) {
		std::cout << "LCP Success: Contact 1 roll, Contact 2 roll" << std::endl;
		CollLCPPointSolution ret_lcp_sol(LCPSolResult::Success);
		ret_lcp_sol.p_sol = psol;
		return ret_lcp_sol;
	} else if (ret_case2a == FeasibleCollLCPPointResult::FRFrictionConeViolation) {
		// check if there is a solution possible due to redundancy
		double violation1 = mu*mu*psol(2)*psol(2) - psol(1)*psol(1);
		double violation2 = mu*mu*psol(5)*psol(5) - psol(4)*psol(4);
		if(violation1 >= 0 && violation2 >= 0) {
			double ximpulse_diff1 = abs(psol(0)) - sqrt(violation1);
			double ximpulse_diff2 = abs(psol(3)) - sqrt(violation2);
			bool redundant_feasible = false;
			if(ximpulse_diff1 > 0 && psol(0) > 0 && abs(psol(3) + ximpulse_diff1) <= sqrt(violation2)) {
				redundant_feasible = true;
				psol(0) -= ximpulse_diff1;
				psol(3) += ximpulse_diff1;
			} else if (ximpulse_diff1 > 0 && psol(0) < 0 && abs(psol(3) - ximpulse_diff1) <= sqrt(violation2)) {
				redundant_feasible = true;
				psol(0) += ximpulse_diff1;
				psol(3) -= ximpulse_diff1;
			} else if (ximpulse_diff2 > 0 && psol(3) > 0 && abs(psol(0) + ximpulse_diff2) <= sqrt(violation1)) { // untested code path
				redundant_feasible = true;
				psol(0) += ximpulse_diff2;
				psol(3) -= ximpulse_diff2;
			} else if (ximpulse_diff2 > 0 && psol(3) < 0 && abs(psol(0) - ximpulse_diff2) <= sqrt(violation1)) { // untested code path
				redundant_feasible = true;
				psol(0) -= ximpulse_diff2;
				psol(3) += ximpulse_diff2;
			}
			if (redundant_feasible) {
				std::cout << "LCP Success: Contact 1 roll, Contact 2 roll with redundancy adjustment" << std::endl;
				CollLCPPointSolution ret_lcp_sol(LCPSolResult::Success);
				ret_lcp_sol.p_sol = psol;
				return ret_lcp_sol;
			}
		}
	}
	// --------- case 2b: contact 1 rolling, contact 2 sliding
	uint c2brind = 0;
	uint c2bsind = 1;
	psol.setZero();
	vsol.setZero();
	vsol(0*3 + 2) = -epsilon * pre_v(0*3 + 2);
	vsol(1*3 + 2) = -epsilon * pre_v(1*3 + 2);
	slip1_direction = pre_v.segment(c2bsind*3, 2);
	if(slip1_direction.norm() < 1e-5) {
		sliding_with_no_preslip = true;
	} else {
		slip1_direction /= slip1_direction.norm();
		A_slide_roll.block(0,0,3,3) = A.block(c2brind*3, c2brind*3, 3, 3);
		A_slide_roll.block(3,0,1,3) = A.block(c2bsind*3+2, c2brind*3, 1, 3);
		A_slide_roll.block(0,3,3,1) = A.block(c2brind*3, c2bsind*3+2, 3, 1) - mu*A.block(c2brind*3, c2bsind*3, 3, 2)*slip1_direction;
		A_slide_roll(3,3) = A(c2bsind*3+2, c2bsind*3+2) - mu*A.row(c2bsind*3+2).segment(c2bsind*3, 2).dot(slip1_direction);
		lhs_slide_roll.segment(0,3) = vsol.segment(c2brind*3, 3) - b.segment(c2brind*3, 3);
		lhs_slide_roll(3) = vsol(c2bsind*3+2) - b(c2bsind*3+2);
		// Note that A_slide_roll is not symmetric. So we use QR decomposition to solve.
		psol_slide_roll = A_slide_roll.householderQr().solve(lhs_slide_roll);
		psol.segment(c2brind*3, 3) = psol_slide_roll.segment(0, 3);
		psol.segment(c2bsind*3, 3) << -mu*slip1_direction, 1.0;
		psol.segment(c2bsind*3, 3) *= psol_slide_roll(3);
		FeasibleCollLCPPointResult ret_case2b = testFeasibleCollLCPPoint(2,
			0,
			pre_v, //unused
			psol,
			vsol, //unused
			epsilon,
			mu
		);
		if(ret_case2b == FeasibleCollLCPPointResult::FRSuccess) {
			std::cout << "LCP Success: Contact 1 roll, Contact 2 slide" << std::endl;
			CollLCPPointSolution ret_lcp_sol(LCPSolResult::Success);
			ret_lcp_sol.p_sol = psol;
			return ret_lcp_sol;
		}
	}
	// --------- case 2c: contact 2 rolling, contact 1 sliding
	uint c2crind = 1;
	uint c2csind = 0;
	psol.setZero();
	vsol.setZero();
	vsol(0*3 + 2) = -epsilon * pre_v(0*3 + 2);
	vsol(1*3 + 2) = -epsilon * pre_v(1*3 + 2);
	slip1_direction = pre_v.segment(c2csind*3, 2);
	if(slip1_direction.norm() < 1e-5) {
		sliding_with_no_preslip = true;
	} else {
		slip1_direction /= slip1_direction.norm();
		A_slide_roll.block(0,0,3,3) = A.block(c2crind*3, c2crind*3, 3, 3);
		A_slide_roll.block(3,0,1,3) = A.block(c2csind*3+2, c2crind*3, 1, 3);
		A_slide_roll.block(0,3,3,1) = A.block(c2crind*3, c2csind*3+2, 3, 1) - mu*A.block(c2crind*3, c2csind*3, 3, 2)*slip1_direction;
		A_slide_roll(3,3) = A(c2csind*3+2, c2csind*3+2) - mu*A.row(c2csind*3+2).segment(c2csind*3, 2).dot(slip1_direction);
		lhs_slide_roll.segment(0,3) = vsol.segment(c2crind*3, 3) - b.segment(c2crind*3, 3);
		lhs_slide_roll(3) = vsol(c2csind*3+2) - b(c2csind*3+2);
		// Note that A_slide_roll is not symmetric. So we use QR decomposition to solve.
		psol_slide_roll = A_slide_roll.householderQr().solve(lhs_slide_roll);
		psol.segment(c2crind*3, 3) = psol_slide_roll.segment(0, 3);
		psol.segment(c2csind*3, 3) << -mu*slip1_direction, 1.0;
		psol.segment(c2csind*3, 3) *= psol_slide_roll(3);
		FeasibleCollLCPPointResult ret_case2c = testFeasibleCollLCPPoint(2,
			0,
			pre_v, //unused
			psol,
			vsol, //unused
			epsilon,
			mu
		);
		if(ret_case2c == FeasibleCollLCPPointResult::FRSuccess) {
			std::cout << "LCP Success: Contact 1 slide, Contact 2 roll" << std::endl;
			CollLCPPointSolution ret_lcp_sol(LCPSolResult::Success);
			ret_lcp_sol.p_sol = psol;
			return ret_lcp_sol;
		}
	}
	// --------- case 2d: both contacts sliding
	psol.setZero();
	vsol.setZero();
	vsol(0*3 + 2) = -epsilon * pre_v(0*3 + 2);
	vsol(1*3 + 2) = -epsilon * pre_v(1*3 + 2);
	slip1_direction = pre_v.segment(0*3, 2);
	slip2_direction = pre_v.segment(1*3, 2);
	if(slip1_direction.norm() < 1e-5 || slip2_direction.norm() < 1e-5) {
		sliding_with_no_preslip = true;
	} else {
		slip1_direction /= slip1_direction.norm();
		slip2_direction /= slip2_direction.norm();
		A_slide_slide(0,0) = A(2, 2) - mu*A.row(0*3+2).segment(0*3, 2).dot(slip1_direction);
		A_slide_slide(0,1) = A(2, 5) - mu*A.row(0*3+2).segment(1*3, 2).dot(slip2_direction);
		A_slide_slide(1,1) = A(5, 5) - mu*A.row(1*3+2).segment(1*3, 2).dot(slip2_direction);
		A_slide_slide(1,0) = A(5, 2) - mu*A.row(1*3+2).segment(0*3, 2).dot(slip1_direction);
		lhs_slide_slide(0) = vsol(0*3 + 2) - b(0*3 + 2);
		lhs_slide_slide(1) = vsol(1*3 + 2) - b(1*3 + 2);
		// Note that A_slide_slide is not symmetric, and small 2x2. So we use inverse to solve.
		psol_slide_slide = A_slide_slide.inverse()*lhs_slide_slide;
		psol.segment(0*3, 3) << -mu*slip1_direction, 1.0;
		psol.segment(0*3, 3) *= psol_slide_slide(0);
		psol.segment(1*3, 3) << -mu*slip2_direction, 1.0;
		psol.segment(1*3, 3) *= psol_slide_slide(1);
		FeasibleCollLCPPointResult ret_case2d = testFeasibleCollLCPPoint(2,
			0,
			pre_v, //unused
			psol,
			vsol, //unused
			epsilon,
			mu
		);
		if(ret_case2d == FeasibleCollLCPPointResult::FRSuccess) {
			std::cout << "LCP Success: Contact 1 slide, Contact 2 slide" << std::endl;
			CollLCPPointSolution ret_lcp_sol(LCPSolResult::Success);
			ret_lcp_sol.p_sol = psol;
			return ret_lcp_sol;
		}
	}
	// --------- case3: contact 2 only ---------
	cind = 1;
	nind = 0;
	// --------- case 3a: rolling contact ---------
	psol.setZero();
	vsol.setZero();
	vsol(cind*3 + 2) = -epsilon * pre_v(cind*3 + 2);
	//TODO: ^ use composed A_contact for active contacts
	psol.segment(cind*3, 3) = A.block(cind*3, cind*3, 3, 3).ldlt().solve(
		vsol.segment(cind*3, 3) - b.segment(cind*3, 3)
	);
	// compute non active velocity
	vsol.segment(nind*3, 3) = A.block(nind*3, cind*3, 3, 3) * psol.segment(cind*3, 3) + b.segment(nind*3, 3);
	FeasibleCollLCPPointResult ret_case3a = testFeasibleCollLCPPoint(1,
		1,
		pre_v.segment(nind*3,3),
		psol.segment(cind*3, 3),
		vsol.segment(nind*3, 3),
		epsilon,
		mu
	);
	if(ret_case3a == FeasibleCollLCPPointResult::FRSuccess) {
		std::cout << "LCP Success: No contact 1, Contact 2 roll" << std::endl;
		CollLCPPointSolution ret_lcp_sol(LCPSolResult::Success);
		ret_lcp_sol.p_sol = psol;
		return ret_lcp_sol;
	}
	// --------- case 1b: sliding contact ---------
	vsol.setZero();
	slip1_direction = pre_v.segment(cind*3, 2);
	if(slip1_direction.norm() < 1e-5) {
		sliding_with_no_preslip = true;
	} else {
		slip1_direction /= slip1_direction.norm();
		psol.segment(cind*3, 3) << -mu*slip1_direction, 1.0;
		double normal_impulse = 1/(A.block(cind*3, cind*3, 3, 3).row(2).dot(psol.segment(cind*3, 3)))
							* (-epsilon * pre_v(cind*3 + 2) - b(cind*3 + 2));
		psol.segment(cind*3, 3) *= normal_impulse;
		// compute non active velocity
		vsol.segment(nind*3, 3) = A.block(nind*3, cind*3, 3, 3) * psol.segment(cind*3, 3) + b.segment(nind*3, 3);
		FeasibleCollLCPPointResult ret_case3b = testFeasibleCollLCPPoint(1,
			1,
			pre_v.segment(nind*3,3),
			psol.segment(cind*3, 3),
			vsol.segment(nind*3, 3),
			epsilon,
			mu
		);
		if(ret_case3b == FeasibleCollLCPPointResult::FRSuccess) {
			std::cout << "LCP Success: No contact 1, Contact 2 slide" << std::endl;
			CollLCPPointSolution ret_lcp_sol(LCPSolResult::Success);
			ret_lcp_sol.p_sol = psol;
			return ret_lcp_sol;
		}
	}

	// When all fails:
	if (sliding_with_no_preslip) {
		return CollLCPPointSolution(LCPSolResult::UnimplementedCase);
	}
	return CollLCPPointSolution(LCPSolResult::NoSolution);
}


/* subset of above 2-pt LCP solver. This is just a one point solver. */
CollLCPPointSolution solveCollLCPOnePoint (const Eigen::Matrix3d& A,
	const Eigen::Vector3d& b,
	const Eigen::Vector3d& pre_v,
	double epsilon,
	double mu
) {
	Eigen::Vector3d psol;
	Eigen::Vector3d vsol;
	Eigen::Vector2d slip_direction;

	// --------- case 1a: sliding contact ---------
	bool sliding_with_no_preslip = false;
	psol.setZero();
	vsol.setZero();
	vsol(2) = -epsilon * pre_v(2);
	psol = A.ldlt().solve(vsol - b);
	FeasibleCollLCPPointResult ret_case1a = testFeasibleCollLCPPoint(1,
		0,
		pre_v, // unused
		psol,
		vsol, // unused
		epsilon,
		mu
	);
	if(ret_case1a == FeasibleCollLCPPointResult::FRSuccess) {
		std::cout << "LCP Success: Contact 1 roll, No contact 2" << std::endl;
		CollLCPPointSolution ret_lcp_sol(LCPSolResult::Success);
		ret_lcp_sol.p_sol = psol;
		return ret_lcp_sol;
	}
	// --------- case 1b: sliding contact ---------
	vsol.setZero();
	slip_direction = pre_v.segment(0, 2);
	if(slip_direction.norm() < 1e-5) {
		sliding_with_no_preslip = true;
	} else {
		slip_direction /= slip_direction.norm();
		psol << -mu*slip_direction, 1.0;
		double normal_impulse = 1/(A.row(2).dot(psol))
							* (-epsilon * pre_v(2) - b(2));
		psol *= normal_impulse;
		FeasibleCollLCPPointResult ret_case1a = testFeasibleCollLCPPoint(1,
			0,
			pre_v, // unused
			psol,
			vsol, // unused
			epsilon,
			mu
		);
		if(ret_case1b == FeasibleCollLCPPointResult::FRSuccess) {
			std::cout << "LCP Success: Contact 1 slide, No contact 2" << std::endl;
			CollLCPPointSolution ret_lcp_sol(LCPSolResult::Success);
			ret_lcp_sol.p_sol = psol;
			return ret_lcp_sol;
		}
	}
	// When all fails:
	if (sliding_with_no_preslip) {
		return CollLCPPointSolution(LCPSolResult::UnimplementedCase);
	}
	return CollLCPPointSolution(LCPSolResult::NoSolution);
}

#endif //LCPSOLVER_H