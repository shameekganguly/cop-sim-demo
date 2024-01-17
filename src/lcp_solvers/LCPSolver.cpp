// LCPSolver.cpp

#include <iostream>
#include "LCPSolver.h"
#include "LCPSolverInternal.h"

using namespace Eigen;

namespace Sai2LCPSolver {
namespace {

class Composer {
public:
void composeMatrices(
	// TA and Trhs are assumed to be instantiated and resized outside this function
	Eigen::MatrixXd& TA,
	Eigen::VectorXd& Trhs,
	uint& TA_size,
	const std::vector<PointState>& pt_states,
	const std::vector<Eigen::Vector2d>& pt_chosen_sliding_directions,
	const std::vector<RollingFrictionRedundancyDir>& pt_rolling_redundancy_directions,
	const std::vector<MomentConstraints>& pt_moment_constraints
) {
	// assemble TA and Trhs. TODO: think about how we can do this incrementally.
	// maybe save the mapping from point id to TA row id?

	const uint num_points = pt_states.size();
	left_mask.fill(0);
	right_mask.fill(0);
	Trhs.fill(0);

	TA_size = 0;
	for(uint i = 0; i < num_points; i++) {
		const uint i_ind_start = i*contact_size;
		if(pt_states[i] == PointState::NoContact ||
		   pt_states[i] == PointState::NoContactAgain) {
			continue;
		}

		// Force components
		else if(pt_states[i] == PointState::Frictionless ||
		        (pt_states[i] == PointState::Rolling &&
		        pt_rolling_redundancy_directions[i] ==
		        	RollingFrictionRedundancyDir::DirXandY)) {
			right_mask(i_ind_start+2, TA_size) = 1;
			left_mask(TA_size, i_ind_start+2) = 1;
			TA_size += 1;
		}

		else if(pt_states[i] == PointState::Rolling &&
			    pt_rolling_redundancy_directions[i] ==
			    	RollingFrictionRedundancyDir::None) {
			right_mask.block(i_ind_start, TA_size, 3, 3).setIdentity();
			left_mask.block(TA_size, i_ind_start, 3, 3).setIdentity();
			TA_size += 3;
		}

		else if(pt_states[i] == PointState::Rolling &&
				pt_rolling_redundancy_directions[i] ==
					RollingFrictionRedundancyDir::DirXOnly) {
			right_mask.block(i_ind_start+1, TA_size, 2, 2).setIdentity();
			left_mask.block(TA_size, i_ind_start+1, 2, 2).setIdentity();
			TA_size += 2;
		}

		else if(pt_states[i] == PointState::Rolling &&
				pt_rolling_redundancy_directions[i] ==
					RollingFrictionRedundancyDir::DirYOnly) {
			right_mask(i_ind_start, TA_size) = 1;
			left_mask(TA_size, i_ind_start) = 1;

			right_mask(i_ind_start+2, TA_size+1) = 1;
			left_mask(TA_size+1, i_ind_start+2) = 1;

			TA_size += 2;
		}

		else if(pt_states[i] == PointState::Sliding) {
			Vector2d slip_dir = pt_chosen_sliding_directions[i];
			right_mask.col(TA_size).segment<2>(i_ind_start) = -mu*slip_dir;
			right_mask(i_ind_start+2, TA_size) = 1;

			left_mask(TA_size, i_ind_start+2) = 1;
			TA_size += 1;
		}

		// Moment components
		if(contact_size == 3) {
			continue;
		}

		if(pt_moment_constraints[i] == MomentConstraints::NoMoment) {
			// nothing to do
			continue;
		} else if(pt_moment_constraints[i] == MomentConstraints::XMoment) {
			right_mask(i_ind_start+3, TA_size) = 1;
			left_mask(TA_size, i_ind_start+3) = 1;
			TA_size += 1;
		} else if(pt_moment_constraints[i] == MomentConstraints::YMoment) {
			right_mask(i_ind_start+4, TA_size) = 1;
			left_mask(TA_size, i_ind_start+4) = 1;
			TA_size += 1;
		} else if(pt_moment_constraints[i] == MomentConstraints::XandYMoments) {
			right_mask.block(i_ind_start+3, TA_size, 2, 2).setIdentity();
			left_mask.block(TA_size, i_ind_start+3, 2, 2).setIdentity();
			TA_size += 2;
		}
	}

	TA.block(0, 0, TA_size, TA_size) = left_mask * A * right_mask;
	Trhs.head(TA_size) = left_mask * (post_v_min - b);
}

void solveReducedMatrices(
	const Eigen::MatrixXd& red_A,
	const Eigen::VectorXd& red_rhs,
    // We need this because only red_A.block(0, 0, red_A_size, red_A_size) is usable.
    // Rest is garbage.
	const uint red_A_size,
	const std::vector<PointState>& pt_states,
	const std::vector<Eigen::Vector2d>& pt_chosen_sliding_directions,
	const std::vector<RollingFrictionRedundancyDir>& pt_rolling_redundancy_directions,
	const std::vector<MomentConstraints>& pt_moment_constraints,
	Eigen::VectorXd& full_p_sol
) {
	const uint num_points = pt_states.size();

	// In debug mode, assert that red_A is full rank before solving.
	// std::cout << "red A determinant: "
	// 		  << abs(red_A.block(0,0,red_A_size,red_A_size).determinant()) << std::endl;
	assert(abs(red_A.block(0,0,red_A_size,red_A_size).determinant()) > 1e-15);

	// Solve the reduced equation red_A*p_sol = red_rhs that we assembled in composeMatrices
	VectorXd p_sol = red_A.block(0,0,red_A_size,red_A_size)
						  .partialPivLu()
						  .solve(red_rhs.segment(0,red_A_size));

	// copy impulses from p_sol to full_p_sol
	full_p_sol.setZero();
	uint red_A_ind = 0;
	for(uint i = 0; i < num_points; i++) {
		const size_t full_p_sol_ind = i*contact_size;
		if(pt_states[i] == PointState::NoContact || pt_states[i] == PointState::NoContactAgain) {
			// nothing to do
			continue;
		}
		// Force component
		if(pt_states[i] == PointState::Frictionless ||
			(pt_states[i] == PointState::Rolling && pt_rolling_redundancy_directions[i] == RollingFrictionRedundancyDir::DirXandY)
		) {
			// Only copy z impulse
			full_p_sol(full_p_sol_ind + 2) = p_sol(red_A_ind);
			red_A_ind += 1;
		} else if(pt_states[i] == PointState::Rolling && pt_rolling_redundancy_directions[i] == RollingFrictionRedundancyDir::None) {
			// Copy x, y and z impulse
			full_p_sol.segment<3>(full_p_sol_ind) = p_sol.segment<3>(red_A_ind);
			red_A_ind += 3;
		} else if(pt_states[i] == PointState::Rolling && pt_rolling_redundancy_directions[i] == RollingFrictionRedundancyDir::DirXOnly) {
			// Copy y and z impulse
			full_p_sol.segment<2>(full_p_sol_ind+1) = p_sol.segment<2>(red_A_ind);
			red_A_ind += 2;
		} else if(pt_states[i] == PointState::Rolling && pt_rolling_redundancy_directions[i] == RollingFrictionRedundancyDir::DirYOnly) {
			// Copy x and z impulse
			full_p_sol(full_p_sol_ind) = p_sol(red_A_ind);
			full_p_sol(full_p_sol_ind+2) = p_sol(red_A_ind+1);
			red_A_ind += 2;
		} else if(pt_states[i] == PointState::Sliding) {
			// z impulse = solution impulse
			full_p_sol(full_p_sol_ind+2) = p_sol(red_A_ind);

			// [x impulse, y impulse] = -mu * solution impulse * slip direction vector
			Vector2d slip_dir = pt_chosen_sliding_directions[i];
			full_p_sol.segment<2>(full_p_sol_ind) = -mu*slip_dir*p_sol(red_A_ind);
			red_A_ind += 1;
		}

		// Moment component
		if(contact_size == 3) {
			continue;
		}

		if(pt_moment_constraints[i] == MomentConstraints::NoMoment) {
			// nothing to do
		} else if(pt_moment_constraints[i] == MomentConstraints::XMoment) {
			full_p_sol(full_p_sol_ind+3) = p_sol(red_A_ind);
			red_A_ind += 1;
		} else if(pt_moment_constraints[i] == MomentConstraints::YMoment) {
			full_p_sol(full_p_sol_ind+4) = p_sol(red_A_ind);
			red_A_ind += 1;
		} else if(pt_moment_constraints[i] == MomentConstraints::XandYMoments) {
			full_p_sol.segment<2>(full_p_sol_ind+3) = p_sol.segment<2>(red_A_ind);
			red_A_ind += 2;
		}
	}
	assert(red_A_ind == red_A_size);
}

Composer(const Eigen::MatrixXd& p_A,
		 const Eigen::VectorXd& p_b,
		 const Eigen::VectorXd& p_pre_v,
		 const Eigen::VectorXd& p_post_v_min,
		 uint p_contact_size,  // either 3 for pts or 5 for line/ surface contacts
		 const double p_mu)
: A(p_A), b(p_b), pre_v(p_pre_v), post_v_min(p_post_v_min), contact_size(p_contact_size),
  mu(p_mu) {
  	assert(p_contact_size == 3 || p_contact_size == 5);
	left_mask = MatrixXd::Zero(A.rows(), A.cols());
	right_mask = MatrixXd::Zero(A.rows(), A.cols());
}

public:
	MatrixXd left_mask;
	MatrixXd right_mask;
	const MatrixXd& A;
	const VectorXd& b;
	const VectorXd& pre_v;
	const VectorXd& post_v_min;
	const uint contact_size;
	const double mu;
};

}  // namespace

CollLCPPointSolution LCPSolver::solveWithMoments(
		const Eigen::MatrixXd& A,
		const Eigen::VectorXd& b,
		const Eigen::VectorXd& pre_v,
		const std::vector<MomentConstraints>& initial_moment_constraints,
		const double mu
) {
	return solveInternal(/*contact_size=*/ 5, A, b, pre_v, initial_moment_constraints,
		                 /*epsilon=*/ 0, mu, /*force_sliding_if_pre_slip=*/ true);
}

CollLCPPointSolution LCPSolver::solve(
		const Eigen::MatrixXd& A,
		const Eigen::VectorXd& b,
		const Eigen::VectorXd& pre_v,
		const double epsilon,
		const double mu,
		const bool force_sliding_if_pre_slip
) {
	num_points = A.rows()/3;
	std::vector<MomentConstraints> dummy_moment_constraints(
									num_points,
									MomentConstraints::NoMoment);
	return solveInternal(/*contact_size=*/ 3, A, b, pre_v, dummy_moment_constraints,
		                 epsilon, mu, force_sliding_if_pre_slip);
}

CollLCPPointSolution LCPSolver::solveInternal(
		const uint contact_size,
		const Eigen::MatrixXd& A,
		const Eigen::VectorXd& b,
		const Eigen::VectorXd& pre_v,
		const std::vector<MomentConstraints>& initial_moment_constraints,
		const double epsilon,
		const double mu,
		const bool force_sliding_if_pre_slip
) {
	this->contact_size = contact_size;
	num_points = A.rows()/contact_size;

	if (contact_size == 3 && num_points == 1) {
		// TODO: extend to contact size > 3
		// Use optimized 1 pt solution
		return solveCollLCPOnePoint(A, b, pre_v, epsilon, mu);
	}

	CollLCPPointSolution ret_sol;

	// setup
	TA = A;
	TA_size = 0;
	Trhs = -b;
	VectorXd post_v_min;

	// initially set all points to be noContact and noMoment
	states = std::vector<PointState>(num_points, PointState::NoContact);
	frictionless_sliding_directions = std::vector<Vector2d>(num_points, Vector2d::Zero());
	rolling_sliding_directions = std::vector<Vector2d>(num_points, Vector2d::Zero());
	chosen_sliding_directions = std::vector<Vector2d>(num_points, Vector2d::Zero());
	rolling_redundancy_directions =
		std::vector<RollingFrictionRedundancyDir>(num_points, RollingFrictionRedundancyDir::None);

	post_v_min.setZero(b.size());
	for(uint i = 0; i < num_points; i++) {
		assert(pre_v(contact_size*i + 2) <= 2e-4);
		post_v_min(contact_size*i + 2) = -epsilon*pre_v(contact_size*i + 2);
	}

	// Initialize moment constraints with the given constraints. But note that if a point
	// is not in contact, moment constraints will be ignored.
	moment_constraints = initial_moment_constraints;

	Composer composer(A, b, pre_v, post_v_min, contact_size, mu);

	uint last_point_index = 0;
	VectorXd full_v_sol = b;
	VectorXd full_p_sol = VectorXd::Zero(b.size());
	int max_iters = num_points * 10;
	// std::cout << max_iters << std::endl;
	int iters = 0;
	while(iters < max_iters) {
		iters++;

		// TODO: think more about the state transitions

		// points go from being inactive to frictionless
		// if a point that was previously noContact has penetration, then it is made frictionless
		// - this can cause previously rolling solutions to lose friction cone, this is
		// handled further below
		solver_state = SolverState::DeterminingActiveFrictionlessContacts;
		if(solver_state == SolverState::DeterminingActiveFrictionlessContacts) {
			bool any_pt_penetrating = false;
			for(uint i = 0; i < num_points; i++) {
				const uint i_ind_start = contact_size*i;
				if(full_v_sol(i_ind_start + 2) < (-1e-10 + post_v_min(i_ind_start + 2))) {
					// we ignore points where pre_v is positive
					// TODO: think more about this. currently, it often causes failure
					// to find a solution.

					// std::cout << "Pt " << i << " penetrating with vel "
						// << full_v_sol(i_ind_start + 2) << std::endl;

					if(epsilon > 0 && pre_v(i_ind_start+2) > 1e-10) {
						if(LCP_LOG_DEBUG) {
							std::cout << "Pre-V positive " << pre_v(i_ind_start+2)
									  << " at pt " << i << ". Ignoring penetration.\n";
						}
						continue;
					}

					// - once a point becomes NoContactAgain, it is not made active again
					// std::cout <<"Full v sol: " << full_v_sol.transpose() << std::endl;
					if(states[i] == PointState::NoContactAgain) {
						// TODO: should this be no solution? hard to say because we are
						// not doing an exhaustive search.

						// force frictionless for now
						// We don't use PointState::Frictionless because the point can be
						// set to PointState::Rolling state with friction again and keep
						// looping between these states.
						// This heuristic circumvents Painleve's problem by sacrificing
						// frictional correctness to preserve non-penetration.

						if(LCP_LOG_DEBUG) {
							std::cout << "Force NoContactAgain to frictionless contact " << i << std::endl;
						}
						any_pt_penetrating = true;
						states[i] = PointState::Rolling;
						rolling_redundancy_directions[i] = RollingFrictionRedundancyDir::DirXandY;
						continue;
					}

					// If pt is already set to contact, but still penetrating, this is
					// likely due to numerical error. so we check again with relaxed
					// penetration velocity constraint
					if(states[i] != PointState::NoContact) {
						if (full_v_sol(i_ind_start + 2) > (-1e-6 + post_v_min(i_ind_start + 2))) {
							if(LCP_LOG_DEBUG) {
								std::cout << "penetration for contact pt " << i
										  << ". No penetration after relaxing constraint." << std::endl;
							}
							continue;
						} else {
							// We shouldn't be here! This is a solver failure.
							std::cout << "LCPSolver failed because contact point " << i
									  << " was found in state " << states[i]
									  << " to be in penetration with velocity"
									  << full_v_sol(i_ind_start + 2)
									  << ". Needed at least " << post_v_min(i_ind_start + 2)
									  << std::endl;
							solver_state = SolverState::Failed;
							break;
						}
					}

					// Point is in NoContact state and penetrating.
					any_pt_penetrating = true;

					// If force_sliding_if_pre_slip is false, start in frictionless state.
					// If force_sliding_if_pre_slip is true but there is no slip at the
					// point, also start in frictionless state.
					if(!force_sliding_if_pre_slip ||
					   pre_v.segment<2>(i_ind_start).norm() < 1e-8) {
						enableFrictionlessContact(i);
					} else {
						// If force_sliding_if_pre_slip is true and point is slipping,
						// skip directly to sliding state.
						if(LCP_LOG_DEBUG) std::cout << "force_sliding_if_pre_slip " << i << std::endl;
						enableSlidingFriction(i, pre_v);
					}
					break;
				}
			}  // end looping over points

			if(solver_state == SolverState::Failed) {
				break;
			} else if(!any_pt_penetrating) {
				if (mu > 0) {
					solver_state = SolverState::EnforcingRollingFriction;
				} else {
					solver_state = SolverState::EnforcingNonNegativeNormalForce;
				}
			}
		}

		// then if no new contact points are being added, points are checked for frictional
		// correctness
		// if slip is occuring, try enforcing rolling friction at the point
		if(solver_state == SolverState::EnforcingRollingFriction) {
			bool any_frictionless_pt_slipping = false;
			for(uint i = 0; i < num_points; i++) {
				const uint i_ind_start = contact_size*i;
				if(states[i] != PointState::Frictionless) {
					// If the point is not in contact or is already set to be in rolling
					// or sliding contact, there's nothing to do.
					continue;
				}

				const Vector2d curr_slip_speed = full_v_sol.segment<2>(i_ind_start);
				if(curr_slip_speed.norm() > 1e-8) {
					// save current slip direction for this point
					frictionless_sliding_directions[i] = curr_slip_speed/curr_slip_speed.norm();

					enableRollingFriction(i, &composer);
					// this checks for redundancy directions with existing contacts
					// TODO: if rolling_redundancy_directions[i] == RollingFrictionRedundancyDir::DirXandY,
					// then rolling_friction will be zero. So we will never turn on sliding friction at this
					// point. Perhaps we can turn on sliding friction directly?

					any_frictionless_pt_slipping = true;
					break;
				}
			}  // end looping over points

			if(solver_state == SolverState::Failed) {
				break;
			} else if(!any_frictionless_pt_slipping) {
				solver_state = SolverState::EnforcingFrictionCone;
			}
		}
		// if friction cone violation occurs, clamp friction force direction to either
		// the initial slip direction if available, or to its current slip direction if
		// available, or to the rolling friction force direction.
		// Also change state to sliding friction.
		// - once a point is sliding, it will remain in sliding contact.
		if(solver_state == SolverState::EnforcingFrictionCone) {
			bool any_pt_friction_cone_violation = false;
			for(uint i = 0; i < num_points; i++) {
				const uint i_ind_start = contact_size*i;
				if(states[i] != PointState::Rolling) {
					// Either point is not in contact or is frictionless or
					// is sliding and the slip direction is already known.
					continue;
				}

				const Vector2d rolling_friction = full_p_sol.segment<2>(i_ind_start);
				if(rolling_friction.norm() - mu*abs(full_p_sol(i_ind_start + 2)) > 1e-15) {
					// save current rolling direction for this point
					if(rolling_friction.norm() > 1e-8) {
						rolling_sliding_directions[i] = -rolling_friction/rolling_friction.norm();
					}
					enableSlidingFriction(i, pre_v); // this also handles impending slip
					any_pt_friction_cone_violation = true;
					break;
				}
			}  // end looping over points

			if(solver_state == SolverState::Failed) {
				break;
			} else if(!any_pt_friction_cone_violation) {
				solver_state = SolverState::EnforcingNonNegativeNormalForce;
			}
		}

		if (solver_state == SolverState::EnforcingNonNegativeNormalForce) {
			// check for normal force non-negative violation
			// if a point violates its non-negative normal force status, it is made noContactAgain
			bool did_disable_any_contacts = false;
			bool did_ignore_rolling_contacts = false;
			for(uint i = 0; i < num_points; i++) {
				const uint i_ind_start = contact_size*i;
				if(states[i] == PointState::Rolling) {
					// for rolling points, if friction cone is violated, then the point will
					// eventually be sliding. so ignore the non-negative force check for now.
					// We can get into this state if another point previously in contact was
					// disabled in the last solver iteration.
					Vector2d rolling_friction = full_p_sol.segment<2>(i_ind_start);
					if(rolling_friction.norm() - mu*abs(full_p_sol(i_ind_start + 2)) > 1e-15) {
						did_ignore_rolling_contacts = true;
						continue;
					}
				}
				if(full_p_sol(i_ind_start+2) < -1e-8) {
					// std::cout << "TRhs: " << Trhs.segment(0, TA_size).transpose() << std::endl;
					// std::cout << "Disable contact" << i << std::endl;
					disableContact(i);
					did_disable_any_contacts = true;
					break;
				}
			}  // end looping over points

			if(solver_state == SolverState::Failed) {
				break;
			} else if(!did_ignore_rolling_contacts && !did_disable_any_contacts) {
				// We are done!
				if(LCP_LOG_DEBUG) {
					std::cout << "Solved LCP in " << iters << " iterations\n";
				}
				ret_sol.p_sol = full_p_sol;
				ret_sol.result = LCPSolResult::Success;
				return ret_sol;
			}
		}

		// - Alternatively, it can cause some previously sliding points to become energy injecting
		// as they can roll. We ignore this case. Ideally, should switch those points to rolling.
		// if a point that was previously noContactAgain has penetration, then fail for now

		// solve reduced equation
		// for(uint i = 0; i < num_points; i++) {
		// 	std::cout << i << " State: " << states[i] << ", ";
		// }
		// std::cout << std::endl;

		composer.composeMatrices(TA, Trhs, TA_size,
								 states, chosen_sliding_directions,
								 rolling_redundancy_directions, moment_constraints);
		// std::cout << "TA" << std::endl;
		// std::cout << TA.block(0,0,TA_size,TA_size) << std::endl;
		// std::cout << "Trhs: " <<  Trhs.segment(0, TA_size).transpose() << std::endl;
		composer.solveReducedMatrices(TA, Trhs, TA_size,
									  states, chosen_sliding_directions,
									  rolling_redundancy_directions,
									  moment_constraints, full_p_sol);
		// std::cout << "TP sol: " << Tp_sol.transpose() << std::endl;

		// solve full equation for full_v_sol
		full_v_sol = A*full_p_sol + b;
		// std::cout << "Full p sol: " << full_p_sol.transpose() << std::endl;
		// std::cout << "Full v sol: " << full_v_sol.transpose() << std::endl;
	}

	// if we are here, that means the resolution failed
	// Eigen::IOFormat logVecFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ");
	std::cerr << "LCPSolver(2) failed to converge" << std::endl;
	// std::cerr << "A:" << std::endl;
	// std::cerr << A.format(logVecFmt) << std::endl;
	// std::cerr << "b: " << b.transpose().format(logVecFmt) << std::endl;
	// std::cerr << "pre_v: " << pre_v.transpose().format(logVecFmt) << std::endl;
	ret_sol.result = LCPSolResult::NoSolution;
	return ret_sol;
}

void LCPSolver::enableFrictionlessContact(uint i) {
	if(LCP_LOG_DEBUG) std::cout << "Enabling contact " << i << std::endl;
	states[i] = PointState::Frictionless;
}

void LCPSolver::disableContact(uint i) {
	if(LCP_LOG_DEBUG) std::cout << "Disabling contact " << i << std::endl;
	states[i] = PointState::NoContactAgain;
}

void LCPSolver::enableRollingFriction(uint i, Composer* composer) {
	states[i] = PointState::Rolling;
	if(LCP_LOG_DEBUG) std::cout << "Enable rolling " << i << std::endl;

	// TODO: accomodate redundancy hints from geometry/ rigid body kinematics that
	// can be computed at a slow model update rate
	// try enabling all 3 axes
	rolling_redundancy_directions[i] = RollingFrictionRedundancyDir::None;
	composer->composeMatrices(TA, Trhs, TA_size,
							  states, chosen_sliding_directions,
							  rolling_redundancy_directions, moment_constraints);
	if(abs(TA.block(0,0,TA_size,TA_size).determinant()) > 1e-5) {
		// TODO: think of a better way of checking
		if(LCP_LOG_DEBUG) std::cout << "Rolling redundancy none " << i << std::endl;
		if(LCP_LOG_DEBUG) std::cout << "Det: " << TA.block(0,0,TA_size,TA_size).determinant() << std::endl;
		return;
	}
	// try enabling Y only
	rolling_redundancy_directions[i] = RollingFrictionRedundancyDir::DirXOnly;
	composer->composeMatrices(TA, Trhs, TA_size,
							  states, chosen_sliding_directions,
							  rolling_redundancy_directions, moment_constraints);
	if(abs(TA.block(0,0,TA_size,TA_size).determinant()) > 1e-5) {
		if(LCP_LOG_DEBUG) std::cout << "Rolling redundancy DirXOnly " << i << std::endl;
		return;
	}
	// try enabling X only
	rolling_redundancy_directions[i] = RollingFrictionRedundancyDir::DirYOnly;
	composer->composeMatrices(TA, Trhs, TA_size,
							  states, chosen_sliding_directions,
							  rolling_redundancy_directions, moment_constraints);
	if(abs(TA.block(0,0,TA_size,TA_size).determinant()) > 1e-5) {
		if(LCP_LOG_DEBUG) std::cout << "Rolling redundancy DirYOnly " << i << std::endl;
		return;
	}
	if(LCP_LOG_DEBUG) std::cout << "Rolling redundancy DirXandY " << i << std::endl;
	rolling_redundancy_directions[i] = RollingFrictionRedundancyDir::DirXandY;
	return;
}

void LCPSolver::enableSlidingFriction(uint i, const Eigen::VectorXd& pre_v) {
	states[i] = PointState::Sliding;
	if(LCP_LOG_DEBUG) std::cout << "Enable sliding " << i << std::endl;
	Vector2d pre_slip = pre_v.segment<2>(i*contact_size);
	if(pre_slip.norm() > 1e-8) {
		chosen_sliding_directions[i] = pre_slip/pre_slip.norm();
		if(LCP_LOG_DEBUG) std::cout << "Use pre_v slip direction " << i << std::endl;
		return;
	}
	Vector2d anti_roll_slip = rolling_sliding_directions[i];
	if(anti_roll_slip.norm() > 1e-8) {
		// TODO: blend rolling and frictionless slip directions
		chosen_sliding_directions[i] = anti_roll_slip/anti_roll_slip.norm();
		if(LCP_LOG_DEBUG) std::cout << "Use anti-roll slip direction " << i << std::endl;
		return;
	}
	Vector2d frictionless_slip = frictionless_sliding_directions[i];
	if(frictionless_slip.norm() > 1e-8) {
		// TODO: blend rolling and frictionless slip directions
		chosen_sliding_directions[i] = frictionless_slip/frictionless_slip.norm();
		if(LCP_LOG_DEBUG) std::cout << "Use frictionless slip direction (slip" <<
			"with no friction force applied, but due to forces at other pts ) " << i << std::endl;
		return;
	}
	// if all fails, simply set the point to be forced frictionless
	states[i] = PointState::Rolling;
	if(LCP_LOG_DEBUG) std::cout << "Force frictionless " << i << std::endl;
	rolling_redundancy_directions[i] = RollingFrictionRedundancyDir::DirXandY;
}

};
