// LCPSolver2.cpp

#include <iostream>
#include "LCPSolver2.h"

using namespace Eigen;

namespace Sai2LCPSolver {

CollLCPPointSolution LCPSolver::solve(
		const Eigen::MatrixXd& A,
		const Eigen::VectorXd& b,
		const Eigen::VectorXd& pre_v,
		double epsilon,
		double mu,
		bool force_sliding_if_pre_slip
) {
	CollLCPPointSolution ret_sol;

	// setup
	states.clear();
	num_points = A.rows()/3;
	TA = A;
	TA_size = 0;
	Trhs = -b;


	// initially set all points to be noContact
	for(uint i = 0; i < num_points; i++) {
		states.push_back(PointState::NoContact);
		frictionless_sliding_directions.push_back(Vector2d::Zero());
		rolling_sliding_directions.push_back(Vector2d::Zero());
		chosen_sliding_directions.push_back(Vector2d::Zero());
		rolling_redundancy_directions.push_back(RollingFrictionRedundancyDir::None);
	}

	uint last_point_index = 0;
	solver_state = SolverState::DeterminingActiveFrictionlessContacts;
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
				if(full_v_sol(3*i + 2) < (-1e-10 + -epsilon*pre_v(3*i + 2))) {
					// we ignore points where pre_v is positive
					// TODO: think more about this. currently, it often causes failure to find a solution
					if(pre_v(3*i+2) > 1e-10) {
						continue;
					}
					// - once a point becomes NoContactAgain, it is not made active again
					// std::cout <<"Full v sol: " << full_v_sol.transpose() << std::endl;
					if(states[i] == PointState::NoContactAgain) {
						// ret_sol.result = LCPSolResult::UnimplementedCase;
						// std::cout << "A: " << std::endl;
						// std::cout << A << std::endl;
						// std::cout << "b: " << b.transpose() <<  std::endl;
						// std::cout << "pre_v: " << pre_v.transpose() << std::endl;
						// std::cout << "epsilon: " << epsilon << std::endl;
						// TODO: ^should this be no solution? hard to say because we are not doing an exhaustive search
						// return ret_sol;

						// force frictionless for now
						states[i] = PointState::Rolling;
						// std::cout << "Force NoContactAgain to frictionless contact " << i << std::endl;
						rolling_redundancy_directions[i] = RollingFrictionRedundancyDir::DirXandY;
					} else {
						// assert(states[i] == PointState::NoContact); // shouldn't be any other state really
						// std::cout << "Enable contact " << i << std::endl;
						if(states[i] != PointState::NoContact) {
							// indicates numerical error. so we relax the penetration velocity constraint
							if(full_v_sol(3*i + 2) > (-1e-6 + -epsilon*pre_v(3*i + 2))) {
								continue;
							}
						} else {
							if(force_sliding_if_pre_slip) {
								Vector2d pre_slip = pre_v.segment<2>(i*3);
								if(pre_slip.norm() > 1e-8) {
									if(LCP_LOG_DEBUG) std::cout << "force_sliding_if_pre_slip " << i << std::endl;
									states[i] = PointState::Sliding;
									chosen_sliding_directions[i] = pre_slip/pre_slip.norm();
								} else {
									enableContact(i);
								}
							} else {
								enableContact(i);
							}
							any_pt_penetrating = true;
						}
					}
					break;
				}
			}
			if(!any_pt_penetrating) {
				solver_state = SolverState::EnforcingRollingFriction;
			}
		}

		// then if no new contact points are being added, points are checked for frictional
		// correctness
		// if slip is occuring, try enforcing rolling friction at the point
		if(solver_state == SolverState::EnforcingRollingFriction) {
			bool any_frictionless_pt_slipping = false;
			for(uint i = 0; i < num_points; i++) {
				if(states[i] != PointState::Frictionless) continue;
				Vector2d curr_slip_speed = full_v_sol.segment<2>(3*i);

				if(curr_slip_speed.norm() > 1e-8) {
					// save current slip direction for this point
					frictionless_sliding_directions[i] = curr_slip_speed/curr_slip_speed.norm();

					enableRollingFriction(i, A, b, pre_v, epsilon, mu);
					// this checks for redundancy directions with existing contacts
					// TODO: if rolling_redundancy_directions[i] == RollingFrictionRedundancyDir::DirXandY,
					// then rolling_friction will be zero. So we will never turn on sliding friction at this
					// point. Perhaps we can turn on sliding friction directly?

					any_frictionless_pt_slipping = true;
					break;
				}
			}
			if(!any_frictionless_pt_slipping) {
				solver_state = SolverState::EnforcingSlidingFriction;
			}
		}
		// if cone violation occurs, clamp its direction to either its initial slip direction
		// if available, or to its current slip direction if available, or to the rolling
		// friction direction. change state to sliding friction.
		// - once a point is sliding, it will remain in sliding contact
		if(solver_state == SolverState::EnforcingSlidingFriction) {
			bool any_pt_friction_cone_violation = false;
			for(uint i = 0; i < num_points; i++) {
				if(states[i] != PointState::Rolling) continue;
				Vector2d rolling_friction = full_p_sol.segment<2>(3*i);

				if(rolling_friction.norm() - mu*abs(full_p_sol(3*i + 2)) > 1e-15) {
					// save current rolling direction for this point
					// std::cout << "Full p sol " << full_p_sol.transpose() << std::endl;
					if(rolling_friction.norm() > 1e-8) {
						rolling_sliding_directions[i] = -rolling_friction/rolling_friction.norm();
					}
					// std::cout << "Enable sliding " << i << std::endl;
					enableSlidingFriction(i, pre_v); // this also handles impending slip
					any_pt_friction_cone_violation = true;
					break;
				}
			}
			if(!any_pt_friction_cone_violation) {
				solver_state = SolverState::EnforcingNonNegativeNormalForce;
			}
		}

		if (solver_state == SolverState::EnforcingNonNegativeNormalForce) {
			// check for normal force non-negative violation
			// if a point violates its non-negative normal force status, it is made noContactAgain
			bool did_disable_any_contacts = false;
			for(uint i = 0; i < num_points; i++) {
				if(states[i] == PointState::Rolling) {
					// for rolling points, if friction cone is violated, then the point will
					// eventually be sliding. so ignore the non-negative force check for now
					Vector2d rolling_friction = full_p_sol.segment<2>(3*i);
					if(rolling_friction.norm() - mu*abs(full_p_sol(3*i + 2)) > 1e-15) {
						continue;
					}
				}
				if(full_p_sol(3*i+2) < -1e-8) {
					// std::cout << "TRhs: " << Trhs.segment(0, TA_size).transpose() << std::endl;
					// std::cout << "Disable contact" << i << std::endl;
					disableContact(i);
					did_disable_any_contacts = true;
					break;
				}
			}
			if(!did_disable_any_contacts) {
				// we are done!
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

		composeMatrices(A, b, pre_v, epsilon, mu);
		// std::cout << "TA" << std::endl;
		// std::cout << TA.block(0,0,TA_size,TA_size) << std::endl;
		// std::cout << "Trhs: " <<  Trhs.segment(0, TA_size).transpose() << std::endl;
		solveReducedMatrices(mu, full_p_sol);
		// std::cout << "TP sol: " << Tp_sol.transpose() << std::endl;

		// solve full equation for full_v_sol
		full_v_sol = A*full_p_sol + b;
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

void LCPSolver::enableContact(uint i) {
	if(LCP_LOG_DEBUG) std::cout << "Enabling contact " << i << std::endl;
	states[i] = PointState::Frictionless;
}

void LCPSolver::disableContact(uint i) {
	if(LCP_LOG_DEBUG) std::cout << "Disabling contact " << i << std::endl;
	states[i] = PointState::NoContactAgain;
}

void LCPSolver::enableRollingFriction(uint i,
	const Eigen::MatrixXd& A,
	const Eigen::VectorXd& b,
	const Eigen::VectorXd& pre_v,
	double epsilon,
	double mu
) {
	states[i] = PointState::Rolling;
	if(LCP_LOG_DEBUG) std::cout << "Enable rolling " << i << std::endl;
	// try enabling all 3 axes
	rolling_redundancy_directions[i] = RollingFrictionRedundancyDir::None;
	composeMatrices(A, b, pre_v, epsilon, mu);
	if(abs(TA.block(0,0,TA_size,TA_size).determinant()) > 1e-5) {
		// TODO: think of a better way of checking
		if(LCP_LOG_DEBUG) std::cout << "Rolling redundancy none " << i << std::endl;
		if(LCP_LOG_DEBUG) std::cout << "Det: " << TA.block(0,0,TA_size,TA_size).determinant() << std::endl;
		return;
	}
	// try enabling Y only
	rolling_redundancy_directions[i] = RollingFrictionRedundancyDir::DirXOnly;
	composeMatrices(A, b, pre_v, epsilon, mu);
	if(abs(TA.block(0,0,TA_size,TA_size).determinant()) > 1e-5) {
		if(LCP_LOG_DEBUG) std::cout << "Rolling redundancy DirXOnly " << i << std::endl;
		return;
	}
	// try enabling X only
	rolling_redundancy_directions[i] = RollingFrictionRedundancyDir::DirYOnly;
	composeMatrices(A, b, pre_v, epsilon, mu);
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
	Vector2d pre_slip = pre_v.segment<2>(i*3);
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

void LCPSolver::composeMatrices(
	const Eigen::MatrixXd& A,
	const Eigen::VectorXd& b,
	const Eigen::VectorXd& pre_v,
	double epsilon,
	double mu
) {
	// assemble TA and Trhs. TODO: think about how we can do this incrementally.
	// maybe save the mapping from point id to TA row id?
	TA_size = 0;
	for(uint i = 0; i < num_points; i++) {
		if(states[i] == PointState::NoContact || states[i] == PointState::NoContactAgain) {
			continue;
		}

		else if(states[i] == PointState::Frictionless ||
				(states[i] == PointState::Rolling && rolling_redundancy_directions[i] == RollingFrictionRedundancyDir::DirXandY)
		) {
			TA(TA_size, TA_size) = A(i*3+2, i*3+2);
			Trhs(TA_size) = -b(i*3+2) - epsilon*pre_v(i*3+2);
			// fill in non-diagonal blocks of TA
			uint r = 0;
			for(uint j = 0; j < i; j++) {
				if(states[j] == PointState::NoContact || states[j] == PointState::NoContactAgain) {
					continue;
				}
				if(states[j] == PointState::Frictionless) {
					TA(r, TA_size) = A(j*3+2, i*3+2);
					TA(TA_size, r) = TA(r, TA_size);
					r += 1;
				}
				if(states[j] == PointState::Rolling) {
					if(rolling_redundancy_directions[j] == RollingFrictionRedundancyDir::None) {
						TA.block<3, 1>(r, TA_size) = A.block<3, 1>(j*3, i*3+2);
						TA.block<1, 3>(TA_size, r) = TA.block<3, 1>(r, TA_size).transpose();
						r += 3;
					} else if(rolling_redundancy_directions[j] == RollingFrictionRedundancyDir::DirXOnly) {
						TA.block<2, 1>(r, TA_size) = A.block<2, 1>(j*3+1, i*3+2);
						TA.block<1, 2>(TA_size, r) = TA.block<2, 1>(r, TA_size).transpose();
						r += 2;
					} else if(rolling_redundancy_directions[j] == RollingFrictionRedundancyDir::DirYOnly) {
						TA(r, TA_size) = A(j*3, i*3+2);
						TA(r+1, TA_size) = A(j*3+2, i*3+2);
						TA(TA_size, r) = TA(r, TA_size);
						TA.block<1, 2>(TA_size, r) = TA.block<2, 1>(r, TA_size).transpose();
						r += 2;
					} else if(rolling_redundancy_directions[j] == RollingFrictionRedundancyDir::DirXandY) {
						TA(r, TA_size) = A(j*3+2, i*3+2);
						TA(TA_size, r) = TA(r, TA_size);
						r += 1;
					}
				}
				if(states[j] == PointState::Sliding) {
					Vector2d slip_dir = chosen_sliding_directions[j];
					TA(r, TA_size) = A(j*3+2, i*3+2);
					TA(TA_size, r) = A(i*3+2, j*3+2) - mu*A.block<1,2>(i*3+2, j*3)*slip_dir;
					r += 1;
				}
			}
			TA_size += 1;
		}

		else if(states[i] == PointState::Rolling && rolling_redundancy_directions[i] == RollingFrictionRedundancyDir::None) {
			TA.block<3,3>(TA_size, TA_size) = A.block<3,3>(i*3, i*3);
			Trhs.segment<2>(TA_size) = -b.segment<2>(i*3);
			Trhs(TA_size + 2) = -b(i*3+2) - epsilon*pre_v(i*3+2);
			// fill in non-diagonal blocks of TA
			uint r = 0;
			for(uint j = 0; j < i; j++) {
				if(states[j] == PointState::NoContact || states[j] == PointState::NoContactAgain) {
					continue;
				}
				if(states[j] == PointState::Frictionless) {
					TA.block<1,3>(r, TA_size) = A.block<1,3>(j*3+2, i*3);
					TA.block<3,1>(TA_size, r) = TA.block<1,3>(r, TA_size).transpose();
					r += 1;
				}
				if(states[j] == PointState::Rolling) {
					if(rolling_redundancy_directions[j] == RollingFrictionRedundancyDir::None) {
						TA.block<3, 3>(r, TA_size) = A.block<3, 3>(j*3, i*3);
						TA.block<3, 3>(TA_size, r) = TA.block<3, 3>(r, TA_size).transpose();
						r += 3;
					} else if(rolling_redundancy_directions[j] == RollingFrictionRedundancyDir::DirXOnly) {
						TA.block<2, 3>(r, TA_size) = A.block<2, 3>(j*3+1, i*3);
						TA.block<3, 2>(TA_size, r) = TA.block<2, 3>(r, TA_size).transpose();
						r += 2;
					} else if(rolling_redundancy_directions[j] == RollingFrictionRedundancyDir::DirYOnly) {
						TA.block<1,3>(r, TA_size) = A.block<1,3>(j*3, i*3);
						TA.block<1,3>(r+1, TA_size) = A.block<1,3>(j*3+2, i*3);
						TA.block<3, 2>(TA_size, r) = TA.block<2, 3>(r, TA_size).transpose();
						r += 2;
					} else if(rolling_redundancy_directions[j] == RollingFrictionRedundancyDir::DirXandY) {
						TA.block<1,3>(r, TA_size) = A.block<1,3>(j*3+2, i*3);
						TA.block<3,1>(TA_size, r) = TA.block<1,3>(r, TA_size).transpose();
						r += 1;
					}
				}
				if(states[j] == PointState::Sliding) {
					Vector2d slip_dir = chosen_sliding_directions[j];
					TA.block<1,3>(r, TA_size) = A.block<1,3>(j*3+2, i*3);
					TA.block<3,1>(TA_size, r) = A.block<3,1>(i*3, j*3+2) - mu*A.block<3,2>(i*3, j*3)*slip_dir;
					r += 1;
				}
			}
			TA_size += 3;
		}

		else if(states[i] == PointState::Rolling && rolling_redundancy_directions[i] == RollingFrictionRedundancyDir::DirXOnly) {
			TA.block<2,2>(TA_size, TA_size) = A.block<2,2>(i*3+1, i*3+1);
			Trhs(TA_size) = -b(i*3+1);
			Trhs(TA_size + 1) = -b(i*3+2) - epsilon*pre_v(i*3+2);
			// fill in non-diagonal blocks of TA
			uint r = 0;
			for(uint j = 0; j < i; j++) {
				if(states[j] == PointState::NoContact || states[j] == PointState::NoContactAgain) {
					continue;
				}
				if(states[j] == PointState::Frictionless) {
					TA.block<1,2>(r, TA_size) = A.block<1,2>(j*3+2, i*3+1);
					TA.block<2,1>(TA_size, r) = TA.block<1,2>(r, TA_size).transpose();
					r += 1;
				}
				if(states[j] == PointState::Rolling) {
					if(rolling_redundancy_directions[j] == RollingFrictionRedundancyDir::None) {
						TA.block<3, 2>(r, TA_size) = A.block<3, 2>(j*3, i*3+1);
						TA.block<2, 3>(TA_size, r) = TA.block<3, 2>(r, TA_size).transpose();
						r += 3;
					} else if(rolling_redundancy_directions[j] == RollingFrictionRedundancyDir::DirXOnly) {
						TA.block<2, 2>(r, TA_size) = A.block<2, 2>(j*3+1, i*3+1);
						TA.block<2, 2>(TA_size, r) = TA.block<2, 2>(r, TA_size).transpose();
						r += 2;
					} else if(rolling_redundancy_directions[j] == RollingFrictionRedundancyDir::DirYOnly) {
						TA.block<1,2>(r, TA_size) = A.block<1,2>(j*3, i*3+1);
						TA.block<1,2>(r+1, TA_size) = A.block<1,2>(j*3+2, i*3+1);
						TA.block<2, 2>(TA_size, r) = TA.block<2, 2>(r, TA_size).transpose();
						r += 2;
					} else if(rolling_redundancy_directions[j] == RollingFrictionRedundancyDir::DirXandY) {
						TA.block<1,2>(r, TA_size) = A.block<1,2>(j*3+2, i*3+1);
						TA.block<2,1>(TA_size, r) = TA.block<1,2>(r, TA_size).transpose();
						r += 1;
					}
				}
				if(states[j] == PointState::Sliding) {
					Vector2d slip_dir = chosen_sliding_directions[j];
					TA.block<1,2>(r, TA_size) = A.block<1,2>(j*3+2, i*3+1);
					TA.block<2,1>(TA_size, r) = A.block<2,1>(i*3+1, j*3+2) - mu*A.block<2,2>(i*3+1, j*3)*slip_dir;
					r += 1;
				}
			}
			TA_size += 2;
		}

		else if(states[i] == PointState::Rolling && rolling_redundancy_directions[i] == RollingFrictionRedundancyDir::DirYOnly) {
			TA(TA_size, TA_size) = A(i*3, i*3);
			TA(TA_size, TA_size+1) = A(i*3, i*3+2);
			TA(TA_size+1, TA_size) = A(i*3+2, i*3);
			TA(TA_size+1, TA_size+1) = A(i*3+2, i*3+2);
			Trhs(TA_size) = -b(i*3);
			Trhs(TA_size + 1) = -b(i*3+2) - epsilon*pre_v(i*3+2);
			// fill in non-diagonal blocks of TA
			uint r = 0;
			for(uint j = 0; j < i; j++) {
				if(states[j] == PointState::NoContact || states[j] == PointState::NoContactAgain) {
					continue;
				}
				if(states[j] == PointState::Frictionless) {
					TA(r, TA_size) = A(j*3+2, i*3);
					TA(r, TA_size+1) = A(j*3+2, i*3+2);
					TA.block<2,1>(TA_size, r) = TA.block<1,2>(r, TA_size).transpose();
					r += 1;
				}
				if(states[j] == PointState::Rolling) {
					if(rolling_redundancy_directions[j] == RollingFrictionRedundancyDir::None) {
						TA.block<3, 1>(r, TA_size) = A.block<3, 1>(j*3, i*3);
						TA.block<3, 1>(r, TA_size+1) = A.block<3, 1>(j*3, i*3+2);
						TA.block<2, 3>(TA_size, r) = TA.block<3, 2>(r, TA_size).transpose();
						r += 3;
					} else if(rolling_redundancy_directions[j] == RollingFrictionRedundancyDir::DirXOnly) {
						TA.block<2, 1>(r, TA_size) = A.block<2, 1>(j*3+1, i*3);
						TA.block<2, 1>(r, TA_size+1) = A.block<2, 1>(j*3+1, i*3+2);
						TA.block<2, 2>(TA_size, r) = TA.block<2, 2>(r, TA_size).transpose();
						r += 2;
					} else if(rolling_redundancy_directions[j] == RollingFrictionRedundancyDir::DirYOnly) {
						TA(r, TA_size) = A(j*3, i*3);
						TA(r+1, TA_size) = A(j*3+2, i*3);
						TA(r, TA_size+1) = A(j*3, i*3+2);
						TA(r+1, TA_size+1) = A(j*3+2, i*3+2);
						TA.block<2, 2>(TA_size, r) = TA.block<2, 2>(r, TA_size).transpose();
						r += 2;
					} else if(rolling_redundancy_directions[j] == RollingFrictionRedundancyDir::DirXandY) {
						TA(r, TA_size) = A(j*3+2, i*3);
						TA(r, TA_size+1) = A(j*3+2, i*3+2);
						TA.block<2,1>(TA_size, r) = TA.block<1,2>(r, TA_size).transpose();
						r += 1;
					}
				}
				if(states[j] == PointState::Sliding) {
					Vector2d slip_dir = chosen_sliding_directions[j];
					TA(r, TA_size) = A(j*3+2, i*3);
					TA(r, TA_size+1) = A(j*3+2, i*3+2);
					TA(TA_size, r) = A(i*3, j*3+2) - mu*A.block<1,2>(i*3, j*3)*slip_dir;
					TA(TA_size+1, r) = A(i*3+2, j*3+2) - mu*A.block<1,2>(i*3+2, j*3)*slip_dir;
					r += 1;
				}
			}
			TA_size += 2;
		}

		else if(states[i] == PointState::Sliding) {
			Vector2d this_slip_dir = chosen_sliding_directions[i];
			TA(TA_size, TA_size) = A(i*3+2, i*3+2) - mu*A.block<1,2>(i*3+2, i*3)*this_slip_dir;
			Trhs(TA_size) = -b(i*3+2) - epsilon*pre_v(i*3+2);
			// fill in non-diagonal blocks of TA
			uint r = 0;
			for(uint j = 0; j < i; j++) {
				if(states[j] == PointState::NoContact || states[j] == PointState::NoContactAgain) {
					continue;
				}
				if(states[j] == PointState::Frictionless) {
					TA(r, TA_size) = A(j*3+2, i*3+2) - mu*A.block<1,2>(j*3+2, i*3)*this_slip_dir;
					TA(TA_size, r) = TA(r, TA_size);
					r += 1;
				}
				if(states[j] == PointState::Rolling) {
					if(rolling_redundancy_directions[j] == RollingFrictionRedundancyDir::None) {
						TA.block<3, 1>(r, TA_size) = A.block<3, 1>(j*3, i*3+2) - mu*A.block<3, 2>(j*3, i*3)*this_slip_dir;
						TA.block<1, 3>(TA_size, r) = A.block<1, 3>(i*3+2, j*3);
						r += 3;
					} else if(rolling_redundancy_directions[j] == RollingFrictionRedundancyDir::DirXOnly) {
						TA.block<2, 1>(r, TA_size) = A.block<2, 1>(j*3+1, i*3+2) - mu*A.block<2, 2>(j*3+1, i*3)*this_slip_dir;
						TA.block<1, 2>(TA_size, r) = A.block<1, 2>(i*3+2, j*3+1);
						r += 2;
					} else if(rolling_redundancy_directions[j] == RollingFrictionRedundancyDir::DirYOnly) {
						TA(r, TA_size) = A(j*3, i*3+2) - mu*A.block<1, 2>(j*3, i*3)*this_slip_dir;
						TA(r+1, TA_size) = A(j*3+2, i*3+2) - mu*A.block<1, 2>(j*3+2, i*3)*this_slip_dir;
						TA(TA_size, r) = A(i*3+2, j*3);
						TA(TA_size, r+1) = A(i*3+2, j*3+2);
						r += 2;
					} else if(rolling_redundancy_directions[j] == RollingFrictionRedundancyDir::DirXandY) {
						TA(r, TA_size) = A(j*3+2, i*3+2) - mu*A.block<1,2>(j*3+2, i*3)*this_slip_dir;
						TA(TA_size, r) = TA(r, TA_size);
						r += 1;
					}
				}
				if(states[j] == PointState::Sliding) {
					Vector2d slip_dir = chosen_sliding_directions[j];
					TA(r, TA_size) = A(j*3+2, i*3+2) - mu*A.block<1, 2>(j*3+2, i*3)*this_slip_dir;;
					TA(TA_size, r) = A(i*3+2, j*3+2) - mu*A.block<1, 2>(i*3+2, j*3)*slip_dir;
					r += 1;
				}
			}
			TA_size += 1;
		}
	}
}

void LCPSolver::solveReducedMatrices(
	double mu,
	Eigen::VectorXd& full_p_sol
) {
	// At any point, in debug mode, if determinant goes to zero, then abort for now
	// TODO: accomodate redundancy hints from geometry/ rigid body kinematics that
	// can be computed at a slow model update rate
	assert(abs(TA.block(0,0,TA_size,TA_size).determinant()) > 1e-15);

	Tp_sol = TA.block(0,0,TA_size,TA_size).partialPivLu().solve(Trhs.segment(0,TA_size));

	// reassemble full_p_sol
	full_p_sol.setZero();
	uint TA_ind = 0;
	for(uint i = 0; i < num_points; i++) {
		if(states[i] == PointState::NoContact || states[i] == PointState::NoContactAgain) {
			// nothing to do
		} else if(states[i] == PointState::Frictionless ||
			(states[i] == PointState::Rolling && rolling_redundancy_directions[i] == RollingFrictionRedundancyDir::DirXandY)
		) {
			full_p_sol(i*3 + 2) = Tp_sol(TA_ind);
			TA_ind += 1;
		} else if(states[i] == PointState::Rolling && rolling_redundancy_directions[i] == RollingFrictionRedundancyDir::None) {
			full_p_sol.segment<3>(i*3) = Tp_sol.segment<3>(TA_ind);
			TA_ind += 3;
		} else if(states[i] == PointState::Rolling && rolling_redundancy_directions[i] == RollingFrictionRedundancyDir::DirXOnly) {
			full_p_sol.segment<2>(i*3+1) = Tp_sol.segment<2>(TA_ind);
			TA_ind += 2;
		} else if(states[i] == PointState::Rolling && rolling_redundancy_directions[i] == RollingFrictionRedundancyDir::DirYOnly) {
			full_p_sol(i*3) = Tp_sol(TA_ind);
			full_p_sol(i*3+2) = Tp_sol(TA_ind+1);
			TA_ind += 2;
		} else if(states[i] == PointState::Sliding) {
			Vector2d slip_dir = chosen_sliding_directions[i];
			full_p_sol.segment<2>(i*3) = -mu*slip_dir*Tp_sol(TA_ind);
			full_p_sol(i*3+2) = Tp_sol(TA_ind);
			TA_ind += 1;
		}
	}
}

};
