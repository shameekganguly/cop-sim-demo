// LCPSolver2.h

#ifndef LCP_SOLVER2_H
#define LCP_SOLVER2_H

#include <vector>

#include "LCPSolver.h"

namespace Sai2LCPSolver {

enum PointState {
	NoContact,
	NoContactAgain,
	Frictionless,
	Rolling,
	// Impending, // TODO: do we need this?
	Sliding
};

enum MomentConstraints {
	NoMoment,
	XMoment,
	YMoment,
	XandYMoments
};

enum RollingFrictionRedundancyDir {
	None,
	DirXOnly,
	DirYOnly,
	DirXandY
};

namespace {
// Internal helper class
class Composer;
}

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

// Pivoting solver. TODO: create a main interface then subclass?
class LCPSolver {
public:
	LCPSolver()
	:solver_state(NoContacts)
	{
		// nothing to do
	}

	// Special case for all pt contacts without moment constraints
	CollLCPPointSolution solve(
		const Eigen::MatrixXd& A,
		const Eigen::VectorXd& b,
		const Eigen::VectorXd& pre_v,
		const double epsilon,
		const double mu,
		// This is set to true for computing steady contact forces instead of collision
		// impulses.
		bool force_sliding_if_pre_slip = false
	);

public: // internal functions
	void enableFrictionlessContact(uint i);

	void disableContact(uint i);

	void enableRollingFriction(uint i, Composer* composer);

	// computes and sets the sliding friction constraint
	void enableSlidingFriction(uint i, const Eigen::VectorXd& pre_v);

public: //internal variables
	enum SolverState {
		NoContacts,
		DeterminingActiveFrictionlessContacts,
		EnforcingRollingFriction,
		EnforcingFrictionCone,
		EnforcingNonNegativeNormalForce,
		// EnforcingImpendingSlip,
		Failed
		// Succeeded
	};

	uint num_points;
	SolverState solver_state;
	Eigen::MatrixXd TA;
	Eigen::VectorXd Trhs;
	uint TA_size;
	std::vector<PointState> states;
	// std::list<uint> no_contact_points;
	// std::list<uint> no_contact_again_points;
	// std::list<uint> frictionless_points;
	// std::list<uint> rolling_points;
	// std::list<uint> impending_points;
	// std::list<uint> sliding_points;
	std::vector<Eigen::Vector2d> frictionless_sliding_directions;
	std::vector<Eigen::Vector2d> rolling_sliding_directions;
	std::vector<Eigen::Vector2d> chosen_sliding_directions;
	std::vector<RollingFrictionRedundancyDir> rolling_redundancy_directions;
	std::vector<MomentConstraints> initial_moment_constraints;
	std::vector<MomentConstraints> moment_constraints;
};

}

#endif //LCP_SOLVER2_H