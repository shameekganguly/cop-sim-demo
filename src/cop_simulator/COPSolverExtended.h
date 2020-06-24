// COPSolverExtended.h

//TODO: remove the old stub class and replace with this class. rename this to COPSolver

#ifndef COP_SOLVER_EXTENDED_H
#define COP_SOLVER_EXTENDED_H

#include "COPSolver.h"
#include "geometry/Primitive.h"

namespace Sai2COPSim {

// solver class will be instantiated at each call
// TODO: consider whether to track the COP position in body or in world frame
// for sliding, it makes sense to track in body frame
// But for slip rotations, it might be better to track in world frame?
// e.g. cylinder rolling on its side at one place, on a frictionless surface
// NOTE: between two calls to the COP solver, the contact geometry doesn't change.
// it only changes when the distances are recomputed.
// the last cop solution will still be valid, in case there is no new collision to resolve

// NOTE: this solver still assumes a polygonal contact patch geometry
// TODO: extend to more general boundary surfaces set
// The reference frame is assumed to be point0. All the quantities are defined at point0.
class COPSolver {
public:
	COPSolver(){ }

	// solve simultaneously for n contact patches to obtain COP position and forces
	// NOTE: contact coordinate is assumed to be Body B displacement - Body A displacement
	ContactCOPSolution solveStartWithLastCOP(
		double friction_coeff,
		const Eigen::MatrixXd& A_constraint, // defined at point0
		const Eigen::VectorXd& rhs_constraint, // defined at point 0. Includes Jdot_qdot
		const std::vector<std::vector<Eigen::Vector3d>>& boundary_points, // points are in the COP frame with point0 being origin
		const std::vector<uint>& patch_indices, // ith entry gives starting row of ith contact patch in A_constraint, rhs_constraint
		const std::vector<ContactType>& contact_types, // ith entry gives contact type of ith contact patch
		const std::vector<Eigen::Vector3d>& omega_bodyA, // vector of body A angular velocities, in COP frame, one entry per contact patch
		const std::vector<Eigen::Vector3d>& omega_bodyB, // vector of body B angular velocities, in COP frame, one entry per contact patch
		const std::vector<Eigen::Vector3d>& linear_contact_velocity, // 3 dof relative translation velocity at point 0, in COP frame, one entry per contact patch
		//^ expected to be zero in the z direction, but we don't explicitly check
		const std::vector<Eigen::Vector3d>& last_cop_point, // position vector from point0 to last_COP_sol point, in COP frame, one entry per contact patch
		//^ TODO: is this required? this is already in last_COP_sol
		ContactCOPSolution last_COP_sol // last COP solution, in the current COP frame
	);

	// solve simultaneously for n contact patches to obtain COP position and forces
	// NOTE: contact coordinate is assumed to be Body B displacement - Body A displacement
	// TODO: cache cenroids? But they would change with change in active pairs, active points
	// TODO: in the general case where we don't have a polyhedral geometry, consider using
	// centroid of the collision points
	// TODO: in general case, how do we estimate the rolling friction coefficient? we 
	// will need a distance estimator to the patch boundary
	ContactCOPSolution solveStartWithPatchCentroid(
		double friction_coeff,
		const Eigen::MatrixXd& A_constraint, // defined at point0
		const Eigen::VectorXd& rhs_constraint, // defined at point 0. Includes Jdot_qdot
		const std::vector<std::vector<Eigen::Vector3d>>& boundary_points, // points are in the COP frame with point0 being origin
		const std::vector<uint>& patch_indices, // ith entry gives starting row of ith contact patch in A_constraint, rhs_constraint
		const std::vector<ContactType>& contact_types, // ith entry gives contact type of ith contact patch
		const std::vector<Eigen::Vector3d>& omega_bodyA, // vector of body A angular velocities, in COP frame, one entry per contact patch
		const std::vector<Eigen::Vector3d>& omega_bodyB, // vector of body B angular velocities, in COP frame, one entry per contact patch
		const std::vector<Eigen::Vector3d>& linear_contact_velocity // 3 dof relative translation velocity at point 0, in COP frame, one entry per contact patch
		//^ expected to be zero in the z direction, but we don't explicitly check
	);

public:
	// saved initialization variables

public:
	// internal variables
};

}


#endif //