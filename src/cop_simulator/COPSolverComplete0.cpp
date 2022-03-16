// COPSolverComplete0.cpp

#include "COPSolverExtended.h"
#include <iostream>

/*  THIS IS WIP CODE */

using namespace Eigen;

namespace Sai2COPSim {

static void getCOPLineMultipleContactsDisplacedMatrices(
	Eigen::MatrixXd& A_disp,
	Eigen::VectorXd& rhs_disp,
	Eigen::Vector3d& lin_vel_disp,
	double signed_dist,
	const Eigen::MatrixXd& A,
	const Eigen::VectorXd& rhs,
	const Eigen::Vector3d& line_omega_bodyA,
	const Eigen::Vector3d& line_omega_bodyB,
	const Eigen::Vector3d& line_lin_vel,
	uint line_start_row_id
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
	A_disp = A; //TODO: since we are using delta changes in the signed distance, we might be able to remove this to save some time
	// col blocks
	A_disp.block(0, line_start_row_id, A.rows(), 5) = A.block(0, line_start_row_id, A.rows(), 5) * cross_mat.transpose();
	// row blocks
	A_disp.block(line_start_row_id, 0, 5, A.cols()) = cross_mat * A_disp.block(line_start_row_id, 0, 5, A.cols());

	// compute RHS at new position
	rhs_disp = rhs; //TODO: since we are using delta changes in the signed distance, we might be able to remove this to save some time
	rhs_disp.segment<5>(line_start_row_id) = cross_mat*rhs.segment<5>(line_start_row_id);

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


MultipleContactCOPSolution COPSolver::solveLinesandPointsStartWithPatchCentroid(
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
) {
	const uint total_dof = rhs_constraint.size();
	const uint num_contacts = patch_indices.size();

	// TODO: include contact velocity-active check here, automatically force to point/line contact
	// TODO: handle redundancy

	// reset state
	contacts_translational_friction_states.clear(); // TODO: reuse to save computation time
	contacts_rotational_friction_states.clear();
	contacts_translation_rolling_dir.clear();
	contacts_translation_slip_dir.clear();
	contacts_rotational_rolling_dir.clear();
	contacts_rotational_slip_dir.clear();
	contacts_active.clear();

	Tstart_rowids.resize(num_contacts);
	Tnum_rows.resize(num_contacts);
	TA = A_constraint;
	Trhs = -rhs_constraint;
	TA_size = 0;

	vector<double> signed_distances;
	vector<double> max_distances;
	vector<double> min_distances;

	MultipleContactCOPSolution ret_sol;

	// --- initializations ---
	MatrixXd A_disp = A_constraint;
	VectorXd rhs_disp = rhs_constraint;
	vector<Vector3d> line_vel_disp;
	for(uint i = 0; i < num_contacts; i++) {
		// set contact to inactive initially
		contacts_active.push_back(false);

		// add contact cop sol to ret_sol
		ret_sol.contact_solutions.push_back(ContactCOPSolution());

		// initial computation of displaced line velocity
		line_vel_disp.push_back(linear_contact_velocity[i]);

		if(contact_types[i] == ContactType::LINE) {
			// line contact: compute line centroids
			Vector3d patch_centroid_point = boundary_points[i][1]/2.0;
			double signed_distance_from_point0 = patch_centroid_point.dot(Vector3d(1, 0, 0));
			signed_distances.push_back(signed_distance_from_point0);

			// line contact: compute min and max signed distance from point 0
			double max_signed_dist = 0.0;
			double min_signed_dist = 0.0;
			if(signed_distance_from_point0 > 0) {
				max_signed_dist = 2*signed_distance_from_point0;
			} else {
				min_signed_dist = 2*signed_distance_from_point0;
			}
			max_distances.push_back(max_signed_dist);
			min_distances.push_back(min_signed_dist);

			// line contact: get first displaced matrix
			// NOTE: we use relative displacements since we need to accomodate changing the
			// displacements for each line contact, one at a time
			getCOPLineMultipleContactsDisplacedMatrices(
				A_disp,
				rhs_disp,
				line_vel_disp[i],
				signed_distance_from_point0,
				A_disp,
				rhs_disp,
				omega_bodyA[i],
				omega_bodyB[i],
				line_vel_disp[i],
				patch_indices[i]
			);

			// line contact: initialize return solution cop type
			ret_sol.contact_solutions[i].cop_type = COPContactType::LineCenter;

			// line contact: compute rotation friction state
			contacts_rotational_rolling_dir.push_back(0);
			contacts_rotational_slip_dir.push_back(omega_bodyB[i](2) - omega_bodyA[i](2));
			if(abs(contacts_rotational_slip_dir[i]) > 1e-8) {
				contacts_rotational_friction_states.push_back(FrictionState::Sliding);
				contacts_rotational_slip_dir[i] /= abs(contacts_rotational_slip_dir[i]);
			} else {
				// initialize to rolling
				contacts_rotational_friction_states.push_back(FrictionState::Rolling);
			}

			// add placeholder for translation state and slip
			// gets overwritten at the beginning of the solver loop
			contacts_translational_friction_states.push_back(FrictionState::Frictionless);
			contacts_translation_rolling_dir.push_back(Vector2d::Zero());
			contacts_translation_slip_dir.push_back(Vector2d::Zero());

			// initialize return COP acceleration with rhs values
			ret_sol.contact_solutions[i].acc_sol = rhs_disp.segment<5>(patch_indices[i]);

		} else if(contact_types[i] == ContactType::POINT) {
			// point contact: compute translation friction state
			contacts_translation_slip_dir.push_back(linear_contact_velocity[i].segment<2>(0));
			if (contacts_translation_slip_dir[i].norm() > 1e-6) {
				contacts_translational_friction_states.push_back(FrictionState::Sliding);
				contacts_translation_slip_dir[i] /= contacts_translation_slip_dir[i].norm();
			} else {
				// initialize to rolling
				contacts_translational_friction_states.push_back(FrictionState::Rolling);
			}

			// add vector fillers
			contacts_rotational_friction_states.push_back(FrictionState::Frictionless);
			contacts_rotational_rolling_dir.push_back(0);
			contacts_rotational_slip_dir.push_back(0);
			signed_distances.push_back(0);
			max_distances.push_back(0);
			min_distances.push_back(0);

			// initialize return COP acceleration with rhs values
			ret_sol.contact_solutions[i].acc_sol = rhs_disp.segment<3>(patch_indices[i]);
		}
		// TODO: else throw exception
	}

	// --- solver loop ---
	// loop variables
	VectorXd full_a_sol(total_dof), full_f_sol(total_dof);
	vector<bool> fForceIgnoreSlip(num_contacts, false);
	vector<bool> fForceLineCenter(num_contacts, false);
	vector<bool> fdidUpdateCOPPos(num_contacts, false);
	const int max_iters = 25*num_contacts;
	int iter_ind = 0;
	double mu_rot;
	const double mu = friction_coeff;
	COPSolResult ret_result;

	while(iter_ind < max_iters) {
		iter_ind++;
		// line contacts: set translation friction state and dir
		for(uint i = 0; i < num_contacts; i++) {
			if(fdidUpdateCOPPos[i]) { // only ever set to true for line contacts
				fdidUpdateCOPPos[i] = false;

				// update cop type
				double dist1 = abs(signed_distances[i]);
				double dist2 = abs(signed_distances[i] - boundary_points[i][1](0));
				// check whether we are on the boundary or inside the contact segment
				if(fForceLineCenter[i]) {
					ret_sol.contact_solutions[i].cop_type = COPContactType::LineCenter;
				} else {
					if(dist1 < 1e-15 || dist2 < 1e-15) {
						ret_sol.contact_solutions[i].cop_type = COPContactType::LineEnd;
					} else {
						ret_sol.contact_solutions[i].cop_type = COPContactType::LineCenter;
					}
				}

				// update slip state
				contacts_translation_slip_dir[i] = line_vel_disp[i].segment<2>(0);
				if (!fForceIgnoreSlip[i] && contacts_translation_slip_dir[i].norm() > 1e-6) {
					contacts_translational_friction_states[i] = FrictionState::Sliding;
					contacts_translation_slip_dir[i] /= contacts_translation_slip_dir[i].norm();
				} else {
					contacts_translational_friction_states[i] = FrictionState::Rolling;
					// if rolling at line-COP, set fForceIgnoreSlip true
					fForceIgnoreSlip[i] = true;
				}
			}
		}
		// assemble matrices
		computeMatricesForMultipleContacts(
			A_disp,
			rhs_disp,
			mu,
			mu_rot,
			ret_sol.contact_solutions,
			patch_indices,
			contact_types
		);

		// solve
		solveInternalForMultipleContacts(
			mu,
			mu_rot,
			A_disp,
			rhs_disp,
			full_f_sol,
			full_a_sol,
			ret_sol.contact_solutions,
			patch_indices,
			contact_types
		);

		// check for constraint violations
		for(uint i = 0; i < num_contacts; i++) {


			// line contacts: check for zero moment constraint violations
			if(contact_types[i] == ContactType::LINE) {
				if(ret_sol.contact_solutions[i].cop_type == COPContactType::LineCenter) {
					if(abs(full_f_sol[patch_indices[i]+3]) > 1e-3) {
						// std::cout << "Full F sol line: " << full_f_sol.segment<5>(line_start_row_id).transpose() << std::endl;
						// - - if so, compute distance to new line-COP.
						double dist_move_cop = -full_f_sol[patch_indices[i]+3]/fmax(full_f_sol[patch_indices[i]+2], 1e-5);
						// - - if distance to new line-COP is very small, continue to check other violations
						double last_signed_dist = signed_distances[i];
						if(abs(dist_move_cop) > 0.02) {// TODO: integrate bisection search and lower this threshold
							signed_distances[i] += dist_move_cop;
							signed_distances[i] = fmin(signed_distances[i], max_distances[i]);
							signed_distances[i] = fmax(signed_distances[i], min_distances[i]);
							fdidUpdateCOPPos[i] = true;
							// std::cout << "New cop pos " << signed_distance_from_point0 << std::endl;
							// std::cout << "Update dist " << dist_move_cop << std::endl;
							// std::cout << "Max dist " << max_signed_dist << std::endl;
							// std::cout << "Min dist " << min_signed_dist << std::endl;
							//TODO: if we have fForceLineCenter set to true, consider setting it false

							// - get delta change - displaced matrices
							dist_move_cop = signed_distance_from_point0 - last_signed_dist;
							getCOPLineMultipleContactsDisplacedMatrices(
								A_disp,
								rhs_disp,
								line_vel_disp[i],
								dist_move_cop,
								A_disp,
								rhs_disp,
								omega_bodyA[i],
								omega_bodyB[i],
								line_vel_disp[i],
								patch_indices[i]
							);
							continue;
						}
						// - - if distance is small, check first if we are currently doing a binary search
						// - - if so, continue binary search with point state fixed.
						// - - if not, check last line-COP moment. if we switch moment directions, start binary 
						// - - search with point state fixed.
					}
				}
			}
		}

		// if still here, we were successful in solving
		ret_result = COPSolResult::Success;
		ret_sol.full_force_sol = full_f_sol;
		break;
	}

	if(iter_ind >= max_iters) {
		std::cerr << "solveLinesandPointsStartWithPatchCentroid failed to converge" << std::endl;
		ret_result = COPSolResult::NoSolution;
	}
	// when returning result, set result type to all 
	ret_sol.result = ret_result;
	for(uint i = 0; i < num_contacts; i++) {
		ret_sol.contact_solutions[i].result = ret_sol.result;
	}
	return ret_sol;
}

inline void getPointPointMats(
	MatrixXd& ptA_row_ptB_col_block, // note that the right size is expected to already be set
	MatrixXd& ptA_col_ptB_row_block, // note that the right size is expected to already be set
	double mu,
	FrictionState ptA_state,
	const Vector2d& ptA_slip_dir, // ignored if pt is rolling
	FrictionState ptB_state,
	const Vector2d& ptB_slip_dir, // ignored if pt is rolling
	const MatrixXd& A_block_ptA_row_ptB_col,
	const MatrixXd& A_block_ptA_col_ptB_row
) {
	if(ptA_state == FrictionState::Sliding ||
		ptA_state == FrictionState::Impending
	) {
		if(ptB_state == FrictionState::Sliding ||
			ptB_state == FrictionState::Impending
		) {
			ptA_row_ptB_col_block = A_block_ptA_row_ptB_col(2,2)
									- mu*A_block_ptA_row_ptB_col.block<1,2>(2,0)*ptB_slip_dir;

			ptA_col_ptB_row_block = A_block_ptA_col_ptB_row(2,2)
									- mu*A_block_ptA_col_ptB_row.block<1,2>(2,0)*ptA_slip_dir;
		} else {
			ptA_row_ptB_col_block = A_block_ptA_row_ptB_col.block<1,3>(2,0);

			ptA_col_ptB_row_block = A_block_ptA_col_ptB_row.block<3,1>(0,2)
									- mu*A_block_ptA_col_ptB_row.block<3,2>(0,0)*ptA_slip_dir;
		}
	} else {
		if(ptB_state == FrictionState::Sliding ||
			ptB_state == FrictionState::Impending
		) {
			ptA_row_ptB_col_block = A_block_ptA_row_ptB_col.block<3,1>(0,2)
									- mu*A_block_ptA_row_ptB_col.block<3,2>(0,0)*ptB_slip_dir;

			ptA_col_ptB_row_block = A_block_ptA_col_ptB_row.block<1,3>(2,0);
		} else {
			ptA_row_ptB_col_block = A_block_ptA_row_ptB_col;
			ptA_col_ptB_row_block = A_block_ptA_col_ptB_row;
		}
	}
}

inline void getPointLineCenterMats( // Line assumed to be with COP at center, else call pt pt
	MatrixXd& pt_row_line_col_block, // note that the right size is expected to already be set
	MatrixXd& pt_col_line_row_block, // note that the right size is expected to already be set
	double mu,
	double mu_rot,
	FrictionState pt_state,
	const Vector2d& pt_slip_dir, // ignored if pt is rolling
	FrictionState line_translation_state,
	const Vector2d& line_trans_slip_dir, // ignored if line translation is rolling
	FrictionState line_rotation_state,
	double line_rot_slip_dir, // ignored if line rotation is rolling
	const MatrixXd& A_block_pt_row_line_col,
	const MatrixXd& A_block_pt_col_line_row
) {
	if(pt_state == FrictionState::Sliding ||
		pt_state == FrictionState::Impending
	) {
		if(line_translation_state == FrictionState::Sliding ||
			line_translation_state == FrictionState::Impending
		) { // line translation sliding
			if(line_rotation_state == FrictionState::Sliding ||
				line_rotation_state == FrictionState::Impending
			) {
				pt_row_line_col_block(0,0) = A_block_pt_row_line_col(2,2)
												- mu*A_block_pt_row_line_col.block<1,2>(2,0)*line_trans_slip_dir
												- mu_rot*A_block_pt_row_line_col(2,4)*line_rot_slip_dir;
				pt_row_line_col_block(0,1) = A_block_pt_row_line_col(2,3);

				pt_col_line_row_block = A_block_pt_col_line_row.block<2,1>(2,2)
												- mu*A_block_pt_col_line_row.block<2,2>(2,0)*pt_slip_dir;
			} else { // line rotation rolling
				pt_row_line_col_block(0,0) = A_block_pt_row_line_col(2,2)
												- mu*A_block_pt_row_line_col.block<1,2>(2,0)*line_trans_slip_dir;
				pt_row_line_col_block.block<1,2>(0,1) = A_block_pt_row_line_col.block<1,2>(2,3);

				pt_col_line_row_block = A_block_pt_col_line_row.block<3,1>(2,2)
												- mu*A_block_pt_col_line_row.block<3,2>(2,0)*pt_slip_dir;
			}
		} else { // line translation rolling
			if(line_rotation_state == FrictionState::Sliding ||
				line_rotation_state == FrictionState::Impending
			) {
				pt_row_line_col_block.block<1,2>(0,0) = A_block_pt_row_line_col.block<1,2>(2,0);
				pt_row_line_col_block(0,2) = A_block_pt_row_line_col(2,2)
												- mu_rot*A_block_pt_row_line_col(2,4)*line_rot_slip_dir;
				pt_row_line_col_block(0,3) = A_block_pt_row_line_col(2,3);

				pt_col_line_row_block = A_block_pt_col_line_row.block<4,1>(0,2)
												- mu*A_block_pt_col_line_row.block<4,2>(0,0)*pt_slip_dir;
			} else { // line rotation rolling
				pt_row_line_col_block = A_block_pt_row_line_col.block<1,5>(2,0);

				pt_col_line_row_block = A_block_pt_col_line_row.block<5,1>(0,2)
												- mu*A_block_pt_col_line_row.block<5,2>(0,0)*pt_slip_dir;
			}
		}
	} else { // pt state rolling
		if(line_translation_state == FrictionState::Sliding ||
			line_translation_state == FrictionState::Impending
		) { // line translation sliding
			if(line_rotation_state == FrictionState::Sliding ||
				line_rotation_state == FrictionState::Impending
			) {
				pt_row_line_col_block.block<3,1>(0,0) = A_block_pt_row_line_col.block<3,1>(0,2)
												- mu*A_block_pt_row_line_col.block<3,2>(0,0)*line_trans_slip_dir
												- mu_rot*A_block_pt_row_line_col.block<3,1>(0,4)*line_rot_slip_dir;
				pt_row_line_col_block.block<3,1>(0,1) = A_block_pt_row_line_col.block<3,1>(0,3);

				pt_col_line_row_block = A_block_pt_col_line_row.block<2,3>(2,0);
			} else { // line rotation rolling
				pt_row_line_col_block.block<3,1>(0,0) = A_block_pt_row_line_col.block<3,1>(0,2)
												- mu*A_block_pt_row_line_col.block<3,2>(0,0)*line_trans_slip_dir;
				pt_row_line_col_block.block<3,2>(0,1) = A_block_pt_row_line_col.block<3,2>(0,3);

				pt_col_line_row_block = A_block_pt_col_line_row.block<3,3>(2,0);
			}
		} else {// line translation rolling
			if(line_rotation_state == FrictionState::Sliding ||
				line_rotation_state == FrictionState::Impending
			) {
				pt_row_line_col_block.block<3,2>(0,0) = A_block_pt_row_line_col.block<3,2>(0,0);
				pt_row_line_col_block.block<3,1>(0,2) = A_block_pt_row_line_col.block<3,1>(0,2)
												- mu_rot*A_block_pt_row_line_col.block<3,1>(0,4)*line_rot_slip_dir;
				pt_row_line_col_block.block<3,1>(0,3) = A_block_pt_row_line_col.block<3,1>(0,3);

				pt_col_line_row_block = A_block_pt_col_line_row.block<4,3>(0,0);
			} else { // line rotation rolling
				pt_row_line_col_block = A_block_pt_row_line_col;

				pt_col_line_row_block = A_block_pt_col_line_row;
			}
		}
	}
}

inline void getLineCenterLineCenterMats( // Line assumed to be with COP at center, else call pt pt
	MatrixXd& lineA_row_lineB_col_block, // note that the right size is expected to already be set
	MatrixXd& lineA_col_lineB_row_block, // note that the right size is expected to already be set
	double mu,
	double mu_rot,
	FrictionState lineA_translation_state,
	const Vector2d& lineA_trans_slip_dir, // ignored if line translation is rolling
	FrictionState lineA_rotation_state,
	double lineA_rot_slip_dir, // ignored if line rotation is rolling
	FrictionState lineB_translation_state,
	const Vector2d& lineB_trans_slip_dir, // ignored if line translation is rolling
	FrictionState lineB_rotation_state,
	double lineB_rot_slip_dir, // ignored if line rotation is rolling
	const MatrixXd& A_block_lineA_row_lineB_col,
	const MatrixXd& A_block_lineA_col_lineB_row
) {
	if(lineA_translation_state == FrictionState::Sliding ||
		lineA_translation_state == FrictionState::Impending
	) { // lineA translation sliding
		if(lineA_rotation_state == FrictionState::Sliding ||
			lineA_rotation_state == FrictionState::Impending
		) { // lineA rotation sliding
			if(lineB_translation_state == FrictionState::Sliding ||
				lineB_translation_state == FrictionState::Impending
			) { // lineB translation sliding
				if(lineB_rotation_state == FrictionState::Sliding ||
					lineB_rotation_state == FrictionState::Impending
				) { // lineB rotation sliding
					lineA_col_lineB_row_block.block<2,1>(0,0) = A_block_lineA_col_lineB_row.block<2,1>(2,2)
												- mu*A_block_lineA_col_lineB_row.block<2,2>(2,0)*lineA_trans_slip_dir
												- mu_rot*A_block_lineA_col_lineB_row.block<2,1>(2,4)*lineA_rot_slip_dir;
					lineA_col_lineB_row_block.block<2,1>(0,1) = A_block_lineA_col_lineB_row.block<2,1>(2,3);

					lineA_row_lineB_col_block.block<2,1>(0,0) = A_block_lineA_row_lineB_col.block<2,1>(2,2)
													- mu*A_block_lineA_row_lineB_col.block<2,2>(2,0)*lineB_trans_slip_dir
													- mu_rot*A_block_lineA_row_lineB_col.block<2,1>(2,4)*lineB_rot_slip_dir;
					lineA_row_lineB_col_block.block<2,1>(0,1) = A_block_lineA_row_lineB_col.block<2,1>(2,3);
				} else { // lineB rotation rolling
					lineA_col_lineB_row_block.block<3,1>(0,0) = A_block_lineA_col_lineB_row.block<3,1>(2,2)
												- mu*A_block_lineA_col_lineB_row.block<3,2>(2,0)*lineA_trans_slip_dir
												- mu_rot*A_block_lineA_col_lineB_row.block<3,1>(2,4)*lineA_rot_slip_dir;
					lineA_col_lineB_row_block.block<3,1>(0,1) = A_block_lineA_col_lineB_row.block<3,1>(2,3);

					lineA_row_lineB_col_block.block<2,1>(0,0) = A_block_lineA_row_lineB_col.block<2,1>(2,2)
													- mu*A_block_lineA_row_lineB_col.block<2,2>(2,0)*lineB_trans_slip_dir;
					lineA_row_lineB_col_block.block<2,2>(0,1) = A_block_lineA_row_lineB_col.block<2,2>(2,3);
				}
			} else { // lineB translation rolling
				if(lineA_rotation_state == FrictionState::Sliding ||
					lineA_rotation_state == FrictionState::Impending
				) { // lineB rotation sliding
					lineA_col_lineB_row_block.block<4,1>(0,0) = A_block_lineA_col_lineB_row.block<4,1>(0,2)
							- mu*A_block_lineA_col_lineB_row.block<4,2>(0,0)*lineA_trans_slip_dir
							- mu_rot*A_block_lineA_col_lineB_row.block<4,1>(0,4)*lineA_rot_slip_dir;
					lineA_col_lineB_row_block.block<4,1>(0,1) = A_block_lineA_col_lineB_row.block<4,1>(0,3);

					lineA_row_lineB_col_block.block<2,2>(0,0) = A_block_lineA_row_lineB_col.block<2,2>(2,0);
					lineA_row_lineB_col_block.block<2,1>(0,2) = A_block_lineA_row_lineB_col.block<2,1>(2,2)
													- mu_rot*A_block_lineA_row_lineB_col.block<2,1>(2,4)*lineB_rot_slip_dir;
					lineA_row_lineB_col_block.block<2,1>(0,3) = A_block_lineA_row_lineB_col.block<2,1>(2,3);
				} else { // lineB rotation rolling
					lineA_col_lineB_row_block.block<5,1>(0,0) = A_block_lineA_col_lineB_row.block<5,1>(0,2)
							- mu*A_block_lineA_col_lineB_row.block<5,2>(0,0)*lineA_trans_slip_dir
							- mu_rot*A_block_lineA_col_lineB_row.block<5,1>(0,4)*lineA_rot_slip_dir;
					lineA_col_lineB_row_block.block<5,1>(0,1) = A_block_lineA_col_lineB_row.block<5,1>(0,3);

					lineA_row_lineB_col_block = A_block_lineA_row_lineB_col.block<2,5>(2,0);
				}
			}
		} else { // lineA rotation rolling
			if(lineB_translation_state == FrictionState::Sliding ||
				lineB_translation_state == FrictionState::Impending
			) { // lineB translation sliding
				if(lineB_rotation_state == FrictionState::Sliding ||
					lineB_rotation_state == FrictionState::Impending
				) { // lineB rotation sliding
					lineA_col_lineB_row_block.block<2,1>(0,0) = A_block_lineA_col_lineB_row.block<2,1>(2,2)
												- mu*A_block_lineA_col_lineB_row.block<2,2>(2,0)*lineA_trans_slip_dir;
					lineA_col_lineB_row_block.block<2,2>(0,1) = A_block_lineA_col_lineB_row.block<2,2>(2,3);

					lineA_row_lineB_col_block.block<3,1>(0,0) = A_block_lineA_row_lineB_col.block<3,1>(2,2)
													- mu*A_block_lineA_row_lineB_col.block<3,2>(2,0)*lineB_trans_slip_dir
													- mu_rot*A_block_lineA_row_lineB_col.block<3,1>(2,4)*lineB_rot_slip_dir;
					lineA_row_lineB_col_block.block<3,1>(0,1) = A_block_lineA_row_lineB_col.block<3,1>(2,3);
				} else { // lineB rotation rolling
					lineA_col_lineB_row_block.block<3,1>(0,0) = A_block_lineA_col_lineB_row.block<3,1>(2,2)
												- mu*A_block_lineA_col_lineB_row.block<3,2>(2,0)*lineA_trans_slip_dir;
					lineA_col_lineB_row_block.block<3,2>(0,1) = A_block_lineA_col_lineB_row.block<3,2>(2,3);

					lineA_row_lineB_col_block.block<3,1>(0,0) = A_block_lineA_row_lineB_col.block<3,1>(2,2)
													- mu*A_block_lineA_row_lineB_col.block<3,2>(2,0)*lineB_trans_slip_dir;
					lineA_row_lineB_col_block.block<3,2>(0,1) = A_block_lineA_row_lineB_col.block<3,2>(2,3);
				}
			} else { // lineB translation rolling
				if(lineB_rotation_state == FrictionState::Sliding ||
					lineB_rotation_state == FrictionState::Impending
				) { // lineB rotation sliding
					lineA_col_lineB_row_block.block<4,1>(0,0) = A_block_lineA_col_lineB_row.block<4,1>(0,2)
												- mu*A_block_lineA_col_lineB_row.block<4,2>(0,0)*lineA_trans_slip_dir;
					lineA_col_lineB_row_block.block<4,2>(0,1) = A_block_lineA_col_lineB_row.block<4,2>(0,3);

					lineA_row_lineB_col_block.block<3,2>(0,0) = A_block_lineA_row_lineB_col.block<3,2>(2,0);
					lineA_row_lineB_col_block.block<3,1>(0,2) = A_block_lineA_row_lineB_col.block<3,1>(2,2)
													- mu_rot*A_block_lineA_row_lineB_col.block<3,1>(2,4)*lineB_rot_slip_dir;
					lineA_row_lineB_col_block.block<3,1>(0,3) = A_block_lineA_row_lineB_col.block<3,1>(2,3);
				} else { // lineB rotation rolling
					lineA_col_lineB_row_block.block<5,1>(0,0) = A_block_lineA_col_lineB_row.block<5,1>(0,2)
												- mu*A_block_lineA_col_lineB_row.block<5,2>(0,0)*lineA_trans_slip_dir;
					lineA_col_lineB_row_block.block<5,2>(0,1) = A_block_lineA_col_lineB_row.block<5,2>(0,3);

					lineA_row_lineB_col_block = A_block_lineA_row_lineB_col.block<3,5>(2,0);
				}
			}
		}
	} else { // lineA translation rolling
		if(lineA_rotation_state == FrictionState::Sliding ||
			lineA_rotation_state == FrictionState::Impending
		) { // lineA rotation sliding
			if(lineB_translation_state == FrictionState::Sliding ||
				lineB_translation_state == FrictionState::Impending
			) { // lineB translation sliding
				if(lineB_rotation_state == FrictionState::Sliding ||
					lineB_rotation_state == FrictionState::Impending
				) { // lineB rotation sliding
					lineA_col_lineB_row_block.block<2,2>(0,0) = A_block_lineA_col_lineB_row.block<2,2>(2,0);
					lineA_col_lineB_row_block.block<2,1>(0,2) = A_block_lineA_col_lineB_row.block<2,1>(2,2)
												- mu_rot*A_block_lineA_col_lineB_row.block<2,1>(2,4)*lineA_rot_slip_dir;
					lineA_col_lineB_row_block.block<2,1>(0,3) = A_block_lineA_col_lineB_row.block<2,1>(2,3);

					lineA_row_lineB_col_block.block<4,1>(0,0) = A_block_lineA_row_lineB_col.block<4,1>(0,2)
													- mu*A_block_lineA_row_lineB_col.block<4,2>(0,0)*lineB_trans_slip_dir
													- mu_rot*A_block_lineA_row_lineB_col.block<4,1>(0,4)*lineB_rot_slip_dir;
					lineA_row_lineB_col_block.block<4,1>(0,1) = A_block_lineA_row_lineB_col.block<4,1>(0,3);
				} else { // lineB rotation rolling
					lineA_col_lineB_row_block.block<3,2>(0,0) = A_block_lineA_col_lineB_row.block<3,2>(2,0);
					lineA_col_lineB_row_block.block<3,1>(0,2) = A_block_lineA_col_lineB_row.block<3,1>(2,2)
												- mu_rot*A_block_lineA_col_lineB_row.block<3,1>(2,4)*lineA_rot_slip_dir;
					lineA_col_lineB_row_block.block<3,1>(0,3) = A_block_lineA_col_lineB_row.block<3,1>(2,3);

					lineA_row_lineB_col_block.block<4,1>(0,0) = A_block_lineA_row_lineB_col.block<4,1>(0,2)
													- mu*A_block_lineA_row_lineB_col.block<4,2>(0,0)*lineB_trans_slip_dir;
					lineA_row_lineB_col_block.block<4,2>(0,1) = A_block_lineA_row_lineB_col.block<4,2>(0,3);
				}
			} else { // lineB translation rolling
				if(lineB_rotation_state == FrictionState::Sliding ||
					lineB_rotation_state == FrictionState::Impending
				) { // lineB rotation sliding
					lineA_col_lineB_row_block.block<4,2>(0,0) = A_block_lineA_col_lineB_row.block<4,2>(0,0);
					lineA_col_lineB_row_block.block<4,1>(0,2) = A_block_lineA_col_lineB_row.block<4,1>(0,2)
												- mu_rot*A_block_lineA_col_lineB_row.block<4,1>(0,4)*lineA_rot_slip_dir;
					lineA_col_lineB_row_block.block<4,1>(0,3) = A_block_lineA_col_lineB_row.block<4,1>(0,3);

					lineA_row_lineB_col_block.block<4,2>(0,0) = A_block_lineA_row_lineB_col.block<4,2>(0,0);
					lineA_row_lineB_col_block.block<4,1>(0,2) = A_block_lineA_row_lineB_col.block<4,1>(0,2)
													- mu_rot*A_block_lineA_row_lineB_col.block<4,1>(0,4)*lineB_rot_slip_dir;
					lineA_row_lineB_col_block.block<4,1>(0,3) = A_block_lineA_row_lineB_col.block<4,1>(0,3);
				} else { // lineB rotation rolling
					lineA_col_lineB_row_block.block<5,2>(0,0) = A_block_lineA_col_lineB_row.block<5,2>(0,0);
					lineA_col_lineB_row_block.block<5,1>(0,2) = A_block_lineA_col_lineB_row.block<5,1>(0,2)
												- mu_rot*A_block_lineA_col_lineB_row.block<5,1>(0,4)*lineA_rot_slip_dir;
					lineA_col_lineB_row_block.block<5,1>(0,3) = A_block_lineA_col_lineB_row.block<5,1>(0,3);

					lineA_row_lineB_col_block = A_block_lineA_row_lineB_col.block<4,5>(0,0);
				}
			}
		} else { // lineA rotation rolling
			if(lineB_translation_state == FrictionState::Sliding ||
				lineB_translation_state == FrictionState::Impending
			) { // lineB translation sliding
				if(lineB_rotation_state == FrictionState::Sliding ||
					lineB_rotation_state == FrictionState::Impending
				) { // lineB rotation sliding
					lineA_col_lineB_row_block.block<2,5>(0,0) = A_block_lineA_col_lineB_row.block<2,5>(2,0);

					lineA_row_lineB_col_block.block<5,1>(0,0) = A_block_lineA_row_lineB_col.block<5,1>(0,2)
													- mu*A_block_lineA_row_lineB_col.block<5,2>(0,0)*lineB_trans_slip_dir
													- mu_rot*A_block_lineA_row_lineB_col.block<5,1>(0,4)*lineB_rot_slip_dir;
					lineA_row_lineB_col_block.block<5,1>(0,1) = A_block_lineA_row_lineB_col.block<5,1>(0,3);
				} else { // lineB rotation rolling
					lineA_col_lineB_row_block.block<3,5>(0,0) = A_block_lineA_col_lineB_row.block<3,5>(2,0);

					lineA_row_lineB_col_block.block<5,1>(0,0) = A_block_lineA_row_lineB_col.block<5,1>(0,2)
													- mu*A_block_lineA_row_lineB_col.block<5,2>(0,0)*lineB_trans_slip_dir;
					lineA_row_lineB_col_block.block<5,2>(0,1) = A_block_lineA_row_lineB_col.block<5,2>(0,3);
				}
			} else { // lineB translation rolling
				if(lineB_rotation_state == FrictionState::Sliding ||
					lineB_rotation_state == FrictionState::Impending
				) { // lineB rotation sliding
					lineA_col_lineB_row_block.block<4,5>(0,0) = A_block_lineA_col_lineB_row.block<4,5>(0,0);

					lineA_row_lineB_col_block.block<5,2>(0,0) = A_block_lineA_row_lineB_col.block<5,2>(0,0);
					lineA_row_lineB_col_block.block<5,1>(0,2) = A_block_lineA_row_lineB_col.block<5,1>(0,2)
													- mu_rot*A_block_lineA_row_lineB_col.block<5,1>(0,4)*lineB_rot_slip_dir;
					lineA_row_lineB_col_block.block<5,1>(0,3) = A_block_lineA_row_lineB_col.block<5,1>(0,3);
				} else { // lineB rotation rolling
					lineA_col_lineB_row_block = A_block_lineA_col_lineB_row.block<5,5>(0,0);
					lineA_row_lineB_col_block = A_block_lineA_row_lineB_col.block<5,5>(0,0);
				}
			}
		}
	}
}

void COPSolver::computeMatricesForMultipleContacts(
	const Eigen::MatrixXd& A_disp,
	const Eigen::VectorXd& rhs_disp,
	double mu,
	double mu_rot,
	const std::vector<ContactCOPSolution>& ret_cop_sols, //NOTE: the contact type must be set here
	const std::vector<uint>& patch_indices,
	const std::vector<ContactType>& contact_types
) {
	TA_size = 0;

	for(uint i = 0; i < patch_indices.size(); i++) {
		// skip inactive contacts
		if(!contacts_active[i]) {
			Tstart_rowids[i] = -1;
			continue;
		}

		uint trid = TA_size;
		Tstart_rowids[i] = trid;
		uint tnum_rows = 0;
		uint full_rid = patch_indices[i];

		// get translational and rotational slip dirs
		Vector2d tslip_dir;
		double rslip_dir;

		// TODO: properly handle impending slip
		if(contacts_translational_friction_states[i] == FrictionState::Sliding) {
			tslip_dir = contacts_translation_slip_dir[i];
		} else {
			tslip_dir = contacts_translation_rolling_dir[i];
		}
		if(contacts_rotational_friction_states[i] == FrictionState::Sliding) {
			rslip_dir = contacts_rotational_slip_dir[i];
		} else {
			rslip_dir = contacts_rotational_rolling_dir[i];
		}

		// contact inertia diagonal block and rhs
		switch(contact_types[i]) {
			case ContactType::POINT:
				if(contacts_translational_friction_states[i] == FrictionState::Sliding 
					|| contacts_translational_friction_states[i] == FrictionState::Impending
				) {
					Trhs(trid) = -rhs_disp(full_rid+2);
					TA(trid,trid) = A_disp(full_rid+2, full_rid+2)
								- mu*A_disp.block<1,2>(full_rid+2,full_rid+0)*(tslip_dir);
					TA_size += 1;
				} else { // rolling
					Trhs.segment<3>(trid) = -rhs_disp.segment<3>(full_rid+0);
					TA.block<3,3>(trid,trid) = A_disp.block<3,3>(full_rid+0, full_rid+0);
					TA_size += 3;
				}
				break;


			case ContactType::LINE:
				if(ret_cop_sols[i].cop_type == COPContactType::LineCenter) {
					if(contacts_translational_friction_states[i] == FrictionState::Sliding ||
						contacts_translational_friction_states[i] == FrictionState::Impending
					) {
						if(line_rotational_state == FrictionState::Sliding ||
							line_rotational_state == FrictionState::Impending
						) {
							Trhs.segment<2>(trid) = -rhs_disp.segment<2>(full_rid+2);
							TA.block<2,1>(trid,trid) = A_disp.block<2,1>(full_rid+2,full_rid+2) 
										- mu*A_disp.block<2,2>(full_rid+2, full_rid+0)*(tslip_dir)
										- mu_rot*rslip_dir*A_disp.block<2,1>(full_rid+2,full_rid+4);
							TA.block<2,1>(trid,trid+1) = A_disp.block<2,1>(full_rid+2,full_rid+3);
							TA_size += 2;
						} else { // rotation rolling
							Trhs.segment<3>(trid) = -rhs_disp.segment<3>(full_rid+2);
							TA.block<3,1>(trid,trid) = A_disp.block<3,1>(full_rid+2,full_rid+2)
										- mu*A_disp.block<3,2>(full_rid+2, full_rid+0)*(tslip_dir);
							TA.block<3,2>(trid,trid+1) = A_disp.block<3,2>(full_rid+2,full_rid+3);
							TA_size += 3;
						}
					} else { // translation rolling
						if(line_rotational_state == FrictionState::Sliding ||
							line_rotational_state == FrictionState::Impending
						) {
							Trhs.segment<4>(trid) = -rhs_disp.segment<4>(full_rid);
							TA.block<4,2>(trid,trid) = A_disp.block<4,2>(full_rid+0,full_rid+0);
							TA.block<4,1>(trid,trid+2) = A_disp.block<4,1>(full_rid+0,full_rid+2)
												- mu_rot*rslip_dir*A_disp.block<4,1>(full_rid+0,full_rid+4);
							TA.block<4,1>(trid,trid+3) = A_disp.block<4,1>(full_rid+0,full_rid+3);
							TA_size += 4;
						} else { // rotation rolling
							Trhs.segment<5>(trid) = -rhs_disp.segment<5>(full_rid);
							TA.block<5,5>(trid,trid) = A_disp.block<5,5>(full_rid+0,full_rid+0);
							TA_size += 5;
						}
					}
				} else { // Line end 
					if(contacts_translational_friction_states[i] == FrictionState::Sliding ||
						contacts_translational_friction_states[i] == FrictionState::Impending
					) {
						Trhs(trid) = -rhs_disp(full_rid+2);
						TA(trid,trid) = A_disp(full_rid+2, full_rid+2)
									- mu*A_disp.block<1,2>(full_rid+2,full_rid+0)*(tslip_dir);
						TA_size += 1;
					} else { // rolling
						Trhs.segment<3>(trid) = -rhs_disp.segment<3>(full_rid+0);
						TA.block<3,3>(trid,trid) = A_disp.block<3,3>(full_rid+0, full_rid+0);
						TA_size += 3;
					}
				}
				break;


			case ContactType::SURFACE:				
				break;
			default:
				break;
		}
		tnum_rows = TA_size - trid;
		Tnum_rows[i] = tnum_rows;

		// contact inertia coupling blocks
		for(uint j = 0; j < i; i++) {
			// skip inactive contacts
			if(!contacts_active[j]) {
				continue;
			}
			uint ofull_rid = patch_indices[j];
			uint otrid = Tstart_rowids[j];
			uint otnum_rows = Tnum_rows[j];

			// get translational and rotational slip dirs
			Vector2d otslip_dir;
			double orslip_dir;

			// TODO: properly handle impending slip
			if(contacts_translational_friction_states[j] == FrictionState::Sliding) {
				otslip_dir = contacts_translation_slip_dir[j];
			} else {
				otslip_dir = contacts_translation_rolling_dir[j];
			}
			if(contacts_rotational_friction_states[j] == FrictionState::Sliding) {
				orslip_dir = contacts_rotational_slip_dir[j];
			} else {
				orslip_dir = contacts_rotational_rolling_dir[j];
			}

			switch(contact_types[i]) {
				case ContactType::POINT:
					switch(contact_types[j]) {
						case ContactType::POINT:
							getPointPointMats(
								TA.block(trid, otrid, tnum_rows, otnum_rows),
								TA.block(otrid, trid, otnum_rows, tnum_rows),
								mu,
								contacts_translational_friction_states[i],
								tslip_dir,
								contacts_translational_friction_states[j],
								otslip_dir,
								A_disp.block<3,3>(full_rid, ofull_rid),
								A_disp.block<3,3>(ofull_rid, full_rid)
							);
							break;
						case ContactType::LINE:
							if(ret_cop_sols[j].cop_type == COPContactType::LineCenter) {
								getPointLineCenterMats(
									TA.block(trid, otrid, tnum_rows, otnum_rows),
									TA.block(otrid, trid, otnum_rows, tnum_rows),
									double mu,
									double mu_rot,
									contacts_translational_friction_states[i],
									tslip_dir,
									contacts_translational_friction_states[j],
									otslip_dir,
									contacts_rotational_friction_states[j],
									orslip_dir,
									A_disp.block<3,5>(full_rid, ofull_rid), //A_block_pt_row_line_col
									A_disp.block<5,3>(ofull_rid, full_rid)
								);
							} else { // j is LineEnd
								getPointPointMats(
									TA.block(trid, otrid, tnum_rows, otnum_rows),
									TA.block(otrid, trid, otnum_rows, tnum_rows),
									mu,
									contacts_translational_friction_states[i],
									tslip_dir,
									contacts_translational_friction_states[j],
									otslip_dir,
									A_disp.block<3,3>(full_rid, ofull_rid),
									A_disp.block<3,3>(ofull_rid, full_rid)
								);
							}
							break;
						case ContactType::SURFACE:
							break;
						default:
							break;
					}
					break;
				case ContactType::LINE:
					switch(contact_types[j]) {
						case ContactType::POINT:
							if(ret_cop_sols[i].cop_type == COPContactType::LineCenter) {
								getPointLineCenterMats(
									TA.block(otrid, trid, otnum_rows, tnum_rows),
									TA.block(trid, otrid, tnum_rows, otnum_rows),
									double mu,
									double mu_rot,
									contacts_translational_friction_states[j],
									otslip_dir,
									contacts_translational_friction_states[i],
									tslip_dir,
									contacts_rotational_friction_states[i],
									rslip_dir,
									A_disp.block<3,5>(ofull_rid, full_rid), //A_block_pt_row_line_col
									A_disp.block<5,3>(full_rid, ofull_rid)
								);
							} else { // i LineEnd
								getPointPointMats(
									TA.block(trid, otrid, tnum_rows, otnum_rows),
									TA.block(otrid, trid, otnum_rows, tnum_rows),
									mu,
									contacts_translational_friction_states[i],
									tslip_dir,
									contacts_translational_friction_states[j],
									otslip_dir,
									A_disp.block<3,3>(full_rid, ofull_rid),
									A_disp.block<3,3>(ofull_rid, full_rid)
								);
							}
							break;
						case ContactType::LINE:
							if(ret_cop_sols[i].cop_type == COPContactType::LineCenter) {
								if(ret_cop_sols[j].cop_type == COPContactType::LineCenter) {
									getLineCenterLineCenterMats( // Line assumed to be with COP at center, else call pt pt
										TA.block(trid, otrid, tnum_rows, otnum_rows),
										TA.block(otrid, trid, otnum_rows, tnum_rows),
										mu,
										mu_rot,
										contacts_translational_friction_states[i],
										tslip_dir,
										contacts_rotational_friction_states[i],
										rslip_dir,
										contacts_translational_friction_states[j],
										otslip_dir,
										contacts_rotational_friction_states[j],
										orslip_dir,
										A_disp.block<5,5>(full_rid, ofull_rid),
										A_disp.block<5,5>(ofull_rid, full_rid)
									);
								} else { // j LineEnd
									getPointLineCenterMats(
										TA.block(otrid, trid, otnum_rows, tnum_rows),
										TA.block(trid, otrid, tnum_rows, otnum_rows),
										double mu,
										double mu_rot,
										contacts_translational_friction_states[j],
										otslip_dir,
										contacts_translational_friction_states[i],
										tslip_dir,
										contacts_rotational_friction_states[i],
										rslip_dir,
										A_disp.block<3,5>(ofull_rid, full_rid), //A_block_pt_row_line_col
										A_disp.block<5,3>(full_rid, ofull_rid)
									);
								}
							} else { // i LineEnd
								if(ret_cop_sols[j].cop_type == COPContactType::LineCenter) {
									getPointLineCenterMats(
										TA.block(trid, otrid, tnum_rows, otnum_rows),
										TA.block(otrid, trid, otnum_rows, tnum_rows),
										double mu,
										double mu_rot,
										contacts_translational_friction_states[i],
										tslip_dir,
										contacts_translational_friction_states[j],
										otslip_dir,
										contacts_rotational_friction_states[j],
										orslip_dir,
										A_disp.block<3,5>(full_rid, ofull_rid), //A_block_pt_row_line_col
										A_disp.block<5,3>(ofull_rid, full_rid)
									);
								} else { // j LineEnd
									getPointPointMats(
										TA.block(trid, otrid, tnum_rows, otnum_rows),
										TA.block(otrid, trid, otnum_rows, tnum_rows),
										mu,
										contacts_translational_friction_states[i],
										tslip_dir,
										contacts_translational_friction_states[j],
										otslip_dir,
										A_disp.block<3,3>(full_rid, ofull_rid),
										A_disp.block<3,3>(ofull_rid, full_rid)
									);
								}
							}
							break;
						case ContactType::SURFACE:
							break;
						default:
							break;
					}
					break;
				case ContactType::SURFACE:
					break;
				default:
					break;
			}
		}
	}
}

void COPSolver::solveInternalForMultipleContacts(
	double mu,
	double mu_rot,
	const Eigen::MatrixXd& A_disp,
	const Eigen::VectorXd& rhs_disp,
	Eigen::VectorXd& full_f_sol,
	Eigen::VectorXd& full_a_sol,
	std::vector<ContactCOPSolution>& ret_cop_sols, //NOTE: the contact type must be set here
	const std::vector<uint>& patch_indices,
	const std::vector<ContactType>& contact_types
) {
	// also populates the per primitive solution vectors, both acceleration and forces
	if(TA_size > 0) {
		if(TA_size > 1) {
			Tf_sol = TA.block(0,0,TA_size,TA_size).partialPivLu().solve(Trhs.segment(0,TA_size));
		} else {
			Tf_sol = VectorXd::Zero(1);
			Tf_sol(0) = Trhs(0)/TA(0,0);
		}
	}

	// compute full_f_sol
	for(uint i = 0; i < patch_indices.size(); i++) {
		// save the friction dirs
		Vector2d tslip_dir;
		double rslip_dir;

		// TODO: properly handle impending slip
		if(contacts_translational_friction_states[i] == FrictionState::Sliding) {
			tslip_dir = contacts_translation_slip_dir[i];
		} else {
			tslip_dir = contacts_translation_rolling_dir[i];
		}
		if(contacts_rotational_friction_states[i] == FrictionState::Sliding) {
			rslip_dir = contacts_rotational_slip_dir[i];
		} else {
			rslip_dir = contacts_rotational_rolling_dir[i];
		}

		uint full_rid = patch_indices[i];
		uint trid = Tstart_rowids[i];

		switch(contact_types[i]) {
			case ContactType::POINT:
				if(!contacts_active[i]) {
					full_f_sol.segment<3>(full_rid).setZero();
				} else if(contacts_translational_friction_states[i] == FrictionState::Sliding 
					|| contacts_translational_friction_states[i] == FrictionState::Impending
				) {
					full_f_sol.segment<2>(full_rid+0) = -mu*Tf_sol(trid+0)*tslip_dir;
					full_f_sol(full_rid+2) = Tf_sol(trid+0);
				} else { // rolling solution
					full_f_sol.segment<3>(full_rid+0) = Tf_sol.segment<3>(trid+0);
				}
				break;


			case ContactType::LINE:
				if(!contacts_active[i]) {
					full_f_sol.segment<5>(full_rid).setZero();
				} else if(ret_cop_sols[i].cop_type == COPContactType::LineCenter) {
					if(contacts_translational_friction_states[i] == FrictionState::Sliding ||
						contacts_translational_friction_states[i] == FrictionState::Impending
					) {
						full_f_sol.segment<2>(full_rid+0) = -mu*Tf_sol(trid+0)*tslip_dir;
						full_f_sol.segment<2>(full_rid+2) = Tf_sol.segment<2>(trid+0);
						if(line_rotational_state == FrictionState::Sliding ||
							line_rotational_state == FrictionState::Impending
						) {
							full_f_sol(full_rid+4) = -mu_rot*Tf_sol(trid+0)*rslip_dir;
						} else { // rotation rolling
							full_f_sol(full_rid+4) = Tf_sol(trid+2);
						}
					} else { // translation rolling
						full_f_sol.segment<4>(full_rid+0) = Tf_sol.segment<4>(trid+0);
						if(contacts_rotational_friction_states[i] == FrictionState::Sliding ||
							contacts_rotational_friction_states[i] == FrictionState::Impending
						) {
							full_f_sol(full_rid+4) = -mu_rot*Tf_sol(trid+2)*rslip_dir;
						} else { // rotation rolling
							full_f_sol(full_rid+4) = Tf_sol(trid+4);
						}
					}
				} else { //if(ret_cop_sols[i].cop_type == COPContactType::LineEnd) { // check is redundant
					full_f_sol.segment<2>(full_rid+3) = Vector2d::Zero();
					if(line_translation_state == FrictionState::Sliding ||
						line_translation_state == FrictionState::Impending
					) {
						full_f_sol.segment<2>(full_rid+0) = -mu*Tf_sol(trid+0)*tslip_dir;
						full_f_sol(full_rid+2) = Tf_sol(trid+0);
					} else { // translation rolling
						full_f_sol.segment<3>(full_rid+0) = Tf_sol.segment<3>(trid+0);
					}
				}
				break;


			case ContactType::SURFACE:
				break;
			default:
				break;
		}
	}

	// compute full_a_sol
	full_a_sol = A_disp * full_f_sol + rhs_disp;

	// compute acc_sol for each contact
	for(uint i = 0; i < patch_indices.size(); i++) {
		uint full_rid = patch_indices[i];
		switch(contact_types[i]) {
			case ContactType::POINT:
				ret_cop_sols[i].force_sol = full_f_sol.segment<3>(full_rid);
				ret_cop_sols[i].acc_sol = full_a_sol.segment<3>(full_rid);
				break;
			case ContactType::LINE:
				ret_cop_sols[i].force_sol = full_f_sol.segment<5>(full_rid);
				ret_cop_sols[i].acc_sol = full_a_sol.segment<5>(full_rid);
				break;
			case ContactType::SURFACE:
				ret_cop_sols[i].force_sol = full_f_sol.segment<6>(full_rid);
				ret_cop_sols[i].acc_sol = full_a_sol.segment<6>(full_rid);
				break;
			default:
				break;
		}
	}
}

}