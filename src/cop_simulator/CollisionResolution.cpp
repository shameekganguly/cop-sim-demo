// CollisionResolution.cpp

#include <iostream>
#include <algorithm>
#include "ContactSpaceModel.h"
#include "lcp_solvers/LCPSolver.h"

using namespace Eigen;

namespace Sai2COPSim {

bool ContactIslandModel::resolveCollisions(double friction_coeff, double restitution_coeff) {
	bool ret_coll_flag = false;
	// if(numContactPoints() > 10) {
	// 	throw(std::runtime_error("Unimplemented collision resolution case."));
	// }
	// NOTE: currently we use a pt based solver for collisions

	bool did_update_active_contacts = false;

	// rhs collision vector for active contacts
	VectorXd pt_contact_rhs_coll_active;
	std::vector<uint> jind_map; // TODO: think about using this

	// start while loop for resolving collisions.
	// while loop is required because we only simultaneously resolve the active collisions.
	// such a resolution can result in an inactive collision becoming active.
	// the resolution continues until all collisions have been resolved.
	// number of resolutions required can be large, roughly proportional to the number of
	// prim-pairs.
	// TODO: think about a proper bound
	long tryCollSolveCounter = 0;
	// std::cout << _pair_state[0]._geom_prim_pair->primA->_name << std::endl;
	// for(auto it: _arb_index_map) {
			// auto arb_model = _arb_manager->getBody(it.first)->_model;
			// std::cout << it.first << std::endl;
			// std::cout << arb_model->_dq.transpose() << std::endl;
		// }
	while(true) {
		tryCollSolveCounter++;

		// update RHS for collisions
		updatePtContactRHSCollVector();

		bool is_colliding = false;
		bool did_update_coll_active_contacts = false;
		double max_collision_speed = 0;

		uint max_num_pts_any_prim = 0; // TODO: get rid of this

		// std::cout << "Coll vel: " << _pt_contact_rhs_coll.transpose() << std::endl;
		// std::cout << tryCollSolveCounter << std::endl;

		// - loop over primitive pairs
		for(auto& p: _pair_state) {
			uint rstart = _pt_contact_Jacobian_prim_start_ind[p._id];
			bool any_pt_active = false;
			// - - check which points are active
			// loop over all points in this primitive pair
			for(uint ptid = 0; ptid < p._num_contact_pts; ptid++) {
				// pre collision velocity at contact points is the same as _pt_contact_rhs_coll
				double pre_coll_separation_speed = _pt_contact_rhs_coll[rstart + ptid*3 + 2];
				max_collision_speed = fmax(max_collision_speed, -pre_coll_separation_speed);

				// check for collision
				if(tryCollSolveCounter > 1) {
					// If we have already solved for collisions once, we use a threshold
					// to numerically quantify if a body is still colliding.
					// This is to accomodate for numerical precision of the solver
					is_colliding = is_colliding || (pre_coll_separation_speed <
														-COPAlgorithmicConstants::SEPARATION_SPEED_NUMERICAL_ZERO_THRESHOLD);
				} else {
					// else, if it is the first time, we use absolute 0
					is_colliding = is_colliding || (pre_coll_separation_speed < 0);
				}

				// check if active.
				if(pre_coll_separation_speed >
					COPAlgorithmicConstants::MAX_SEPARATION_SPEED_FOR_ACTIVE_CONTACT
				) {
					did_update_coll_active_contacts = p.deactivateContactPoint(ptid) || did_update_coll_active_contacts;
					// std::cout << "Deactivate pt " << ptid << std::endl;
				} else {
					did_update_coll_active_contacts = p.activateContactPoint(ptid) || did_update_coll_active_contacts;
					any_pt_active = true;
					// std::cout << "Activate pt " << ptid << " # " << p._active_points.size() << std::endl;
				}
			}
			// if no points are active for this primitive pair, deactivate the pair completely
			// TODO: can a primitive be active even though none of the contact points are active?
			// 	Think particularly about non-polygonal surface contact patches
			if(any_pt_active) {
				// std::cout << "Activate contact pair" << std::endl;
				did_update_coll_active_contacts = activateContactPair(p._id) || did_update_coll_active_contacts;
			} else {
				did_update_coll_active_contacts = deactivateContactPair(p._id) || did_update_coll_active_contacts;
				// if(did_update_coll_active_contacts) std::cout << "Deactivate contact pair "  << std::endl;
			}
			max_num_pts_any_prim = std::max(max_num_pts_any_prim, static_cast<uint>(p._active_points.size()));
		}
		// track whether we ever changed the active contacts during the iterations
		did_update_active_contacts = did_update_active_contacts || did_update_coll_active_contacts;


		if(did_update_coll_active_contacts) {
			// update all of the active contact matrices, both for collisions and for steady contacts
			// NOTE: we update these matrices even if we are not colliding, because the
			// active contact state might remain the same for the next simulation iteration.
			// so we want to preserve them
			// std::cout << "Updated active contacts " << _pair_state[0]._active_points.size() << std::endl;
			getActivePtContactCollisionMatrices(
				_pt_contact_Jacobian_active,
				_pt_contact_Lambda_inv_active,
				pt_contact_rhs_coll_active,
				jind_map
			);
		} else if(is_colliding) {
			// just update the active RHS vector, if we are colliding
			getActivePtContactCollisionRHSVector(pt_contact_rhs_coll_active);
		}

		// if not colliding, we are done resolving collisions
		if(!is_colliding) break;
		ret_coll_flag = true;

		// - compute epsilon based on max collision speed. set to 0 if it is very small
		double adjusted_restitution_coeff = restitution_coeff;
		if (max_collision_speed < COPAlgorithmicConstants::MIN_COLLISION_SPEED_FOR_STEADY_CONTACT) {
			// We use the fastest collision to determine the coefficient of restitution
			// std::cout << "Force stick" << std::endl;
			adjusted_restitution_coeff = 0.0; // force inelastic collision to bring to steady contact
		} else {
			// TODO: we need another condition here to enable steady contact
			// if the prim pair is in collision for multiple iterations where
			// contact detection has been updated. Otherwise we can get limit
			// cycles where the object keeps bouncing indefinitely since it
			// separates each time with a small speed > MIN_COLLISION_SPEED_FOR_STEADY_CONTACT
			// and then collides again due to gravity with the same speed each
			// time.

			// std::cout << max_collision_speed << std::endl;
		}

		// - call LCP solver.
		const uint num_active_contact_pts = pt_contact_rhs_coll_active.size()/3;
		if (num_active_contact_pts > 0) {
			auto solver = Sai2LCPSolver::LCPSolver();
			// std::cout << _pt_contact_Lambda_inv_active << std::endl;
			// std::cout << pt_contact_rhs_coll_active.transpose() << std::endl;
			// std::cout << adjusted_restitution_coeff <<std::endl;
			_last_coll_lcp_sol = solver.solve(
				_pt_contact_Lambda_inv_active,
				pt_contact_rhs_coll_active,
				pt_contact_rhs_coll_active,
				adjusted_restitution_coeff,
				friction_coeff
			);
			if (_last_coll_lcp_sol.result == Sai2LCPSolver::LCPSolResult::Success) {
				// std::cout << "LCP Impulse: " << _last_coll_lcp_sol.p_sol.transpose() << std::endl;
				// std::cout << "Num contacts " << num_active_contact_pts << std::endl;
				// std::cout << _pt_contact_Jacobian_active << std::endl;
			} else {
				std::cerr << "LCP failed with type: " << static_cast<int>(_last_coll_lcp_sol.result) << std::endl;
				// cout << "Num contacts " << contact_model._activeContacts.size() << endl;
				// cout << "pre coll vel: " << pre_collision_contact_vel.transpose() << endl;
				// cout << "lambda inv " << contact_lambda_inv << endl;
				// cout << "restitution " << coll_restitution << endl;
				// cout << "q " << model->_q.transpose() << endl;
				// cout << "contact Jacobian " << contact_jacobian << endl;
				throw(std::out_of_range("Collision LCP failure."));
			}
		}

		// - loop over ARBs in this island and compute dq+
		VectorXd Jtranspose_P = _pt_contact_Jacobian_active.transpose()*_last_coll_lcp_sol.p_sol;
		for(auto it: _arb_index_map) {
			uint arb_Jind = it.second;
			auto arb_model = _arb_manager->getBody(it.first)->_model;
			// std::cout << it.first << std::endl;
			// std::cout << arb_model->_dq.transpose() << std::endl;
			// std::cout << _pt_contact_Jacobian_active << std::endl;
			arb_model->_dq += arb_model->_M_inv*Jtranspose_P.segment(arb_Jind, arb_model->dof());
			// std::cout << "T: " << Jtranspose_P.segment(arb_Jind, arb_model->dof()).transpose() << std::endl;
			// std::cout << "Dq: " << arb_model->_dq.transpose() << std::endl;
		}

		// check if too many iters
		if(tryCollSolveCounter == COPAlgorithmicConstants::COLLISION_RESOLUTION_SIMULTANEOUS_MAX_ITERATIONS) {

			std::cerr << "Post coll vel: " << _pt_contact_rhs_coll.transpose() << std::endl;
			std::cerr << "LCP impulse: " << _last_coll_lcp_sol.p_sol.transpose() << std::endl;
			std::cerr << "A" << "\n" << _pt_contact_Lambda_inv_active << std::endl;
			std::cerr << "rhs: " << pt_contact_rhs_coll_active.transpose() << std::endl;
			std::cerr << "vels: " << pt_contact_rhs_coll_active.transpose() << std::endl;
			std::cerr << "restitution " <<  adjusted_restitution_coeff << std::endl;
			throw(std::runtime_error("Too many collision retries"));
		}
	}

	// if a collision occurred, velocity dependent terms will change. so update the RHS vectors
	if(tryCollSolveCounter > 1) {
		// update nonlin acceleration terms for bodies
		updateBodyNonlinAccelerations();
		// update only RHS vectors.
		updateRHSVectors(); //TODO: optimize by removing the collision RHS vectors from here
	}

	// check if active contacts changed since beginning. If so, update the active contact
	// model matrices
	if(did_update_active_contacts) {
		// update all active matrices for steady contact
		// TODO: optimize by removing the contact RHS vectors which are unused
		VectorXd cop_full_RHS_act, cop_const_RHS_act;
		std::vector<uint> temp;
		getActiveFullCOPMatrices(_cop_full6_Jacobian_active, _cop_full6_Lambda_inv_active, cop_full_RHS_act, temp);
		getActiveConstraintCOPMatrices(_cop_constraint_Jacobian_active, _cop_constraint_Lambda_inv_active, cop_const_RHS_act, temp);
		//TODO: what about the last COP solution?
	}

	// // initialize COP from inelastic collision result
	// if (is_post_steady_contact) {
	// 	if(contact_model._activeContacts.size() == 2) {
	// 		if(abs(lcp_sol.p_sol[2]) > 1e-8 && abs(lcp_sol.p_sol[5]) > 1e-8) {
	// 			if(LOG_DEBUG) cout << "2 pt impulse too small. Not updating COP." << endl;
	// 		} else {
	// 			// update cop sol from collision lcp sol
	// 			last_cop_sol = getLCPForInelasticCollResult(lcp_sol.p_sol, contact_model._contactPointsLocal);
	// 			Vector3d global_cop_pos;
	// 			model->position(global_cop_pos, object_link_name, last_cop_sol.local_cop_pos);
	// 			if(LOG_DEBUG) cout << "Coll deduced COP in global: " << (object_in_world*global_cop_pos).transpose() << endl;
	// 		}
	// 	} else {
	// 		if(LOG_DEBUG) cout << "Rolling with 1 point contact. No COP deduced." << endl;
	// 	}
	// }
	return ret_coll_flag;
}

}
