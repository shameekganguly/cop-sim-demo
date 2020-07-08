// SteadyContactResolution.cpp

#include <iostream>
#include "ContactSpaceModel.h"

using namespace Eigen;

namespace Sai2COPSim {

bool ContactIslandModel::resolveSteadyContacts(double friction_coeff, double restitution_coeff) {
	if(_active_contacts.size() == 0) {
		// std::cout<< "No active contacts" << std::endl;
		return false;
	}

	// assemble solver matrices
	std::vector<std::vector<Vector3d>> boundary_points;
	std::vector<ContactType> contact_types;
	std::vector<Eigen::Vector3d> omegaAs, omegaBs;
	std::vector<Eigen::Vector3d> linear_contact_velocity;
	for(uint prim_id: _active_contacts) {
		auto& prim = _pair_state[prim_id];
		// std::cout << "Prim id " << prim_id << std::endl;
		// std::cout << prim._geom_prim_pair->primA->_name << ", " << prim._geom_prim_pair->primB->_name << std::endl;

		// compute body omega's
		//TODO: move to ContactModel::updateVelocityTerms
		Vector3d omegaA, omegaB;
		omegaA.setZero();
		omegaB.setZero();
		auto primA = prim._geom_prim_pair->primA;
		auto primB = prim._geom_prim_pair->primB;
		ArticulatedRigidBody* arbA = NULL; // TODO: consider moving these to the contact_pair_state
		ArticulatedRigidBody* arbB = NULL;
		if(primA->_is_static) {
			std::string arb_name = primB->_articulated_body_name;
			uint arb_ind = _arb_index_map[arb_name];
			arbB = _arb_manager->getBody(arb_name);
			auto arb_model = _arb_manager->getBody(arb_name)->_model;
			omegaB = _cop_full6_Jacobian.block(prim_id*6+3, arb_ind, 3, arb_model->dof()) * arb_model->_dq;
		} else if (primB->_is_static) {
			std::string arb_name = primA->_articulated_body_name;
			uint arb_ind = _arb_index_map[arb_name];
			arbB = _arb_manager->getBody(arb_name);
			auto arb_model = _arb_manager->getBody(arb_name)->_model;
			omegaB = _cop_full6_Jacobian.block(prim_id*6+3, arb_ind, 3, arb_model->dof()) * arb_model->_dq;
		} else {
			std::string arbA_name = primA->_articulated_body_name;
			std::string arbB_name = primB->_articulated_body_name;
			arbA = _arb_manager->getBody(arbA_name);
			arbB = _arb_manager->getBody(arbB_name);
			uint arbA_ind = _arb_index_map[arbA_name];
			uint arbB_ind = _arb_index_map[arbB_name];
			auto arbA_model = _arb_manager->getBody(arbA_name)->_model;
			auto arbB_model = _arb_manager->getBody(arbB_name)->_model;
			omegaA = -_cop_full6_Jacobian.block(prim_id*6+3, arbA_ind, 3, arbA_model->dof()) * arbA_model->_dq;
			omegaB = _cop_full6_Jacobian.block(prim_id*6+3, arbB_ind, 3, arbB_model->dof()) * arbB_model->_dq;
		}
		omegaAs.push_back(omegaA);
		omegaBs.push_back(omegaB);

		// add boundary points in COP frame //TODO: move this to ContactModel get active matrices
		boundary_points.push_back(std::vector<Vector3d>());
		uint bpsize = boundary_points.size();
		Vector3d pt0;
		if(prim._geom_prim_pair->info->type == ContactType::SURFACE) {
			pt0 = prim._geom_prim_pair->info->contact_patch._interior_point;
		} else {
			pt0 = prim._geom_prim_pair->info->contact_points[0];
		}
		for(auto pt: prim._geom_prim_pair->info->contact_points) { //TODO: handle some points not being active
			boundary_points[bpsize - 1].push_back(prim._rot_contact_frame_to_world.transpose()*(pt - pt0));
		}

		// add contact types
		contact_types.push_back(prim._geom_prim_pair->info->type);

		// linear_contact_velocity
		Vector3d constraint_lin_vel = Vector3d::Zero();
		if(arbB != NULL) {
			uint arb_ind = _arb_index_map[arbB->_name];
			constraint_lin_vel += _cop_full6_Jacobian.block(prim_id*6, arb_ind, 3, arbB->_model->dof())*arbB->_model->_dq;
		}
		if(arbA != NULL) {
			uint arb_ind = _arb_index_map[arbA->_name];
			constraint_lin_vel += _cop_full6_Jacobian.block(prim_id*6, arb_ind, 3, arbA->_model->dof())*arbA->_model->_dq;
		}
		linear_contact_velocity.push_back(constraint_lin_vel);
	}

	// get RHS vectors. no need to update the complete vectors as that would have happened either 
	// in model update or after collision resolution
	VectorXd full_cop_rhs;
	VectorXd constraint_cop_rhs;
	std::vector<uint> full_Jrow_ind_to_contact_pair_map;
	std::vector<uint> constraint_Jrow_ind_to_contact_pair_map;
	getActiveFullCOPRHSVector(full_cop_rhs, full_Jrow_ind_to_contact_pair_map);
	getActiveConstraintCOPRHSVector(constraint_cop_rhs, constraint_Jrow_ind_to_contact_pair_map);


	// TODO: send restitution coeff in case a shock condition occurs
	if(_active_contacts.size() == 1) {
		auto& prim = _pair_state[_active_contacts.front()];
		// const auto& contact_points = prim._geom_prim_pair->info->contact_points;
		if(prim._active_points.size() == 1) {
			// --------------------- LCP solver solution ---------------------
			// // if point contact:
			// // - get displaced matrices from the point0 constraint matrix to the active point
			// Matrix3d active_pt_lambda_inv;
			// Vector3d active_pt_rhs, active_pt_lin_vel, angular_vel;

			// Vector3d r_point0_to_active_pt = contact_points[prim._active_points.front()] - contact_points[0]; // in world frame
			// r_point0_to_active_pt = prim._rot_contact_frame_to_world.transpose() * r_point0_to_active_pt;

			// VectorXd pt_contact_active_lin_vel;
			// getActivePtContactCollisionRHSVector(pt_contact_active_lin_vel); // 1 active pair, 1 active pt contact. so this is just size 3
			// active_pt_lin_vel = pt_contact_active_lin_vel;

			// active_pt_rhs = full_cop_rhs.segment<3>(0) + full_cop_rhs.segment<3>(3).cross(r_point0_to_active_pt);
			// active_pt_rhs += omegaB.cross(omegaB.cross(r_point0_to_active_pt)) - omegaA.cross(omegaA.cross(r_point0_to_active_pt));
			// active_pt_lambda_inv = _pt_contact_Lambda_inv_active;

			// // - call the single point LCP contact solver
			// auto lcp_sol = solveCollLCPOnePoint (
			// 	active_pt_lambda_inv,
			// 	active_pt_rhs,
			// 	active_pt_lin_vel,
			// 	0.0,
			// 	friction_coeff
			// );
			// if (lcp_sol.result == LCPSolResult::Success) {
			// 	// std::cout << "q " << (arbB->_model->_T_world_robot.linear()*arbB->_model->_q.segment<3>(0)).transpose() << std::endl;
			// 	// std::cout << "Contact LCP Contact force: " << lcp_sol.p_sol.transpose() << std::endl;
			// } else {
			// 	std::cerr << "Contact LCP failed with type: " << static_cast<int>(lcp_sol.result) << std::endl;
			// // 	cout << "Num contacts " << contact_model._activeContacts.size() << endl;
			// // 	cout << "dq: " << model->_dq.transpose() << endl;
			// // 	cout << "Last post contact vel: " << post_collision_contact_vel.transpose() << endl;
			// // 	cout << "Contact lambda inv: " << contact_lambda_inv << endl;
			// // 	cout << "Contact rhs: " << rhs_contact.transpose() << endl;
			// // 	break;
			// }
			// // TODO: build cop sol from the point sol
			// Vector3d contact_force = lcp_sol.p_sol;
			// if(arbB != NULL) {
			// 	uint arb_ind = _arb_index_map[arbB->_name];
			// 	arbB->jtau_contact = _pt_contact_Jacobian_active.block(0, arb_ind, 3, arbB->_model->dof()).transpose()*contact_force;
			// 	// std::cout << arbB->jtau_contact.transpose() << std::endl;
			// 	// std::cout << "J in contact: " << std::endl;
			// 	// std::cout << _pt_contact_Jacobian_active << std::endl;
			// 	// std::cout << "Complete j pt:" << std::endl;
			// 	// std::cout << _pt_contact_Jacobian << std::endl;
			// }
			// if(arbA != NULL) {
			// 	uint arb_ind = _arb_index_map[arbA->_name];
			// 	arbA->jtau_contact = _pt_contact_Jacobian_active.block(0, arb_ind, 3, arbA->_model->dof()).transpose()*contact_force;
			// }

			// --------------------- COP solver solution ---------------------
			// Check for case where this is actually a line contact with one pt disabled
			Matrix3d A_constraint;
			Vector3d rhs_constraint;
			Vector3d pt0_to_pt_pos;
			Vector3d lin_vel_pt;
			if(prim._geom_prim_pair->info->type == ContactType::LINE) {//
				// TODO: add surface
				// TODO: this should be handled in getActiveMatrices
				if(prim._active_points.front() == 0) {
					A_constraint = _cop_constraint_Lambda_inv_active.block<3,3>(0,0);
					rhs_constraint = constraint_cop_rhs.segment<3>(0);
					lin_vel_pt = linear_contact_velocity[0];
					pt0_to_pt_pos.setZero();
				} else {
					MatrixXd tempA;
					VectorXd temprhs;
					COPSolver::getCOPLineContactDisplacedMatricesExtended(
						tempA,
						temprhs,
						lin_vel_pt,
						boundary_points[0][1](0),
						_cop_constraint_Lambda_inv_active,
						constraint_cop_rhs,
						omegaAs[0],
						omegaBs[0],
						linear_contact_velocity[0]
					);
					A_constraint = tempA.block<3,3>(0,0);
					rhs_constraint = temprhs.segment<3>(0);
					pt0_to_pt_pos = boundary_points[0][1];
				}
			} else if(prim._geom_prim_pair->info->type == ContactType::POINT) {
				A_constraint = _cop_constraint_Lambda_inv_active;
				rhs_constraint = constraint_cop_rhs;
				lin_vel_pt = linear_contact_velocity[0];
			} else {
				// check last cop solution. if it is success and lies on patch end,
				// then solve for a pt contact at the last cop position
				if(prim._last_cop_sol.result == COPSolResult::Success
					&& prim._last_cop_sol.cop_type != COPContactType::PatchCenter
				) {
					MatrixXd tempA;
					VectorXd temprhs;
					COPSolver::getCOPSurfaceContactDisplacedMatrices(
						tempA,
						temprhs,
						lin_vel_pt,
						prim._last_cop_sol.local_cop_pos, //TODO: should we use the active point?
						_cop_constraint_Lambda_inv_active,
						constraint_cop_rhs,
						omegaAs[0],
						omegaBs[0],
						linear_contact_velocity[0]
					);
					A_constraint = tempA.block<3,3>(0,0);
					rhs_constraint = temprhs.segment<3>(0);
					pt0_to_pt_pos = prim._last_cop_sol.local_cop_pos;
				} else {
					// std::cout << "Last cop pos: " << prim._last_cop_sol.local_cop_pos.transpose() << std::endl;
					// std::cout << "Last cop force: " << prim._last_cop_sol.force_sol.transpose() << std::endl;
					// throw(std::runtime_error("Unimplemented Surface contact case"));
					MatrixXd tempA;
					VectorXd temprhs;
					COPSolver::getCOPSurfaceContactDisplacedMatrices(
						tempA,
						temprhs,
						lin_vel_pt,
						boundary_points[0][prim._active_points.front()], //TODO: should we use the active point?
						_cop_constraint_Lambda_inv_active,
						constraint_cop_rhs,
						omegaAs[0],
						omegaBs[0],
						linear_contact_velocity[0]
					);
					A_constraint = tempA.block<3,3>(0,0);
					rhs_constraint = temprhs.segment<3>(0);
					pt0_to_pt_pos = boundary_points[0][prim._active_points.front()];
				}
			}
			COPSolver solver;
			prim._last_cop_sol = solver.solvePtOnly(
				friction_coeff,
				A_constraint,
				rhs_constraint,
				lin_vel_pt
			);
			if(prim._last_cop_sol.result != COPSolResult::Success) {
				// std::cerr << "COP solution failed" << std::endl;
				throw(std::runtime_error("COP solution failed"));
			} else {
				// translate pt solution to line solution if necessary
				if(prim._geom_prim_pair->info->type == ContactType::LINE) {
					Vector3d force_pt = prim._last_cop_sol.force_sol;
					prim._last_cop_sol.force_sol.setZero(5);
					prim._last_cop_sol.force_sol.segment<3>(0) = force_pt;
					prim._last_cop_sol.local_cop_pos = pt0_to_pt_pos;
					prim._last_cop_sol.cop_type = COPContactType::LineEnd;
				} else if(prim._geom_prim_pair->info->type == ContactType::SURFACE) {
					Vector3d force_pt = prim._last_cop_sol.force_sol;
					prim._last_cop_sol.force_sol.setZero(6);
					prim._last_cop_sol.force_sol.segment<3>(0) = force_pt;
					prim._last_cop_sol.local_cop_pos = pt0_to_pt_pos;
					prim._last_cop_sol.cop_type = COPContactType::PatchCurvePoint;
				}
			}
		} else {
			// else if line contact:
			if(contact_types[0] == ContactType::LINE) {
				COPSolver solver;
				prim._last_cop_sol = solver.solveStartWithPatchCentroid(
					friction_coeff,
					_cop_constraint_Lambda_inv_active,
					constraint_cop_rhs,
					boundary_points,
					constraint_Jrow_ind_to_contact_pair_map,
					contact_types,
					omegaAs,
					omegaBs,
					linear_contact_velocity
				);
				if(prim._last_cop_sol.result != COPSolResult::Success) {
					// std::cerr << "COP solution failed" << std::endl;
					throw(std::runtime_error("COP solution failed"));
				}
			} else if(contact_types[0] == ContactType::SURFACE) {
				COPSolver solver;
				prim._last_cop_sol = solver.solveSurfaceContact(
					friction_coeff,
					_cop_constraint_Lambda_inv_active,
					constraint_cop_rhs,
					prim._geom_prim_pair->info->contact_patch,
					constraint_Jrow_ind_to_contact_pair_map,
					contact_types,
					omegaAs,
					omegaBs,
					linear_contact_velocity
				);
				// std::cout << "Solver call Lin vel: " << linear_contact_velocity[0].transpose() << std::endl;
				// std::cout << "Solver cop pos: " << prim._last_cop_sol.local_cop_pos.transpose() << std::endl;
				if(prim._last_cop_sol.result != COPSolResult::Success) {
					std::cout << (uint)(prim._last_cop_sol.result) << std::endl;
					// std::cerr << "COP solution failed" << std::endl;
					throw(std::runtime_error("COP solution failed"));
				}
			}
		}
	} else if(_active_contacts.size() == 2) {
		if( (contact_types[0] == ContactType::LINE && contact_types[1] == ContactType::POINT) ||
			(contact_types[0] == ContactType::POINT && contact_types[1] == ContactType::LINE)
		) {
			COPSolver solver;
			ContactCOPSolution two_sol = solver.solveStartWithPatchCentroidOnePointOneLine(
				friction_coeff,
				_cop_constraint_Lambda_inv_active,
				constraint_cop_rhs,
				boundary_points,
				constraint_Jrow_ind_to_contact_pair_map,
				contact_types,
				omegaAs,
				omegaBs,
				linear_contact_velocity
			);
			if(two_sol.result != COPSolResult::Success) {
				// std::cerr << "COP solution failed" << std::endl;
				throw(std::runtime_error("COP solution failed"));
			}
			uint ind = 0;
			uint rind = 0;
			for(uint prim_id: _active_contacts) {
				auto& prim = _pair_state[prim_id];
				if(contact_types[ind] == ContactType::LINE && prim._active_points.size() != 2) {
					std::cout << contact_types[0] << " " << contact_types[1] << std::endl;
					std::cout << "Num active pts: " << prim._active_points.size() << std::endl;
					std::cout << _pt_contact_rhs_coll.transpose() << std::endl;
					// throw(std::runtime_error("Unimplemented steady contact resolution case. Line and pt contact, but line both points not active."));
				}
				if(contact_types[ind] == ContactType::LINE) {
					prim._last_cop_sol.result = two_sol.result;
					prim._last_cop_sol.force_sol = two_sol.force_sol.segment<5>(rind);
					prim._last_cop_sol.local_cop_pos = two_sol.local_cop_pos;
					prim._last_cop_sol.cop_type = two_sol.cop_type;
					rind += 5;
				} else {
					prim._last_cop_sol.result = two_sol.result;
					prim._last_cop_sol.force_sol = two_sol.force_sol.segment<3>(rind);
					prim._last_cop_sol.local_cop_pos.setZero();
					prim._last_cop_sol.cop_type = COPContactType::Unknown;
					rind += 3;
				}
				ind++;
			}
		} else {
			throw(std::runtime_error("Unimplemented steady contact resolution case. Two primitives, not line and pt."));
		}
	} else {
		throw(std::runtime_error("Unimplemented steady contact resolution case. Num primitives."));
	}

	// set torques on bodies
	for(uint prim_id: _active_contacts) {
		auto& prim = _pair_state[prim_id];

		if(prim._last_cop_sol.result != COPSolResult::Success) {
			// throws above for now
			continue;
		}

		auto primA = prim._geom_prim_pair->primA;
		auto primB = prim._geom_prim_pair->primB;
		auto ctype = prim._geom_prim_pair->info->type;
		ArticulatedRigidBody* arbA = NULL; // TODO: consider moving these to the contact_pair_state
		ArticulatedRigidBody* arbB = NULL;
		if(primA->_is_static) {
			std::string arb_name = primB->_articulated_body_name;
			arbB = _arb_manager->getBody(arb_name);
		} else if (primB->_is_static) {
			std::string arb_name = primA->_articulated_body_name;
			arbB = _arb_manager->getBody(arb_name);
		} else {
			std::string arbA_name = primA->_articulated_body_name;
			std::string arbB_name = primB->_articulated_body_name;
			arbA = _arb_manager->getBody(arbA_name);
			arbB = _arb_manager->getBody(arbB_name);
		}
		VectorXd cop_point0_force = VectorXd::Zero(6);
		switch (ctype) {
			case ContactType::POINT:
				if(arbB != NULL) {
					uint arb_ind = _arb_index_map[arbB->_name];
					arbB->jtau_contact = _cop_full6_Jacobian.block(prim_id*6, arb_ind, 3, arbB->_model->dof()).transpose()*prim._last_cop_sol.force_sol;
				}
				if(arbA != NULL) {
					uint arb_ind = _arb_index_map[arbA->_name];
					arbA->jtau_contact = _cop_full6_Jacobian.block(prim_id*6, arb_ind, 3, arbA->_model->dof()).transpose()*prim._last_cop_sol.force_sol;
				}
				// update COP world pos
				prim._cop_pos = prim._reference_point;
				break;
			case ContactType::LINE:
				cop_point0_force.segment<3>(0) = prim._last_cop_sol.force_sol.segment<3>(0);
				cop_point0_force.segment<2>(4) = prim._last_cop_sol.force_sol.segment<2>(3);
				// transform force from local COP point to point0
				cop_point0_force.segment<3>(3) += prim._last_cop_sol.local_cop_pos.cross(cop_point0_force.segment<3>(0));
				if(arbB != NULL) {
					uint arb_ind = _arb_index_map[arbB->_name];
					arbB->jtau_contact = _cop_full6_Jacobian.block(prim_id*6, arb_ind, 6, arbB->_model->dof()).transpose()*cop_point0_force;
				}
				if(arbA != NULL) {
					uint arb_ind = _arb_index_map[arbA->_name];
					arbA->jtau_contact = _cop_full6_Jacobian.block(prim_id*6, arb_ind, 6, arbA->_model->dof()).transpose()*cop_point0_force;
				}
				prim._cop_pos = prim._reference_point + prim._rot_contact_frame_to_world*prim._last_cop_sol.local_cop_pos;
				break;
			case ContactType::SURFACE:
				cop_point0_force = prim._last_cop_sol.force_sol;
				// transform force from local COP point to contact patch interior point, which 
				// is "point 0" for surface contact
				cop_point0_force.segment<3>(3) += prim._last_cop_sol.local_cop_pos.cross(cop_point0_force.segment<3>(0));
				if(arbB != NULL) {
					uint arb_ind = _arb_index_map[arbB->_name];
					arbB->jtau_contact = _cop_full6_Jacobian.block(prim_id*6, arb_ind, 6, arbB->_model->dof()).transpose()*cop_point0_force;
				}
				if(arbA != NULL) {
					uint arb_ind = _arb_index_map[arbA->_name];
					arbA->jtau_contact = _cop_full6_Jacobian.block(prim_id*6, arb_ind, 6, arbA->_model->dof()).transpose()*cop_point0_force;
				}
				// update COP world pos
				prim._cop_pos = prim._geom_prim_pair->info->contact_patch._interior_point
									+ prim._rot_contact_frame_to_world*prim._last_cop_sol.local_cop_pos;
				break;
			default:
				throw(std::runtime_error("Unimplemented steady contact resolution case. Prim type."));
				break;
		}
	}

	return true;
}

}