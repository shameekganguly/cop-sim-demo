// SteadyContactResolution.cpp

#include <iostream>
#include "ContactSpaceModel.h"

using namespace Eigen;

namespace Sai2COPSim {

void ContactIslandModel::resolveSteadyContacts(double friction_coeff, double restitution_coeff) {
	if(numContactPoints() > 2) {
		throw(std::runtime_error("Unimplemented steady contact resolution case. Num contact points."));
	}

	if(_active_contacts.size() > 1) {
		throw(std::runtime_error("Unimplemented steady contact resolution case. Num primitives."));
	}


	if(_active_contacts.size() == 0) {
		// std::cout<< "No active contacts" << std::endl;
		return;
	}

	auto& prim = _pair_state[_active_contacts.front()];

	// get RHS vectors. no need to update the complete vectors as that would have happened either 
	// in model update or after collision resolution
	VectorXd full_cop_rhs;
	VectorXd constraint_cop_rhs;
	std::vector<uint> full_Jrow_ind_to_contact_pair_map;
	std::vector<uint> constraint_Jrow_ind_to_contact_pair_map;
	getActiveFullCOPRHSVector(full_cop_rhs, full_Jrow_ind_to_contact_pair_map);
	getActiveConstraintCOPRHSVector(constraint_cop_rhs, constraint_Jrow_ind_to_contact_pair_map);

	// compute body omega's . TODO: extend to multiple contact pairs
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
		VectorXd pt0_full_vel = _cop_full6_Jacobian_active.block(0, arb_ind, 6, arb_model->dof()) * arb_model->_dq;
		omegaB = pt0_full_vel.segment<3>(3);
	} else if (primB->_is_static) {
		std::string arb_name = primA->_articulated_body_name;
		uint arb_ind = _arb_index_map[arb_name];
		arbB = _arb_manager->getBody(arb_name);
		auto arb_model = _arb_manager->getBody(arb_name)->_model;
		VectorXd pt0_full_vel = _cop_full6_Jacobian_active.block(0, arb_ind, 6, arb_model->dof()) * arb_model->_dq;
		omegaB = pt0_full_vel.segment<3>(3);
	} else {
		std::string arbA_name = primA->_articulated_body_name;
		std::string arbB_name = primB->_articulated_body_name;
		arbA = _arb_manager->getBody(arbA_name);
		arbB = _arb_manager->getBody(arbB_name);
		uint arbA_ind = _arb_index_map[arbA_name];
		uint arbB_ind = _arb_index_map[arbB_name];
		auto arbA_model = _arb_manager->getBody(arbA_name)->_model;
		auto arbB_model = _arb_manager->getBody(arbB_name)->_model;
		VectorXd pt0_full_vel_bdA = -_cop_full6_Jacobian_active.block(0, arbA_ind, 6, arbA_model->dof()) * arbA_model->_dq;
		VectorXd pt0_full_vel_bdB = _cop_full6_Jacobian_active.block(0, arbB_ind, 6, arbB_model->dof()) * arbB_model->_dq;
		omegaA = pt0_full_vel_bdA.segment<3>(3);
		omegaB = pt0_full_vel_bdB.segment<3>(3);
	}

	// TODO: send restitution coeff in case a shock condition occurs
	const auto& contact_points = prim._geom_prim_pair->info.contact_points;
	if(prim._active_points.size() == 1) {
		// if point contact:
		// - get displaced matrices from the point0 constraint matrix to the active point
		Matrix3d active_pt_lambda_inv;
		Vector3d active_pt_rhs, active_pt_lin_vel, angular_vel;

		Vector3d r_point0_to_active_pt = contact_points[prim._active_points.front()] - contact_points[0]; // in world frame
		r_point0_to_active_pt = prim._rot_contact_frame_to_world.transpose() * r_point0_to_active_pt;

		VectorXd pt_contact_active_lin_vel;
		getActivePtContactCollisionRHSVector(pt_contact_active_lin_vel); // 1 active pair, 1 active pt contact. so this is just size 3
		active_pt_lin_vel = pt_contact_active_lin_vel;

		active_pt_rhs = full_cop_rhs.segment<3>(0) + full_cop_rhs.segment<3>(3).cross(r_point0_to_active_pt);
		active_pt_rhs += omegaB.cross(omegaB.cross(r_point0_to_active_pt)) - omegaA.cross(omegaA.cross(r_point0_to_active_pt));
		active_pt_lambda_inv = _pt_contact_Lambda_inv_active;

		// - call the single point LCP contact solver
		auto lcp_sol = solveCollLCPOnePoint (
			active_pt_lambda_inv,
			active_pt_rhs,
			active_pt_lin_vel,
			0.0,
			friction_coeff
		);
		if (lcp_sol.result == LCPSolResult::Success) {
			// std::cout << "q " << (arbB->_model->_T_world_robot.linear()*arbB->_model->_q.segment<3>(0)).transpose() << std::endl;
			// std::cout << "Contact LCP Contact force: " << lcp_sol.p_sol.transpose() << std::endl;
		} else {
			std::cerr << "Contact LCP failed with type: " << static_cast<int>(lcp_sol.result) << std::endl;
		// 	cout << "Num contacts " << contact_model._activeContacts.size() << endl;
		// 	cout << "dq: " << model->_dq.transpose() << endl;
		// 	cout << "Last post contact vel: " << post_collision_contact_vel.transpose() << endl;
		// 	cout << "Contact lambda inv: " << contact_lambda_inv << endl;
		// 	cout << "Contact rhs: " << rhs_contact.transpose() << endl;
		// 	break;
		}
		Vector3d contact_force = lcp_sol.p_sol;
		if(arbB != NULL) {
			uint arb_ind = _arb_index_map[arbB->_name];
			arbB->jtau_contact = _pt_contact_Jacobian_active.block(0, arb_ind, 3, arbB->_model->dof()).transpose()*contact_force;
			// std::cout << arbB->jtau_contact.transpose() << std::endl;
			// std::cout << "J in contact: " << std::endl;
			// std::cout << _pt_contact_Jacobian_active << std::endl;
			// std::cout << "Complete j pt:" << std::endl;
			// std::cout << _pt_contact_Jacobian << std::endl;
		}
		if(arbA != NULL) {
			uint arb_ind = _arb_index_map[arbA->_name];
			arbA->jtau_contact = _pt_contact_Jacobian_active.block(0, arb_ind, 3, arbA->_model->dof()).transpose()*contact_force;
		}
	} else {
		// else if line contact:
		// - assemble matrices
		std::vector<std::vector<Vector3d>> boundary_points;
		boundary_points.push_back(std::vector<Vector3d>());
		boundary_points[0].push_back(prim._rot_contact_frame_to_world.transpose()*(contact_points[0] - contact_points[0]));
		boundary_points[0].push_back(prim._rot_contact_frame_to_world.transpose()*(contact_points[1] - contact_points[0]));

		std::vector<ContactType> contact_types;
		contact_types.push_back(ContactType::LINE);

		std::vector<Eigen::Vector3d> omegaAs, omegaBs;
		omegaAs.push_back(omegaA);
		omegaBs.push_back(omegaB);

		std::vector<Eigen::Vector3d> linear_contact_velocity;
		VectorXd constraint_lin_vel = VectorXd::Zero(6);
		if(arbB != NULL) {
			uint arb_ind = _arb_index_map[arbB->_name];
			constraint_lin_vel += _cop_full6_Jacobian_active.block(0, arb_ind, 6, arbB->_model->dof())*arbB->_model->_dq;
		}
		if(arbA != NULL) {
			uint arb_ind = _arb_index_map[arbA->_name];
			constraint_lin_vel += _cop_full6_Jacobian_active.block(0, arb_ind, 6, arbA->_model->dof())*arbA->_model->_dq;
		}
		// std::cout << "dq " << (arbB->_model->_T_world_robot.linear()*arbB->_model->_dq).transpose() << std::endl;
		// std::cout << "q " << (arbB->_model->_T_world_robot.linear()*arbB->_model->_q.segment<3>(0)).transpose() << std::endl;
		linear_contact_velocity.push_back(constraint_lin_vel.segment<3>(0));

		COPSolver solver;
		auto cop_sol = solver.solveStartWithPatchCentroid(
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
		if(cop_sol.result == COPSolResult::Success) {
			// set COP force
			// std::cout << "v:" << linear_contact_velocity[0].transpose() << std::endl;
			// std::cout << cop_sol.force_sol.transpose() << std::endl;
			VectorXd cop_point0_force = VectorXd::Zero(6);
			cop_point0_force.segment<3>(0) = cop_sol.force_sol.segment<3>(0);
			cop_point0_force.segment<2>(4) = cop_sol.force_sol.segment<2>(3);
			// transform force from local COP point to point0
			cop_point0_force.segment<3>(3) += cop_sol.local_cop_pos.cross(cop_point0_force.segment<3>(0));
			if(arbB != NULL) {
				uint arb_ind = _arb_index_map[arbB->_name];
				arbB->jtau_contact = _cop_full6_Jacobian_active.block(0, arb_ind, 6, arbB->_model->dof()).transpose()*cop_point0_force;
			}
			if(arbA != NULL) {
				uint arb_ind = _arb_index_map[arbA->_name];
				arbA->jtau_contact = _cop_full6_Jacobian_active.block(0, arb_ind, 6, arbA->_model->dof()).transpose()*cop_point0_force;
			}
		} else {
			std::cerr << "COP solution failed" << std::endl;
		}
	}
}

}