// ContactSpaceModel.cpp

#include <iostream>
#include "ContactSpaceModel.h"

using namespace Eigen;

namespace Sai2COPSim {

// GENERAL CODING NOTE FOR THIS FILE: THE CLASSES ARE BUILT MANY TIME OVER.
// SO NO EXCEPTIONS ARE THROWN. THERE ARE A FEW ASSERTS THOUGH

/* ------ Contact Island Model ------- */
ContactIslandModel::ContactIslandModel(const ContactIsland* geom_island, ArticulatedRigidBodyManager* arb_manager)
: _geom_island(geom_island), _arb_manager(arb_manager)
{
	assert(geom_island != NULL);
	// build indexed contact pair list
	uint prim_pair_id = 0;
	for(auto it = geom_island->_contact_prim_pairs.begin();
		it != geom_island->_contact_prim_pairs.end();
		it++, prim_pair_id++
	) {
		// add a new contact pair state, initially invalid
		_pair_state.push_back(ContactPairState(prim_pair_id, &(*it)));

		// save iterator with index
		_index_to_geom_island_contact_list.push_back(it);

		// mark as an active contact initially
		_active_contacts.push_back(prim_pair_id);
	}

	// create the Contact Jacobian and Lambda Inv matrices
	createContactJacobianAndLambdaInv();

	// save active matrices to full matrices, because we start with all points active
	// TODO: optimize by already determining which points are active
	_cop_full6_Jacobian_active = _cop_full6_Jacobian;
	_cop_constraint_Jacobian_active = _cop_constraint_Jacobian;
	_pt_contact_Jacobian_active = _pt_contact_Jacobian;
	_cop_full6_Lambda_inv_active = _cop_full6_Lambda_inv;
	_cop_constraint_Lambda_inv_active = _cop_constraint_Lambda_inv;
	_pt_contact_Lambda_inv_active = _pt_contact_Lambda_inv;

	// update the RHS vectors
	updateRHSVectors();
}

ContactIslandModel::~ContactIslandModel() {

}

void ContactIslandModel::createContactJacobianAndLambdaInv() {
	// Note: we create three versions of the Jacobian and LambdaInv matrices
	// version 1: 6 dof per primitive pair including Jv and Jw
	//		- this is used in the COP solver to compute projected dynamics from the last
	//      - COP solution to the current contact geometry.
	//		- This is also used for rotating the contact Jacobian from one frame to another.
	//		- TODO: think really hard about this. Do we really need it? The whole contact
	//		- model is torn down and rebuilt when the contact geometry changes currently.
	// version 2: cop constraint dof per primitive pair
	//		- 3 dof for point, 5 dof for line, 6 dof for surface
	// version 3: pt contact, 3dof per point. used for collisions

	// compute contact dof
	uint cop_constraint_contact_dof = 0; // based on COP constraints at each location
	uint full6_contact_dof = 0; // for COP solver input, we need to send the full Jacobians
	uint pt_contact_dof = 0; // based on number of points at 
	uint dq_dof = 0;
	for (auto& p: _pair_state) {
		auto p_info = p._geom_prim_pair->info;
		full6_contact_dof += 6;
		pt_contact_dof += p_info.contact_points.size()*3;
		switch(p_info.type) {
			case ContactType::POINT:
				cop_constraint_contact_dof += 3;
				break;
			case ContactType::LINE:
				cop_constraint_contact_dof += 5;
				break;
			case ContactType::SURFACE:
				cop_constraint_contact_dof += 6;
				break;
			case ContactType::UNDEFINED:
				std::cerr << "Undefined contact type" << std::endl;
				throw(std::runtime_error("Bad contact"));
		}
	}

	// compute dq dof
	for(auto bname: _geom_island->_articulated_bodies) {
		auto arb = _arb_manager->getBody(bname);
		dq_dof += arb->_model->dof();
	}

	// size Jacobian matrices
	_cop_full6_Jacobian.setZero(full6_contact_dof, dq_dof);
	_cop_constraint_Jacobian.setZero(cop_constraint_contact_dof, dq_dof);
	_pt_contact_Jacobian.setZero(pt_contact_dof, dq_dof);

	// fill in the Jacobians
	uint p_constr_J_start_ind = 0;
	uint p_ptct_J_start_ind = 0;
	uint arb_col_J_start_ind = 0;
	for (auto& p: _pair_state) {
		auto primA = p._geom_prim_pair->primA;
		auto primB = p._geom_prim_pair->primB;
		auto contact_points = p._geom_prim_pair->info.contact_points;
		if(primA->_is_static) {
			// get ARB for p->primB
			auto arb = _arb_manager->getBody(primB->_articulated_body_name);
			auto arb_model = arb->_model;
			if(_arb_index_map.find(primB->_articulated_body_name) != _arb_index_map.end()) {
				// do not extend arb_col_J_start_ind
			} else {
				_arb_index_map[primB->_articulated_body_name] = arb_col_J_start_ind;
				arb_col_J_start_ind += arb_model->dof();
			}
			uint this_arb_col_ind = _arb_index_map[primB->_articulated_body_name];

			MatrixXd Jv(3, arb_model->dof());
			MatrixXd Jw(3, arb_model->dof());
			// get body point for point 0
			// this is used for the cop Jacobian
			// TODO: should we consider using the centroid of the patch instead?
			Vector3d body_point0 = arb->worldToBodyPosition(primB->_link_name, contact_points[0]);

			arb_model->Jv(Jv, primB->_link_name, body_point0);
			arb_model->Jw(Jw, primB->_link_name);
			// rotate:
			Jv = p._rot_contact_frame_to_world.transpose()* arb_model->_T_world_robot.linear() * Jv;
			Jw = p._rot_contact_frame_to_world.transpose()* arb_model->_T_world_robot.linear() * Jw;

			// fill in the Jacobians
			_cop_full6_Jacobian.block(p._id*6, this_arb_col_ind, 3, arb_model->dof()) = Jv;
			_cop_full6_Jacobian.block(p._id*6 + 3, this_arb_col_ind, 3, arb_model->dof()) = Jw;
			
			_cop_constraint_Jacobian_prim_start_ind.push_back(p_constr_J_start_ind);
			_cop_constraint_Jacobian.block(p_constr_J_start_ind, this_arb_col_ind, 3, arb_model->dof()) = Jv;
			if(p._geom_prim_pair->info.type == ContactType::LINE) {
				_cop_constraint_Jacobian.block(p_constr_J_start_ind + 3, this_arb_col_ind, 2, arb_model->dof()) = Jw.block(1, 0, 2, arb_model->dof());
				p_constr_J_start_ind += 5;
			} else if(p._geom_prim_pair->info.type == ContactType::SURFACE) {
				_cop_constraint_Jacobian.block(p_constr_J_start_ind + 3, this_arb_col_ind, 3, arb_model->dof()) = Jw;
				p_constr_J_start_ind += 6;
			} else {
				p_constr_J_start_ind += 3;
			}

			_pt_contact_Jacobian_prim_start_ind.push_back(p_ptct_J_start_ind);
			for(auto& pt: contact_points) {
				Vector3d body_point = arb->worldToBodyPosition(primB->_link_name, pt);
				arb_model->Jv(Jv, primB->_link_name, body_point);
				Jv = p._rot_contact_frame_to_world.transpose()* arb_model->_T_world_robot.linear() * Jv;
				//TODO: do we need to save the start index for each primitive in the pt contact Jacobian?
				_pt_contact_Jacobian.block(p_ptct_J_start_ind, this_arb_col_ind, 3, arb_model->dof()) = Jv;
				p_ptct_J_start_ind += 3;
			}
		} else if(primB->_is_static) {
			// get ARB for p->primA
			auto arb = _arb_manager->getBody(primA->_articulated_body_name);
			auto arb_model = arb->_model;
			if(_arb_index_map.find(primA->_articulated_body_name) != _arb_index_map.end()) {
				// do not extend arb_col_J_start_ind
			} else {
				_arb_index_map[primA->_articulated_body_name] = arb_col_J_start_ind;
				arb_col_J_start_ind += arb_model->dof();
			}
			uint this_arb_col_ind = _arb_index_map[primA->_articulated_body_name];

			MatrixXd Jv(3, arb_model->dof());
			MatrixXd Jw(3, arb_model->dof());
			// get body point for point 0
			// this is used for the cop Jacobian
			Vector3d body_point0 = arb->worldToBodyPosition(primA->_link_name, contact_points[0]);
			// std::cout << "W point0 " << contact_points[0].transpose() << std::endl;
			// std::cout << "B point0 " << body_point0.transpose() << std::endl;

			arb_model->Jv(Jv, primA->_link_name, body_point0);
			// std::cout <<"Jv " << Jv <<std::endl;
			arb_model->Jw(Jw, primA->_link_name);
			// rotate:
			Jv = p._rot_contact_frame_to_world.transpose()* arb_model->_T_world_robot.linear() * Jv;
			Jw = p._rot_contact_frame_to_world.transpose()* arb_model->_T_world_robot.linear() * Jw;

			// fill in the Jacobians
			_cop_full6_Jacobian.block(p._id*6, this_arb_col_ind, 3, arb_model->dof()) = Jv;
			_cop_full6_Jacobian.block(p._id*6 + 3, this_arb_col_ind, 3, arb_model->dof()) = Jw;
			
			_cop_constraint_Jacobian_prim_start_ind.push_back(p_constr_J_start_ind);
			_cop_constraint_Jacobian.block(p_constr_J_start_ind, this_arb_col_ind, 3, arb_model->dof()) = Jv;
			if(p._geom_prim_pair->info.type == ContactType::LINE) {
				_cop_constraint_Jacobian.block(p_constr_J_start_ind + 3, this_arb_col_ind, 2, arb_model->dof()) = Jw.block(1, 0, 2, arb_model->dof());
				p_constr_J_start_ind += 5;
			} else if(p._geom_prim_pair->info.type == ContactType::SURFACE) {
				_cop_constraint_Jacobian.block(p_constr_J_start_ind + 3, this_arb_col_ind, 3, arb_model->dof()) = Jw;
				p_constr_J_start_ind += 6;
			} else {
				p_constr_J_start_ind += 3;
			}

			_pt_contact_Jacobian_prim_start_ind.push_back(p_ptct_J_start_ind);
			for(auto& pt: contact_points) {
				Vector3d body_point = arb->worldToBodyPosition(primA->_link_name, pt);
				arb_model->Jv(Jv, primA->_link_name, body_point);
				// std::cout <<"pt Jv " << Jv <<std::endl;
				Jv = p._rot_contact_frame_to_world.transpose()* arb_model->_T_world_robot.linear() * Jv;
				//TODO: do we need to save the start index for each primitive in the pt contact Jacobian?
				_pt_contact_Jacobian.block(p_ptct_J_start_ind, this_arb_col_ind, 3, arb_model->dof()) = Jv;
				p_ptct_J_start_ind += 3;
			}
		} else {
			// relative contact space coordinates directed from A to B
			// get ARBs for p->primA and p->primB
			auto arbA = _arb_manager->getBody(primA->_articulated_body_name);
			auto arbA_model = arbA->_model;
			auto arbB = _arb_manager->getBody(primB->_articulated_body_name);
			auto arbB_model = arbB->_model;
			if(_arb_index_map.find(primA->_articulated_body_name) != _arb_index_map.end()) {
				// do not extend arb_col_J_start_ind
			} else {
				_arb_index_map[primA->_articulated_body_name] = arb_col_J_start_ind;
				arb_col_J_start_ind += arbA_model->dof();
			}
			if(_arb_index_map.find(primB->_articulated_body_name) != _arb_index_map.end()) {
				// do not extend arb_col_J_start_ind
			} else {
				_arb_index_map[primB->_articulated_body_name] = arb_col_J_start_ind;
				arb_col_J_start_ind += arbB_model->dof();
			}
			uint arbA_col_ind = _arb_index_map[primA->_articulated_body_name];
			uint arbB_col_ind = _arb_index_map[primB->_articulated_body_name];

			MatrixXd JvA(3, arbA_model->dof());
			MatrixXd JwA(3, arbA_model->dof());
			MatrixXd JvB(3, arbB_model->dof());
			MatrixXd JwB(3, arbB_model->dof());
			// get body point for point 0
			// this is used for the cop Jacobian
			Vector3d bodyA_point0 = arbA->worldToBodyPosition(primA->_link_name, contact_points[0]);
			Vector3d bodyB_point0 = arbB->worldToBodyPosition(primB->_link_name, contact_points[0]);

			arbA_model->Jv(JvA, primA->_link_name, bodyA_point0);
			arbA_model->Jw(JwA, primA->_link_name);
			arbB_model->Jv(JvB, primB->_link_name, bodyB_point0);
			arbB_model->Jw(JwB, primB->_link_name);
			// rotate:
			JvA = p._rot_contact_frame_to_world.transpose()* arbA_model->_T_world_robot.linear() * JvA;
			JwA = p._rot_contact_frame_to_world.transpose()* arbA_model->_T_world_robot.linear() * JwA;
			JvB = p._rot_contact_frame_to_world.transpose()* arbB_model->_T_world_robot.linear() * JvB;
			JwB = p._rot_contact_frame_to_world.transpose()* arbB_model->_T_world_robot.linear() * JwB;

			// fill in the Jacobians
			// contact normal is directed from arb A to B
			_cop_full6_Jacobian.block(p._id*6, arbA_col_ind, 3, arbA_model->dof()) = -JvA;
			_cop_full6_Jacobian.block(p._id*6 + 3, arbA_col_ind, 3, arbA_model->dof()) = -JwA;
			_cop_full6_Jacobian.block(p._id*6, arbB_col_ind, 3, arbB_model->dof()) = JvB;
			_cop_full6_Jacobian.block(p._id*6 + 3, arbB_col_ind, 3, arbB_model->dof()) = JwB;
			
			_cop_constraint_Jacobian_prim_start_ind.push_back(p_constr_J_start_ind);
			_cop_constraint_Jacobian.block(p_constr_J_start_ind, arbA_col_ind, 3, arbA_model->dof()) = -JvA;
			_cop_constraint_Jacobian.block(p_constr_J_start_ind, arbB_col_ind, 3, arbB_model->dof()) = JvB;
			if(p._geom_prim_pair->info.type == ContactType::LINE) {
				_cop_constraint_Jacobian.block(p_constr_J_start_ind + 3, arbA_col_ind, 2, arbA_model->dof()) = -JwA.block(1, 0, 2, arbA_model->dof());
				_cop_constraint_Jacobian.block(p_constr_J_start_ind + 3, arbB_col_ind, 2, arbB_model->dof()) = JwB.block(1, 0, 2, arbB_model->dof());
				p_constr_J_start_ind += 5;
			} else if(p._geom_prim_pair->info.type == ContactType::SURFACE) {
				_cop_constraint_Jacobian.block(p_constr_J_start_ind + 3, arbA_col_ind, 3, arbA_model->dof()) = -JwA;
				_cop_constraint_Jacobian.block(p_constr_J_start_ind + 3, arbB_col_ind, 3, arbB_model->dof()) = JwB;
				p_constr_J_start_ind += 6;
			} else {
				p_constr_J_start_ind += 3;
			}

			_pt_contact_Jacobian_prim_start_ind.push_back(p_ptct_J_start_ind);
			for(auto& pt: contact_points) {
				Vector3d bodyA_point = arbA->worldToBodyPosition(primA->_link_name, pt);
				Vector3d bodyB_point = arbB->worldToBodyPosition(primB->_link_name, pt);
				arbA_model->Jv(JvA, primA->_link_name, bodyA_point);
				arbB_model->Jv(JvB, primB->_link_name, bodyB_point);
				JvA = p._rot_contact_frame_to_world.transpose()* arbA_model->_T_world_robot.linear() * JvA;
				JvB = p._rot_contact_frame_to_world.transpose()* arbB_model->_T_world_robot.linear() * JvB;
				//TODO: do we need to save the start index for each primitive in the pt contact Jacobian?
				_pt_contact_Jacobian.block(p_ptct_J_start_ind, arbA_col_ind, 3, arbA_model->dof()) = -JvA;
				_pt_contact_Jacobian.block(p_ptct_J_start_ind, arbB_col_ind, 3, arbB_model->dof()) = JvB;
				p_ptct_J_start_ind += 3;
			}
		}
	}

	// create lambda inv matrices.
	_cop_full6_Lambda_inv.setZero(full6_contact_dof, full6_contact_dof);
	_cop_constraint_Lambda_inv.setZero(cop_constraint_contact_dof, cop_constraint_contact_dof);
	_pt_contact_Lambda_inv.setZero(pt_contact_dof, pt_contact_dof);
	for(auto it: _arb_index_map) {
		uint arb_Jind = it.second;
		auto arb_model = _arb_manager->getBody(it.first)->_model;
		uint dq_dof = arb_model->dof();
		MatrixXd Jseg;

		// cop_full6
		Jseg = _cop_full6_Jacobian.block(0, arb_Jind, full6_contact_dof, dq_dof);
		_cop_full6_Lambda_inv += Jseg * arb_model->_M_inv * Jseg.transpose();

		// cop_constraint
		Jseg = _cop_constraint_Jacobian.block(0, arb_Jind, cop_constraint_contact_dof, dq_dof);
		_cop_constraint_Lambda_inv += Jseg * arb_model->_M_inv * Jseg.transpose();

		// pt contact
		Jseg = _pt_contact_Jacobian.block(0, arb_Jind, pt_contact_dof, dq_dof);
		_pt_contact_Lambda_inv += Jseg * arb_model->_M_inv * Jseg.transpose();
	}
}

void ContactIslandModel::updateRHSVectors() {
	// compute (- J A^-1 * non_linear_torques)
	uint full6_contact_dof = _cop_full6_Jacobian.rows();
	_cop_full6_rhs_contact.setZero(full6_contact_dof);

	uint cop_constraint_contact_dof = _cop_constraint_Jacobian.rows();
	_cop_constraint_rhs_contact.setZero(cop_constraint_contact_dof);

	uint pt_contact_dof = _pt_contact_Jacobian.rows();
	_pt_contact_rhs_coll.setZero(pt_contact_dof);

	// we assume here that model->nonLinearEffects has been called already after updating
	// the dq
	for(auto it: _arb_index_map) {
		uint arb_Jind = it.second;
		auto arb = _arb_manager->getBody(it.first);
		auto arb_model = arb->_model;

		uint dq_dof = arb_model->dof();
		MatrixXd Jseg;

		// cop_full6
		Jseg = _cop_full6_Jacobian.block(0, arb_Jind, full6_contact_dof, dq_dof);
		_cop_full6_rhs_contact += Jseg * (arb->jacc_nonlinear + arb_model->_M_inv*arb->jtau_act);
		// std::cout << _cop_full6_rhs_contact.transpose() << std::endl;

		// cop_constraint
		Jseg = _cop_constraint_Jacobian.block(0, arb_Jind, cop_constraint_contact_dof, dq_dof);
		_cop_constraint_rhs_contact += Jseg * (arb->jacc_nonlinear + arb_model->_M_inv*arb->jtau_act);

		// pt_ct for collision
		Jseg = _pt_contact_Jacobian.block(0, arb_Jind, pt_contact_dof, dq_dof);
		_pt_contact_rhs_coll += Jseg * arb_model->_dq; // collision RHS is simply dx-
	}

	// compute dot(J)*dq and add to RHS
	// TODO: cache the local body coordinates in the ContactPairState
	// Use _model->linearAccelerationInWorld at point i on Body B - _model->linearAcceleration at point j on Body A
	// TODO: J_contact = R_contact_frame * J_world
	// Consider the effects of ignoring dot(R_contact_frame) in dot(J_contact)
	for (auto& p: _pair_state) {
		auto primA = p._geom_prim_pair->primA;
		auto primB = p._geom_prim_pair->primB;
		auto contact_point0 = p._geom_prim_pair->info.contact_points[0];
		uint p_ind_start = _cop_constraint_Jacobian_prim_start_ind[p._id];
		if(primA->_is_static) {
			// get body location for primB
			auto arb = _arb_manager->getBody(primB->_articulated_body_name);
			auto arb_model = arb->_model;

			// get body point for point 0
			// TODO: cache the local body position for point in primitive
			Vector3d body_point0 = arb->worldToBodyPosition(primB->_link_name, contact_point0);
			// TODO: move update kinematics call elsewhere to avoid duplication
			// for the same articulated body multiple 
			
			// Note: updateKinematics must be called between 
			// model->nonLinearEffects() call and model->JvdotTimesQdotInWorld() call
			// because nonLinearEffects computes the internal body accelerations with
			// the gravity effect. we don't want that when computing dot(J)*dq
			arb_model->updateKinematics();
			Vector3d djv_times_dq, djw_times_dq;
			arb_model->JvdotTimesQdotInWorld(djv_times_dq, primB->_link_name, body_point0);
			arb_model->JwdotTimesQdotInWorld(djw_times_dq, primB->_link_name);
			// std::cout << "djvdq " << djv_times_dq.transpose() << std::endl;
			// std::cout << "djwdq " << djw_times_dq.transpose() << std::endl;

			// update full cop
			_cop_full6_rhs_contact.segment<3>(p._id*6) += p._rot_contact_frame_to_world.transpose() * djv_times_dq;
			_cop_full6_rhs_contact.segment<3>(p._id*6 + 3) += p._rot_contact_frame_to_world.transpose() * djw_times_dq;

			// update constraint cop
			if(p._geom_prim_pair->info.type == ContactType::LINE) {
				_cop_constraint_rhs_contact.segment<3>(p_ind_start) += p._rot_contact_frame_to_world.transpose() * djv_times_dq;
				_cop_constraint_rhs_contact.segment<2>(p_ind_start + 3) += (p._rot_contact_frame_to_world.transpose() * djw_times_dq).segment<2>(1);
			} else if(p._geom_prim_pair->info.type == ContactType::SURFACE) {
				_cop_constraint_rhs_contact.segment<3>(p_ind_start) += p._rot_contact_frame_to_world.transpose() * djv_times_dq;
				_cop_constraint_rhs_contact.segment<3>(p_ind_start + 3) += p._rot_contact_frame_to_world.transpose() * djw_times_dq;
			} else { // ContactType::Point
				_cop_constraint_rhs_contact.segment<3>(p_ind_start) += p._rot_contact_frame_to_world.transpose() * djv_times_dq;
			}
		} else if(primB->_is_static) {
			// get body location for primA
			auto arb = _arb_manager->getBody(primA->_articulated_body_name);
			auto arb_model = arb->_model;

			// get body point for point 0
			// TODO: cache the local body position for point in primitive
			Vector3d body_point0 = arb->worldToBodyPosition(primA->_link_name, contact_point0);
			// TODO: move update kinematics call elsewhere to avoid duplication
			// for the same articulated body multiple 
			
			// Note: updateKinematics must be called between 
			// model->nonLinearEffects() call and model->JvdotTimesQdotInWorld() call
			// because nonLinearEffects computes the internal body accelerations with
			// the gravity effect. we don't want that when computing dot(J)*dq
			arb_model->updateKinematics();
			Vector3d djv_times_dq, djw_times_dq;
			arb_model->JvdotTimesQdotInWorld(djv_times_dq, primA->_link_name, body_point0);
			arb_model->JwdotTimesQdotInWorld(djw_times_dq, primA->_link_name);
			// std::cout << "djvdq " << djv_times_dq.transpose() << std::endl;
			// std::cout << "djwdq " << djw_times_dq.transpose() << std::endl;

			// update full cop
			_cop_full6_rhs_contact.segment<3>(p._id*6) += p._rot_contact_frame_to_world.transpose() * djv_times_dq;
			_cop_full6_rhs_contact.segment<3>(p._id*6 + 3) += p._rot_contact_frame_to_world.transpose() * djw_times_dq;

			// update constraint cop
			if(p._geom_prim_pair->info.type == ContactType::LINE) {
				_cop_constraint_rhs_contact.segment<3>(p_ind_start) += p._rot_contact_frame_to_world.transpose() * djv_times_dq;
				_cop_constraint_rhs_contact.segment<2>(p_ind_start + 3) += (p._rot_contact_frame_to_world.transpose() * djw_times_dq).segment<2>(1);
			} else if(p._geom_prim_pair->info.type == ContactType::SURFACE) {
				_cop_constraint_rhs_contact.segment<3>(p_ind_start) += p._rot_contact_frame_to_world.transpose() * djv_times_dq;
				_cop_constraint_rhs_contact.segment<3>(p_ind_start + 3) += p._rot_contact_frame_to_world.transpose() * djw_times_dq;
			} else { // ContactType::Point
				_cop_constraint_rhs_contact.segment<3>(p_ind_start) += p._rot_contact_frame_to_world.transpose() * djv_times_dq;
			}
		} else { // !primA->_is_static && !primB->_is_static
			auto arbA = _arb_manager->getBody(primA->_articulated_body_name);
			auto arbA_model = arbA->_model;
			auto arbB = _arb_manager->getBody(primB->_articulated_body_name);
			auto arbB_model = arbB->_model;

			// get body point for point 0
			// TODO: cache the local body position for point in primitive
			Vector3d bodyA_point0 = arbA->worldToBodyPosition(primA->_link_name, contact_point0);
			Vector3d bodyB_point0 = arbB->worldToBodyPosition(primB->_link_name, contact_point0);
			// TODO: move update kinematics call elsewhere to avoid duplication
			// for the same articulated body multiple 
			
			// Note: updateKinematics must be called between 
			// model->nonLinearEffects() call and model->JvdotTimesQdotInWorld() call
			// because nonLinearEffects computes the internal body accelerations with
			// the gravity effect. we don't want that when computing dot(J)*dq
			Vector3d djv_times_dq_A, djw_times_dq_A, djv_times_dq_B, djw_times_dq_B;
			arbA_model->updateKinematics();
			arbA_model->JvdotTimesQdotInWorld(djv_times_dq_A, primA->_link_name, bodyA_point0);
			arbA_model->JwdotTimesQdotInWorld(djw_times_dq_A, primA->_link_name);
			arbB_model->updateKinematics();
			arbB_model->JvdotTimesQdotInWorld(djv_times_dq_B, primB->_link_name, bodyB_point0);
			arbB_model->JwdotTimesQdotInWorld(djw_times_dq_B, primB->_link_name);
			// std::cout << "djvdqA " << djv_times_dq_A.transpose() << std::endl;
			// std::cout << "djwdqA " << djw_times_dq_A.transpose() << std::endl;
			// std::cout << "djvdqB " << djv_times_dq_B.transpose() << std::endl;
			// std::cout << "djwdqB " << djw_times_dq_B.transpose() << std::endl;
			

			// update full cop RHS
			_cop_full6_rhs_contact.segment<3>(p._id*6) += p._rot_contact_frame_to_world.transpose() * (djv_times_dq_B - djv_times_dq_A);
			_cop_full6_rhs_contact.segment<3>(p._id*6 + 3) += p._rot_contact_frame_to_world.transpose() * (djw_times_dq_B - djw_times_dq_A);

			// update constraint cop RHS
			if(p._geom_prim_pair->info.type == ContactType::LINE) {
				_cop_constraint_rhs_contact.segment<3>(p_ind_start) += p._rot_contact_frame_to_world.transpose() * (djv_times_dq_B - djv_times_dq_A);
				_cop_constraint_rhs_contact.segment<2>(p_ind_start + 3) += (p._rot_contact_frame_to_world.transpose() * (djw_times_dq_B - djw_times_dq_A)).segment<2>(1);
			} else if(p._geom_prim_pair->info.type == ContactType::SURFACE) {
				_cop_constraint_rhs_contact.segment<3>(p_ind_start) += p._rot_contact_frame_to_world.transpose() * (djv_times_dq_B - djv_times_dq_A);
				_cop_constraint_rhs_contact.segment<3>(p_ind_start + 3) += p._rot_contact_frame_to_world.transpose() * (djw_times_dq_B - djw_times_dq_A);
			} else { // ContactType::Point
				_cop_constraint_rhs_contact.segment<3>(p_ind_start) += p._rot_contact_frame_to_world.transpose() * (djv_times_dq_B - djv_times_dq_A);
			}
		}
	}
}

void ContactIslandModel::updatePtContactRHSCollVector() {
	// collision RHS is simply dx-
	uint pt_contact_dof = _pt_contact_Jacobian.rows();
	_pt_contact_rhs_coll.setZero(pt_contact_dof);

	MatrixXd Jseg;
	for(auto it: _arb_index_map) {
		uint arb_Jind = it.second;
		auto arb = _arb_manager->getBody(it.first);
		auto arb_model = arb->_model;

		uint dq_dof = arb_model->dof();

		// pt_ct for collision
		Jseg = _pt_contact_Jacobian.block(0, arb_Jind, pt_contact_dof, dq_dof);
		_pt_contact_rhs_coll += Jseg * arb_model->_dq; // collision RHS is simply dx-
	}
}

void ContactIslandModel::updateBodyNonlinAccelerations() {
	for(auto it: _arb_index_map) {
		auto arb = _arb_manager->getBody(it.first);
		arb->updateNonLinearJAcc();
	}
}

void ContactIslandModel::getActiveFullCOPMatrices(
		Eigen::MatrixXd& J_full_cop,
		Eigen::MatrixXd& Lambda_inv_full_cop,
		Eigen::VectorXd& rhs_full_cop,
		std::vector<uint>& Jrow_ind_to_contact_pair_map 
		//TODO: ^consider making this a map from ContactPair.id -> row start ind in active Jacobian
) const {
	if(_active_contacts.size() == 0) return;

	J_full_cop.setZero(6*_active_contacts.size(), _cop_full6_Jacobian.cols());
	Lambda_inv_full_cop.setZero(6*_active_contacts.size(), 6*_active_contacts.size());
	rhs_full_cop.setZero(6*_active_contacts.size());
	int active_contacts_ind = 0;
	Jrow_ind_to_contact_pair_map.clear();
	
	// fill in Jacobian and rhs
	for(uint cid: _active_contacts) {
		J_full_cop.block(active_contacts_ind*6, 0, 6, _cop_full6_Jacobian.cols()) = _cop_full6_Jacobian.block(cid*6, 0, 6, _cop_full6_Jacobian.cols());
		rhs_full_cop.segment(active_contacts_ind*6, 6) = _cop_full6_rhs_contact.segment(cid*6, 6);
		Jrow_ind_to_contact_pair_map.push_back(active_contacts_ind*6);
		active_contacts_ind++;
	}
	// fill in lambda_inv
	uint i = 0;
	for(auto cind1: _active_contacts) {
		uint row_start = Jrow_ind_to_contact_pair_map[i];
		i++;
		uint j = 0;
		for(auto cind2: _active_contacts) {
			uint col_start = Jrow_ind_to_contact_pair_map[j];
			j++;
			Lambda_inv_full_cop.block(row_start, col_start, 6, 6) = _cop_full6_Lambda_inv.block(cind1*6, cind2*6, 6, 6);
		}
	}
}

void ContactIslandModel::getActiveConstraintCOPMatrices(
		Eigen::MatrixXd& J_constraint_cop,
		Eigen::MatrixXd& Lambda_inv_constraint_cop,
		Eigen::VectorXd& rhs_constraint_cop,
		std::vector<uint>& Jrow_ind_to_contact_pair_map
) const {
	if(_active_contacts.size() == 0) return;

	uint active_J_constraint_cop_dof = 0;
	std::vector<uint> dof_count_map;
	for(uint cid: _active_contacts) {
		if(_pair_state[cid]._geom_prim_pair->info.type == ContactType::POINT) {
			active_J_constraint_cop_dof += 3;
			dof_count_map.push_back(3);
		} else if(_pair_state[cid]._geom_prim_pair->info.type == ContactType::LINE) {
			active_J_constraint_cop_dof += 5;
			dof_count_map.push_back(5);
		} else if(_pair_state[cid]._geom_prim_pair->info.type == ContactType::SURFACE) {
			active_J_constraint_cop_dof += 6;
			dof_count_map.push_back(6);
		}
	}

	J_constraint_cop.setZero(active_J_constraint_cop_dof, _cop_constraint_Jacobian.cols());
	Lambda_inv_constraint_cop.setZero(active_J_constraint_cop_dof, active_J_constraint_cop_dof);
	rhs_constraint_cop.setZero(active_J_constraint_cop_dof);
	int active_contacts_ind = 0;
	Jrow_ind_to_contact_pair_map.clear(); // this is to the start of the constraint rows

	// fill in Jacobian and rhs
	uint active_prim_start_Jrow_ind = 0;
	for(uint cid: _active_contacts) {
		uint J_start_ind = _cop_constraint_Jacobian_prim_start_ind[cid];
		uint dof_count = dof_count_map[active_contacts_ind];
		J_constraint_cop.block(active_prim_start_Jrow_ind, 0, dof_count, _cop_constraint_Jacobian.cols()) = _cop_constraint_Jacobian.block(J_start_ind, 0, dof_count, _cop_constraint_Jacobian.cols());
		rhs_constraint_cop.segment(active_prim_start_Jrow_ind, dof_count) = _cop_constraint_rhs_contact.segment(J_start_ind, dof_count);
		Jrow_ind_to_contact_pair_map.push_back(active_prim_start_Jrow_ind);
		active_prim_start_Jrow_ind += dof_count;
		active_contacts_ind++;
	}
	// fill in lambda_inv
	uint i = 0;
	for(auto cind1: _active_contacts) {
		uint Jstart1 = _cop_constraint_Jacobian_prim_start_ind[cind1];
		uint dof_count1 = dof_count_map[i];
		uint Jactive_start1 = Jrow_ind_to_contact_pair_map[i];
		i++;
		uint j = 0;
		for(auto cind2: _active_contacts) {
			uint Jstart2 = _cop_constraint_Jacobian_prim_start_ind[cind2];
			uint dof_count2 = dof_count_map[j];
			uint Jactive_start2 = Jrow_ind_to_contact_pair_map[j];
			j++;
			Lambda_inv_constraint_cop.block(Jactive_start1, Jactive_start2, dof_count1, dof_count2) = _cop_constraint_Lambda_inv.block(Jstart1, Jstart2, dof_count1, dof_count2);
		}
	}
}

void ContactIslandModel::getActivePtContactCollisionMatrices(
		Eigen::MatrixXd& J_pt_contacts,
		Eigen::MatrixXd& Lambda_inv_pt_contacts,
		Eigen::VectorXd& rhs_pt_contacts_collision,
		std::vector<uint>& Jrow_ind_to_contact_pair_map
) const {
	if(_active_contacts.size() == 0) return;

	// number of active points per contact pair
	// std::vector<uint> contact_num_active_pts;
	uint pt_contact_active_dof = 0;
	for(uint cind: _active_contacts) {
		// contact_num_active_pts.push_back(_active_points.size());
		pt_contact_active_dof += 3*_pair_state[cind]._active_points.size();
	}
	// Take into account which points are active
	J_pt_contacts.setZero(pt_contact_active_dof, _pt_contact_Jacobian.cols());
	Lambda_inv_pt_contacts.setZero(pt_contact_active_dof, pt_contact_active_dof);
	rhs_pt_contacts_collision.setZero(pt_contact_active_dof);

	Jrow_ind_to_contact_pair_map.clear();
	
	// fill in Jacobian and rhs
	uint row_ind = 0;
	for(uint cid: _active_contacts) {
		Jrow_ind_to_contact_pair_map.push_back(row_ind);
		uint Jstart_ind = _pt_contact_Jacobian_prim_start_ind[cid];
		for(uint pid: _pair_state[cid]._active_points) {
			J_pt_contacts.block(row_ind, 0, 3, _pt_contact_Jacobian.cols()) = _pt_contact_Jacobian.block(Jstart_ind + pid*3, 0, 3, _pt_contact_Jacobian.cols());
			rhs_pt_contacts_collision.segment(row_ind, 3) = _pt_contact_rhs_coll.segment(Jstart_ind + pid*3, 3);
			row_ind += 3;
		}
	}
	// fill in lambda_inv
	row_ind = 0;
	for(auto cind1: _active_contacts) {
		uint Jstart1_ind = _pt_contact_Jacobian_prim_start_ind[cind1];
		for(auto pind1: _pair_state[cind1]._active_points) {
			uint col_ind = 0;
			for(auto cind2: _active_contacts) {
				uint Jstart2_ind = _pt_contact_Jacobian_prim_start_ind[cind2];
				for(auto pind2: _pair_state[cind2]._active_points) {
					Lambda_inv_pt_contacts.block(row_ind, col_ind, 3, 3) = _pt_contact_Lambda_inv.block(Jstart1_ind + pind1*3, Jstart2_ind + pind2*3, 3, 3);
					col_ind += 3;
				}
			}
			row_ind += 3;
		}
	}
}

void ContactIslandModel::getActivePtContactCollisionRHSVector(
		Eigen::VectorXd& rhs_pt_contacts_collision
) const {
	if(_active_contacts.size() == 0) return;

	// NOTE: We assume that the active contacts have not changed since the last call to getActivePtContactCollisionMatrices
	// Therefore, we assume that rhs_pt_contacts_collision is already the correct size
	rhs_pt_contacts_collision.setZero(_pt_contact_Jacobian_active.rows());

	uint row_ind = 0;
	// std::cout << _pt_contact_rhs_coll.size() <<" "<< rhs_pt_contacts_collision.size() <<std::endl;
	for(uint cid: _active_contacts) {
		uint Jstart_ind = _pt_contact_Jacobian_prim_start_ind[cid];
		// std::cout << cid << " " << Jstart_ind << std::endl;
		for(uint pid: _pair_state[cid]._active_points) {
			rhs_pt_contacts_collision.segment(row_ind, 3) = _pt_contact_rhs_coll.segment(Jstart_ind + pid*3, 3);
			row_ind += 3;
		}
	}
}

void ContactIslandModel::getActiveFullCOPRHSVector(
	Eigen::VectorXd& rhs_full_cop,
	std::vector<uint>& Jrow_ind_to_contact_pair_map
) const {
	// we assume that the correct size has been set already for Jacobian
	rhs_full_cop.setZero(_cop_full6_Jacobian_active.rows());

	int active_contacts_ind = 0;
	Jrow_ind_to_contact_pair_map.clear();

	// fill in rhs
	for(uint cid: _active_contacts) {
		rhs_full_cop.segment(active_contacts_ind*6, 6) = _cop_full6_rhs_contact.segment(cid*6, 6);
		Jrow_ind_to_contact_pair_map.push_back(active_contacts_ind*6);
		active_contacts_ind++;
	}
}

void ContactIslandModel::getActiveConstraintCOPRHSVector(
	Eigen::VectorXd& rhs_constraint_cop,
	std::vector<uint>& Jrow_ind_to_contact_pair_map
) const {
	// we assume that the correct size has been set already for Jacobian
	rhs_constraint_cop.setZero(_cop_constraint_Jacobian_active.rows());

	int active_contacts_ind = 0;
	Jrow_ind_to_contact_pair_map.clear(); // this is to the start of the constraint rows
	
	// fill in Jacobian and rhs
	uint active_prim_start_Jrow_ind = 0;
	for(uint cid: _active_contacts) {
		uint J_start_ind = _cop_constraint_Jacobian_prim_start_ind[cid];
		uint dof_count = 0;
		if(_pair_state[cid]._geom_prim_pair->info.type == ContactType::POINT) {
			dof_count = 3;
		} else if(_pair_state[cid]._geom_prim_pair->info.type == ContactType::LINE) {
			dof_count = 5;
		} else if(_pair_state[cid]._geom_prim_pair->info.type == ContactType::SURFACE) {
			dof_count = 6;
		}
		rhs_constraint_cop.segment(active_prim_start_Jrow_ind, dof_count) = _cop_constraint_rhs_contact.segment(J_start_ind, dof_count);
		Jrow_ind_to_contact_pair_map.push_back(active_prim_start_Jrow_ind);
		active_prim_start_Jrow_ind += dof_count;
		active_contacts_ind++;
	}
}

bool ContactIslandModel::activateContactPair(uint id) {
	assert(id <= _pair_state.size());
	for(auto cind: _active_contacts) {
		if(cind == id) {
			return false;
		}
	}
	_active_contacts.push_back(id);
	return true;
}

bool ContactIslandModel::deactivateContactPair(uint id) {
	assert(id <= _pair_state.size());
	uint pre_size = _active_contacts.size();
	_active_contacts.remove(id);
	return pre_size != _active_contacts.size();
}

/* ------ Contact Space Model ------- */
ContactSpaceModel::ContactSpaceModel(ArticulatedRigidBodyManager* arb_manager)
: _arb_manager(arb_manager)
{
	assert(arb_manager != NULL);
}

ContactSpaceModel::~ContactSpaceModel() {

}

void ContactSpaceModel::build(const WorldContactMap* geom_map) {
	clear();
	// loop over islands in the geometric contact map
	for(auto geom_island_it = geom_map->_islands.begin(); geom_island_it != geom_map->_islands.end(); geom_island_it++) {
		// create a ContactIslandModel
		_contact_island_models.push_back(
				ContactIslandModel(&(*geom_island_it), _arb_manager)
		);
	}
}

void ContactSpaceModel::updateVelocityTerms() {
	// update body non-linear accelerations
	for(auto arb_it: _arb_manager->_articulated_bodies) {
		auto arb = arb_it.second;
		arb->updateNonLinearJAcc();
	}

	// update RHS vectors for each island
	for(auto& island: _contact_island_models) {
		island.updateRHSVectors();
	}
}

void ContactSpaceModel::resolveCollisions(double friction_coeff, double restitution_coeff) {
	// resolve collision for each island
	for(auto& island: _contact_island_models) {
		//TODO: implement island sleep
		island.resolveCollisions(friction_coeff, restitution_coeff);
	}
}

void ContactSpaceModel::resolveSteadyContacts(double friction_coeff, double restitution_coeff) {
	// resolve steady contact for each island
	for(auto& island: _contact_island_models) {
		//TODO: implement island sleep
		island.resolveSteadyContacts(friction_coeff, restitution_coeff);
	}
}

/* ------ Contact Pair State ------- */
ContactPairState::ContactPairState(uint id, const ContactPrimitivePair* geom_prim_pair)
: _id(id), _geom_prim_pair(geom_prim_pair)
{
	assert(geom_prim_pair != NULL);

	// compute the rotation matrix from the directions in the geom_prim_pair
	_rot_contact_frame_to_world.col(0) = geom_prim_pair->info.constraint_dir1;
	_rot_contact_frame_to_world.col(1) = geom_prim_pair->info.constraint_dir2;
	//TODO: the above is dangerous. For a point contact, the primitive check might
	// accidentally not set a proper direction. then the rotation matrix can be completely
	// wrong. consider computing the rotation matrix directly in ContactPairState.
	// TODO: is this robust? does the rotation matrix flip wildly sometimes? e.g. capsule
	// becomes vertical and falls?. We should consider caching and using the rotation matrix
	// from the last time if possible, especially for point contact
	_rot_contact_frame_to_world.col(2) = geom_prim_pair->info.normal_dir;
	assert(abs(1.0 - _rot_contact_frame_to_world.determinant()) < 1e-5);

	// save number of points
	_num_contact_pts = geom_prim_pair->info.contact_points.size();

	// add initial active points
	for(uint i = 0; i < _num_contact_pts; i++) {
		_active_points.push_back(i);
	}
}

ContactPairState::~ContactPairState() {

}

bool ContactPairState::activateContactPoint(uint id) {
	assert(id <= _geom_prim_pair->info.contact_points.size());
	for(auto cind: _active_points) {
		if(cind == id) {
			return false;
		}
	}
	_active_points.push_back(id);
	return true;
}

bool ContactPairState::deactivateContactPoint(uint id) {
	assert(id <= _geom_prim_pair->info.contact_points.size());
	uint pre_size = _active_points.size();
	_active_points.remove(id);
	return pre_size != _active_points.size();
}

}
