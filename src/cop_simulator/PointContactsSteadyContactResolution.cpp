// PointContactSteadyContacts.cpp

#include <iostream>
#include "ContactSpaceModel.h"
#include "lcp_solvers/LCPSolver2.h"

using namespace Eigen;

namespace Sai2COPSim {

bool ContactIslandModel::resolveSteadyContactsWithPointContacts(
	double friction_coeff,
	double restitution_coeff
) {
	if(_active_contacts.size() == 0) {
		// std::cout<< "No active contacts" << std::endl;
		return false;
	}
	if(!_f_support_pt_contact_steady_contact) {
		return false;
	}

	VectorXd pt_contacts_rhs;
	getActivePtContactSteadyContactRHSVector(pt_contacts_rhs);

	VectorXd pre_V;
	pre_V.setZero(_pt_contact_Lambda_inv_active.rows());
	MatrixXd Jseg;

	uint row_ind = 0;
	uint prim_row_start_ind = 0;
	for(uint prim_id: _active_contacts) {
		prim_row_start_ind = row_ind;
		auto& prim = _pair_state[prim_id];
		auto primA = prim._geom_prim_pair->primA;
		auto primB = prim._geom_prim_pair->primB;
		std::string arbA_name = "";
		std::string arbB_name = "";
		if(primA->_is_static) {
			arbB_name = primB->_articulated_body_name;
		} else if (primB->_is_static) {
			arbB_name = primA->_articulated_body_name;
		} else {
			arbA_name = primA->_articulated_body_name;
			arbB_name = primB->_articulated_body_name;
		}
		uint Jstart_ind = _pt_contact_Jacobian_prim_start_ind[prim_id];

		{ // ArbB
			row_ind = prim_row_start_ind;
			uint arbB_ind = _arb_index_map[arbB_name];
			ArticulatedRigidBody* arbB = _arb_manager->getBody(arbB_name);
			uint dqB_dof = arbB->_model->dof();
			for(uint pid: prim._active_points) {
				Jseg = _pt_contact_Jacobian.block(Jstart_ind + pid*3, arbB_ind, 3, dqB_dof);
				pre_V.segment(row_ind, 3) += Jseg*arbB->_model->_dq;
				row_ind += 3;
			}
		}

		if(!arbA_name.empty()) {
			row_ind = prim_row_start_ind;
			uint arbA_ind = _arb_index_map[arbA_name];
			ArticulatedRigidBody* arbA = _arb_manager->getBody(arbA_name);
			uint dqA_dof = arbA->_model->dof();
			for(uint pid: prim._active_points) {
				Jseg = _pt_contact_Jacobian.block(Jstart_ind + pid*3, arbA_ind, 3, dqA_dof);
				pre_V.segment(row_ind, 3) += Jseg*arbA->_model->_dq;
				row_ind += 3;
			}
		}
	}

	auto solver = Sai2LCPSolver::LCPSolver();
	// std::cout << _pt_contact_Lambda_inv_active << std::endl;
	// std::cout << pt_contact_rhs_coll_active.transpose() << std::endl;
	// std::cout << adjusted_restitution_coeff <<std::endl;
	auto lcp_sol = solver.solve(
		_pt_contact_Lambda_inv_active,
		pt_contacts_rhs,
		pre_V,
		0, /* restitution*/
		friction_coeff,
		true /* force_sliding_if_pre_slip */
	);
	if(lcp_sol.result != LCPSolResult::Success) {
		throw(std::runtime_error("Pt contacts Steady Contact LCP failed"));
	}

	// NOTE: It is assumed that arb->jtau_contact is set to zero before
	// resolveSteadyContactsWithPointContacts is called

	for(uint prim_id: _active_contacts) {
		prim_row_start_ind = row_ind;
		auto& prim = _pair_state[prim_id];
		auto primA = prim._geom_prim_pair->primA;
		auto primB = prim._geom_prim_pair->primB;
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
		uint Jstart_ind = _pt_contact_Jacobian_prim_start_ind[prim_id];

		if(arbB != NULL) {
			row_ind = prim_row_start_ind;
			uint arbB_ind = _arb_index_map[arbB->_name];
			uint dqB_dof = arbB->_model->dof();
			for(uint pid: prim._active_points) {
				Jseg = _pt_contact_Jacobian.block(Jstart_ind + pid*3, arbB_ind, 3, dqB_dof);
				arbB->jtau_contact += Jseg.transpose()*lcp_sol.p_sol;
				row_ind += 3;
			}
		}

		if(arbA != NULL) {
			row_ind = prim_row_start_ind;
			uint arbA_ind = _arb_index_map[arbA->_name];
			uint dqA_dof = arbA->_model->dof();
			for(uint pid: prim._active_points) {
				Jseg = _pt_contact_Jacobian.block(Jstart_ind + pid*3, arbA_ind, 3, dqA_dof);
				arbA->jtau_contact += Jseg.transpose()*lcp_sol.p_sol;
				row_ind += 3;
			}
		}

		// TODO: update primitive last COP sol force and COP
	}

	return true;
}

}