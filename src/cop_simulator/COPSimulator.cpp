// COPSimulator.cpp
#include <iostream>
#include "COPSimulator.h"
#include <rbdl/rbdl.h>

using namespace Eigen;

namespace Sai2COPSim {
// constructor
COPSimulator::COPSimulator(double friction_coeff, double restitution_coeff)
:_friction_coeff(friction_coeff), _resitution_coeff(restitution_coeff), _iterations(0)
{
	_contact_model = new ContactSpaceModel(&_arb_manager);
}

// destructor
COPSimulator::~COPSimulator() {
}

// add an inifinite plane object
void COPSimulator::addPlane(const std::string& name, const Eigen::Vector3d& planeNormal, const Eigen::Vector3d& planePoint) {
	// check if planeNormal is aligned with +Z. Fail otherwise
	if(_static_objects.find(name) != _static_objects.end()) {
		std::cerr << name << std::endl;
		throw(std::runtime_error("Plane with name already exists"));
	}
	_static_objects[name] = new PlanePrimitive(name, planeNormal, planePoint);
	_static_objects[name]->_is_static = true;
	_geom_manager.addPrimitive(_static_objects[name]);

	//TODO: check if intersecting with articulated bodies, if so throw
	//NOTE: It is ok if there is intersection with other static objects
	//TODO: force update model
}

// add a capsule object
// object->_T_world_robot is assumed to have been set already
// capsule is assumed to be aligned with X axis. local frame at the center of the capsule.
void COPSimulator::addCapsuleObject(const std::string& articulated_body_name,
					const std::string& link_name, // name for the link on which capsule primitive will be attached
					const std::string& primitive_name,
					Sai2Model::Sai2Model* object,
					double radius,
					double length
) {
	auto arb = new ArticulatedRigidBody(articulated_body_name, object);
	_arb_manager.addBody(arb);

	// create a new capsule primitive
	auto capsule = new CapsulePrimitive(primitive_name, radius, length);
	capsule->_is_static = false;
	capsule->_articulated_body_name = articulated_body_name;
	capsule->_link_name = link_name;
	// TODO: set primitive _transform_in_link

	_geom_manager.addPrimitive(capsule);

	arb->addPrimitive(link_name, capsule);
	//TODO: check if intersecting with anything, if so throw
	//TODO: force update model
}

void COPSimulator::addCylinderObject(const std::string& articulated_body_name,
					const std::string& link_name,
					const std::string& primitive_name,
					Sai2Model::Sai2Model* object,
					double radius,
					double length,
					uint num_points
) {
	auto arb = new ArticulatedRigidBody(articulated_body_name, object);
	_arb_manager.addBody(arb);

	// create a new cylinder primitive
	auto cylinder = new CylinderPrimitive(primitive_name, radius, length, num_points);
	cylinder->_is_static = false;
	cylinder->_articulated_body_name = articulated_body_name;
	cylinder->_link_name = link_name;
	cylinder->_transform_in_link.translation() << 0, 0, -length/2; 
	// to ensure that link center is at the center of the cylinder

	_geom_manager.addPrimitive(cylinder);

	arb->addPrimitive(link_name, cylinder);
	//TODO: check if intersecting with anything, if so throw
	//TODO: force update model
}

// automatically sets q and dq for the objects
void COPSimulator::integrate(double dt) {
	bool f_force_update_dynamics = false;
	//TODO: this might be an internal property of the class, in case we need to 
	// set the flag from further below, such as after resolveCollisions or resolveSteadyContacts
	// check for collisions
	if(_iterations % COPAlgorithmicConstants::NUM_ITERS_BEFORE_COLLISION_CHECK_UPDATE == 0) {
		for(auto arb_it: _arb_manager._articulated_bodies) {
			auto arb = arb_it.second;
			arb->_model->updateKinematics();
		}
		computeWorldContactMap();

		f_force_update_dynamics = true;
	} else {
		// TODO: delta updates in geometry positions
	}
	// update contact dynamics
	if(_iterations % COPAlgorithmicConstants::NUM_ITERS_BEFORE_INERTIA_AND_CONTACT_MODEL_UPDATE == 0 || f_force_update_dynamics) {
		for(auto arb_it: _arb_manager._articulated_bodies) {
			auto arb = arb_it.second;
			if(f_force_update_dynamics) {
				// kinematics has just been updated, so just update dynamics
				arb->_model->updateDynamics();
			} else {
				// update everything
				arb->_model->updateModel();
			}
			arb->updateNonLinearJAcc();
		}
		// rebuild the contact model
		_contact_model->build(&_contact_map);
		// TODO: consider persisting the COP solution? If so, rotate and displace it to the current contact patch
	}

	// update nonlinear accelerations
	if(_iterations % COPAlgorithmicConstants::NUM_ITERS_BEFORE_NONLINEAR_ACCERLATION_UPDATE == 0 && !f_force_update_dynamics) {
		// unnecessary if f_force_update_dynamics was called
		_contact_model->updateVelocityTerms();
	}

	// resolve collisions
	_contact_model->resolveCollisions(_friction_coeff, _resitution_coeff);

	// TODO: integrate with substepping, including recomputation of the contact forces
	// resolve steady contacts
	_contact_model->resolveSteadyContacts(_friction_coeff, _resitution_coeff);

	for(auto arb_it: _arb_manager._articulated_bodies) {
		auto arb = arb_it.second;
		auto model = arb->_model;
		auto rbdl = model->_rbdl_model;

		model->_ddq = model->_M_inv * (arb->jtau_contact + arb->jtau_act) + arb->jacc_nonlinear;

		Eigen::VectorXd half_dq_update(model->dof());
		half_dq_update = model->_dq + 0.5*model->_ddq * dt;

		Eigen::VectorXd last_q = model->_q;
		
		// for each joint
		for(uint body_ind = 0; body_ind < rbdl->mJoints.size(); ++body_ind) {
			// NOTE: the iteration starts at body_ind because the first joint is a root joint that connects
			// the base body
			const auto& joint = rbdl->mJoints[body_ind];
			if(joint.mJointType == RigidBodyDynamics::JointTypeUndefined) {
				continue;
			} else if(joint.mJointType == RigidBodyDynamics::JointTypeSpherical) {
				// - if it is spherical,
				if (half_dq_update.segment<3>(joint.q_index).norm() > 1e-15) {
					// - - get current quaternion
					auto jquat = rbdl->GetQuaternion(body_ind, last_q);
					// - - get body orientation for parent link in arb base
					Matrix3d body_orientation_in_base;
					body_orientation_in_base = CalcBodyWorldOrientation(*rbdl, last_q, body_ind, false).transpose();
					// - - update quaternion
					jquat = jquat.timeStep(body_orientation_in_base * half_dq_update.segment<3>(joint.q_index), dt);
					// - - set back to _q
					rbdl->SetQuaternion(body_ind, jquat, model->_q);
				} else {
					// TODO: some first order approximation?
				}
	    	} else {
	    		// - if it is not spherical, just integrate
	    		model->_q.segment(joint.q_index, joint.mDoFCount) += half_dq_update.segment(joint.q_index, joint.mDoFCount) * dt;
	    	}
		}
		model->_dq += model->_ddq * dt;

		// clear the contact torques for this body
		arb->jtau_contact.setZero(model->dof());

		// std::cout << arb_it.first << " q:" << model->_q.transpose() << std::endl;
	}
	_iterations++;	
}

}