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

	_time_total = 0.0;
	_time_update_geometry = 0.0;
	_time_update_dynamics = 0.0;
	_time_update_velocity_terms = 0.0;
	_time_resolve_collisions = 0.0;
	_time_resolve_steady_contact = 0.0;
	_time_integrate = 0.0;
	_time_update_dynamics_rbdl = 0.0;
	_time_update_dynamics_contactmap = 0.0;
	_time_bcm_compute_matrices = 0.0;
	_time_bcm_copy_matrices = 0.0;
	_time_bcm_compute_vectors = 0.0;
	_time_bcm_ci_total = 0.0;
	_time_bcm_build = 0.0;
	_iters_geom_update = 0;
	_iters_dynamics_update = 0;
	_iters_velterms_update = 0;
	_iters_resolve_colls = 0;
	_iters_resolve_steady_contact = 0;
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
	// time point 1
	auto time_pt1 = std::chrono::high_resolution_clock::now();
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
		_iters_geom_update++;
	} else {
		// TODO: delta updates in geometry positions
	}
	auto time_pt2 = std::chrono::high_resolution_clock::now();
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
		auto time_pt2a = std::chrono::high_resolution_clock::now();
		_time_update_dynamics_rbdl += std::chrono::duration<double>(time_pt2a - time_pt2).count();

		// rebuild the contact model
		_contact_model->build(&_contact_map);
		// TODO: consider persisting the COP solution? If so, rotate and displace it to the current contact patch
		_iters_dynamics_update++;
		auto time_pt2b = std::chrono::high_resolution_clock::now();
		_time_update_dynamics_contactmap += std::chrono::duration<double>(time_pt2b - time_pt2a).count();
		_time_bcm_compute_matrices += _contact_model->_time_compute_matrices;
		_time_bcm_copy_matrices += _contact_model->_time_copy_matrices;
		_time_bcm_compute_vectors += _contact_model->_time_compute_vectors;
		_time_bcm_build += _contact_model->_time_build;
		_time_bcm_ci_total += _contact_model->_time_ci_total;
	}
	auto time_pt3 = std::chrono::high_resolution_clock::now();

	// update nonlinear accelerations
	if(_iterations % COPAlgorithmicConstants::NUM_ITERS_BEFORE_NONLINEAR_ACCERLATION_UPDATE == 0 && !f_force_update_dynamics) {
		// unnecessary if f_force_update_dynamics was called
		_contact_model->updateVelocityTerms();
		_iters_velterms_update++;
	}
	auto time_pt4 = std::chrono::high_resolution_clock::now();

	// resolve collisions
	bool was_colliding = _contact_model->resolveCollisions(_friction_coeff, _resitution_coeff);
	if(was_colliding) _iters_resolve_colls++;

	auto time_pt5 = std::chrono::high_resolution_clock::now();

	// TODO: integrate with substepping, including recomputation of the contact forces
	// resolve steady contacts
	bool was_in_contact = _contact_model->resolveSteadyContacts(_friction_coeff, _resitution_coeff);
	if(was_in_contact) _iters_resolve_steady_contact++;
	auto time_pt6 = std::chrono::high_resolution_clock::now();

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
	auto time_pt7 = std::chrono::high_resolution_clock::now();
	_iterations++;

	// update timings
	_time_total += std::chrono::duration<double>(time_pt7 - time_pt1).count();
	_time_update_geometry += std::chrono::duration<double>(time_pt2 - time_pt1).count();
	_time_update_dynamics += std::chrono::duration<double>(time_pt3 - time_pt2).count();
	_time_update_velocity_terms += std::chrono::duration<double>(time_pt4 - time_pt3).count();
	_time_resolve_collisions += std::chrono::duration<double>(time_pt5 - time_pt4).count();
	_time_resolve_steady_contact += std::chrono::duration<double>(time_pt6 - time_pt5).count();
	_time_integrate += std::chrono::duration<double>(time_pt7 - time_pt6).count();
}

void COPSimulator::printTimeAnalytics() {
	using namespace std;
	cout << "Total sim time " << _time_total << endl;
    cout << "Average sim time (ms/1000 calls) " << (_time_total*1e6)/_iterations << endl;
    cout << "Total geometry update time " << _time_update_geometry << endl;
    cout << "Average geom update time (ms/1000 calls) " << (_time_update_geometry*1e6)/_iters_geom_update << endl;
    cout << "Total collision resolution time " << _time_resolve_collisions << endl;
    cout << "Average collision resolution time (ms/1000 calls)" << (_time_resolve_collisions*1e6)/_iters_resolve_colls << endl;
    cout << "Total contact resolution time " << _time_resolve_steady_contact << endl;
    cout << "Average contact resolution time (ms/1000 calls)" << (_time_resolve_steady_contact*1e6)/_iters_resolve_steady_contact << endl;
    cout << "Total update dynamics time " << _time_update_dynamics << endl;
    cout << "Average update dynamics time (ms/1000 calls)" << (_time_update_dynamics*1e6)/_iters_dynamics_update << endl;
    cout << "Total update dynamics (rbdl) time " << _time_update_dynamics_rbdl << endl;
    cout << "Average update dynamics (rbdl) time (ms/1000 calls)" << (_time_update_dynamics_rbdl*1e6)/_iters_dynamics_update << endl;
    cout << "Total update dynamics (contact map) time " << _time_update_dynamics_contactmap << endl;
    cout << "Average update dynamics (contact map) time (ms/1000 calls)" << (_time_update_dynamics_contactmap*1e6)/_iters_dynamics_update << endl;
    cout << "Total integrate time " << _time_integrate << endl;
    cout << "Average integrate time (ms/1000 calls)" << (_time_integrate*1e6)/_iterations << endl;
    cout << "Total vel terms update time " << _time_update_velocity_terms << endl;
    cout << "Average vel terms update time (ms/1000 calls)" << (_time_update_velocity_terms*1e6)/_iters_velterms_update << endl;

    cout << "--- Update dynamics break down --- " << endl;
    cout << "Compute matrices time " << _time_bcm_compute_matrices << endl;
    cout << "Average compute matrices time (ms/1000 calls)" << (_time_bcm_compute_matrices*1e6)/_iters_dynamics_update << endl;
    cout << "Copy matrices time " << _time_bcm_copy_matrices << endl;
    cout << "Average copy matrices time (ms/1000 calls)" << (_time_bcm_copy_matrices*1e6)/_iters_dynamics_update << endl;
    cout << "Compute rhs vectors time " << _time_bcm_compute_vectors << endl;
    cout << "Average rhs vectors time (ms/1000 calls)" << (_time_bcm_compute_vectors*1e6)/_iters_dynamics_update << endl;
    cout << "Compute build total time " << _time_bcm_build << endl;
    cout << "Average build total time (ms/1000 calls)" << (_time_bcm_build*1e6)/_iters_dynamics_update << endl;
    cout << "Compute island construct total time " << _time_bcm_ci_total << endl;
    cout << "Average island construct total time (ms/1000 calls)" << (_time_bcm_ci_total*1e6)/_iters_dynamics_update << endl;
}

}