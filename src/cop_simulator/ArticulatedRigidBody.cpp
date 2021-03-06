// ArticulatedRigidBody.cpp

#include "ArticulatedRigidBody.h"
#include <rbdl/rbdl.h>

namespace Sai2COPSim {
	// constructor
	ArticulatedRigidBody::ArticulatedRigidBody(const std::string& name, Sai2Model::Sai2Model* model) {
		if(name.empty()) throw(std::runtime_error("Name cannot be empty"));
		if(model == NULL) throw(std::runtime_error("Model is NULL"));
		_name = name;
		_model = model;

		// set joint torques to zero initially. TODO: add option for auto gravity compensate
		jtau_act.setZero(_model->dof());

		// set contact torques to zero initially
		jtau_contact.setZero(_model->dof());

		// get initial non-linear torques from the model
		updateNonLinearJAcc();
	}

	// destructor
	ArticulatedRigidBody::~ArticulatedRigidBody() {
		// nothing to do
	}

	// add primitive. Note: Primitive origin and orientation is assumed to be the same
	// as the local frame of the rigid link to which it is connected
	// NOTE: primitive ownership is assumed to be with ArticulatedRigidBody
	// TODO: allow primitives to be directly instantiated by ArticulatedRigidBody
	void ArticulatedRigidBody::addPrimitive(const std::string& link_name, Primitive* primitive) {
		if(primitive == NULL) throw(std::runtime_error("primitive is NULL"));

		int model_linkID = _model->linkId(link_name);
		//^ if link is not present in model, an exception is thrown

		if(_link_primitives.find(link_name) == _link_primitives.end()) {
			_link_primitives[link_name] = std::unordered_map<std::string /* primitive name*/, Primitive*>();
		}

		if(_link_primitives[link_name].find(primitive->_name) != _link_primitives[link_name].end()) {
			std::cerr << primitive->_name << std::endl;
			throw(std::runtime_error("primitive already exists in the link."));
		}

		_link_primitives[link_name][primitive->_name] = primitive;
	}

	void ArticulatedRigidBody::updateNonLinearJAcc() {
		Eigen::VectorXd jtau_nonlin(_model->dof());
		_model->coriolisPlusGravity(jtau_nonlin);
		// reset the link accelerations in rbdl model to ensure they do not include
		// acceleration due to gravity and joint acceleration
		_model->updateKinematicsCustom(false /*q*/, false /*dq*/, true /*ddx*/, false /*use ddq*/);
		jacc_nonlinear = -_model->_M_inv*jtau_nonlin; //NOTE THAT THIS IS RHS ACCELERATION
	}

	Eigen::Vector3d ArticulatedRigidBody::worldToBodyPosition(const std::string& link_name, const Eigen::Vector3d& world_point) {
		Eigen::Affine3d T_link_world;
		_model->transformInWorld(T_link_world, link_name);
		// cout << T_object_base.linear() << endl;
		// cout << T_object_base.translation().transpose() << endl;
		return T_link_world.inverse() * world_point;
	}

	void ArticulatedRigidBody::JvdotTimesQdotInWorld(Eigen::Vector3d& jvdotTimesQdot,
		const std::string& link_name,
		const Eigen::Vector3d& pos_in_link)
	{
		jvdotTimesQdot = CalcPointAcceleration(*(_model->_rbdl_model),_model->_q,_model->_dq,Eigen::VectorXd::Zero(_model->_dof),_model->linkId(link_name),pos_in_link,false);
		jvdotTimesQdot = _model->_T_world_robot.linear()*jvdotTimesQdot;
	}


	void ArticulatedRigidBody::JwdotTimesQdotInWorld(Eigen::Vector3d& jwdotTimesQdot,
		const std::string& link_name)
	{
		Eigen::VectorXd a_temp = CalcPointAcceleration6D(*(_model->_rbdl_model),_model->_q,_model->_dq,Eigen::VectorXd::Zero(_model->_dof),_model->linkId(link_name),Eigen::Vector3d::Zero(),false);
		jwdotTimesQdot = _model->_T_world_robot.linear()*a_temp.head(3);
	}
}