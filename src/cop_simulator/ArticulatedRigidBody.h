// ArticulatedRigidBody.h

#ifndef ARTICULATED_RIGID_BODY_H
#define ARTICULATED_RIGID_BODY_H

#include <string>
#include <unordered_map>

#include <Eigen/Dense>
#include "Sai2Model.h"

#include "geometry/Primitive.h"

namespace Sai2COPSim {

// TODO: think about merging this with SAI2Model directly, so that it can be
// used both for simulation and for control
class ArticulatedRigidBody {
public:
	// constructor
	// NOTE: Model ownership is assumed to be with client application
	ArticulatedRigidBody(const std::string& name, Sai2Model::Sai2Model* model);

	// destructor
	~ArticulatedRigidBody();

	// add primitive. Note: Primitive origin and orientation is assumed to be the same
	// as the local frame of the rigid link to which it is connected
	// NOTE: primitive ownership is with COP GeometryManager. Therefore, it is possible
	// that the returned pointer is NULL
	// TODO: allow primitives to be directly instantiated by ArticulatedRigidBody
	// TODO: allow arbitrary local transform to primitive frame. we can then set
	// primitive->_transform_in_link
	void addPrimitive(const std::string& link_name, Primitive* primitive);

	// TODO: add access functions for q, dq, jtau

	// joint motor torques. Set by the controller.
	// TODO: extend to arbitrary actuation model such as muscles 
	// or some actuator group with a particular J_actuator^T
	Eigen::VectorXd jtau_act;

	// update joint acceleration from non-linear torques
	void updateNonLinearJAcc();

	// get body position from world position
	// utility function
	Eigen::Vector3d worldToBodyPosition(const std::string& link_name, const Eigen::Vector3d& world_point);

	// -- Jdot*dq functions --
	/**
	* @brief      Product of time derivative of Jv of a link, and dq. Output is with
	*              respect to world frame
	* @note       We dont ask RBDL to update kinematics. This assumes no dynamics call
	*              has been made between the last updateKinematics() call and this call
	*              dynamics calls include coriolisForce and coriolisForcePlusGravity in
	*              Sai2Model, as well as any of the functions in RBDL/Dynamics.h
	*				It is ok to call the updateNonLinearJAcc() function however, as it
	*				resets the link accelerations to not include gravity and ddq terms.
	*
	* @param      jvdotTimesQdot  Vector to which the result is written
	* @param      link_name        name of the link for which to compute the angular
	*                              acceleration
	* @param      pos_in_link      the position of the point in the link, in local
	*                              link frame
	*/
	void JvdotTimesQdotInWorld(Eigen::Vector3d& jvdotTimesQdot,
								const std::string& link_name,
								const Eigen::Vector3d& pos_in_link = Eigen::Vector3d::Zero());
	/* @param     jwdotTimesQdot     Vector to which the result is written
	* @param      link_name  name of the link for which to compute the angular
	*                        acceleration
	*/
	void JwdotTimesQdotInWorld(Eigen::Vector3d& jwdotTimesQdot,
								const std::string& link_name);


public:
	// string name
	std::string _name;

	// RBDL model
	Sai2Model::Sai2Model* _model;

	// geometry primitives. NOTE: this membership is really for introspection and logical membership
	// all primitives in the world are stored in a GeometryManager repository under the
	// COPSimulator in order to optimize contact detection and updates
	// NOTE: primitive ownership is with COP GeometryManager. Therefore, it is possible
	// that the returned pointer is NULL
	std::unordered_map< std::string /* link name*/,
				std::unordered_map<std::string /* primitive name*/, Primitive*>> _link_primitives;
	// TODO: add boundary relationship for primitives. This is required for computing the contact 
	// patch in surface-surface contacts

	// contact torques in joint space. From the contact space model.
	// TODO: separate out to ARBSim class in future, when we merge this to Sai2Model
	Eigen::VectorXd jtau_contact; // RHS vector of contact torques.

	// cached acceleration from non-linear torques in joint space.
	// TODO: separate out to ARBSim class in future, when we merge this to Sai2Model
	Eigen::VectorXd jacc_nonlinear; // RHS vector of -A^{-1}(g+b) accelerations.
};

}

#endif //ARTICULATED_RIGID_BODY_H