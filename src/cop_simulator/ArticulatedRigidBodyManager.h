// ArticulatedRigidBodyManager.h

#ifndef ARTICULATED_RIGID_BODY_MANAGER_H
#define ARTICULATED_RIGID_BODY_MANAGER_H

#include <string>
#include <unordered_map>

#include "ArticulatedRigidBody.h"

namespace Sai2COPSim {

// shared repository of articulated bodies in this simulation between different components
// owns and manages all the articulated bodies
class ArticulatedRigidBodyManager {
public:
	// ctor
	ArticulatedRigidBodyManager();

	// dtor
	~ArticulatedRigidBodyManager();

	// add ARB. this also takes ownership
	void addBody(ArticulatedRigidBody* body);

	// TODO: remove body

	// get body non-const interface
	ArticulatedRigidBody* getBody(const std::string& name);

public:
	// all articulated bodies in this simulator.
	// Assumed to be single rigid bodies for now in the contact model.
	std::unordered_map<std::string, ArticulatedRigidBody*> _articulated_bodies;
};

}

#endif //ARTICULATED_RIGID_BODY_MANAGER_H