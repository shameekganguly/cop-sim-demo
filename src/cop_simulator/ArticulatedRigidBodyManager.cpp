// ArticulatedRigidBodyManager.cpp

#include "ArticulatedRigidBodyManager.h"

namespace Sai2COPSim {

ArticulatedRigidBodyManager::ArticulatedRigidBodyManager() {
	// nothing to do
}

ArticulatedRigidBodyManager::~ArticulatedRigidBodyManager() {
	// destroy managed articulated bodies
	for(auto it = _articulated_bodies.begin(); it != _articulated_bodies.end(); it++) {
		delete it->second;
		it->second = NULL;
	}
}

void ArticulatedRigidBodyManager::addBody(ArticulatedRigidBody* body) {
	assert(body != NULL);
	assert(!body->_name.empty());
	if(_articulated_bodies.find(body->_name) != _articulated_bodies.end()) {
		std::cerr << body->_name << std::endl;
		throw(std::runtime_error("Body with name already exists in manager"));
	}
	_articulated_bodies[body->_name] = body;
}

ArticulatedRigidBody* ArticulatedRigidBodyManager::getBody(const std::string& name) {
	assert(_articulated_bodies.find(name) != _articulated_bodies.end());
	return _articulated_bodies[name];
}

}