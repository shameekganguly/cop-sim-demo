// WorldContactMap.cpp

#include "WorldContactMap.h"

namespace Sai2COPSim {

void ContactIsland::merge(ContactIsland& other) {
	assert(!other._is_merged);
	// add all ContactPrimitivePairs
	_contact_prim_pairs.splice(_contact_prim_pairs.end(), other._contact_prim_pairs);
	// add names of all articulated bodies
	for(auto name: other._articulated_bodies) {
		_articulated_bodies.insert(name);
	}
	// set merged flag
	other._is_merged = true;
}

}