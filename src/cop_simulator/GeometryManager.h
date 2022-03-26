// GeometryManager.h

#ifndef GEOMETRY_MANAGER_H
#define GEOMETRY_MANAGER_H

#include <vector>

#include "geometry/PrimPrimContactInfo.h"

namespace Sai2COPSim {

class Primitive;

class GeometryManager {
public:
	// ctor
	GeometryManager() {
		// nothing to do right now
	}

	// dtor
	~GeometryManager() {
		// destroy allocated static primitives
		for(auto it = _primitives.begin(); it != _primitives.end(); it++) {
			delete *it;
			*it = NULL;
		}
		for(auto it_a = _prim_prim_distances.begin(); it_a != _prim_prim_distances.end(); it_a++) {
			for(auto it_b = (*it_a).begin(); it_b != (*it_a).end(); it_b++) {
				delete *it_b;
				*it_b = NULL;
			}
		}
	}

	// add primitive. This takes ownership as well
	// NOTE: we don't check if this primitive already exists
	void addPrimitive(Primitive* prim) {
		_primitives.push_back(prim);
		// initialize distance to other primitives
		_prim_prim_distances.push_back(std::vector<PrimPrimContactInfo*>());
		auto it_last = _prim_prim_distances.rbegin();
		for(uint i = 0; i < _primitives.size() - 1; i++) {
			(*it_last).push_back(new PrimPrimContactInfo());
		}
	}

public:
	// TODO: add sphere hierarchy

	// all geometric primitives in this simulation.
	std::vector<Primitive*> _primitives;

	// distance information between primitives
	// NOTE: this uses a vector of vector. We preallocate the distance 
	// infos to save runtime
	// NOTE: i-jth entry is distance from ith primitive to jth primitive where i < j
	// i and j denote the indices for the primitive in the _primitives vector
	std::vector<std::vector<PrimPrimContactInfo*>> _prim_prim_distances; 
};

}

#endif //GEOMETRY_MANAGER_H