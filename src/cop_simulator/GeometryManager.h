// GeometryManager.h

#ifndef GEOMETRY_MANAGER_H
#define GEOMETRY_MANAGER_H

#include <vector>

#include "geometry/Primitive.h"

namespace Sai2COPSim {

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
	}

public:
	// TODO: add sphere hierarchy

	// all geometric primitives in this simulation.
	std::vector<Primitive*> _primitives;
};

}

#endif //GEOMETRY_MANAGER_H