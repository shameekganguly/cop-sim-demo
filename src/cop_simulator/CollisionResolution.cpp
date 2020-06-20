// CollisionResolution.cpp

#include <iostream>
#include "ContactSpaceModel.h"

using namespace Eigen;

namespace Sai2COPSim {

void ContactIslandModel::resolveCollisions(double friction_coeff, double restitution_coeff) {
	if(numContactPoints() > 2) {
		throw(std::runtime_error("Unimplemented collision resolution case."));
	}
}

}
