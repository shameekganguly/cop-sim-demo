// CopSimulator.h
#ifndef COPSIMULATOR_H
#define COPSIMULATOR_H

#include <string>
#include <unordered_map>
#include <Eigen/Dense>

#include "geometry/Primitive.h"
#include "ArticulatedRigidBodyManager.h"
#include "WorldContactMap.h"
#include "GeometryManager.h"
#include "ContactSpaceModel.h"

namespace Sai2COPSim {

// Constants used
namespace COPAlgorithmicConstants {
	const double GEOMETRIC_CONTACT_DISTANCE_THRESHOLD = 0.003; //m
	const long NUM_ITERS_BEFORE_COLLISION_CHECK_UPDATE = 20; //count
	const long NUM_ITERS_BEFORE_INERTIA_AND_CONTACT_MODEL_UPDATE = 10; //count
	const long NUM_ITERS_BEFORE_NONLINEAR_ACCERLATION_UPDATE = 4; //count
};

class COPSimulator {
public:

// constructor
COPSimulator(double friction_coeff, double restitution_coeff);

// destructor
~COPSimulator();

// add an inifinite plane object
void addPlane(const std::string& name, const Eigen::Vector3d& planeNormal, const Eigen::Vector3d& planePoint);

// add a capsule object
// object->_T_world_robot is assumed to have been set already
// capsule is assumed to be aligned with X axis. local frame at the center of the capsule.
void addCapsuleObject(const std::string& articulated_body_name,
		const std::string& link_name, // name for the link on which capsule primitive will be attached
		const std::string& primitive_name,
		Sai2Model::Sai2Model* object,
		double radius,
		double length);

// add a cylinder object
// object->_T_world_robot is assumed to have been set already
// cylinder is assumed to be aligned with +Z axis. local frame at the center of one end
// face of the cylinder.
void addCylinderObject(const std::string& articulated_body_name,
		const std::string& link_name, // name for the link on which capsule primitive will be attached
		const std::string& primitive_name,
		Sai2Model::Sai2Model* object,
		double radius,
		double length,
		uint num_points);

// automatically sets q and dq for the objects
void integrate(double dt);

// TODO: function to add/delete ARB

public:
	// Coulomb friction coefficient
	// TODO: this should be property of a particular primitive once we have blending
	// of material properties in the contact solver
	double _friction_coeff;

	// Restitution coefficient
	// TODO: this should be property of a particular primitive once we have blending
	// of material properties in the contact solver
	// TODO: currently, this is a Newtonian coefficient. need to switch to Poisson later
	double _resitution_coeff;

	// map of named fixed primitives
	std::unordered_map<std::string, Primitive*> _static_objects;

	// repository of articulated bodies in this sim
	ArticulatedRigidBodyManager _arb_manager;

public: //internal functions of the COPSimulator
	void computeWorldContactMap();

public: //internal members of the COPSimulator
	// contact map of the world
	WorldContactMap _contact_map;

	// geometry manager for all articulated bodies in the world
	GeometryManager _geom_manager;

	// contact space model
	// this has to be a pointer because it cannot be initialized by default
	ContactSpaceModel* _contact_model;

	// number of sim iterations since instantiation
	unsigned long long _iterations; 
	//TODO: consider using an internal timer for adaptive stepping and updating 
	// collision meshes at a slower rate.
	// TODO: Later we need to have a 5 rate simulation with the following:
	// Rate 1: Integrator sub-stepping (for higher order methods)
	// Rate 2: Integrator gross dt, collision resolution, contact forces
	// Rate 3: Update model kinematics, and non-linear torques (update kinematics is required for calculating dot_J*dq, which is required for the contact force evaluation)
	// Rate 4: Update model dynamics: ABA algo, calculate M_inv, Contact model
	// Rate 5: Compute contact geometry

	// Rate 1 is the fastest, Rate 5 is the slowest
};

}

#endif //COPSIMULATOR_H