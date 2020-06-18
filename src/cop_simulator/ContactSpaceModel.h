// ContactSpaceModel.h

#ifndef CONTACT_SPACE_MODEL_H
#define CONTACT_SPACE_MODEL_H

#include <string>
#include <list>
#include <vector>
#include <Eigen/Dense>
#include "WorldContactMap.h"
#include "ArticulatedRigidBodyManager.h"
#include "COPSolver.h"
#include "lcp_solvers/LCPSolver.h"

namespace Sai2COPSim {

// COP contact model for a single prim pair
class ContactPairState {
public:
	ContactPairState(uint id, const ContactPrimitivePair* geom_prim_pair);

	~ContactPairState();

	bool isValid() {
		return (_last_cop_sol.result != COPSolResult::None);
	}

public:
	// id corresponding to ContactIslandModel::_index_to_geom_island_contact_list
	uint _id;

	// position of the last computed COP, in global frame
	Eigen::Vector3d _cop_pos;

	// computed contact force
	Eigen::VectorXd _cop_force; // dimension = COP constraints dimension

	// Rotation matrix for the contact Jacobian
	// designed such that the axes have the following notations:
	// axis x & y: tangential directions. For line contacts, axis x is aligned with the constraint line
	// axis z: normal direction
	Eigen::Matrix3d _rot_contact_frame_to_world; //TODO: consider moving to prim-prim contact info

	// contact description. for introspection mostly
	// TODO: this can be a set of prim-prim pairs because of a shared COP in future
	const ContactPrimitivePair* _geom_prim_pair;

	// list of points in this contact that are active
	std::list<uint> _active_points;

public: // internal
	// last returned COP result object from the solver
	ContactCOPSolution _last_cop_sol;

public: // internal
	// add active Contact Point by id
	// returns whether the active points list changed
	bool activateContactPoint(uint id);

	// remove active Contact Point by id
	// returns whether the active points list changed
	bool deactivateContactPoint(uint id);

};

// contact model for a single ContactIsland in WorldContactMap.
class ContactIslandModel {
public:
	// default constructor does nothing
	ContactIslandModel() {
		// nothing
	}

	ContactIslandModel(const ContactIsland* geom_island, ArticulatedRigidBodyManager* arb_manager);

	~ContactIslandModel();

	// resolve collisions for this island
	// changes the values of ARB::_model::dq to non-colliding
	// assumes that the kinematic and dynamic models of the ARBs have been updated 
	// already
	void resolveCollisions();

	// resolve steady contacts for this island
	// changes the values of ARB::jtau_contact to achieve non-penetrating accelerations
	// assumes that the kinematic and dynamic models of the ARBs have been updated 
	// already
	void resolveSteadyContacts();

	// add active Contact Pair by id
	// returns whether the active contact list changed
	bool activateContactPair(uint id);

	// remove active Contact Pair by id
	// returns whether the active contact list changed
	bool deactivateContactPair(uint id);

	// get active matrices
	void getActiveFullCOPMatrices(
		Eigen::MatrixXd& J_full_cop,
		Eigen::MatrixXd& Lambda_inv_full_cop,
		Eigen::VectorXd& rhs_full_cop,
		std::vector<uint>& Jrow_ind_to_contact_pair_map
	);

	void getActiveConstraintCOPMatrices(
		Eigen::MatrixXd& J_constraint_cop,
		Eigen::MatrixXd& Lambda_inv_constraint_cop,
		Eigen::VectorXd& rhs_constraint_cop,
		std::vector<uint>& Jrow_ind_to_contact_pair_map
	);

	void getActivePtContactCollisionMatrices(
		Eigen::MatrixXd& J_pt_contacts,
		Eigen::MatrixXd& Lambda_inv_pt_contacts,
		Eigen::VectorXd& rhs_pt_contacts_collision,
		std::vector<uint>& Jrow_ind_to_contact_pair_map
	);

public:
	// Pointer to the ContactIsland for which this model holds
	const ContactIsland* _geom_island;

	// articulated body manager
	ArticulatedRigidBodyManager* _arb_manager;

	// start index for each ARB dq in the Jacobian
	std::unordered_map<std::string, uint> _arb_index_map;

	// Contact Jacobian for this island
	// Note: row groups correspond to contact pairs with indices given by 
	// _index_to_geom_island_contact_list
	// column groups correspond to particular ARBs with indices given by _arbs
	// Note: COP Jacobians are all with respect to point 0 in the contact point list.
	Eigen::MatrixXd _cop_full6_Jacobian;
	Eigen::MatrixXd _cop_constraint_Jacobian;
	Eigen::MatrixXd _pt_contact_Jacobian;
	
	// TODO: add constansts for COP constraint dof count (POINT = 3, LINE = 5, SURFACE = 6)

	// the starting row index in the _cop_constraint_Jacobian for each primitive
	std::vector<uint> _cop_constraint_Jacobian_prim_start_ind;

	// the starting row index in the _pt_contact_Jacobian for each primitive
	std::vector<uint> _pt_contact_Jacobian_prim_start_ind;

	// Contact LambdaInv for this island
	Eigen::MatrixXd _cop_full6_Lambda_inv;
	Eigen::MatrixXd _cop_constraint_Lambda_inv;
	Eigen::MatrixXd _pt_contact_Lambda_inv;

	// RHS vector for collisions
	Eigen::VectorXd _pt_contact_rhs_coll;

	// RHS vector for stead contacts
	Eigen::VectorXd _cop_full6_rhs_contact;
	Eigen::VectorXd _cop_constraint_rhs_contact;

	// last collision LCP solution
	CollLCPPointSolution _last_coll_lcp_sol;

	// list of active contact indices. This corresponds to particular primitive pairs in
	// the _geom_island->_contact_prim_pairs which are either colliding or at steady contact
	std::list<uint> _active_contacts;
	std::vector<ContactIsland::ContactList::const_iterator> _index_to_geom_island_contact_list;
	//TODO: in future, we might consider putting primitives together into one single higher
	// pair.

	// vector of contact pair states
	std::vector<ContactPairState> _pair_state;

public: //internal functions
	// creates the full contact space Jacobian from the contact primitive pairs
	void createContactJacobianAndLambdaInv();

	// computes the RHS vectors based on last non-linear torques on the bodies
	void updateRHSVectors();
};

class ContactSpaceModel {
public:
	// ctor
	ContactSpaceModel(ArticulatedRigidBodyManager* arb_manager);

	// dtor
	~ContactSpaceModel();

	// clear this model
	void clear() { _contact_island_models.clear(); }

	// clear and rebuild this model
	void build(const WorldContactMap* geom_map);

	// resolve collisions
	// changes the values of ARB::_model::dq to non-colliding
	// assumes that the kinematic and dynamic models of the ARBs have been updated 
	// already
	void resolveCollisions();

	// resolve steady contacts
	// changes the values of ARB::jtau_contact to achieve non-penetrating accelerations
	// assumes that the kinematic and dynamic models of the ARBs have been updated 
	// already
	void resolveSteadyContacts();

public:
	// reference to sim articulated body manager
	ArticulatedRigidBodyManager* _arb_manager;

	// contact islands
	std::vector<ContactIslandModel> _contact_island_models;
};

}

#endif //CONTACT_SPACE_MODEL_H