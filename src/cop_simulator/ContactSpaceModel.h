// ContactSpaceModel.h

#ifndef CONTACT_SPACE_MODEL_H
#define CONTACT_SPACE_MODEL_H

#include <string>
#include <list>
#include <vector>
#include <Eigen/Dense>
#include "WorldContactMap.h"
#include "ArticulatedRigidBodyManager.h"
#include "COPSolverExtended.h"
#include "lcp_solvers/LCPSolver.h"

namespace Sai2COPSim {

// algorithmic constants
namespace COPAlgorithmicConstants {
	const long COLLISION_RESOLUTION_SIMULTANEOUS_MAX_ITERATIONS = 500; //count
	const double MAX_SEPARATION_SPEED_FOR_ACTIVE_CONTACT = 1e-4; //m/s
	const double SEPARATION_SPEED_NUMERICAL_ZERO_THRESHOLD = 1e-10; //m/s, absolute value
	const double MIN_COLLISION_SPEED_FOR_STEADY_CONTACT = 1e-3; //m/s
};

// COP contact model for a single prim pair
class ContactPairState {
public:
	ContactPairState(uint id, const ContactPrimitivePair* geom_prim_pair);

	~ContactPairState();

	bool isValid() const {
		return (_last_cop_sol.result != COPSolResult::None);
	}

	double getAccelerationDueToCurvedContact(
		const Eigen::Vector3d& omegaA,
		const Eigen::Vector3d& omegaB,
		const Eigen::Vector3d& linear_contact_velocity
	) const;

public:
	// id corresponding to ContactIslandModel::_index_to_geom_island_contact_list
	uint _id;

	// number of contact points
	uint _num_contact_pts;

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

	// reference point in world frame at which cop solution is computed
	// NOTE: this is not the COP point.
	// rather, COP_world_pos = _rot_contact_frame_to_world*_last_cop_sol.local_cop_pos
	//							+ _reference_point;
	Eigen::Vector3d _reference_point;

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

	// optimization function to allow reconstruction without actually reallocating memory
	void build(const ContactIsland* geom_island, ArticulatedRigidBodyManager* arb_manager);

	~ContactIslandModel();

	// resolve collisions for this island
	// changes the values of ARB::_model::dq to non-colliding
	// assumes that the kinematic and dynamic models of the ARBs have been updated
	// already
	// returns True if there was a colliding contact that was resolved, else false
	bool resolveCollisions(double friction_coeff, double restitution_coeff);

	// resolve steady contacts for this island
	// changes the values of ARB::jtau_contact to achieve non-penetrating accelerations
	// assumes that the kinematic and dynamic models of the ARBs have been updated
	// already
	// returns True if a steady contact was resolved else false
	bool resolveSteadyContacts(double friction_coeff, double restitution_coeff);

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
	) const;

	void getActiveConstraintCOPMatrices(
		Eigen::MatrixXd& J_constraint_cop,
		Eigen::MatrixXd& Lambda_inv_constraint_cop,
		Eigen::VectorXd& rhs_constraint_cop,
		std::vector<uint>& Jrow_ind_to_contact_pair_map
	) const;

	void getActivePtContactCollisionMatrices(
		Eigen::MatrixXd& J_pt_contacts,
		Eigen::MatrixXd& Lambda_inv_pt_contacts,
		Eigen::VectorXd& rhs_pt_contacts_collision,
		std::vector<uint>& Jrow_ind_to_contact_pair_map
	) const;

	// this is an optimization required as we update the RHS vector for collisions many
	// times during a single collision resolution cycle
	// NOTE: if the active contact changes, then getActivePtContactCollisionMatrices must
	// be called once first
	void getActivePtContactCollisionRHSVector(
		Eigen::VectorXd& rhs_pt_contacts_collision
	) const;

	// functions to get just the COP RHS vectors
	void getActiveFullCOPRHSVector(
		Eigen::VectorXd& rhs_full_cop,
		std::vector<uint>& Jrow_ind_to_contact_pair_map
	) const;

	void getActiveConstraintCOPRHSVector(
		Eigen::VectorXd& rhs_constraint_cop,
		std::vector<uint>& Jrow_ind_to_contact_pair_map
	) const;

	uint numContactPoints() const {
		return _pt_contact_Jacobian.rows() / 3;
	}

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
	Eigen::MatrixXd _cop_full6_Jacobian_active;
	Eigen::MatrixXd _cop_constraint_Jacobian_active;
	Eigen::MatrixXd _pt_contact_Jacobian_active;

	// TODO: add constansts for COP constraint dof count (POINT = 3, LINE = 5, SURFACE = 6)

	// the starting row index in the _cop_constraint_Jacobian for each primitive
	std::vector<uint> _cop_constraint_Jacobian_prim_start_ind;

	// the starting row index in the _pt_contact_Jacobian for each primitive
	std::vector<uint> _pt_contact_Jacobian_prim_start_ind;

	// Contact LambdaInv for this island
	Eigen::MatrixXd _cop_full6_Lambda_inv;
	Eigen::MatrixXd _cop_constraint_Lambda_inv;
	Eigen::MatrixXd _pt_contact_Lambda_inv;
	Eigen::MatrixXd _cop_full6_Lambda_inv_active;
	Eigen::MatrixXd _cop_constraint_Lambda_inv_active;
	Eigen::MatrixXd _pt_contact_Lambda_inv_active;

	// RHS vector for collisions
	Eigen::VectorXd _pt_contact_rhs_coll;

	// RHS vector for steady contacts
	Eigen::VectorXd _cop_full6_rhs_contact;
	Eigen::VectorXd _cop_constraint_rhs_contact;

	// last collision LCP solution
	CollLCPPointSolution _last_coll_lcp_sol;

	// list of active contact indices. This corresponds to particular primitive pairs in
	// the _geom_island->_contact_prim_pairs which are either colliding or at steady contact
	std::list<uint> _active_contacts;
	// std::vector<ContactIsland::ContactList::const_iterator> _index_to_geom_island_contact_list;
	//TODO: in future, we might consider putting primitives together into one single higher
	// pair.

	// vector of contact pair states
	std::vector<ContactPairState> _pair_state;

public:
	bool _f_support_pt_contact_steady_contact;
	// - RHS vector for steady contacts with pt-contacts
	Eigen::VectorXd _pt_contact_rhs_steady_contact;

	void getActivePtContactSteadyContactRHSVector(
		Eigen::VectorXd& rhs_pt_contacts_steady_contact
	) const;

	// resolve steady contacts with multi-point contact rather than COP
	// changes the values of ARB::jtau_contact to achieve non-penetrating accelerations
	// assumes that the kinematic and dynamic models of the ARBs have been updated
	// already
	// returns True if a steady contact was resolved else false
	bool resolveSteadyContactsWithPointContacts(double friction_coeff, double restitution_coeff);

public: //internal functions
	// creates the full contact space Jacobian from the contact primitive pairs
	void createContactJacobianAndLambdaInv();

	// updates the RHS acceleration terms for each arb in this island
	void updateBodyNonlinAccelerations();

	// computes the RHS vectors based on last non-linear torques on the bodies
	void updateRHSVectors();

	// computes the RHS only for collisions based on updated _dq.
	// this is an optimization over calling updateRHSVectors, required for multiple
	// calls from resolveCollisions()
	void updatePtContactRHSCollVector();

public:
	// build time analytics
	double _time_compute_matrices;
	double _time_copy_matrices;
	double _time_compute_vectors;
	double _time_total;
};

class ContactSpaceModel {
public:
	// ctor
	ContactSpaceModel(ArticulatedRigidBodyManager* arb_manager);

	// dtor
	~ContactSpaceModel();

	// clear this model
	void clear() {
		_contact_island_models_size = 0;
		_time_compute_matrices = 0;
		_time_copy_matrices = 0;
		_time_compute_vectors = 0;
		_time_ci_total = 0;
	}

	// clear and rebuild this model
	void build(const WorldContactMap* geom_map);

	// update terms dependent on body velocities
	void updateVelocityTerms();

	// resolve collisions
	// changes the values of ARB::_model::dq to non-colliding
	// assumes that the kinematic and dynamic models of the ARBs have been updated
	// already
	// returns True if there was a colliding contact that was resolved, else false
	bool resolveCollisions(double friction_coeff, double restitution_coeff);

	// resolve steady contacts
	// changes the values of ARB::jtau_contact to achieve non-penetrating accelerations
	// assumes that the kinematic and dynamic models of the ARBs have been updated
	// already
	// returns True if a steady contact was resolved else false
	bool resolveSteadyContacts(double friction_coeff, double restitution_coeff);

public:
	// reference to sim articulated body manager
	ArticulatedRigidBodyManager* _arb_manager;

	// contact islands
	// NOTE: we preallocate space to save time
	// therefore we track an index
	std::vector<ContactIslandModel> _contact_island_models;
	uint _contact_island_models_size;

private:
	bool _f_support_pt_contact_steady_contact;
public:
	void setSupportPtContactSteadyContact(bool value);
	bool getSupportPtContactSteadyContact() const {
		return _f_support_pt_contact_steady_contact;
	}

public: // build time analytics
	double _time_compute_matrices;
	double _time_copy_matrices;
	double _time_compute_vectors;
	double _time_build;
	double _time_ci_total;
};

}

#endif //CONTACT_SPACE_MODEL_H