/*  01-single-contact - main.cpp
This demo shows collision and contact resolution between a single cylinder and a plane.

Author: Shameek Ganguly shameekg@stanford.edu
Date: 11/8/2017
*/
#include <iostream>
#include <string>
#include <thread>
#include <chrono>
#include <math.h>
#include <vector>

#include "Sai2Model.h"
#include "Sai2Graphics.h"
#include "Sai2Simulation.h"
#include "dynamics3d.h"
#include "rbdl/rbdl.h"
#include "chai_extension/Capsule.h"

#include <Eigen/Dense>
using namespace Eigen;

#include "timer/LoopTimer.h"

#include <GLFW/glfw3.h> //must be loaded after loading opengl/glew as part of graphicsinterface

using namespace std;

const string world_fname = "resources/01-single-contact/world.urdf";
const string object_fname = "resources/01-single-contact/sphere_object.urdf";
const string object_name = "Sphere";
const string object_link_name = "object";
const string camera_name = "camera_front";

// geometry enum
enum GeometryType {SPHERE, CAPSULE, BOX, CYLINDER};

const GeometryType geom_type = GeometryType::CAPSULE;

// sphere constants
const double radius = 0.15; //TODO: remove this hardcode.
const uint sphere_height_joint_id = 2;

// capsule constants
const double cap_radius = 0.1;
const double cap_length = 0.2;
const Vector3d end1_local(-cap_length/2, 0.0, 0.0); // Capsule is centered
const Vector3d end2_local(cap_length/2, 0.0, 0.0);

// thresholds
const double FREE_TO_COLL_DIST_THRES = 0.003; // mm
const double RESTIT_EPS = 0.3;
const double FRICTION_COEFF = 0.01;

Affine3d object_in_world;

enum ContactType {
	UNDEFINED,
	POINT,
	LINE,
	SURFACE
};

struct ContactInfo {
	double min_distance;
	vector<Vector3d> contact_points;
	ContactType type;
	ContactInfo(): min_distance(0), type(ContactType::UNDEFINED) { }
	ContactInfo(double adist, ContactType atype): min_distance(adist), type(atype) { }
};

ContactInfo getDistanceToSurface(Sai2Model::Sai2Model* object) {
	ContactInfo ret_info;
	if ( geom_type == GeometryType::SPHERE ) {
		double object_height_offset = object_in_world.translation()[2] - radius;
		auto object_joints = object->_q;
		ret_info.min_distance = object_joints[sphere_height_joint_id] + object_height_offset;
		ret_info.contact_points.push_back(Vector3d(object_joints[0], object_joints[1], ret_info.min_distance));
		ret_info.type = ContactType::POINT;
	} else if (geom_type == GeometryType::CAPSULE) {
		Vector3d end1_base, end2_base, end1_world, end2_world;
		double end1_distance, end2_distance;
		object->position(end1_base, object_link_name, end1_local);
		object->position(end2_base, object_link_name, end2_local);
		end1_world = object_in_world * end1_base;
		end2_world = object_in_world * end2_base;
		end1_distance = end1_world[2] - cap_radius; // surface normal is in +z direction
		end2_distance = end2_world[2] - cap_radius;
		if(abs(end1_distance - end2_distance) < 1e-5) {
			ret_type.type = ContactType::LINE;
			ret_info.contact_points.push_back(Vector3d(end1_world[0], end1_world[1], end1_distance));
			ret_info.contact_points.push_back(Vector3d(end2_world[0], end2_world[1], end2_distance));
		} else if(end1_distance < end2_distance) {
			ret_info.type = ContactType::POINT;
			ret_info.contact_points.push_back(Vector3d(end1_world[0], end1_world[1], end1_distance));
		} else {
			ret_info.type = ContactType::POINT;
			ret_info.contact_points.push_back(Vector3d(end2_world[0], end2_world[1], end2_distance));
		}
		ret_type.min_distance = min(end1_distance, end2_distance);
	}
	return ret_info;
}

void getContactPointLocalFrame(Sai2Model::Sai2Model* sphere, const Vector3d& global_point, Vector3d& local_point) {
	static Eigen::Affine3d T_object_base;
	sphere->transformInWorld(T_object_base, object_link_name);
	// cout << T_object_base.linear() << endl;
	// cout << T_object_base.translation().transpose() << endl;
	local_point = (object_in_world*T_object_base).inverse() * global_point;
	// cout << local_point << endl;
}

class ContactSpaceModel {
public:
	ContactInfo _contact_description; //TODO: should this be array?
	MatrixXd _contactJacobian;
	MatrixXd _contactLambdaInv;
	VectorXd _rhs_coll;
	VectorXd _rhs_contact;
	std::list<uint> _activeContacts;

	ContactSpaceModel(Sai2Model::Sai2Model* object, const VectorXd& nonlinear_torques, ContactInfo contact_description) { //TODO: extend to multi object
		uint Ndof = object->dof();
		N_contacts = contact_description.contact_points.size();
		_contactJacobian.setZero(3*N_contacts, Ndof);
		Vector3d local_contact_point;
		uint i = 0;
		for(auto point: contact_description.contact_points) {
			getContactPointLocalFrame(model, point, local_contact_point);
			// TODO: reorient local frame such that first two components are tangent and third is normal
			model->Jv(_contactJacobian.block(i*3, 0, 3, Ndof), object_link_name, local_contact_point);
			_activeContacts.push_back(i);
			i++;
		}
		_contactLambdaInv = _contactJacobian * object->_M_inv * _contactJacobian.transpose();
		_rhs_coll = _contactJacobian * model->_dq; // post_v = Lambda_inv*p_coll + pre_v
		_rhs_contact = -_contactJacobian * object->_M_inv * nonlinear_torques;
	}

	void getActiveContactSpaceMatrices(MatrixXd& JC, MatrixXd& LambdaCInv, VectorXd& rhsC_coll, VectorXd& rhsC_contact) {
		uint Ndof = _contactJacobian.cols();
		JC.resize(_activeContacts.size()*3, Ndof);
		LambdaCInv.resize(_activeContacts.size()*3, _activeContacts.size()*3);
		rhsC_coll.resize(_activeContacts.size()*3);
		rhsC_contact.resize(_activeContacts.size()*3);
		uint cind = 0;
		for(auto index: _activeContacts) {
			JC.block(3*cind, 0, 3, Ndof) = _contactJacobian(3*index, 0, 3, Ndof);
			rhsC_coll.segment(3*cind, 3) = _rhs_coll.segment(3*index, 3);
			rhsC_contact.segment(3*cind, 3) = _rhs_contact.segment(3*index, 3);
			cind++;
		}
		uint cind1 = 0;
		uint cind2 = 0;
		for(auto index1: _activeContacts) {
			for(auto index2: _activeContacts) {
				LambdaCInv.block(3*cind1, 3*cind2, 3, 3) = _contactLambdaInv.block(3*index1, 3*index2, 3, 3);
				cind2++;
			}
			cind1++;
		}
	}

	bool addActiveContact(uint contact_index) {
		//TODO: assert contact_index < N_contacts.size()
		for(auto cind: _activeContacts) {
			if(cind == contact_index) {
				return false;
			}
		}
		_activeContacts.push_back(contact_index);
		return true;
	}

	bool removeActiveContact(uint contact_index) {
		//TODO: assert contact_index < N_contacts.size()
		for(auto cind: _activeContacts) {
			if(cind == contact_index) {
				_activeContacts.remove(contact_index);
				return true;
			}
		}
		return false;
	}

private:
	uint N_contacts;
};

// simulation loop
bool fSimulationRunning = false;
void simulation(Simulation::Sai2Simulation* sim, Sai2Model::Sai2Model* model);

// initialize window manager
GLFWwindow* glfwInitialize();

// callback to print glfw errors
void glfwError(int error, const char* description);

// callback when a key is pressed
void keySelect(GLFWwindow* window, int key, int scancode, int action, int mods);

int main (int argc, char** argv) {
	cout << "Loading URDF world model file: " << world_fname << endl;

	// load simulation world
	auto sim = new Simulation::Sai2Simulation(world_fname, false);
	// set co-efficient of restition to zero to avoid bounce
    // see issue: https://github.com/manips-sai/sai2-simulation/issues/1
    sim->setCollisionRestitution(0.0);
    sim->setCoeffFrictionStatic(0.5);
    sim->setCoeffFrictionDynamic(0.5);
    object_in_world.translation() = sim->_world->getBaseNode("Sphere")->getLocalPos().eigen();
    object_in_world.linear() = sim->_world->getBaseNode("Sphere")->getLocalRot().eigen();

    Vector3d grav_vector = sim->_world->getGravity().eigen();

	// load graphics scene
	auto graphics = new Sai2Graphics::Sai2Graphics(world_fname, false);

	// set object graphics to wireframe, and show frame for last link
	graphics->showLinkFrame(true, object_name, object_link_name);

	// load object
	auto coobject = new Sai2Model::Sai2Model(object_fname, false, Eigen::Affine3d::Identity(), grav_vector);
	// coobject->_dq[1] = 0.1;
	// coobject->_dq[5] = 0.1;

	// initialize GLFW window
	GLFWwindow* window = glfwInitialize();

    // set callbacks
	glfwSetKeyCallback(window, keySelect);

	// start the simulation thread first
    fSimulationRunning = true;
	thread sim_thread(simulation, sim, coobject);

    // while window is open:
    while (!glfwWindowShouldClose(window)) {
    	// update model for object
    	coobject->updateModel();

		// update graphics. this automatically waits for the correct amount of time
		int width, height;
		glfwGetFramebufferSize(window, &width, &height);
		graphics->updateGraphics(object_name, coobject);
		graphics->render(camera_name, width, height);
		glfwSwapBuffers(window);
		glFinish();

	    // poll for events
	    glfwPollEvents();
	}

	// stop simulation
	fSimulationRunning = false;
	sim_thread.join();

    // destroy context
    glfwDestroyWindow(window);

    // terminate
    glfwTerminate();

	return 0;
}

//------------------------------------------------------------------------------
void simulation(Simulation::Sai2Simulation* sim, Sai2Model::Sai2Model* model) {
	fSimulationRunning = true;

	// create a timer
	LoopTimer timer;
	timer.initializeTimer();
	timer.setLoopFrequency(1500); //1500Hz timer

	// sleep for a few moments to let graphics start
	// std::this_thread::sleep_for(std::chrono::seconds(1));

	// enum SphereContactState {
	// 	NoContact,
	// 	Colliding,
	// 	Rolling,
	// 	Sticking
	// } state;
	// state = SphereContactState::NoContact;

	enum CapsuleContactState {
		NoContact,
		Colliding,
		Rolling,
		Sticking
	} state;
	state = CapsuleContactState::NoContact;
	
	// state vars
	const uint dof = model->dof();
	VectorXd contact_torques(dof); contact_torques.setZero();
	VectorXd nonlinear_torques(dof); nonlinear_torques.setZero();
	Vector3d contact_force; contact_force.setZero();
	double contact_moment = 0.0; // for sphere, moment is 1-dof
	double normal_impulse = 0.0;
	double normal_force;
	Vector2d tangent_impulse(0, 0);
	Vector3d contact_point; // in local link frame
	Matrix3d contact_frame = Matrix3d::Identity(); // 1st and 2nd cols are tangent, 3rd is normal
	MatrixXd contact_jacobian;
	MatrixXd contact_lambda_inv;
	ContactSpaceModel contact_model;
	VectorXd rolling_contact_impulse, sliding_contact_impulse;
	VectorXd pre_collision_contact_vel, post_collision_contact_vel;
	VectorXd rhs_coll, rhs_contact;
	Vector2d slip_direction;

	ContactInfo cinfo;

	double last_time = timer.elapsedTime(); //secs
	bool fTimerDidSleep = true;

	while (fSimulationRunning) {
		fTimerDidSleep = timer.waitForNextLoop();
		// if (timer.elapsedCycles() % 10000 == 0) {
		// 	cout << "Simulation loop frequency: " << timer.elapsedCycles()/timer.elapsedTime() << endl; 
		// }

		// integrate forward
		double curr_time = timer.elapsedTime();
		double loop_dt = curr_time - last_time; 
		// sim->integrate(loop_dt);
		
		// check contact state every now and then
		if (timer.elapsedCycles() % 10 == 0) {
			model->updateModel(); // TODO: we shouldn't explicity compute M inverse
			model->coriolisPlusGravity(nonlinear_torques);

			// check if contacts have changed:
			// TODO: implement a change tracker so that we can do delta updates
			// to the contact state model
			cinfo = getDistanceToSurface(model);
			if (cinfo.min_distance < FREE_TO_COLL_DIST_THRES) {
				// compute contact velocities to see if any contact is colliding
				contact_model = ContactSpaceModel(model, nonlinear_torques, cinfo);
				pre_collision_contact_vel = contact_model._contactJacobian * model->_dq;

				bool is_colliding = false;
				for(uint i = 0; i < cinfo.contact_points.size(); i++) {
					is_colliding = is_colliding && (pre_collision_contact_vel[3*i + 2] < 0);
					if (pre_collision_contact_vel[3*i + 2] > 0) {
						contact_model.removeActiveContact(i);
					}
				}
				if (is_colliding) {
					state = CapsuleContactState::Colliding;
					cout << "Collision" << endl;
					// initialize state for CapsuleContactState::Colliding
					contact_model.getActiveContactSpaceMatrices(contact_jacobian, contact_lambda_inv, rhs_coll, rhs_contact);
				}
			}

			switch (state) {
				case CapsuleContactState::NoContact:
					break;
				case CapsuleContactState::Colliding:
					break;
				case CapsuleContactState::Rolling:

				break;
				case CapsuleContactState::Sticking:
				break;
			}
			// cout << model->_q.transpose() << endl;
		}
		// TODO: else, extrapolated model updates

		// resolve collision
		contact_torques.setZero();
		if (state == CapsuleContactState::Colliding) {
			bool is_post_steady_contact = false;
			while(true) {
				pre_collision_contact_vel = contact_jacobian * model->_dq;
				double min_coll_vel = 0.0;
				for (uint i = 0; i < contact_model._activeContacts.size(); i++) {
					min_coll_vel = min(min_coll_vel, pre_collision_contact_vel[3*i + 2]);
				}
				double coll_restitution = RESTIT_EPS;
				if (abs(min_coll_vel) < 1e-3) { // We use the fastest collision to determine the coefficient of restitution
					coll_restitution = 0.0; // force inelastic collision to bring to steady contact
				}
				CollLCPPointSolution lcp_sol;
				if(contact_model._activeContacts.size() == 1) {
					lcp_sol = solveCollLCPOnePoint (
						contact_lambda_inv,
						pre_collision_contact_vel,
						pre_collision_contact_vel,
						coll_restitution,
						FRICTION_COEFF
					);
				} else {
					lcp_sol = solveCollLCPPoint (
						2, //TODO: remove hardcode
						contact_lambda_inv,
						pre_collision_contact_vel,
						pre_collision_contact_vel,
						coll_restitution,
						FRICTION_COEFF
					);
				}
				if (lcp_sol.result == LCPSolResult::Success) {
					cout << "LCP Impulse: " << lcp_sol.p_sol;
				} else {
					cerr << "LCP failed with type: " << lcp_sol.result << endl;
					cout << pre_collision_contact_vel.transpose() << endl;
					break;
				}
				// update joint velocities
				model->_dq += model->_M_inv*(contact_jacobian.transpose()*lcp_sol.p_sol);
				// recheck for collision at ALL contacts:
				post_collision_contact_vel = contact_model._contactJacobian * model->_dq;
				bool did_contact_change = false;
				bool is_post_colliding = false;
				for (uint i = 0; i < cinfo.contact_points.size(); i++) {
					if (post_collision_contact_vel[3*i + 2] <= 0) {
						did_contact_change = did_contact_change || contact_model->addActiveContact(i);
					} else {
						did_contact_change = did_contact_change || contact_model->removeActiveContact(i);
					}
					if(post_collision_contact_vel[3*i + 2] < 0) {
						is_post_colliding = true;
					} else if (post_collision_contact_vel[3*i + 2] <= 1e-4) {
						is_post_steady_contact = true;
					}
				}
				if(did_contact_change) {
					contact_model.getActiveContactSpaceMatrices(contact_jacobian, contact_lambda_inv, rhs_coll, rhs_contact);
				}
				if(!is_post_colliding) break;
			} // end sequential collisions

			// check if steady contact is formed
			if (is_post_steady_contact) {
				state = CapsuleContactState::Rolling;
				cout << "Rolling" << endl;
			} else {
				state = CapsuleContactState::NoContact;
				cout << "No contact" << endl;
			}
		}

		// compute contact torques
		if (state == CapsuleContactState::Rolling) {
			break;

			// recompute active contacts: TODO: this might be optimized when a collision was resolved in this sim cycle
			bool did_contact_change = false;
			post_collision_contact_vel = contact_model._contactJacobian * model->_dq;
			for (uint i = 0; i < cinfo.contact_points.size(); i++) {
				if (post_collision_contact_vel[3*i + 2] <= 1e-4) { // note that we allow some positive contact velocity
					did_contact_change = did_contact_change || contact_model->addActiveContact(i);
				} else {
					did_contact_change = did_contact_change || contact_model->removeActiveContact(i);
				}
			}
			if(did_contact_change) {
				contact_model.getActiveContactSpaceMatrices(contact_jacobian, contact_lambda_inv, rhs_coll, rhs_contact);
			}
			post_collision_contact_vel = contact_jacobian * model->_dq;

			// TODO: switch to COP model
			CollLCPPointSolution lcp_sol;
			if(contact_model._activeContacts.size() == 1) {
				lcp_sol = solveCollLCPOnePoint (
					contact_lambda_inv,
					rhs_contact,
					post_collision_contact_vel,
					0.0,
					FRICTION_COEFF
				);
			} else {
				lcp_sol = solveCollLCPPoint (
					2, //TODO: remove hardcode
					contact_lambda_inv,
					rhs_contact,
					post_collision_contact_vel,
					0.0,
					FRICTION_COEFF
				);
			}
			if (lcp_sol.result == LCPSolResult::Success) {
				cout << "Contact LCP Contact force: " << lcp_sol.p_sol;
			} else {
				cerr << "Contact LCP failed with type: " << lcp_sol.result << endl;
				break;
			}
			contact_force = lcp_sol.p_sol;
			contact_torques = contact_jacobian.transpose() * contact_force;
		}

		// integrate //TODO: switch to Euler-Heun
		model->_ddq = model->_M_inv * (-nonlinear_torques + contact_torques);
		model->_dq += model->_ddq * loop_dt;
		model->_q.segment<3>(0) += model->_dq.segment<3>(0) * loop_dt; //TODO: specifically isolate spherical joints
		if (model->_dq.segment<3>(3).norm() != 0) {
			model->_q.segment<4>(3) = RigidBodyDynamics::Math::Quaternion(model->_q.segment<4>(3)).timeStep(model->_dq.segment<3>(3), loop_dt);
		}

		// if (!fTimerDidSleep) {
		// 	cout << "Warning: timer underflow! dt: " << loop_dt << "\n";
		// }

		// update last time
		last_time = curr_time;
	}
}

//------------------------------------------------------------------------------
GLFWwindow* glfwInitialize() {
		/*------- Set up visualization -------*/
    // set up error callback
    glfwSetErrorCallback(glfwError);

    // initialize GLFW
    glfwInit();

    // retrieve resolution of computer display and position window accordingly
    GLFWmonitor* primary = glfwGetPrimaryMonitor();
    const GLFWvidmode* mode = glfwGetVideoMode(primary);

    // information about computer screen and GLUT display window
	int screenW = mode->width;
    int screenH = mode->height;
    int windowW = 0.8 * screenH;
    int windowH = 0.5 * screenH;
    int windowPosY = (screenH - windowH) / 2;
    int windowPosX = windowPosY;

    // create window and make it current
    glfwWindowHint(GLFW_VISIBLE, 0);
    GLFWwindow* window = glfwCreateWindow(windowW, windowH, "SAI2.0 - CS327a HW4", NULL, NULL);
	glfwSetWindowPos(window, windowPosX, windowPosY);
	glfwShowWindow(window);
    glfwMakeContextCurrent(window);
	glfwSwapInterval(1);

	return window;
}

//------------------------------------------------------------------------------

void glfwError(int error, const char* description) {
	cerr << "GLFW Error: " << description << endl;
	exit(1);
}

//------------------------------------------------------------------------------

void keySelect(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    // option ESC: exit
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    {
        // exit application
         glfwSetWindowShouldClose(window, 1);
    }
}
