/*  01-single-contact - main.cpp
This demo shows collision and contact resolution between a single capsule and a plane.

Author: Shameek Ganguly shameekg@stanford.edu
Date: 11/8/2017
*/

// TODO: Consider deleting this example since it is a subset of 02-two-capsules.

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
#include <chai3d.h>
#include "chai_extension/Capsule.h"

#include <Eigen/Dense>
using namespace Eigen;

#include "timer/LoopTimer.h"
#include "force_sensor/ForceSensorSim.h"
#include "force_sensor/ForceSensorDisplay.h"

#include "lcp_solvers/LCPSolverInternal.h"
#include "cop_simulator/COPSolver.h"

#include <GLFW/glfw3.h> //must be loaded after loading opengl/glew as part of graphicsinterface

using namespace std;
using namespace chai3d;
using namespace Sai2LCPSolver;

const string world_fname = "resources/01-single-contact/world.urdf";
const string object_fname = "resources/01-single-contact/sphere_object.urdf";
const string object_name = "Sphere";
const string object_link_name = "object";
string camera_name = "camera_front";

// flags for scene camera movement
static bool fTransXp = false;
static bool fTransXn = false;
static bool fTransYp = false;
static bool fTransYn = false;
static bool fTransZp = false;
static bool fTransZn = false;

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
const double RESTIT_EPS = 0.1;
const double FRICTION_COEFF = 0.05;

// logging
bool LOG_DEBUG = false;

// which sim to use
bool USE_ORB_SIM = true;

// pause
bool f_pause_sim = false;

Affine3d object_in_world;

chai3d::cShapeSphere pt_contact_display(0.02);
chai3d::cShapeCylinder line_contact_display(0.01, 0.01, 0.2);
ForceSensorDisplay* cop_force_display;
ForceSensorSim* cop_force_sensor;

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
		if(end2_base.segment(0,2).norm() < 1e-10) {
			throw(std::out_of_range("Something wrong with end2base"));
		}
		end1_world = object_in_world * end1_base;
		end2_world = object_in_world * end2_base;
		if(end2_world.segment(0,2).norm() < 1e-10) {
			throw(std::out_of_range("Something wrong with end2_world"));
		}
		end1_distance = end1_world[2] - cap_radius; // surface normal is in +z direction
		end2_distance = end2_world[2] - cap_radius;
		if(abs(end1_distance - end2_distance) < FREE_TO_COLL_DIST_THRES*0.02) {
			ret_info.type = ContactType::LINE;
			// if(LOG_DEBUG) cout << "Line contact " << endl;
			ret_info.contact_points.push_back(Vector3d(end1_world[0], end1_world[1], end1_distance));
			ret_info.contact_points.push_back(Vector3d(end2_world[0], end2_world[1], end2_distance));
		} else if(end1_distance < end2_distance) {
			ret_info.type = ContactType::POINT;
			// if(LOG_DEBUG) cout << "End 2 contact " << endl;
			ret_info.contact_points.push_back(Vector3d(end1_world[0], end1_world[1], end1_distance));
		} else {
			ret_info.type = ContactType::POINT;
			// if(LOG_DEBUG) cout << "End 1 contact " << endl;
			ret_info.contact_points.push_back(Vector3d(end2_world[0], end2_world[1], end2_distance));
		}
		ret_info.min_distance = min(end1_distance, end2_distance);
	}
	return ret_info;
}

void getContactPointLocalFrame(Sai2Model::Sai2Model* sphere, const Vector3d& global_point, Vector3d& local_point) {
	static Eigen::Affine3d T_object_base;
	sphere->transform(T_object_base, object_link_name);
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
	Matrix3d _rotXInterpoint;
	std::list<uint> _activeContacts;
	// for COP model:
	std::vector<Vector3d> _contactPointsLocal;
	Matrix3d _R_body_COP_frame;
	MatrixXd _point0J6;
	MatrixXd _point0_lambdaInv_6;
	VectorXd _point0_rhs_nonlin_6;

	ContactSpaceModel() {
		// nothing to do
	}

	ContactSpaceModel(Sai2Model::Sai2Model* object, const VectorXd& nonlinear_torques, ContactInfo contact_description) { //TODO: extend to multi object
		uint Ndof = object->dof();
		_contact_description = contact_description;
		N_contacts = contact_description.contact_points.size();
		_contactJacobian.setZero(3*N_contacts, Ndof);
		Vector3d local_contact_point;
		uint i = 0;
		_rotXInterpoint = Matrix3d::Identity();
		if(N_contacts == 2) { //TODO: extend to more than 2 points
			Vector3d t = contact_description.contact_points[1] - contact_description.contact_points[0];
			double a = t[0]/sqrt(t[1]*t[1] + t[0]*t[0]);
			double b = t[1]/sqrt(t[1]*t[1] + t[0]*t[0]);
			_rotXInterpoint << a, b, 0,
							-b, a, 0,
							0, 0, 1;
		}
		for(auto point: contact_description.contact_points) {
			getContactPointLocalFrame(object, point, local_contact_point);
			_contactPointsLocal.push_back(local_contact_point);
			if(LOG_DEBUG) cout << "Gobal contact point: " << point.transpose() << endl;
			if(LOG_DEBUG) cout << "Local contact point: " << local_contact_point.transpose() << endl;
			MatrixXd blockJv;
			object->Jv(blockJv, object_link_name, local_contact_point);
			blockJv = _rotXInterpoint*object_in_world.linear()*blockJv;
			// TODO: reorient local frame such that first two components are tangent and third is normal
			_contactJacobian.block(i*3, 0, 3, Ndof) = blockJv;
			_activeContacts.push_back(i);
			i++;
		}
		_contactLambdaInv = _contactJacobian * object->_M_inv * _contactJacobian.transpose();
		_rhs_coll = _contactJacobian * object->_dq; // post_v = Lambda_inv*p_coll + pre_v
		recomputeContactRHS(object, nonlinear_torques);

		// compute cop matrices
		if(N_contacts == 2) {
			_point0J6.setZero(6, Ndof);
			object->J_0(_point0J6, object_link_name, _contactPointsLocal[0]);
			_point0J6.block(0, 0, 3, Ndof) = _rotXInterpoint*object_in_world.linear()*_point0J6.block(0, 0, 3, Ndof);
			_point0J6.block(3, 0, 3, Ndof) = _rotXInterpoint*object_in_world.linear()*_point0J6.block(3, 0, 3, Ndof);
			Matrix3d ori_body_base;
			object->rotation(ori_body_base, object_link_name);
			_R_body_COP_frame = _rotXInterpoint*object_in_world.linear()*ori_body_base;
			_point0_lambdaInv_6 = _point0J6 * object->_M_inv * _point0J6.transpose();
			_point0_rhs_nonlin_6 = -_point0J6 * object->_M_inv * nonlinear_torques;
		}
	}

	void recomputeContactRHS(Sai2Model::Sai2Model* object, const VectorXd& nonlinear_torques) {
		_rhs_contact = -_contactJacobian * object->_M_inv * nonlinear_torques;
		// Add dotJ*dotq
		MatrixXd objectJw;
		object->Jw(objectJw, object_link_name);
		Vector3d omega, linvel;
		omega = object_in_world.linear()* objectJw * object->_dq;
		int i = 0;
		for(auto point: _contact_description.contact_points) {
			Vector3d rworld;
			object->position(rworld, object_link_name, Vector3d::Zero());
			rworld = point - object_in_world * rworld;
			_rhs_contact.segment(3*i, 3) += _rotXInterpoint*(omega.cross(omega.cross(rworld)));
			i++;
		}
		// cop nonlin update
		if(_contact_description.contact_points.size() == 2) {
			_point0_rhs_nonlin_6 = -_point0J6 * object->_M_inv * nonlinear_torques;
		}
	}

	void getActiveContactSpaceMatrices(MatrixXd& JC, MatrixXd& LambdaCInv, VectorXd& rhsC_coll, VectorXd& rhsC_contact) {
		uint Ndof = _contactJacobian.cols();
		JC.resize(_activeContacts.size()*3, Ndof);
		LambdaCInv.resize(_activeContacts.size()*3, _activeContacts.size()*3);
		rhsC_coll.resize(_activeContacts.size()*3);
		rhsC_contact.resize(_activeContacts.size()*3);
		uint cind = 0;
		for(auto index: _activeContacts) {
			JC.block(3*cind, 0, 3, Ndof) = _contactJacobian.block(3*index, 0, 3, Ndof);
			rhsC_coll.segment(3*cind, 3) = _rhs_coll.segment(3*index, 3);
			rhsC_contact.segment(3*cind, 3) = _rhs_contact.segment(3*index, 3);
			cind++;
		}
		uint cind1 = 0;
		for(auto index1: _activeContacts) {
			uint cind2 = 0;
			for(auto index2: _activeContacts) {
				LambdaCInv.block(3*cind1, 3*cind2, 3, 3) = _contactLambdaInv.block(3*index1, 3*index2, 3, 3);
				cind2++;
			}
			cind1++;
		}
	}

	bool addActiveContact(uint contact_index) {
		if(contact_index >= _contact_description.contact_points.size()) {
			throw(std::out_of_range("Added contact index is out of range!"));
		} //TODO: move to assert
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

	Vector3d getOmegaCOPFrame(const VectorXd& dq) {
		return _point0J6.block(3, 0, 3, dq.size()) * dq;
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
	sim->setCollisionRestitution(RESTIT_EPS);
    sim->setCoeffFrictionStatic(FRICTION_COEFF);
    sim->setCoeffFrictionDynamic(FRICTION_COEFF);
    object_in_world.translation() = sim->_world->getBaseNode(object_name)->getLocalPos().eigen();
    object_in_world.linear() = sim->_world->getBaseNode(object_name)->getLocalRot().eigen();
	// VectorXd temp(6);
 //    sim->getJointVelocities(object_name, temp);
 //    cout << temp.transpose() << endl;

    Vector3d grav_vector = sim->_world->getGravity().eigen();

	// load graphics scene
	auto graphics = new Sai2Graphics::Sai2Graphics(world_fname, false);
	auto graphics_capsule_chailink = graphics->findLink(object_name, object_link_name);
	auto cap_mmesh = dynamic_cast<cMultiMesh*>(graphics_capsule_chailink->getChild(0));
	auto cap_mesh = dynamic_cast<cMesh*>(cap_mmesh->getMesh(0));

    // set graphics object as collision for sai2-simulation
    {
    	auto objectsimbase = sim->_world->getBaseNode(object_name);
	    objectsimbase->enableDynamics(false);
	    auto link = objectsimbase->getLink(object_link_name);
	    auto tmp_mmesh = new cMultiMesh();
		tmp_mmesh->m_name = std::string("sai_dyn3d_link_mesh");
		auto tmp_mesh = cap_mesh->copy();
		tmp_mmesh->addMesh(tmp_mesh);
		tmp_mmesh->setLocalPos(cVector3d(-0.1, 0.0, 0.0));
		tmp_mmesh->setLocalRot(Matrix3d::Identity());
		link->setCollisionModel(tmp_mmesh);
		link->buildCollisionHull(0.0001, 0.0001);
		objectsimbase->enableDynamics(true);
    }

	// set object graphics to wireframe, and show frame for last link
	// graphics->showLinkFrame(true, object_name, object_link_name);
	graphics->showWireMeshRender(true, object_name, object_link_name);

	// load object
	auto coobject = new Sai2Model::Sai2Model(object_fname, false, object_in_world, object_in_world.linear().transpose()*grav_vector);
	// cout << coobject->_q.transpose() << endl;
	// coobject->_dq[0] = 0.3;
	// coobject->_dq[1] = 0.1;
	// coobject->_dq[2] = 1.0;
	coobject->_dq[3] = 0.7;
	coobject->_dq[5] = -0.5;

	// force sim/ display
	cop_force_sensor = new ForceSensorSim(object_name, object_link_name, Eigen::Affine3d::Identity(), coobject);
	cop_force_display = new ForceSensorDisplay(cop_force_sensor, graphics);
	cop_force_display->_force_line_scale = 10;
	cop_force_display->_moment_line_scale = 200;
	cop_force_display->_display_line_moment->setLineWidth(8.0);

	// initialize GLFW window
	GLFWwindow* window = glfwInitialize();

    // set callbacks
	glfwSetKeyCallback(window, keySelect);

	// start the simulation thread first
    fSimulationRunning = true;
	thread sim_thread(simulation, sim, coobject);

	// create sphere to see contact point
	pt_contact_display.setShowEnabled(false);
	pt_contact_display.m_material->setBrownMaroon();
	pt_contact_display.m_material->setShininess(100);
	graphics->_world->addChild(&pt_contact_display);
	pt_contact_display.setLocalPos(Vector3d(0.0, 0.0, 0.4));

	line_contact_display.setShowEnabled(false);
	line_contact_display.m_material->setBrownMaroon();
	line_contact_display.m_material->setShininess(100);
	graphics->_world->addChild(&line_contact_display);
	line_contact_display.setLocalPos(Vector3d(0.0, 0.0, 0.4));

    // while window is open:
    Eigen::Vector3d camera_pos, camera_lookat, camera_vertical;
	while (!glfwWindowShouldClose(window)) {
    	// update model for object
    	// coobject->updateModel();

		// update graphics. this automatically waits for the correct amount of time
		int width, height;
		glfwGetFramebufferSize(window, &width, &height);
		graphics->updateGraphics(object_name, coobject);
		cop_force_display->update();
		graphics->render(camera_name, width, height);
		// compute global position of spline and cherry
		// pt_contact_display.computeGlobalPositions();
		glfwSwapBuffers(window);
		glFinish();

	    // poll for events
	    glfwPollEvents();

	    // move scene camera as required
    	Eigen::Vector3d cam_up_axis;
    	cam_up_axis << 0.0, 0.0, 1.0; //TODO: there might be a better way to do this
    	graphics->getCameraPose(camera_name, camera_pos, camera_vertical, camera_lookat);
	    Eigen::Vector3d cam_roll_axis = (camera_lookat - camera_pos).cross(cam_up_axis);
	    Eigen::Vector3d cam_lookat_axis = camera_lookat - camera_pos;
    	cam_roll_axis.normalize();
    	cam_lookat_axis.normalize();
    	if (fTransXp) {
	    	camera_pos = camera_pos + 0.05*cam_roll_axis;
	    	camera_lookat = camera_lookat + 0.05*cam_roll_axis;
	    }
	    if (fTransXn) {
	    	camera_pos = camera_pos - 0.05*cam_roll_axis;
	    	camera_lookat = camera_lookat - 0.05*cam_roll_axis;
	    }
	    if (fTransYp) {
	    	camera_pos = camera_pos + 0.05*cam_up_axis;
	    	camera_lookat = camera_lookat + 0.05*cam_up_axis;
	    }
	    if (fTransYn) {
	    	camera_pos = camera_pos - 0.05*cam_up_axis;
	    	camera_lookat = camera_lookat - 0.05*cam_up_axis;
	    }
	    if (fTransZp) {
	    	camera_pos = camera_pos + 0.05*cam_lookat_axis;
	    	camera_lookat = camera_lookat + 0.05*cam_lookat_axis;
	    }
	    if (fTransZn) {
	    	camera_pos = camera_pos - 0.05*cam_lookat_axis;
	    	camera_lookat = camera_lookat - 0.05*cam_lookat_axis;
	    }
	    graphics->setCameraPose(camera_name, camera_pos, cam_up_axis, camera_lookat);
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
	Matrix3d bodyorient = Matrix3d::Identity();

	ContactInfo cinfo;

	// COP state
	ContactCOPSolution last_cop_sol;
	MatrixXd last_cop_J6(6, dof);
	MatrixXd last_cop_lambdainv_full(6, 6);
	VectorXd last_cop_rhsnonlin_full(6);

	// create a timer
	LoopTimer timer;
	timer.initializeTimer();
	timer.setLoopFrequency(20000); //1500Hz timer
	double last_time = timer.elapsedTime(); //secs
	bool fTimerDidSleep = true;
	model->updateModel();
	// model->coriolisPlusGravity(nonlinear_torques);
	// cout << model->_M_inv << endl;
	// cout << nonlinear_torques.transpose() << endl;
	// cout << model->_world_gravity.transpose() << endl;
	sim->_world->getBaseNode(object_name)->getJoint("jrs")->setVelSpherical(Vector3d(0.7,0.0,0.3));
	while (fSimulationRunning) {
		fTimerDidSleep = timer.waitForNextLoop();
		// if (timer.elapsedCycles() % 10000 == 0) {
		// 	cout << "Simulation loop frequency: " << timer.elapsedCycles()/timer.elapsedTime() << endl; 
		// }

		// integrate forward
		double curr_time = timer.elapsedTime();
		if(f_pause_sim) {
			last_time = curr_time;
			continue;
		}

		double loop_dt = curr_time - last_time; 
		if(!USE_ORB_SIM) {
			sim->integrate(loop_dt);
			sim->getJointPositions(object_name, model->_q);
			sim->getJointVelocities(object_name, model->_dq);
			// cout << model->_dq.transpose() << endl;
			model->updateModel();
			// update last time
			last_time = curr_time;
			continue;
		}
		
		// check contact state every now and then
		if (timer.elapsedCycles() % 10 == 0) {
			model->updateModel(); // TODO: we shouldn't explicity compute M inverse
			model->rotation(bodyorient, object_link_name);
			model->coriolisPlusGravity(nonlinear_torques);
			// cout << model->_q.transpose() << endl;
			// cout << "Nonlinear torques: " << nonlinear_torques[0] << " " << nonlinear_torques[2] << endl;

			// check if contacts have changed:
			// TODO: implement a change tracker so that we can do delta updates
			// to the contact state model
			cinfo = getDistanceToSurface(model);
				
			//if(LOG_DEBUG) cout << "Surface distance " << cinfo.min_distance << endl;
			if (cinfo.min_distance < FREE_TO_COLL_DIST_THRES) {
				// compute contact velocities to see if any contact is colliding
				contact_model = ContactSpaceModel(model, nonlinear_torques, cinfo);
				pre_collision_contact_vel = contact_model._contactJacobian * model->_dq;
				if(LOG_DEBUG) cout << "Model update: pre coll vel: " << pre_collision_contact_vel.transpose() << endl;

				bool is_colliding = false;
				for(uint i = 0; i < cinfo.contact_points.size(); i++) {
					is_colliding = is_colliding || (pre_collision_contact_vel[3*i + 2] < 0);
					if (pre_collision_contact_vel[3*i + 2] > 1e-4) {
						contact_model.removeActiveContact(i);
					}
				}
				if (is_colliding) {
					state = CapsuleContactState::Colliding;
					if(LOG_DEBUG) cout << "Model update: Collision" << endl;
					// initialize state for CapsuleContactState::Colliding
					contact_model.getActiveContactSpaceMatrices(contact_jacobian, contact_lambda_inv, rhs_coll, rhs_contact);
				} else if(contact_model._activeContacts.size() > 0) {
					state = CapsuleContactState::Rolling;
					if(LOG_DEBUG) cout << "Model update: Rolling" << endl;
					contact_model.getActiveContactSpaceMatrices(contact_jacobian, contact_lambda_inv, rhs_coll, rhs_contact);
				} else {
					state = CapsuleContactState::NoContact;
					if(LOG_DEBUG) cout << "Model update: Points close but no contact" << endl;
				}
			} else {
				state = CapsuleContactState::NoContact;
			}
			if(state == CapsuleContactState::Colliding || state == CapsuleContactState::Rolling) {
				// update visual
				if(cinfo.contact_points.size() == 1) {
					pt_contact_display.setLocalPos(Vector3d(cinfo.contact_points[0][0], cinfo.contact_points[0][1], 0));
					pt_contact_display.setShowEnabled(true);
					line_contact_display.setShowEnabled(false);
				} else {
					// code to visual contact line
					Vector3d temp ((cinfo.contact_points[1][0]+cinfo.contact_points[0][0])/2,
						(cinfo.contact_points[1][1]+cinfo.contact_points[0][1])/2, 
						0.0);
					double theta = atan2(cinfo.contact_points[1][1]-cinfo.contact_points[0][1],
										cinfo.contact_points[1][0]-cinfo.contact_points[0][0]);
					// cout << theta << endl;
					Matrix3d rot1, rot2;
					rot1 << 0, 0, 1,
							0, 1, 0,
							-1, 0, 0;
					rot2 << cos(theta), -sin(theta), 0,
							sin(theta), cos(theta), 0,
								0,			0,		1;
					line_contact_display.setLocalRot(rot2*rot1);
					line_contact_display.setLocalPos(rot2*rot1*Vector3d(0, 0, -0.1) + temp);
					// line_contact_display.setShowEnabled(true);
					pt_contact_display.setShowEnabled(true);
				}
			} else {
				pt_contact_display.setShowEnabled(false);
				line_contact_display.setShowEnabled(false);
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
		bool did_update_contact_matrices = false;
		contact_torques.setZero();
		if (state == CapsuleContactState::Colliding) {
			bool is_post_steady_contact = false;
			CollLCPPointSolution lcp_sol;
			try{
				int tryCollSolveCounter = 0;
				while(tryCollSolveCounter < 500) {
					tryCollSolveCounter++;
					pre_collision_contact_vel = contact_jacobian * model->_dq;
					if(LOG_DEBUG) cout << "Pre coll vel: " << pre_collision_contact_vel.transpose() << endl; 
					double min_coll_vel = 0.0;
					for (uint i = 0; i < contact_model._activeContacts.size(); i++) {
						min_coll_vel = min(min_coll_vel, pre_collision_contact_vel[3*i + 2]);
					}
					double coll_restitution = RESTIT_EPS;
					if (abs(min_coll_vel) < 1e-3) { // We use the fastest collision to determine the coefficient of restitution
						coll_restitution = 0.0; // force inelastic collision to bring to steady contact
					}
					if(LOG_DEBUG) cout << "Num colliding contacts " << contact_model._activeContacts.size() << endl;
					if(contact_model._activeContacts.size() == 1) {
						lcp_sol = solveCollLCPOnePoint (
							contact_lambda_inv,
							pre_collision_contact_vel,
							pre_collision_contact_vel,
							coll_restitution,
							FRICTION_COEFF
						);
					} else if (contact_model._activeContacts.size() == 2) {
						// ensure that we do not have an assymetry due to numerical error
						double v_red = pre_collision_contact_vel[0] + pre_collision_contact_vel[3];
						v_red *= 0.5;
						pre_collision_contact_vel[0] = v_red;
						pre_collision_contact_vel[3] = v_red;
						//TODO: generalize above to any redundancy direction
						lcp_sol = solveCollLCPPoint (
							2, //TODO: remove hardcode
							contact_lambda_inv,
							pre_collision_contact_vel,
							pre_collision_contact_vel,
							coll_restitution,
							FRICTION_COEFF
						);
					}
					if(contact_model._activeContacts.size() > 0) {
						if (lcp_sol.result == LCPSolResult::Success) {
							if(LOG_DEBUG) cout << "LCP Impulse: " << lcp_sol.p_sol.transpose() << endl; 
							if(LOG_DEBUG) cout << "Num contacts " << contact_model._activeContacts.size() << endl; 
						} else {
							cerr << "LCP failed with type: " << static_cast<int>(lcp_sol.result) << endl;
							cout << "Num contacts " << contact_model._activeContacts.size() << endl;
							cout << "pre coll vel: " << pre_collision_contact_vel.transpose() << endl;
							cout << "lambda inv " << contact_lambda_inv << endl;
							cout << "restitution " << coll_restitution << endl;
							cout << "q " << model->_q.transpose() << endl;
							cout << "contact Jacobian " << contact_jacobian << endl;
							throw(std::out_of_range("Collision LCP failure."));
						}						
					}
					if(LOG_DEBUG) {
						double max_impact_vel = 0;
						for(int ii = 0; ii < contact_model._activeContacts.size(); ii++) {
							max_impact_vel = max(max_impact_vel, abs(pre_collision_contact_vel[3*ii+2]));
						}
						max_impact_vel = max(max_impact_vel, 1e-5);
						if(lcp_sol.p_sol[2] > 2*max_impact_vel || 
							(contact_model._activeContacts.size() > 1 && lcp_sol.p_sol[5] > 2*max_impact_vel)) {
							cout << "Much higher impact impulse than expected!!" << endl;
							cout << "Num contacts " << contact_model._activeContacts.size() << endl;
							cout << "Pre collision vel: " << pre_collision_contact_vel.transpose() << endl;
							cout << "Collision lambda inv: " << contact_lambda_inv << endl;
							throw(std::out_of_range("Impulse too high"));
						}
					}
					if(isnan(lcp_sol.p_sol[0]) || isnan(lcp_sol.p_sol[1]) || isnan(lcp_sol.p_sol[2])) {
						cout << "Nan encountered in collision LCP: " << lcp_sol.p_sol.transpose() << endl;
						cout << "Num contacts " << contact_model._activeContacts.size() << endl;
						cout << "dq: " << model->_dq.transpose() << endl;
						cout << "Pre collision vel: " << pre_collision_contact_vel.transpose() << endl;
						cout << "Collision lambda inv: " << contact_lambda_inv << endl;
						break;
					}
					// update joint velocities
					model->_dq += model->_M_inv*(contact_jacobian.transpose()*lcp_sol.p_sol);
					// recheck for collision at ALL contacts:
					post_collision_contact_vel = contact_model._contactJacobian * model->_dq;
					if(LOG_DEBUG) cout << "Post coll: " << post_collision_contact_vel.transpose() << endl;
					if(tryCollSolveCounter == 500) {
						cerr << "Post coll vel: " << post_collision_contact_vel.transpose() << endl;
						cerr << "LCP impulse: " << lcp_sol.p_sol.transpose() << endl;
						throw(std::out_of_range("Too many collision retries"));
					}
					bool did_contact_change = false;
					bool is_post_colliding = false;
					for (uint i = 0; i < cinfo.contact_points.size(); i++) {
						if (post_collision_contact_vel[3*i + 2] <= 1e-4) {
							did_contact_change = did_contact_change || contact_model.addActiveContact(i);
						} else {
							did_contact_change = did_contact_change || contact_model.removeActiveContact(i);
						}
						if(post_collision_contact_vel[3*i + 2] < -1e-10) {
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
			} catch (exception& e) {
				cerr << e.what() << endl;
				cerr << "Exception coll, active points: " << contact_model._activeContacts.size() << endl;
				break;
			}

			// check if steady contact is formed
			if (is_post_steady_contact) {
				// update nonlinear torques and contact RHS because they depend on dq which has changed
				model->coriolisPlusGravity(nonlinear_torques);
				contact_model.recomputeContactRHS(model, nonlinear_torques);
				did_update_contact_matrices = true;
				if(contact_model._activeContacts.size() == 2) {
					if(abs(lcp_sol.p_sol[2]) > 1e-8 && abs(lcp_sol.p_sol[5]) > 1e-8) {
						if(LOG_DEBUG) cout << "2 pt impulse too small. Not updating COP." << endl;	
					} else {
						// update cop sol from collision lcp sol
						last_cop_sol = getLCPForInelasticCollResult(lcp_sol.p_sol, contact_model._contactPointsLocal);
						Vector3d global_cop_pos;
						model->position(global_cop_pos, object_link_name, last_cop_sol.local_cop_pos);
						if(LOG_DEBUG) cout << "Coll deduced COP in global: " << (object_in_world*global_cop_pos).transpose() << endl;
					}
				} else {
					if(LOG_DEBUG) cout << "Rolling with 1 point contact. No COP deduced." << endl;
				}
				state = CapsuleContactState::Rolling;
				if(LOG_DEBUG) cout << "Rolling" << endl;
			} else {
				state = CapsuleContactState::NoContact;
				// contact_display.setShowEnabled(false);
				if(LOG_DEBUG) cout << "No contact" << endl;
			}
		}

		// compute contact torques
		if (state == CapsuleContactState::Rolling) {
			// break;

			// recompute active contacts: TODO: this might be optimized when a collision was resolved in this sim cycle
			bool did_contact_change = false;
			bool is_post_colliding = false;
			post_collision_contact_vel = contact_model._contactJacobian * model->_dq;
			for (uint i = 0; i < cinfo.contact_points.size(); i++) {
				if (post_collision_contact_vel[3*i + 2] <= 1e-4) { // note that we allow some positive contact velocity
					did_contact_change = did_contact_change || contact_model.addActiveContact(i);
				} else {
					did_contact_change = did_contact_change || contact_model.removeActiveContact(i);
				}
				if(post_collision_contact_vel[3*i + 2] < -1e-5) {
					is_post_colliding = true;
				}
			}
			if(is_post_colliding) {
				state = CapsuleContactState::Colliding;
				if(LOG_DEBUG) cout << "Rolling to colliding" << endl;
			}
			if(contact_model._activeContacts.size() == 0) {
				state = CapsuleContactState::NoContact;
				// contact_display.setShowEnabled(false);
				if(LOG_DEBUG) cout << "No contact" << endl;
			} else {
				if(did_contact_change || did_update_contact_matrices) {
					contact_model.getActiveContactSpaceMatrices(contact_jacobian, contact_lambda_inv, rhs_coll, rhs_contact);
				}
				post_collision_contact_vel = contact_jacobian * model->_dq;

				// TODO: switch to COP model
				CollLCPPointSolution lcp_sol;
				try {
					if(LOG_DEBUG) cout << "Num contacts " << contact_model._activeContacts.size() << endl;
					if(contact_model._activeContacts.size() == 1) {
						if(rhs_contact[2] > 0) {
							cerr << "RHS contact is positive!" << endl;
							cout << "dq: " << model->_dq.transpose() << endl;
							cerr << "Num contacts " << contact_model._activeContacts.size() << endl;
							cerr << "Contact rhs: " << rhs_contact.transpose() << endl;
							cerr << "Nonlinear torques: " << nonlinear_torques.transpose() << endl;
							break;
						}
						lcp_sol = solveCollLCPOnePoint (
							contact_lambda_inv,
							rhs_contact,
							post_collision_contact_vel,
							0.0,
							FRICTION_COEFF
						);
						if (lcp_sol.result == LCPSolResult::Success) {
							if(LOG_DEBUG) cout << "Contact LCP Contact force: " << lcp_sol.p_sol.transpose() << endl;
						} else {
							cerr << "Contact LCP failed with type: " << static_cast<int>(lcp_sol.result) << endl;
							cout << "Num contacts " << contact_model._activeContacts.size() << endl;
							cout << "dq: " << model->_dq.transpose() << endl;
							cout << "Last post contact vel: " << post_collision_contact_vel.transpose() << endl;
							cout << "Contact lambda inv: " << contact_lambda_inv << endl;
							cout << "Contact rhs: " << rhs_contact.transpose() << endl;
							break;
						}
						contact_force = lcp_sol.p_sol;
						if(isnan(contact_force[0]) || isnan(contact_force[1]) || isnan(contact_force[2])) {
							cout << "Nan encountered " << lcp_sol.p_sol.transpose() << endl;
							cout << "Num contacts " << contact_model._activeContacts.size() << endl;
							cout << "dq: " << model->_dq.transpose() << endl;
							cout << "Last post contact vel: " << post_collision_contact_vel.transpose() << endl;
							cout << "Contact lambda inv: " << contact_lambda_inv << endl;
							cout << "Contact Jacobian: " << contact_jacobian << endl;
							cout << "Mass matrix inv: " << model->_M_inv << endl;
							cout << "Contact rhs: " << rhs_contact.transpose() << endl;
							cout << "cinfo contacts size " << cinfo.contact_points.size() << endl;
							cout << "Contact points, 1: " << contact_model._contact_description.contact_points[0].transpose() <<
							    " 2: " << contact_model._contact_description.contact_points[1].transpose() << endl;
						    Vector3d temp;
						    model->position(temp, object_link_name, end1_local);
							cout << "point 1: " << (object_in_world * temp).transpose() << endl;
							model->position(temp, object_link_name, end2_local);
							cout << "point 2: " << (object_in_world * temp).transpose() << endl;
							break;
						}
						contact_torques = contact_jacobian.transpose() * contact_force;
						// get COP solution from one point contact lcp sol in case we need it
						bool isContact0Active = (*(contact_model._activeContacts.begin()) == 0);
						last_cop_sol = getLCPForOnePtContactResult(lcp_sol.p_sol, contact_model._contactPointsLocal, isContact0Active);
						Vector3d global_cop_pos;
						model->position(global_cop_pos, object_link_name, last_cop_sol.local_cop_pos);
						if(LOG_DEBUG) cout << "1 pt contact deduced COP in global: " << (object_in_world*global_cop_pos).transpose() << endl;
					} else {
						if(rhs_contact[2] > 0 && rhs_contact[5] > 0) {
							cerr << "RHS contact is positive!" << endl;
							cout << "dq: " << model->_dq.transpose() << endl;
							cerr << "Num contacts " << contact_model._activeContacts.size() << endl;
							cerr << "Contact rhs: " << rhs_contact.transpose() << endl;
							cerr << "Nonlinear torques: " << nonlinear_torques.transpose() << endl;
							break;
						}
						// lcp_sol = solveCollLCPPoint (
						// 	2, //TODO: remove hardcode
						// 	contact_lambda_inv,
						// 	rhs_contact,
						// 	post_collision_contact_vel,
						// 	0.0,
						// 	FRICTION_COEFF
						// );
						Vector3d r_point0_to_last_cop_in_cop_frame = contact_model._R_body_COP_frame*
																(last_cop_sol.local_cop_pos - contact_model._contactPointsLocal[0]);
						getCOPJ6FullFromPoint0J6Full(
							contact_model._point0J6,
							r_point0_to_last_cop_in_cop_frame,
							last_cop_J6
						);
						getCOPLambdaInv6FullFromPoint0LambdaInv6Full(
							contact_model._point0_lambdaInv_6,
							r_point0_to_last_cop_in_cop_frame,
							last_cop_lambdainv_full
						);
						getCOPRhsNonlinFullFromPoint0RhsNonlinFull(
							contact_model._point0_rhs_nonlin_6,
							r_point0_to_last_cop_in_cop_frame,
							last_cop_rhsnonlin_full
						);
						if(LOG_DEBUG) cout << "COP Jacobian " << last_cop_J6 << endl;
						ContactCOPSolution cop_sol = resolveCOPLineContactWithLastCOPSol(
							last_cop_lambdainv_full, // 6 x 6 defined at last_COP_sol local point
							last_cop_rhsnonlin_full, // 6 defined at last_COP_sol local point
							contact_model.getOmegaCOPFrame(model->_dq), // defined in global COP frame
							contact_model._R_body_COP_frame * last_cop_sol.local_cop_pos, // position vector from center of object to last_COP_sol local point, defined in global COP frame
							last_cop_J6 * model->_dq, // full 6 dof at last_COP_sol local point
							contact_model._R_body_COP_frame, // rotation from body frame to new global COP frame
							contact_model._contactPointsLocal, // points are in the local body frame, NOT in the global COP frame
							FRICTION_COEFF,
							last_cop_sol
						);
						last_cop_sol = cop_sol;
						if(LOG_DEBUG) {
							cout << "COP result: " << static_cast<int>(cop_sol.result) << endl;
							cout << "COP position: " << cop_sol.local_cop_pos.transpose() << endl;
							cout << "COP force: " << cop_sol.force_sol.transpose() << endl;
						}
						if(last_cop_sol.result == COPSolResult::Success) {
							r_point0_to_last_cop_in_cop_frame = contact_model._R_body_COP_frame*
																(last_cop_sol.local_cop_pos - contact_model._contactPointsLocal[0]);
							getCOPJ6FullFromPoint0J6Full(
								contact_model._point0J6,
								r_point0_to_last_cop_in_cop_frame,
								last_cop_J6
							);
							contact_torques = last_cop_J6.block(0,0,3,dof).transpose() * last_cop_sol.force_sol.segment<3>(0);
							contact_torques += last_cop_J6.block(4,0,2,dof).transpose() * last_cop_sol.force_sol.segment<2>(3);
							// update contact force display
							cop_force_sensor->_data->_transform_in_link.translation() = last_cop_sol.local_cop_pos;
							cop_force_sensor->_data->_force = contact_model._rotXInterpoint.transpose()*last_cop_sol.force_sol.segment<3>(0);
							// cout << last_cop_sol.force_sol[4] << endl;
							cop_force_sensor->_data->_moment << 0.0, 0.0, last_cop_sol.force_sol[4];
							pt_contact_display.setLocalPos(
								cop_force_display->_display_line_force->m_pointA
							);
						} else {
							throw(std::runtime_error("COP sol failed."));
						}
					}
				} catch (exception& e) {
					cerr << e.what() << endl;
					cerr << "Seg fault contact LCP, active contact size: " << contact_model._activeContacts.size() << endl;
					break;
				}
				if(cinfo.min_distance < -FREE_TO_COLL_DIST_THRES) {
					cout << "Num contacts " << contact_model._activeContacts.size() << endl;
					cout << "dq: " << model->_dq.transpose() << endl;
					cout << "Last post contact vel: " << post_collision_contact_vel.transpose() << endl;
					cout << "Contact lambda inv: " << contact_lambda_inv << endl;
					cout << "Contact Jacobian: " << contact_jacobian << endl;
					cout << "Mass matrix inv: " << model->_M_inv << endl;
					cout << "Contact rhs: " << rhs_contact.transpose() << endl;
					cout << "cinfo contacts size " << cinfo.contact_points.size() << endl;
					cout << "Contact points, 1: " << contact_model._contact_description.contact_points[0].transpose() <<
					    " 2: " << contact_model._contact_description.contact_points[1].transpose() << endl;
				    break;
				}
				// contact_torques = contact_jacobian.transpose() * contact_force;
			}
		}

		// integrate //TODO: switch to Euler-Heun
		model->_ddq = model->_M_inv * (-nonlinear_torques + contact_torques);
		VectorXd temp(6), temp2(7), temp3(7);
		temp2 << model->_q; 
		temp = model->_dq + 0.5*model->_ddq * loop_dt;
		// model->_q += temp * loop_dt;
		model->_q.segment<3>(0) += temp.segment<3>(0) * loop_dt;
		if (temp.segment<3>(3).norm() > 1e-15) {
			model->_q.segment<4>(3) = RigidBodyDynamics::Math::Quaternion(model->_q.segment<4>(3)).timeStep(bodyorient*temp.segment<3>(3), loop_dt);
			// model->_q.segment<4>(3) /= model->_q.segment<4>(3).norm();
		}
		temp3 << model->_q;
		// if(LOG_DEBUG && state == CapsuleContactState::Rolling) {
		// 	ContactInfo cinfo2 = getDistanceToSurface(model);
		// 	Vector3d flocal_point;
		// 	model->updateModel();
		// 	getContactPointLocalFrame(model, cinfo2.contact_points[0], flocal_point);
		// 	cout << "Future local point " << flocal_point.transpose() << endl;
		// 	model->_q << temp2;
		// 	model->updateModel();
		// 	MatrixXd Jvfc;
		// 	model->Jv(Jvfc, object_link_name, flocal_point);
		// 	Jvfc = object_in_world.linear() * Jvfc;
		// 	cout << "Future contact velocity " << (Jvfc * model->_dq).transpose() << endl;

		// 	model->_q << temp3;
		// 	Vector3d accel;
		// 	Vector3d local_point;
		// 	model->updateModel();
		// 	getContactPointLocalFrame(model, cinfo.contact_points[0], local_point);
		// 	model->linearAcceleration(accel, object_link_name, local_point);
		// 	cout << "loop dt" << loop_dt << endl;
		// 	cout << "Contact 0 acceleration " << (object_in_world.linear()*accel).transpose() << endl;
		// 	// cout << "J ddq " << (contact_jacobian*model->_ddq).transpose() << endl;
		// 	// MatrixXd objectJw;
		// 	// model->Jw(objectJw, object_link_name);
		// 	// Vector3d omega, linvel;
		// 	// omega = objectJw*model->_dq;
		// 	// linvel = contact_jacobian*model->_dq;
		// 	// Vector3d rworld;
		// 	// model->position(rworld, object_link_name, Vector3d::Zero());
		// 	// rworld = cinfo.contact_points[0] - object_in_world * rworld;
		// 	// cout << rworld.transpose() << endl;
		// 	// cout << ((object_in_world.linear()*omega).cross(linvel)).transpose() << endl;
		// 	// cout << ((object_in_world.linear()*omega).cross((object_in_world.linear()*omega).cross(rworld))).transpose() << endl;
		// 	cout << "Contact 0 velocity " << (contact_model._contactJacobian*model->_dq).transpose() << endl;
		// }
		model->_dq += model->_ddq * loop_dt;

		

		// if (!fTimerDidSleep) {
		// 	cout << "Warning: timer underflow! dt: " << loop_dt << "\n";
		// }

		// update last time
		last_time = curr_time;
	}
	cout << "Simulation loop finished, average loop frequency: "
		<< timer.elapsedCycles()/timer.elapsedTime() << endl;
	cout << "Min distance of object to ground at sim end: " << cinfo.min_distance << endl;
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
    GLFWwindow* window = glfwCreateWindow(windowW, windowH, "SAI2.0 - COP SIM 1", NULL, NULL);
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
    if ((key == 'd' || key == 'D') && action == GLFW_PRESS) {
    	LOG_DEBUG = !LOG_DEBUG;
    }
    if ((key == 'z' || key == 'Z') && action == GLFW_PRESS) {
        // change camera
        camera_name = "camera_zoom";
    }
    if ((key == '1') && action == GLFW_PRESS) {
        // change camera
        camera_name = "camera_front";
    }
    if ((key == '2') && action == GLFW_PRESS) {
        // change camera
        camera_name = "camera_side";
    }
    if ((key == 'p' || key == 'P') && action == GLFW_PRESS) {
    	f_pause_sim = (f_pause_sim)? false: true;
    }
    bool set = (action != GLFW_RELEASE);
    switch(key) {
		case GLFW_KEY_RIGHT:
			fTransXp = set;
			break;
		case GLFW_KEY_LEFT:
			fTransXn = set;
			break;
		case GLFW_KEY_UP:
			fTransYp = set;
			break;
		case GLFW_KEY_DOWN:
			fTransYn = set;
			break;
		case 'i':
		case 'I':
			fTransZp = set;
			break;
		case 'o':
		case 'O':
			fTransZn = set;
			break;
		default:
			break;
    }
}
