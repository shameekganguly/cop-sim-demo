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

#include "Sai2Model.h"
#include "Sai2Graphics.h"
#include "Sai2Simulation.h"
#include "dynamics3d.h"
#include "rbdl/rbdl.h"

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
const double radius = 0.15; //TODO: remove this hardcode.

// thresholds
const double FREE_TO_COLL_DIST_THRES = 0.003; // mm
const double RESTIT_EPS = 0.3;
const double FRICTION_COEFF = 0.1;

Affine3d object_in_world;

const uint sphere_height_joint_id = 2;
double getDistanceToSurface(const VectorXd& object_joints) {
	return object_joints[sphere_height_joint_id] + object_in_world.translation()[2] - radius;
}

void getContactPoint(Sai2Model::Sai2Model* sphere, Vector3d& point) {
	sphere->positionInWorld(point, object_link_name, Vector3d(0, 0, -radius));
}

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
	coobject->_dq[1] = 0.1;

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

	enum SphereContactState {
		NoContact,
		Colliding,
		Rolling,
		Sticking
	} state;
	state = SphereContactState::NoContact;
	
	// state vars
	const uint dof = model->dof();
	VectorXd contact_torques(dof); contact_torques.setZero();
	VectorXd nonlinear_torques(dof); nonlinear_torques.setZero();
	Vector3d contact_force; contact_force.setZero();
	double contact_moment = 0.0; // for sphere, moment is 1-dof
	double normal_impulse = 0.0;
	Vector2d tangent_impulse(0, 0);
	Vector3d contact_point;
	Matrix3d contact_frame = Matrix3d::Identity(); // 1st and 2nd cols are tangent, 3rd is normal
	MatrixXd contact_jacobian(3, model->dof());
	Matrix3d contact_lambda_inv;
	Vector3d rolling_contact_impulse, sliding_contact_impulse;
	Vector3d pre_collision_contact_vel, post_collision_contact_vel;
	Vector2d slip_direction;

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

			switch (state) {
				case SphereContactState::NoContact:
					if (getDistanceToSurface(model->_q) < FREE_TO_COLL_DIST_THRES && model->_dq[sphere_height_joint_id] < 0) {
						state = SphereContactState::Colliding;
						cout << "Collision" << endl;
						// initialize state for SphereContactState::Colliding
						getContactPoint(model, contact_point);
						model->Jv(contact_jacobian, object_link_name, Vector3d(0, 0, -radius));
						// check if rolling
						// TODO: rotate to contact frame
						contact_lambda_inv = contact_jacobian * model->_M_inv * contact_jacobian.transpose();
					}
					break;
				case SphereContactState::Colliding:
					break;
				case SphereContactState::Rolling:

				break;
				case SphereContactState::Sticking:
				break;
			}
			// cout << model->_q.transpose() << endl;
		}

		// compute contact force
		contact_torques.setZero();
		switch (state) {
			case SphereContactState::NoContact:
				// nothing to do
				break;
			case SphereContactState::Colliding:
				pre_collision_contact_vel = contact_jacobian * model->_dq;
				// handle frictional collision
				rolling_contact_impulse = contact_lambda_inv.ldlt().solve(
				(Vector3d(0, 0, -RESTIT_EPS * pre_collision_contact_vel[2]) -
					 pre_collision_contact_vel)
				);
				// check for friction cone
				if (rolling_contact_impulse.segment<2>(0).norm() < rolling_contact_impulse[2]*FRICTION_COEFF) {
					model->_dq += model->_M_inv*(contact_jacobian.transpose()*rolling_contact_impulse);
				} else {
					sliding_contact_impulse = rolling_contact_impulse;
					// TODO: determine slip direction
					slip_direction = pre_collision_contact_vel.segment<2>(0);
					if(slip_direction.norm() < 1e-5) {
						// TODO: force stick
					} else {
						slip_direction = slip_direction/slip_direction.norm();
						sliding_contact_impulse << -FRICTION_COEFF*slip_direction, 1.0;
						normal_impulse = 1/(contact_lambda_inv.row(2).dot(sliding_contact_impulse.transpose())) 
											* -(1 + RESTIT_EPS) * pre_collision_contact_vel[2];
						sliding_contact_impulse *= normal_impulse;
						model->_dq += model->_M_inv*(contact_jacobian.transpose()*sliding_contact_impulse);
					}
				}

				// check if sticking
				if (model->_dq[sphere_height_joint_id] < 1e-3) {
					state = SphereContactState::Rolling;
					model->_dq[sphere_height_joint_id] = 0;
					contact_torques[sphere_height_joint_id] = nonlinear_torques[sphere_height_joint_id];
					cout << "Rolling" << endl;
				} else {
					cout << "No contact" << endl;
					state = SphereContactState::NoContact;
				}
				// no contact torques to apply
				break;
			case SphereContactState::Rolling:
				contact_torques[sphere_height_joint_id] = nonlinear_torques[sphere_height_joint_id];
				break;
			case SphereContactState::Sticking:
				break;
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
