/*  04-two-arm-rolling - main.cpp
This demo shows two robot arms rolling a rolling pin.

Author: Shameek Ganguly shameekg@stanford.edu
Date: 07/08/2020
*/

#include <iostream>
#include <string>
#include <thread>
#include <chrono>
#include <math.h>
#include <vector>
using namespace std;

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

#include "cop_simulator/COPSimulator.h"

#include <GLFW/glfw3.h> //must be loaded after loading opengl/glew as part of graphicsinterface

using namespace chai3d;

const string world_fname = "resources/04-two-arm-rolling/world.urdf";

const string robot_fname = "../resources/kuka_iiwa/kuka_iiwa.urdf";
const string robot_name = "IIWA";
const string robot_ee_name = "link6";
const Eigen::Vector3d arm1_ee_local_pos(0.03, 0, 0.028);
Eigen::VectorXd arm1_home_qpos;
const double robot_ee_capsule_radius = 0.05;
const double robot_ee_capsule_length = 0.1;

const string object_fname = "resources/04-two-arm-rolling/roller_object.urdf";
const string roller_object_name = "Roller";
const double roller_capsule_mid_radius = 0.05;
const double roller_capsule_right_radius = 0.03;
const double roller_capsule_left_radius = 0.03;
const double roller_capsule_mid_length = 0.3;
const double roller_capsule_right_length = 0.15;
const double roller_capsule_left_length = 0.15;
const string object_link_name = "object";
const string box_name = "Box";

string camera_name = "camera_front";

// flags for scene camera movement
static bool fTransXp = false;
static bool fTransXn = false;
static bool fTransYp = false;
static bool fTransYn = false;
static bool fTransZp = false;
static bool fTransZn = false;

// logging
bool LOG_DEBUG = false;

// pause
bool f_pause_sim = false;

// sim state display
chai3d::cLabel* sim_state_label;
chai3d::cFontPtr label_font;

// COP visualization
chai3d::cShapeSphere pt_contact_display1(0.01);
chai3d::cShapeSphere pt_contact_display2(0.01);

// simulation loop
bool fSimulationRunning = false;
void simulation(Sai2COPSim::COPSimulator* sim);
void control_arm1(Sai2Model::Sai2Model* model, Sai2COPSim::COPSimulator* sim);

// initialize window manager
GLFWwindow* glfwInitialize();

// callback to print glfw errors
void glfwError(int error, const char* description);

// callback when a key is pressed
void keySelect(GLFWwindow* window, int key, int scancode, int action, int mods);

int main(int argc, char** argv) {
	const double restitution = 0.4;
    const double friction = 0.2;

    // load sai2 simulation world
    // TODO: this is currently needed in order to load the base transform for each object
    // as well as the transform the base plane and the gravity vector
    // TODO: write a parser for world URDF to the cop sim world
    auto sai2_sim = new Simulation::Sai2Simulation(world_fname, false);

    Affine3d arm_in_world;
    arm_in_world.translation() = sai2_sim->_world->getBaseNode(robot_name)->getLocalPos().eigen();
    arm_in_world.linear() = sai2_sim->_world->getBaseNode(robot_name)->getLocalRot().eigen();

    Affine3d roller_in_world;
    roller_in_world.translation() = sai2_sim->_world->getBaseNode(roller_object_name)->getLocalPos().eigen();
    roller_in_world.linear() = sai2_sim->_world->getBaseNode(roller_object_name)->getLocalRot().eigen();

    Affine3d static_plane_in_world;
    static_plane_in_world.translation() = sai2_sim->_world->getBaseNode(box_name)->getLocalPos().eigen();
    static_plane_in_world.linear() = sai2_sim->_world->getBaseNode(box_name)->getLocalRot().eigen();
    static_plane_in_world.translation() += static_plane_in_world.linear().col(2)*0.1/2.0; // z-height of box
    // cout << static_plane_in_world.translation().transpose() << endl;

    Vector3d grav_vector = sai2_sim->_world->getGravity().eigen();

    // load graphics scene
    auto graphics = new Sai2Graphics::Sai2Graphics(world_fname, false);
    // set object graphics to wireframe, and show frame for last link
    // graphics->showLinkFrame(true, object_name, object_link_name);
    graphics->showWireMeshRender(true, roller_object_name, object_link_name);

    // load objects
    auto roller_object = new Sai2Model::Sai2Model(object_fname, false, roller_in_world, roller_in_world.linear().transpose()*grav_vector);

    // load robots
    auto arm = new Sai2Model::Sai2Model(robot_fname, false, arm_in_world, arm_in_world.linear().transpose()*grav_vector);
    auto arm_control_model = new Sai2Model::Sai2Model(robot_fname, false, arm_in_world, arm_in_world.linear().transpose()*grav_vector);
    arm1_home_qpos.setZero(arm->dof());
    arm1_home_qpos << 90/180.0*M_PI,
                45/180.0*M_PI, //59.2 will cause it to be in contact with the capsule
                0/180.0*M_PI,
                -70.0/180.0*M_PI,
                -10/180.0*M_PI,
                60/180.0*M_PI,
                180/180.0*M_PI;
    arm->_q = arm1_home_qpos;

    // add some initial velocity. TODO: move to parser
    // coobject1->_dq[5] = 0.1;

    // set up COP sim world
    auto cop_sim = new Sai2COPSim::COPSimulator(friction, restitution);

    // add robot arm
    cop_sim->addObject(robot_name, arm);

    // add robot arm primitives
    Affine3d tf_robot_ee_capsule = Affine3d::Identity();
    tf_robot_ee_capsule.translation() = arm1_ee_local_pos;
    cop_sim->addCapsuleToObject(robot_name, robot_ee_name, "capsule_ee", robot_ee_capsule_radius, robot_ee_capsule_length, tf_robot_ee_capsule);
    
    // add roller object
    cop_sim->addObject(roller_object_name, roller_object);
    
    // add primitives to roller object
    cop_sim->addCapsuleToObject(roller_object_name, object_link_name, "capsule_mid", roller_capsule_mid_radius, roller_capsule_mid_length, Affine3d::Identity());
    Affine3d tf_right_capsule = Affine3d::Identity();
    tf_right_capsule.translation() << 0.305, 0, 0;
    cop_sim->addCapsuleToObject(roller_object_name, object_link_name, "capsule_right", roller_capsule_right_radius, roller_capsule_right_length, tf_right_capsule);
    Affine3d tf_left_capsule = Affine3d::Identity();
    tf_left_capsule.translation() << -0.305, 0, 0;
    cop_sim->addCapsuleToObject(roller_object_name, object_link_name, "capsule_left", roller_capsule_left_radius, roller_capsule_left_length, tf_left_capsule);
    
    // add plane
    cop_sim->addPlane(box_name, static_plane_in_world.linear().col(2), static_plane_in_world.translation());

    // TODO: force sim/ display

    // display sim state in a label
    label_font = NEW_CFONTCALIBRI24();
    sim_state_label = new chai3d::cLabel(label_font);
    auto front_camera = graphics->getCamera("camera_front");
    front_camera->m_frontLayer->addChild(sim_state_label);
    auto zoom_camera = graphics->getCamera("camera_zoom");
    zoom_camera->m_frontLayer->addChild(sim_state_label);
    sim_state_label->setText("Initializing");

    // visualize contact
    pt_contact_display1.setShowEnabled(false);
    pt_contact_display1.m_material->setBrownMaroon();
    pt_contact_display1.m_material->setShininess(100);
    graphics->_world->addChild(&pt_contact_display1);
    pt_contact_display1.setLocalPos(Vector3d(0.0, 0.0, 0.4));
    pt_contact_display2.setShowEnabled(false);
    pt_contact_display2.m_material->setBrownMaroon();
    pt_contact_display2.m_material->setShininess(100);
    graphics->_world->addChild(&pt_contact_display2);
    pt_contact_display2.setLocalPos(Vector3d(0.0, 0.0, 0.4));

    // ---- INITIALIZATION COMPLETE. STARTING MULTITHREADED APP ----
    // initialize GLFW window
    GLFWwindow* window = glfwInitialize();

    // set callbacks
    glfwSetKeyCallback(window, keySelect);

    // start the simulation thread first
    fSimulationRunning = true;
    thread sim_thread(simulation, cop_sim);
    sim_state_label->setText("Sim running");

    // start the control thread for arm 1
    thread arm1_thread(control_arm1, arm_control_model, cop_sim);

    // TODO: visualize contact

    // ---- GRAPHICS LOOP ON MAIN THREAD ----
    Eigen::Vector3d camera_pos, camera_lookat, camera_vertical;
    while (!glfwWindowShouldClose(window)) {
        // update model for object. called from simulation thread
        // coobject->updateModel();

        // update graphics. this automatically waits for the correct amount of time
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        graphics->updateGraphics(robot_name, arm);
        graphics->updateGraphics(roller_object_name, roller_object);
        // TODO: update force display
        // cop_force_display->update();

        // render
        graphics->render(camera_name, width, height);
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

    // ---- TEARDOWN MULTITHREADED APP ----
    // stop simulation
    fSimulationRunning = false;
    sim_thread.join();
    arm1_thread.join();

    // destroy context
    glfwDestroyWindow(window);

    // terminate
    glfwTerminate();

	return 0;
}

void simulation(Sai2COPSim::COPSimulator* sim) {
	fSimulationRunning = true;

    auto object = sim->_arb_manager.getBody(roller_object_name);
    auto cmodel = sim->_contact_model;
    auto& geom_manager = sim->_geom_manager;

    // create a timer
    LoopTimer timer;
    timer.initializeTimer();
    timer.setLoopFrequency(20000); //1500Hz timer
    double last_time = timer.elapsedTime(); //secs
    bool fTimerDidSleep = true;
    while (fSimulationRunning) {
        fTimerDidSleep = timer.waitForNextLoop();

        // integrate forward
        double curr_time = timer.elapsedTime();
        if(f_pause_sim) {
            last_time = curr_time;
            continue;
        }
        double loop_dt = curr_time - last_time;

        try {
            sim->integrate(loop_dt);

            if(timer.elapsedCycles() % 500 == 0) {
                // cout << object->jtau_act(1) << endl;

                // visualize cop if present
                if(cmodel != NULL
                    && cmodel->_contact_island_models_size > 0
                    && cmodel->_contact_island_models[0]._active_contacts.size() > 1
                ) {
                    const auto& island = cmodel->_contact_island_models[0];
                    const auto& prim_state1 = island._pair_state[0];
                    if(prim_state1.isValid()) {
                        pt_contact_display1.setLocalPos(prim_state1._cop_pos);
                        pt_contact_display1.setShowEnabled(true);
                    }
                    const auto& prim_state2 = island._pair_state[1];
                    if(prim_state2.isValid()) {
                        pt_contact_display2.setLocalPos(prim_state2._cop_pos);
                        pt_contact_display2.setShowEnabled(true);
                    }
                } else if(geom_manager._prim_prim_distances[4][1]->min_distance > 0.01) {
                    pt_contact_display1.setShowEnabled(false);
                    pt_contact_display2.setShowEnabled(false);
                }
            }

        } catch (exception& e) {
            cerr << e.what() << endl;
            sim_state_label->setText("Simulation failed");
            break;
        }

        // update last time
        last_time = curr_time;
    }
    cout << "Simulation loop finished, average loop frequency: "
        << timer.elapsedCycles()/timer.elapsedTime() << endl;
    cout << "Total simulated time: " << timer.elapsedTime() << endl;
    cout << "Sim speed up over real time "<< timer.elapsedTime()/sim->_time_total << "x" << endl;
    sim->printTimeAnalytics();
}

// TODO: add Sai2Simulation for comparison

void control_arm1(Sai2Model::Sai2Model* model, Sai2COPSim::COPSimulator* sim) {
    // control variables
    VectorXd tau_gravity;
    auto arb = sim->_arb_manager.getBody(robot_name);

    const uint dof = model->dof();

    MatrixXd J0(6, dof);
    MatrixXd Lambda(6, 6);
    MatrixXd NT(dof, dof);

    VectorXd F_control(6); //6DOF in world frame
    const double kp = 10;
    const double kd = 8;
    const double kjd = 8;

    Vector3d des_pos_world(0, -0.01, 0.05+0.05);
    Matrix3d des_ori;
    des_ori << 0, 1, 0,
               1, 0, 0,
               0, 0, -1;

    Vector3d ee_pos;
    Vector3d ori_error;
    Matrix3d ee_ori;

    // create a timer
    LoopTimer timer;
    timer.initializeTimer();
    timer.setLoopFrequency(1000); //1000Hz timer
    while(fSimulationRunning) {
        timer.waitForNextLoop();
        // copy over q and dq from sim
        model->_q = arb->_model->_q;
        model->_dq = arb->_model->_dq;

        // update model
        model->updateModel();
        model->gravityVector(tau_gravity);
        model->J_0(J0, robot_ee_name, arm1_ee_local_pos);
        J0.block(0,0,3,dof) = model->_T_world_robot.linear()*J0.block(0,0,3,dof);
        J0.block(3,0,3,dof) = model->_T_world_robot.linear()*J0.block(3,0,3,dof);
        Lambda = (J0 * model->_M_inv * J0.transpose()).inverse();
        NT = MatrixXd::Identity(dof, dof) - J0.transpose()*Lambda*J0*model->_M_inv;

        model->positionInWorld(ee_pos, robot_ee_name, arm1_ee_local_pos);
        model->rotationInWorld(ee_ori, robot_ee_name);
        Sai2Model::orientationError(ori_error, des_ori, ee_ori);

        // compute control torques
        Vector3d set_pos = des_pos_world;
        set_pos(1) += 0.06*sin(timer.elapsedTime()*1.0);
        F_control.segment<3>(0) = -kp*(ee_pos - set_pos) - kd*J0.block(0,0,3,dof)*model->_dq;
        F_control.segment<3>(3) = -kp*ori_error - kd*J0.block(3,0,3,dof)*model->_dq;


        // set control torques
        arb->jtau_act = tau_gravity + J0.transpose()*Lambda*F_control + NT*(model->_M * -kjd * model->_dq);
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
    GLFWwindow* window = glfwCreateWindow(windowW, windowH, "SAI2.0 - COP SIM 2", NULL, NULL);
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
        sim_state_label->setText((f_pause_sim)? "PAUSED": "UNPAUSED");
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
