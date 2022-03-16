// collisionDetectionDemo.cpp
/*
This demo shows collision detection between two capsules.

Author: Shameek Ganguly shameekg@stanford.edu
Date: 07/10/2020
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

const string world_fname = "resources/02-two-capsules/worldCollDetect.urdf";
const string object_fname = "resources/02-two-capsules/capsule_object.urdf";
const string object1_name = "Capsule1";
const double capsule1_radius = 0.1;
const double capsule2_radius = 0.1;
const double capsule3_radius = 0.1;
const double capsule4_radius = 0.1;
const double capsule5_radius = 0.1;
const double capsule1_length = 0.2;
const double capsule2_length = 0.2;
const double capsule3_length = 0.2;
const double capsule4_length = 0.2;
const double capsule5_length = 0.2;
const string object2_name = "Capsule2";
const string object3_name = "Capsule3";
const string object4_name = "Capsule4";
const string object5_name = "Capsule5";
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

// Distance visualization
chai3d::cShapeSphere ptA_display(0.02);
chai3d::cShapeSphere ptB_display(0.02);
chai3d::cShapeLine min_dist_display;

// simulation loop
bool fSimulationRunning = false;
void simulation(Sai2COPSim::COPSimulator* sim);

// initialize window manager
GLFWwindow* glfwInitialize();

// callback to print glfw errors
void glfwError(int error, const char* description);

// callback when a key is pressed
void keySelect(GLFWwindow* window, int key, int scancode, int action, int mods);

int main(int argc, char** argv) {
	const double restitution = 0.4;
    const double friction = 0.1;

    // load sai2 simulation world
    // TODO: this is currently needed in order to load the base transform for each object
    // as well as the transform the base plane and the gravity vector
    // TODO: write a parser for world URDF to the cop sim world
    auto sai2_sim = new Simulation::Sai2Simulation(world_fname, false);

    Affine3d capsule1_in_world;
    capsule1_in_world.translation() = sai2_sim->_world->getBaseNode(object1_name)->getLocalPos().eigen();
    capsule1_in_world.linear() = sai2_sim->_world->getBaseNode(object1_name)->getLocalRot().eigen();
    // cout << capsule1_in_world.translation().transpose() << endl;

    Affine3d capsule2_in_world;
    capsule2_in_world.translation() = sai2_sim->_world->getBaseNode(object2_name)->getLocalPos().eigen();
    capsule2_in_world.linear() = sai2_sim->_world->getBaseNode(object2_name)->getLocalRot().eigen();
    // cout << capsule2_in_world.translation().transpose() << endl;

    Vector3d grav_vector = sai2_sim->_world->getGravity().eigen();

    // load graphics scene
    auto graphics = new Sai2Graphics::Sai2Graphics(world_fname, false);
    // set object graphics to wireframe, and show frame for last link
    // graphics->showLinkFrame(true, object_name, object_link_name);
    graphics->showWireMeshRender(true, object1_name, object_link_name);
    graphics->showWireMeshRender(true, object2_name, object_link_name);

    // visualize contact
    ptA_display.setShowEnabled(false);
    ptA_display.m_material->setBrownMaroon();
    ptA_display.m_material->setShininess(100);
    graphics->_world->addChild(&ptA_display);
    ptA_display.setLocalPos(Vector3d(0.0, 0.0, 0.4));
    ptB_display.setShowEnabled(false);
    ptB_display.m_material->setBrownMaroon();
    ptB_display.m_material->setShininess(100);
    graphics->_world->addChild(&ptB_display);
    ptB_display.setLocalPos(Vector3d(0.0, 0.0, 0.4));

    min_dist_display.setShowEnabled(false);
    min_dist_display.m_colorPointA.setGreenYellowGreen();
    min_dist_display.m_colorPointB.setGreenYellowGreen();
    graphics->_world->addChild(&min_dist_display);
    min_dist_display.setLineWidth(3.0);

    // load objects
    auto coobject1 = new Sai2Model::Sai2Model(object_fname, false, capsule1_in_world, capsule1_in_world.linear().transpose()*grav_vector);
    auto coobject2 = new Sai2Model::Sai2Model(object_fname, false, capsule2_in_world, capsule2_in_world.linear().transpose()*grav_vector);

    // add some initial velocity. TODO: move to parser
    Vector3d world_object1_tvel(0, 0.1, 0);
    Vector3d world_object2_tvel(0, -0.1, 0);
    Vector3d world_object1_rvel(0, 2.0, 0);
    Vector3d world_object2_rvel(0, 0.4, 1.5);
    coobject1->_dq.segment<3>(0) = coobject1->_T_world_robot.linear().transpose()* world_object1_tvel;
    coobject2->_dq.segment<3>(0) = coobject2->_T_world_robot.linear().transpose()* world_object2_tvel;
    coobject1->_dq.segment<3>(3) = coobject1->_T_world_robot.linear().transpose()* world_object1_rvel;
    coobject2->_dq.segment<3>(3) = coobject2->_T_world_robot.linear().transpose()* world_object2_rvel;

    // set up COP sim world
    auto cop_sim = new Sai2COPSim::COPSimulator(friction, restitution);
    //TODO: test with same primitive name
    cop_sim->addCapsuleObject(object1_name, object_link_name, "capsule1", coobject1, capsule1_radius, capsule1_length);
    cop_sim->addCapsuleObject(object2_name, object_link_name, "capsule2", coobject2, capsule2_radius, capsule2_length);

    // TODO: force sim/ display

    // ---- INITIALIZATION COMPLETE. STARTING MULTITHREADED APP ----
    // initialize GLFW window
    GLFWwindow* window = glfwInitialize();

    // set callbacks
    glfwSetKeyCallback(window, keySelect);

    // start the simulation thread first
    fSimulationRunning = true;
    thread sim_thread(simulation, cop_sim);

    // TODO: visualize contact

    // ---- GRAPHICS LOOP ON MAIN THREAD ----
    Eigen::Vector3d camera_pos, camera_lookat, camera_vertical;
    while (!glfwWindowShouldClose(window)) {
        // update model for object. called from simulation thread
        // coobject->updateModel();

        // update graphics. this automatically waits for the correct amount of time
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        graphics->updateGraphics(object1_name, coobject1);
        graphics->updateGraphics(object2_name, coobject2);
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

    // destroy context
    glfwDestroyWindow(window);

    // terminate
    glfwTerminate();

	return 0;
}

void simulation(Sai2COPSim::COPSimulator* sim) {
	fSimulationRunning = true;

    // create a timer
    LoopTimer timer;
    timer.initializeTimer();
    timer.setLoopFrequency(20000); //1500Hz timer
    double last_time = timer.elapsedTime(); //secs
    bool fTimerDidSleep = true;

    auto& geom_manager = sim->_geom_manager;
    auto prim_dist_ptr = geom_manager._prim_prim_distances[1][0]; //TODO: have a named interface in Geometry Manager
    auto capsule2_model = sim->_arb_manager.getBody(object2_name)->_model;
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

            if(timer.elapsedCycles() % 1000 == 0) {
                // cout << object->jtau_act(1) << endl;

                // visualize cop if present
            	if(prim_dist_ptr->contact_points.size() > 0) {
                    ptA_display.setLocalPos(prim_dist_ptr->contact_points[0]);
                    ptA_display.setShowEnabled(true);
                    ptB_display.setLocalPos(prim_dist_ptr->shadow_contact_points[0]);
                    ptB_display.setShowEnabled(true);
                }

                // // visualize external force
                min_dist_display.setShowEnabled(true);
                min_dist_display.m_pointA = chai3d::cVector3d(prim_dist_ptr->contact_points[0]);
                min_dist_display.m_pointB = chai3d::cVector3d(prim_dist_ptr->shadow_contact_points[0]);
            }
            
        } catch (exception& e) {
            cerr << e.what() << endl;
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