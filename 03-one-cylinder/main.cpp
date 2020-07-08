/*  03-one-cylinder - main.cpp
This demo shows collision and contact resolution between a cylinder and a plane.

Author: Shameek Ganguly shameekg@stanford.edu
Date: 07/02/2020
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

#include <Eigen/Dense>
using namespace Eigen;

#include "timer/LoopTimer.h"
#include "force_sensor/ForceSensorSim.h"
#include "force_sensor/ForceSensorDisplay.h"

#include "cop_simulator/COPSimulator.h"

#include <GLFW/glfw3.h> //must be loaded after loading opengl/glew as part of graphicsinterface

using namespace chai3d;

const string world_fname = "resources/03-one-cylinder/world.urdf";
const string object_fname = "resources/03-one-cylinder/cylinder_object.urdf";
const string object1_name = "Cylinder1";
const double cylinder_radius = 0.1;
const double cylinder_length = 0.2;
const uint cylinder_num_points = 4;
const string object_link_name = "object";
const string box_name = "Box";
Affine3d cylinder1_in_world;

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

// COP visualization
chai3d::cShapeSphere pt_contact_display(0.02);
chai3d::cShapeLine ext_force_display;

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
	const double restitution = 0.3;
    const double friction = 0.1;

    // load sai2 simulation world
    // TODO: this is currently needed in order to load the base transform for each object
    // as well as the transform the base plane and the gravity vector
    // TODO: write a parser for world URDF to the cop sim world
    auto sai2_sim = new Simulation::Sai2Simulation(world_fname, false);

    cylinder1_in_world.translation() = sai2_sim->_world->getBaseNode(object1_name)->getLocalPos().eigen();
    cylinder1_in_world.linear() = sai2_sim->_world->getBaseNode(object1_name)->getLocalRot().eigen();
    // cout << cylinder1_in_world.translation().transpose() << endl;

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
    graphics->showWireMeshRender(true, object1_name, object_link_name);

    // load objects
    auto coobject1 = new Sai2Model::Sai2Model(object_fname, false, cylinder1_in_world, cylinder1_in_world.linear().transpose()*grav_vector);

    // add some initial velocity. TODO: move to parser
    coobject1->_dq[5] = 3.6;

    // set up COP sim world
    auto cop_sim = new Sai2COPSim::COPSimulator(friction, restitution);
    //TODO: test with same primitive name
    cop_sim->addCylinderObject(object1_name, object_link_name, "cylinder1", coobject1, cylinder_radius, cylinder_length, cylinder_num_points);
    // add plane
    cop_sim->addPlane(box_name, static_plane_in_world.linear().col(2), static_plane_in_world.translation());

    // TODO: force sim/ display

    // ---- INITIALIZATION COMPLETE. STARTING MULTITHREADED APP ----
    // initialize GLFW window
    GLFWwindow* window = glfwInitialize();

    // set callbacks
    glfwSetKeyCallback(window, keySelect);

    // start the simulation thread first
    fSimulationRunning = true;
    thread sim_thread(simulation, cop_sim);

    // visualize contact
    pt_contact_display.setShowEnabled(false);
    pt_contact_display.m_material->setBrownMaroon();
    pt_contact_display.m_material->setShininess(100);
    graphics->_world->addChild(&pt_contact_display);
    pt_contact_display.setLocalPos(Vector3d(0.0, 0.0, 0.4));

    ext_force_display.setShowEnabled(false);
    ext_force_display.m_colorPointA.setGreenYellowGreen();
    ext_force_display.m_colorPointB.setGreenYellowGreen();
    graphics->_world->addChild(&ext_force_display);
    ext_force_display.setLineWidth(4.0);

    // ---- GRAPHICS LOOP ON MAIN THREAD ----
    Eigen::Vector3d camera_pos, camera_lookat, camera_vertical;
    while (!glfwWindowShouldClose(window)) {
        // update model for object. called from simulation thread
        // coobject->updateModel();

        // update graphics. this automatically waits for the correct amount of time
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        graphics->updateGraphics(object1_name, coobject1);
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
    auto object = sim->_arb_manager.getBody(object1_name);
    auto cmodel = sim->_contact_model;
    auto& geom_manager = sim->_geom_manager;
    double sim_time_sum = 0.0;
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
            sim_time_sum += loop_dt;
            // object->jtau_act(1) = sim_time_sum*0.08;
            // object->jtau_act(3) = -sim_time_sum*0.08*0.1;
            if(timer.elapsedCycles() % 1000 == 0) {
                // cout << object->jtau_act(1) << endl;

                // visualize cop if present
                if(cmodel != NULL
                    && cmodel->_contact_island_models.size() > 0
                    && cmodel->_contact_island_models[0]._active_contacts.size() > 0
                ) {
                    const auto& island = cmodel->_contact_island_models[0];
                    const auto& prim_state = island._pair_state[island._active_contacts.front()];
                    if(prim_state.isValid()) {
                        pt_contact_display.setLocalPos(prim_state._cop_pos);
                        pt_contact_display.setShowEnabled(true);
                    }
                } else if(geom_manager._prim_prim_distances[1][0]->min_distance > 0.01) {
                    pt_contact_display.setShowEnabled(false);
                }

                // visualize external force
                ext_force_display.setShowEnabled(true);
                ext_force_display.m_pointA = chai3d::cVector3d(
                                                cylinder1_in_world*(object->_model->_q.segment<3>(0)) + Vector3d(0.0, 0.0, 0.1));
                ext_force_display.m_pointB = ext_force_display.m_pointA + chai3d::cVector3d(Vector3d(0.0, 1.0, 0.0)*object->jtau_act(1));
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
    GLFWwindow* window = glfwCreateWindow(windowW, windowH, "SAI2.0 - COP SIM 3", NULL, NULL);
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
