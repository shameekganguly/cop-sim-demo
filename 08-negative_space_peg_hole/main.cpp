/*  08-negative_space_peg_hole - main.cpp
This demo shows a robot inserting a capsule shaped peg into a
hole modeled as a negative capsule shape.

Author: Shameek Ganguly shameekg@alumni.stanford.edu
Date: 03/21/2022
*/

#include <iostream>
#include <string>
#include <thread>
#include <chrono>
#include <math.h>
#include <vector>
#include <execinfo.h>
#include <signal.h>
#include <unistd.h>
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

#include "cop_simulator/geometry/Composite1PkN.h"
#include "cop_simulator/COPSimulator.h"

#include <GLFW/glfw3.h> //must be loaded after loading opengl/glew as part of graphicsinterface

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

using namespace chai3d;

const string world_fname = "resources/08-negative_space_peg_hole/world.urdf";

const string object_fname = "resources/08-negative_space_peg_hole/capsule_object.urdf";
const string capsule_object_name = "Capsule";
const double capsule_radius = 0.01;
const double capsule_length = 0.08;
const string object_link_name = "object";

const double hole_capsule_length = 0.1;
const double hole_capsule_radius = 0.015;

const string block_name = "Block";

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
chai3d::cShapeSphere pt_contact_display1(0.015);
chai3d::cShapeSphere pt_contact_display2(0.015);

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
    signal(SIGSEGV, handler);
	const double restitution = 0;
    const double friction = 0.2;

    // load sai2 simulation world
    // TODO: this is currently needed in order to load the base transform for each object
    // as well as the transform the base plane and the gravity vector
    // TODO: write a parser for world URDF to the cop sim world
    auto sai2_sim = new Simulation::Sai2Simulation(world_fname, false);

    Affine3d capsule_in_world;
    capsule_in_world.translation() = sai2_sim->_world->getBaseNode(capsule_object_name)->getLocalPos().eigen();
    capsule_in_world.linear() = sai2_sim->_world->getBaseNode(capsule_object_name)->getLocalRot().eigen();

    Affine3d block_in_world;
    block_in_world.translation() = sai2_sim->_world->getBaseNode(block_name)->getLocalPos().eigen();
    block_in_world.linear() = sai2_sim->_world->getBaseNode(block_name)->getLocalRot().eigen();
    // cout << block_in_world.translation().transpose() << endl;

    Vector3d grav_vector = sai2_sim->_world->getGravity().eigen();

    // load graphics scene
    auto graphics = new Sai2Graphics::Sai2Graphics(world_fname, false);
    graphics->showWireMeshRender(true, capsule_object_name, object_link_name);

    // load objects
    // -- UNCOMMENT FOR CAPSULE DEMO --
    auto capsule_object = new Sai2Model::Sai2Model(object_fname,
                                                    false,
                                                    capsule_in_world,
                                                    capsule_in_world.linear().transpose()*grav_vector);

    // add some initial velocity. TODO: move to parser
    // coobject1->_dq[5] = 0.1;

    // set up COP sim world
    auto cop_sim = new Sai2COPSim::COPSimulator(friction, restitution);

    // add capsule object
    cop_sim->addObject(capsule_object_name, capsule_object);

    // add primitives to capsule object
    cop_sim->addCapsuleToObject(capsule_object_name, object_link_name, "capsule_mid", capsule_radius, capsule_length, Affine3d::Identity());

    // add block
    auto* composite = cop_sim->addPlaneComposite1PkN(block_name, Vector3d::UnitZ(), Vector3d::Zero());
    composite->_transform_in_link = block_in_world;
    {
        auto* negCap1 = new Sai2COPSim::NegCapsulePrimitive("negCap1", hole_capsule_radius, hole_capsule_length);
        Affine3d negCap1Tf = Eigen::Affine3d::Identity();
        negCap1Tf.linear() = Matrix3d(AngleAxisd(-M_PI/2,  Vector3d::UnitY()));
        composite->addNegativePrimitive(negCap1, negCap1Tf);
    }

    // Allow point-based resolution for CONCAVE contacts
    cop_sim->_contact_model->setSupportPtContactSteadyContact(true);

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
    // pt_contact_display1.setShowEnabled(false);
    // pt_contact_display1.m_material->setBrownMaroon();
    // pt_contact_display1.m_material->setShininess(100);
    // graphics->_world->addChild(&pt_contact_display1);
    // pt_contact_display1.setLocalPos(Vector3d(0.0, 0.0, 0.4));
    // pt_contact_display2.setShowEnabled(false);
    // pt_contact_display2.m_material->setBrownMaroon();
    // pt_contact_display2.m_material->setShininess(100);
    // graphics->_world->addChild(&pt_contact_display2);
    // pt_contact_display2.setLocalPos(Vector3d(0.0, 0.0, 0.4));

    // ---- INITIALIZATION COMPLETE. STARTING MULTITHREADED APP ----
    // initialize GLFW window
    GLFWwindow* window = glfwInitialize();

    // set callbacks
    glfwSetKeyCallback(window, keySelect);

    // start the simulation thread first
    fSimulationRunning = true;
    thread sim_thread(simulation, cop_sim);
    sim_state_label->setText(f_pause_sim? "PAUSED": "Sim running");

    // TODO: visualize contact

    // ---- GRAPHICS LOOP ON MAIN THREAD ----
    Eigen::Vector3d camera_pos, camera_lookat, camera_vertical;
    while (!glfwWindowShouldClose(window)) {
        // update graphics. this automatically waits for the correct amount of time
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);

        // -- UNCOMMENT FOR CAPSULE DEMO --
        graphics->updateGraphics(capsule_object_name, capsule_object);
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

    // ---- Print additional info from sim ----
    Vector3d capsule_pos;
    capsule_object->positionInWorld(capsule_pos, object_link_name);
    cout << capsule_pos.transpose() << endl;

	return 0;
}

void simulation(Sai2COPSim::COPSimulator* sim) {
	fSimulationRunning = true;

    // auto object = sim->_arb_manager.getBody(capsule_object_name);
    // auto cmodel = sim->_contact_model;
    // auto& geom_manager = sim->_geom_manager;

    // create a timer
    constexpr uint kFrequency = 1000;
    LoopTimer timer;
    timer.initializeTimer();
    timer.setLoopFrequency(kFrequency); //1500Hz timer
    double last_time = timer.elapsedTime(); //secs
    bool fTimerDidSleep = true;

    constexpr uint kSimulationTime = 10;

    while (fSimulationRunning && sim->_iterations < kFrequency * kSimulationTime) {
        fTimerDidSleep = timer.waitForNextLoop();

        // integrate forward
        double curr_time = timer.elapsedTime();
        if(f_pause_sim) {
            last_time = curr_time;
            continue;
        }
        // double loop_dt = curr_time - last_time;
        constexpr double loop_dt = 1.0/kFrequency;

        try {
            sim->integrate(loop_dt);

            // if(timer.elapsedCycles() % 500 == 0) {
            //     // cout << object->jtau_act(1) << endl;

            //     // visualize cop if present
            //     if(cmodel != NULL
            //         && cmodel->_contact_island_models_size > 0
            //         && cmodel->_contact_island_models[0]._active_contacts.size() > 1
            //     ) {
            //         const auto& island = cmodel->_contact_island_models[0];
            //         const auto& prim_state1 = island._pair_state[0];
            //         if(prim_state1.isValid()) {
            //             pt_contact_display1.setLocalPos(prim_state1._cop_pos);
            //             pt_contact_display1.setShowEnabled(true);
            //         }
            //         const auto& prim_state2 = island._pair_state[1];
            //         if(prim_state2.isValid()) {
            //             pt_contact_display2.setLocalPos(prim_state2._cop_pos);
            //             pt_contact_display2.setShowEnabled(true);
            //         }
            //     } else if(geom_manager._prim_prim_distances[4][1]->min_distance > 0.01) {
            //         pt_contact_display1.setShowEnabled(false);
            //         pt_contact_display2.setShowEnabled(false);
            //     }
            // }

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
    GLFWwindow* window = glfwCreateWindow(windowW, windowH, "SAI2.0 - Negative space peg in hole insertion", NULL, NULL);
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
