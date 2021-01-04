/*  02-two-capsules - main.cpp
This demo shows collision and contact resolution between two capsules and a plane.

Author: Shameek Ganguly shameekg@stanford.edu
Date: 06/10/2020
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

const string world_fname = "resources/02-two-capsules/world.urdf";
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

    Affine3d capsule3_in_world;
    capsule3_in_world.translation() = sai2_sim->_world->getBaseNode(object3_name)->getLocalPos().eigen();
    capsule3_in_world.linear() = sai2_sim->_world->getBaseNode(object3_name)->getLocalRot().eigen();

    Affine3d capsule4_in_world;
    capsule4_in_world.translation() = sai2_sim->_world->getBaseNode(object4_name)->getLocalPos().eigen();
    capsule4_in_world.linear() = sai2_sim->_world->getBaseNode(object4_name)->getLocalRot().eigen();

    Affine3d capsule5_in_world;
    capsule5_in_world.translation() = sai2_sim->_world->getBaseNode(object5_name)->getLocalPos().eigen();
    capsule5_in_world.linear() = sai2_sim->_world->getBaseNode(object5_name)->getLocalRot().eigen();

    Affine3d static_plane_in_world;
    static_plane_in_world.translation() = sai2_sim->_world->getBaseNode(box_name)->getLocalPos().eigen();
    static_plane_in_world.linear() = sai2_sim->_world->getBaseNode(box_name)->getLocalRot().eigen();
    static_plane_in_world.translation() += static_plane_in_world.linear().col(2)*0.1/2.0; // z-height of box
    // cout << static_plane_in_world.translation().transpose() << endl;


    Vector3d grav_vector = sai2_sim->_world->getGravity().eigen();

    // load graphics scene
    auto graphics = new Sai2Graphics::Sai2Graphics(world_fname, false);
    // set object graphics to wireframe, and show frame for last link
    // graphics->showLinkFrame(true, object1_name, object_link_name);
    // graphics->showLinkFrame(true, object2_name, object_link_name);
    // graphics->showLinkFrame(true, object3_name, object_link_name);
    // graphics->showLinkFrame(true, object4_name, object_link_name);
    // graphics->showLinkFrame(true, object5_name, object_link_name);

    graphics->showWireMeshRender(true, object1_name, object_link_name);
    graphics->showWireMeshRender(true, object2_name, object_link_name);
    graphics->showWireMeshRender(true, object3_name, object_link_name);
    graphics->showWireMeshRender(true, object4_name, object_link_name);
    graphics->showWireMeshRender(true, object5_name, object_link_name);

    // load objects
    auto coobject1 = new Sai2Model::Sai2Model(object_fname, false, capsule1_in_world, capsule1_in_world.linear().transpose()*grav_vector);
    auto coobject2 = new Sai2Model::Sai2Model(object_fname, false, capsule2_in_world, capsule2_in_world.linear().transpose()*grav_vector);
    auto coobject3 = new Sai2Model::Sai2Model(object_fname, false, capsule3_in_world, capsule3_in_world.linear().transpose()*grav_vector);
    auto coobject4 = new Sai2Model::Sai2Model(object_fname, false, capsule4_in_world, capsule4_in_world.linear().transpose()*grav_vector);
    auto coobject5 = new Sai2Model::Sai2Model(object_fname, false, capsule5_in_world, capsule5_in_world.linear().transpose()*grav_vector);

    // add some initial velocity. TODO: move to parser
    // coobject1->_dq[5] = 0.1;

    // set up COP sim world
    auto cop_sim = new Sai2COPSim::COPSimulator(friction, restitution);
    //TODO: test with same primitive name
    cop_sim->addCapsuleObject(object1_name, object_link_name, "capsule1", coobject1, capsule1_radius, capsule1_length);
    cop_sim->addCapsuleObject(object2_name, object_link_name, "capsule2", coobject2, capsule2_radius, capsule2_length);
    cop_sim->addCapsuleObject(object3_name, object_link_name, "capsule3", coobject3, capsule3_radius, capsule3_length);
    cop_sim->addCapsuleObject(object4_name, object_link_name, "capsule4", coobject4, capsule4_radius, capsule4_length);
    cop_sim->addCapsuleObject(object5_name, object_link_name, "capsule5", coobject5, capsule5_radius, capsule5_length);
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
        graphics->updateGraphics(object3_name, coobject3);
        graphics->updateGraphics(object4_name, coobject4);
        graphics->updateGraphics(object5_name, coobject5);
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
            // if(timer.elapsedCycles() % 100 == 0) {
            //     cout << "Update sim" << endl;
            //     cout << sim->_contact_map._islands.size() << endl;
            // }
            // cout << (sim->_contact_map._islands[0]._contact_prim_pairs.front()).primA->_name << endl;
            // cout << (sim->_contact_map._islands[0]._contact_prim_pairs.front()).primB->_name << endl;
            // if(sim->_contact_map._islands.size() > 0) {
                // cout << sim->_contact_model->_contact_island_models.size() << endl;
                // print contact model
                // for(const auto& c_island: sim->_contact_model->_contact_island_models) {
                //     cout << "Island bodies: ";
                //     for(const auto& arb_name: c_island._geom_island->_articulated_bodies) {
                //         cout << arb_name << ", ";
                //     }
                //     cout << endl;
                //     cout << "ARB index map:";
                //     for(auto it: c_island._arb_index_map) {
                //         cout << it.first << ": " << it.second << ", ";
                //     }
                //     cout << endl;
                    // cout << "Num prim pairs: " << c_island._pair_state.size() << endl;
                    // cout << "Num active contacts: " << c_island._active_contacts.size() << endl;
                    // cout << "Pt ct Jacobian: " << endl;
                    // cout << c_island._pt_contact_Jacobian << endl;
                    // cout << "COP Full6 Jacobian: " << endl;
                    // cout << c_island._cop_full6_Jacobian << endl;
                    // cout << "COP Constraint Jacobian: " << endl;
                    // cout << c_island._cop_constraint_Jacobian << endl;
                    // cout << "Pt ct Lambda inv: " << endl;
                    // cout << c_island._pt_contact_Lambda_inv << endl;
                    // cout << "COP Full6 Lambda inv: " << endl;
                    // cout << c_island._cop_full6_Lambda_inv << endl;
                    // cout << "COP Constraint Lambda inv: " << endl;
                    // cout << c_island._cop_constraint_Lambda_inv << endl;
                    // cout << "Pt ct RHS: " << endl;
                    // cout << c_island._pt_contact_rhs_coll.transpose() << endl;
                    // cout << "COP Full6 RHS: " << endl;
                    // cout << c_island._cop_full6_rhs_contact.transpose() << endl;
                    // cout << "COP Constraint RHS: " << endl;
                    // cout << c_island._cop_constraint_rhs_contact.transpose() << endl;
                    // cout << "COP full active matrix" << endl;
                    // Eigen::MatrixXd J_full_cop, Lambda_inv_full_cop;
                    // Eigen::VectorXd rhs_full_cop;
                    // std::vector<uint> Jrow_ind_to_contact_pair_map;
                    // c_island.getActiveFullCOPMatrices(J_full_cop, Lambda_inv_full_cop, rhs_full_cop, Jrow_ind_to_contact_pair_map);
                    // cout << "J" << endl;
                    // cout << J_full_cop << endl;
                    // cout << "L inv" << endl;
                    // cout << Lambda_inv_full_cop << endl;
                    // cout << "rhs" << endl;
                    // cout << rhs_full_cop.transpose() << endl;
                    // cout << "COP constraint active matrix" << endl;
                    // Eigen::MatrixXd J_cons_cop, Lambda_inv_cons_cop;
                    // Eigen::VectorXd rhs_cons_cop;
                    // std::vector<uint> Jrow_ind_to_contact_pair_map;
                    // c_island.getActiveConstraintCOPMatrices(J_cons_cop, Lambda_inv_cons_cop, rhs_cons_cop, Jrow_ind_to_contact_pair_map);
                    // cout << "J" << endl;
                    // cout << J_cons_cop << endl;
                    // cout << "L inv" << endl;
                    // cout << Lambda_inv_cons_cop << endl;
                    // cout << "rhs" << endl;
                    // cout << rhs_cons_cop.transpose() << endl;
                    // cout << "Pt contact active matrix" << endl;
                    // Eigen::MatrixXd J_pt_coll, Lambda_inv_pt_coll;
                    // Eigen::VectorXd rhs_pt_coll;
                    // std::vector<uint> Jrow_ind_to_contact_pair_map;
                    // c_island.getActivePtContactCollisionMatrices(J_pt_coll, Lambda_inv_pt_coll, rhs_pt_coll, Jrow_ind_to_contact_pair_map);
                    // cout << "J" << endl;
                    // cout << J_pt_coll << endl;
                    // cout << "L inv" << endl;
                    // cout << Lambda_inv_pt_coll << endl;
                    // cout << "rhs" << endl;
                    // cout << rhs_pt_coll.transpose() << endl;


                    // for(const auto& pair: c_island._pair_state) {
                    //     cout << pair._geom_prim_pair->primA->_name << endl;
                    //     cout << pair._geom_prim_pair->primB->_name << endl;
                    //     cout << pair._geom_prim_pair->info.contact_points[0].transpose() << endl;
                    // }
            //     }
            // }
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
