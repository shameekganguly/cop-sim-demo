
cmake_minimum_required (VERSION 2.8.12)
project (COP_DEMO)
SET(CMAKE_CXX_FLAGS "-std=c++14")

include_directories(
	/usr/local/include
  /usr/local/include/gte  # required for gte Mathematics
	${PROJECT_SOURCE_DIR}/src
)

# - CHAI3D
find_package(CHAI3D REQUIRED)
include_directories(${CHAI3D_INCLUDE_DIRS})
add_definitions(${CHAI3D_DEFINITIONS})
if (CMAKE_SYSTEM_NAME MATCHES Darwin)
  add_definitions(-D_MACOSX -DMACOSX)
elseif (CMAKE_SYSTEM_NAME MATCHES Linux)
  add_definitions(-D_LINUX -DLINUX)
endif () #NOTE: Windows is currently unsupported
link_directories(${CHAI3D_LIBRARY_DIRS})

# - SAI2-SIMULATION. TODO: this should not be needed. something fishy.
find_package(SAI2-SIMULATION REQUIRED)
include_directories(${SAI2-SIMULATION_INCLUDE_DIRS})
add_definitions(${SAI2-SIMULATION_DEFINITIONS})

# - SAI2-URDF. TODO: this should not be needed. something fishy.
find_package(SAI2-URDF REQUIRED)
include_directories(${SAI2-URDF_INCLUDE_DIRS})
add_definitions(${SAI2-URDF_DEFINITIONS})

# - SAI2-MODEL. TODO: this should not be needed. something fishy.
find_package(SAI2-MODEL REQUIRED)
include_directories(${SAI2-MODEL_INCLUDE_DIRS})
add_definitions(${SAI2-MODEL_DEFINITIONS})

# - SAI2-GRAPHICS. TODO: this should not be needed. something fishy.
find_package(SAI2-GRAPHICS REQUIRED)
include_directories(${SAI2-GRAPHICS_INCLUDE_DIRS})
add_definitions(${SAI2-GRAPHICS_DEFINITIONS})

# - SAI2-PRIMITIVES
find_package(SAI2-PRIMITIVES REQUIRED)
include_directories(${SAI2-PRIMITIVES_INCLUDE_DIRS})
add_definitions(${SAI2-PRIMITIVES_DEFINITIONS})

# - SAI2-COMMON
find_package(SAI2-COMMON REQUIRED)
include_directories(${SAI2-COMMON_INCLUDE_DIRS})
add_definitions(${SAI2-COMMON_DEFINITIONS})

# - glfw3
find_package(glfw3 QUIET)
find_library(GLFW_LIBRARY glfw)

# - rbdl
find_package(RBDL QUIET)
include_directories(${RBDL_INCLUDE_DIR})

set(COP_THIRDPARTY_LIBRARIES
	${CHAI3D_LIBRARIES}
    ${SAI2-URDF_LIBRARIES}
    ${SAI2-SIMULATION_LIBRARIES}
    ${SAI2-MODEL_LIBRARIES}
    ${SAI2-GRAPHICS_LIBRARIES}
    ${SAI2-PRIMITIVES_LIBRARIES}
    ${SAI2-COMMON_LIBRARIES}
	${GLFW_LIBRARY}
)

set(COP_DEMO_COMMON_SOURCE
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/geometry/ContactPatch.cpp
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/geometry/Primitive.cpp
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/geometry/PrimNormal.cpp
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/geometry/PrimPrimContactInfo.cpp
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/geometry/PointPrimDistance.cpp
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/geometry/PrimPrimDistance.cpp
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/geometry/IntersectionEdge.cpp
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/geometry/IntersectionEdgeFactory.cpp
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/geometry/IntersectionEdgePrimDistance.cpp
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/geometry/NegPrimPrimDistance.cpp
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/geometry/GeometryUtils.cpp
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/geometry/Composite1PkN.cpp
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/ArticulatedRigidBody.cpp
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/ArticulatedRigidBodyManager.cpp
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/WorldContactMap.cpp
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/COPSolverExtended.cpp
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/ContactSpaceModel.cpp
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/CollisionResolution.cpp
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/SteadyContactResolution.cpp
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/PointContactsSteadyContactResolution.cpp
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/ContactCurvature.cpp
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/ContactDetection.cpp
  ${PROJECT_SOURCE_DIR}/src/cop_simulator/COPSimulator.cpp
  ${PROJECT_SOURCE_DIR}/src/lcp_solvers/LCPSolver.cpp
  ${PROJECT_SOURCE_DIR}/src/lcp_solvers/LCPSolverInternal.cpp
)

add_library(copsim ${COP_DEMO_COMMON_SOURCE})

set(COP_DEMO_COMMON_LIBRARIES
  copsim
  ${COP_THIRDPARTY_LIBRARIES}
)

# add apps
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY                ${PROJECT_SOURCE_DIR}/bin)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG          ${PROJECT_SOURCE_DIR}/bin)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL     ${PROJECT_SOURCE_DIR}/bin)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE        ${PROJECT_SOURCE_DIR}/bin)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${PROJECT_SOURCE_DIR}/bin)
set (COP_DEMO_BINARY_DIR ${PROJECT_SOURCE_DIR}/bin)

# add_subdirectory(01-single-contact)
# add_subdirectory(02-two-capsules)
# add_subdirectory(03-one-cylinder)
# add_subdirectory(04-two-arm-rolling)
# add_subdirectory(05-one-box)
# add_subdirectory(06-one-pyramid)
# add_subdirectory(07-sphere-curved)
# add_subdirectory(08-negative_space_peg_hole)
# add_subdirectory(09-neg_space_arm_insertion)
add_subdirectory(10-redundant-surface-contacts)

# tests
# add_subdirectory(src/lcp_solvers)
# add_subdirectory(src/cop_simulator_tests)
