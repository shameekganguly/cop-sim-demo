# create an executable
ADD_EXECUTABLE (testCOPGeometry COPGeometryTest.cpp)

ADD_EXECUTABLE (testCOPSolver COPSolverTest.cpp)

# and link the library against the executable
TARGET_LINK_LIBRARIES (testCOPGeometry ${COP_DEMO_COMMON_LIBRARIES})
TARGET_LINK_LIBRARIES (testCOPSolver ${COP_DEMO_COMMON_LIBRARIES})