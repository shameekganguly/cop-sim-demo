# create an executable
ADD_EXECUTABLE (testLCPSolver LCPSolverTest.cpp LCPSolver2.cpp)

# and link the library against the executable
TARGET_LINK_LIBRARIES (testLCPSolver ${COP_DEMO_COMMON_LIBRARIES})