# create an executable
ADD_EXECUTABLE (testCOPGeometry COPGeometryTest.cpp 
	${PROJECT_SOURCE_DIR}/src/cop_simulator/geometry/ContactPatch.cpp
	${PROJECT_SOURCE_DIR}/src/cop_simulator/geometry/Primitive.cpp
	${PROJECT_SOURCE_DIR}/src/cop_simulator/geometry/PrimPrimDistance.cpp
	${PROJECT_SOURCE_DIR}/src/cop_simulator/geometry/GeometryUtils.cpp
)

ADD_EXECUTABLE (testCOPSolver COPSolverTest.cpp
	${PROJECT_SOURCE_DIR}/src/cop_simulator/geometry/ContactPatch.cpp
	${PROJECT_SOURCE_DIR}/src/cop_simulator/geometry/Primitive.cpp
	${PROJECT_SOURCE_DIR}/src/cop_simulator/COPSolverExtended.cpp
)

# and link the library against the executable
TARGET_LINK_LIBRARIES (testCOPGeometry ${COP_DEMO_COMMON_LIBRARIES})
TARGET_LINK_LIBRARIES (testCOPSolver ${COP_DEMO_COMMON_LIBRARIES})