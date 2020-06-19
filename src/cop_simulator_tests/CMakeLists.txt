# create an executable
ADD_EXECUTABLE (testCOPGeometry COPGeometryTest.cpp ${COP_DEMO_COMMON_SOURCE})

# and link the library against the executable
TARGET_LINK_LIBRARIES (testCOPGeometry ${COP_DEMO_COMMON_LIBRARIES})