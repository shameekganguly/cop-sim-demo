# create an executable
ADD_EXECUTABLE (08-negative_space_peg_hole main.cpp)

# and link the library against the executable
TARGET_LINK_LIBRARIES (08-negative_space_peg_hole ${COP_DEMO_COMMON_LIBRARIES})

# export resources such as model files.
# NOTE: this requires an install build
SET(RESOURCE_DIR ${COP_DEMO_BINARY_DIR}/resources)
FILE(MAKE_DIRECTORY ${RESOURCE_DIR})
SET(APP_RESOURCE_DIR ${RESOURCE_DIR}/08-negative_space_peg_hole)
FILE(MAKE_DIRECTORY ${APP_RESOURCE_DIR})
FILE(COPY world.urdf capsule_object.urdf left_half.mtl left_half.obj right_half.mtl right_half.obj DESTINATION ${APP_RESOURCE_DIR})
