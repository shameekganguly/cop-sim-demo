# create an executable
ADD_EXECUTABLE (02-two-capsules main.cpp)

# and link the library against the executable
TARGET_LINK_LIBRARIES (02-two-capsules ${COP_DEMO_COMMON_LIBRARIES})

# export resources such as model files.
# NOTE: this requires an install build
SET(RESOURCE_DIR ${COP_DEMO_BINARY_DIR}/resources)
FILE(MAKE_DIRECTORY ${RESOURCE_DIR})
SET(APP_RESOURCE_DIR ${RESOURCE_DIR}/02-two-capsules)
FILE(MAKE_DIRECTORY ${APP_RESOURCE_DIR})
FILE(COPY world.urdf capsule_object.urdf DESTINATION ${APP_RESOURCE_DIR})
