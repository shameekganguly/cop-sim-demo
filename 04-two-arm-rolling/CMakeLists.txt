# create an executable
ADD_EXECUTABLE (04-two-arm-rolling main.cpp)

# and link the library against the executable
TARGET_LINK_LIBRARIES (04-two-arm-rolling ${COP_DEMO_COMMON_LIBRARIES})

# export resources such as model files.
# NOTE: this requires an install build
SET(RESOURCE_DIR ${COP_DEMO_BINARY_DIR}/resources)
FILE(MAKE_DIRECTORY ${RESOURCE_DIR})
SET(APP_RESOURCE_DIR ${RESOURCE_DIR}/04-two-arm-rolling)
FILE(MAKE_DIRECTORY ${APP_RESOURCE_DIR})
FILE(COPY world.urdf roller_object.urdf DESTINATION ${APP_RESOURCE_DIR})
