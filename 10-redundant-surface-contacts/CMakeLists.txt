# create an executable
ADD_EXECUTABLE (10-redundant-surface-contacts main.cpp)

# and link the library against the executable
TARGET_LINK_LIBRARIES (10-redundant-surface-contacts ${COP_DEMO_COMMON_LIBRARIES})

# export resources such as model files.
# NOTE: this requires an install build
SET(RESOURCE_DIR ${COP_DEMO_BINARY_DIR}/resources)
FILE(MAKE_DIRECTORY ${RESOURCE_DIR})
SET(APP_RESOURCE_DIR ${RESOURCE_DIR}/10-redundant-surface-contacts)
FILE(MAKE_DIRECTORY ${APP_RESOURCE_DIR})
FILE(COPY world.urdf box_object.urdf DESTINATION ${APP_RESOURCE_DIR})
