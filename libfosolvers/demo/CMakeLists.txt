SET(demo_src demo.f90)
ADD_EXECUTABLE(demo ${demo_src})
TARGET_LINK_LIBRARIES(demo fosolvers)
