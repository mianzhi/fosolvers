SET(fons_src fons.f90)
ADD_EXECUTABLE(fons ${fons_src})
TARGET_LINK_LIBRARIES(fons fosolvers)
