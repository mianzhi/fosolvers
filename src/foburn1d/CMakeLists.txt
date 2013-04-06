SET(foburn1d_src foburn1d.f90 miscBurn1D.f90)
ADD_EXECUTABLE(foburn1d ${foburn1d_src})
TARGET_LINK_LIBRARIES(foburn1d fosolvers)
