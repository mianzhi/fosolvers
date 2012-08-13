SET(foeuler_src foeuler.f90)
ADD_EXECUTABLE(foeuler ${foeuler_src})
TARGET_LINK_LIBRARIES(foeuler fosolvers)
