SET(foeuler_src
  foeuler.f90
)

ADD_EXECUTABLE(foeuler ${foeuler_src})

TARGET_LINK_LIBRARIES(foeuler
  fosolvers
  cgns
  ${HDF5_LIBRARIES}
  blas
  lapack
  ${CVODE_LIBRARIES}
  matheval
)

INSTALL(TARGETS foeuler DESTINATION bin)
