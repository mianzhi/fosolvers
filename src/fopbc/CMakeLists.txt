SET(fopbc_src
  fopbc.f90
  modPbc.f90
)

ADD_EXECUTABLE(fopbc ${fopbc_src})

TARGET_LINK_LIBRARIES(fopbc
  fosolvers
  cgns
  ${HDF5_LIBRARIES}
  blas
  lapack
  ${CVODE_LIBRARIES}
  matheval
)

INSTALL(TARGETS fopbc DESTINATION bin)
