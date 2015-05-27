SET(foeulerss_src
  foeulerss.f90
)

ADD_EXECUTABLE(foeulerss ${foeulerss_src})

TARGET_LINK_LIBRARIES(foeulerss
  fosolvers
  cgns
  ${HDF5_LIBRARIES}
  lapack
)
