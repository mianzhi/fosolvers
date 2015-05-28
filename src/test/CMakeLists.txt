SET(tests_src
  test.f90
  otGrid1.f90
  otGrid2.f90
  polyX1.f90
  polyMesh1.f90
  polyGrid1.f90
  polyFvGrid1.f90
  readGTS1.f90
  readVTK1.f90
  writeVTK1.f90
  readCondTab1.f90
  gradient1.f90
  advection1.f90
  advection2.f90
  euler.f90
)

ADD_EXECUTABLE(tests ${tests_src})

TARGET_LINK_LIBRARIES(tests
  fosolvers
  cgns
  ${HDF5_LIBRARIES}
  lapack
)

ADD_CUSTOM_COMMAND(
  TARGET tests
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/runTest.sh
                                   ${PROJECT_BINARY_DIR}/bin/runTest
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data
                                             ${PROJECT_BINARY_DIR}/bin/data
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/sums.md5
                                   ${PROJECT_BINARY_DIR}/bin/sums.md5
)
