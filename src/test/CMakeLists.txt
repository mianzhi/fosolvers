SET(test_src
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
)

ADD_EXECUTABLE(test ${test_src})

TARGET_LINK_LIBRARIES(test fosolvers)

ADD_CUSTOM_COMMAND(
  TARGET test
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/runTest.sh
                                   ${PROJECT_BINARY_DIR}/bin/runTest
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data
                                             ${PROJECT_BINARY_DIR}/bin/data
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/sums.md5
                                   ${PROJECT_BINARY_DIR}/bin/sums.md5
)
