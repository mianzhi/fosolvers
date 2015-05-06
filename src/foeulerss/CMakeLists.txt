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

#ADD_CUSTOM_COMMAND(
#  TARGET foeulerss
#  POST_BUILD
#  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/runTest.sh
#                                   ${PROJECT_BINARY_DIR}/bin/runTest
#  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data
#                                             ${PROJECT_BINARY_DIR}/bin/data
#  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/sums.md5
#                                   ${PROJECT_BINARY_DIR}/bin/sums.md5
#)
