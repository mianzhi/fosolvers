SET(tests_src tests.f90 showTestRst.f90)
ADD_EXECUTABLE(tests ${tests_src})
TARGET_LINK_LIBRARIES(tests fosolvers)

SET(test_src test.f90 showTestRst.f90)
ADD_EXECUTABLE(test ${test_src})
TARGET_LINK_LIBRARIES(test fosolvers)
ADD_CUSTOM_COMMAND(
  TARGET test
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/gridGMSH1.msh
                                   ${PROJECT_BINARY_DIR}/bin/gridGMSH1.msh
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/gridGMSH2.msh
                                   ${PROJECT_BINARY_DIR}/bin/gridGMSH2.msh
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/gridGMSH3.msh
                                   ${PROJECT_BINARY_DIR}/bin/gridGMSH3.msh
)
