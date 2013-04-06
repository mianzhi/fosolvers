SET(foflow_src foflow.f90 miscNS.f90 resFuncNS.f90)
ADD_EXECUTABLE(foflow ${foflow_src})
TARGET_LINK_LIBRARIES(foflow fosolvers)
ADD_CUSTOM_COMMAND(
  TARGET foflow
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/sim
                                   ${PROJECT_BINARY_DIR}/bin/sim
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/cond
                                   ${PROJECT_BINARY_DIR}/bin/cond
)
