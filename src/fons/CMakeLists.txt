SET(fons_src fons.f90 miscNS.f90 resFuncNS.f90)
ADD_EXECUTABLE(fons ${fons_src})
TARGET_LINK_LIBRARIES(fons fosolvers)
ADD_CUSTOM_COMMAND(
  TARGET fons
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/sim
                                   ${PROJECT_BINARY_DIR}/bin/sim
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/cond
                                   ${PROJECT_BINARY_DIR}/bin/cond
)
