SET(tests_src
  test.f90
  otGrid1.f90
  otGrid2.f90
  polyX1.f90
  polyMesh1.f90
  polyGrid1.f90
  polyFeGrid1.f90
  readGTS1.f90
  readVTK1.f90
  writeVTK1.f90
  readCondTab1.f90
  UDF1.f90
  fixPt1.f90
  NewtonKrylov1.f90
  BDFNewtonKrylov1.f90
  ILU1.f90
  multiFront1.f90
  matVecMul.f90
  basicFEM1.f90
  edgeFEM1.f90
)

ADD_EXECUTABLE(tests ${tests_src})

TARGET_LINK_LIBRARIES(tests
  fosolvers
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${SUNDIALS_LIBRARIES}
  ${UMFPACK_LIBRARIES}
  ${MATHEVAL_LIBRARIES}
  ${METIS_LIBRARIES}
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
