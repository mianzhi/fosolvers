SET(foeuler_src
  foeuler.f90
)

ADD_EXECUTABLE(foeuler ${foeuler_src})

TARGET_LINK_LIBRARIES(foeuler
  fosolvers
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${SUNDIALS_LIBRARIES}
  ${UMFPACK_LIBRARIES}
  ${MATHEVAL_LIBRARIES}
  ${METIS_LIBRARIES}
)

INSTALL(TARGETS foeuler DESTINATION bin)
