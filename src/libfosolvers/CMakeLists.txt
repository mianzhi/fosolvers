SET(libfosolvers_src
  BasicDataStruct/moduleBasicDataStruct.f90
)
ADD_LIBRARY(fosolvers STATIC ${libfosolvers_src})
TARGET_LINK_LIBRARIES(fosolvers lapack)
