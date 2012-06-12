SET(libfosolvers_src
  BasicDataStruct/moduleBasicDataStruct.f90
  SimpleGeometry/moduleSimpleGeometry.f90
  SimpleSetLogic/moduleSimpleSetLogic.f90
  Grid/moduleGrid.f90
  GridOperation/moduleGridOperation.f90
  FileIO/moduleFileIO.f90
)
ADD_LIBRARY(fosolvers STATIC ${libfosolvers_src})
TARGET_LINK_LIBRARIES(fosolvers lapack)
