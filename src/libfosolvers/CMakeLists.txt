SET(libfosolvers_src
  otGrid/otGrid.f90
  otGrid/otree.c
  polyX/polyX.f90
  polyGrid/polyGrid.f90
  polyMesh/polyMesh.f90
  fileIO/fileIO.f90
  fileIO/GTS.f90
  fileIO/VTK.f90
  geometry/geometry.f90
)

SET(libfosolvers_lib
  
)

ADD_LIBRARY(fosolvers STATIC ${libfosolvers_src})
TARGET_LINK_LIBRARIES(fosolvers ${libfosolvers_lib})
