SET(libfosolvers_src
  cond.f90
  fileio.f90
  fvm.f90
  grid.f90
  misc.f90
  mp1.f90
  mtl.f90
  simpgeo.f90
  time.f90
  utility.f90
)
ADD_LIBRARY(fosolvers STATIC ${libfosolvers_src})
TARGET_LINK_LIBRARIES(fosolvers lapack)

# subdirectories
ADD_SUBDIRECTORY(demo)
