SET(libfosolvers_src fileio.f90 grid.f90 misc.f90 time.f90 utility.f90)
ADD_LIBRARY(fosolvers STATIC ${libfosolvers_src})
TARGET_LINK_LIBRARIES(fosolvers lapack)

# subdirectories
ADD_SUBDIRECTORY(demo)
