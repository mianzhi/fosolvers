SET(libfosolvers_src
  otGrid/otGrid.f90
  otGrid/otree.c
  polyGrid/polyX.f90
  polyGrid/polyGrid.f90
  polyGrid/polyFvGrid.f90
  polyGrid/polyFeGrid.f90
  polyGrid/polyMesh.f90
  condition/condition.f90
  UDF/UDF.f90
  fileIO/fileIO.f90
  fileIO/GTS.f90
  fileIO/VTK.f90
  fileIO/GMSH.f90
  fileIO/condTab.f90
  fileIO/UDFTab.f90
  geometry/geometry.f90
  FVM/gradient.f90
  FVM/advection.f90
  FVM/diffusion.f90
  FVM/Newtonian.f90
  FVM/pressure.f90
  FVM/Euler.f90
  FEM/SMQ.f90
  FEM/basicFEM.f90
  numerics/SUNDIALS.f90
  numerics/sparse.f90
  contrib/fmetis/metis_interface.f90
  contrib/SPARSKIT2/FORMATS/formats.f
  contrib/SPARSKIT2/FORMATS/unary.f
  contrib/SPARSKIT2/BLASSM/blassm.f
  contrib/SPARSKIT2/BLASSM/matvec.f
  contrib/SPARSKIT2/ITSOL/ilut.f
  misc/sort.f90
)

SET(libfosolvers_lib
  
)

ADD_LIBRARY(fosolvers STATIC ${libfosolvers_src})
TARGET_LINK_LIBRARIES(fosolvers ${libfosolvers_lib})
TARGET_INCLUDE_DIRECTORIES(fosolvers INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
SET_TARGET_PROPERTIES(fosolvers PROPERTIES Fortran_MODULE_DIRECTORY
                      ${CMAKE_CURRENT_BINARY_DIR}/include)
