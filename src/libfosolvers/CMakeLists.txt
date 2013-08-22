SET(libfosolvers_src
  BasicDataStruct/moduleBasicDataStruct.f90
  SimpleGeometry/moduleSimpleGeometry.f90
  SimpleSetLogic/moduleSimpleSetLogic.f90
  Grid/moduleGrid.f90
  Grid/moduleGridInspection.f90
  Grid/moduleGridOperation.f90
  Grid1D/moduleGrid1D.f90
  FileIO/moduleFileIO.f90
  Interpolation/moduleInterpolation.f90
  FVM/moduleFVMGrad.f90
  FVM/moduleFVMDiffus.f90
  FVM/moduleFVMConvect.f90
  Condition/moduleCondition.f90
  AlgebraicSolver/moduleNonlinearSolve.f90
  MultiProcessing/moduleMPIComm.f90
  SpatialHashing/moduleSpatialHashing.f90
  Misc/moduleCLIO.f90
)
SET(libfosolvers_lib
  sundials_fkinsol
  sundials_fnvecserial
  sundials_kinsol
  sundials_nvecserial
  lapack
)
ADD_LIBRARY(fosolvers STATIC ${libfosolvers_src})
TARGET_LINK_LIBRARIES(fosolvers ${libfosolvers_lib})
