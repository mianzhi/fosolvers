!----------------------------------------------------------------------------- best with 100 columns

program test
  use moduleBasicDataStruct
  use moduleSimpleSetLogic
  use moduleFileIO
  use moduleGrid
  use moduleGridInspection
  use moduleGridOperation
  use moduleFVMGrad
  use moduleFVMDiffus
  use moduleFVMConvect
  use moduleInterpolation
  use moduleCondition
  use moduleNonlinearSolve
  use moduleMPIComm
  type(typeGrid)::grid
  double precision,allocatable::u(:),uF(:),gradU(:,:)
  
  call initMPI()
  if(pidMPI==0)then
    open(12,file='bin/gridGMSH5.msh',status='old')
    call readGMSH(12,grid)
    close(12)
    allocate(u(grid%nBlock))
    allocate(uF(grid%nFacet))
    allocate(gradU(DIMS,grid%nBlock))
    call grid%updateBlockPos()
    u(:)=grid%BlockPos(1,:)
    uF(:)=0d0
    gradU=findGrad(u,BIND_BLOCK,grid,FacetVal=uF)
    open(13,file='rstTest.msh',status='replace')
    call writeGMSH(13,grid)
    call writeGMSH(13,gradU,grid,BIND_BLOCK,'u',0,0d0)
    close(13)
  else
  end if
  call finalMPI()
end program
