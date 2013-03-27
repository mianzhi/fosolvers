!----------------------------------------------------------------------------- best with 100 columns

program test
  use moduleBasicDataStruct
  use moduleSimpleSetLogic
  use moduleFileIO
  use moduleGrid
  use moduleGridInspection
  use moduleGridOperation
  use moduleGrid1D
  use moduleFVMGrad
  use moduleFVMDiffus
  use moduleFVMConvect
  use moduleInterpolation
  use moduleCondition
  use moduleNonlinearSolve
  use moduleMPIComm
  type(typeGrid)::grid
  type(typeGrid1D)::grid1d
  double precision,allocatable::u(:),uF(:),gradU(:,:),v(:),gamm(:)
  
  call initMPI()
  if(pidMPI==0)then
    call grid1d%genUniform(1d0,2d0,10)
    write(*,*),grid1d%nNode,grid1d%nCell
    write(*,*),grid1d%NodePos
    write(*,*),grid1d%CellPos
    open(12,file='bin/gridGMSH5.msh',status='old')
    call readGMSH(12,grid)
    close(12)
    allocate(u(grid%nBlock))
    allocate(uF(grid%nFacet))
    allocate(gradU(DIMS,grid%nBlock))
    allocate(v(grid%nBlock))
    allocate(gamm(grid%nBlock))
    call grid%updateBlockPos()
    u(:)=0d0
    uF(:)=1d0
    gradU=findGrad(u,BIND_BLOCK,grid,FacetVal=uF)
    gamm=1d0
    v=findDiffus(gamm,BIND_BLOCK,u,grid,gradU,FacetVal=uF)
    open(13,file='rstTest.msh',status='replace')
    call writeGMSH(13,grid)
    call writeGMSH(13,v,grid,BIND_BLOCK,'v',0,0d0)
    close(13)
  else
  end if
  call finalMPI()
end program
