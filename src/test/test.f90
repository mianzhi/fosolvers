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
  type(typeGrid1D)::grid
  double precision,allocatable::r(:),v(:)
  
  call initMPI()
  if(pidMPI==0)then
    call grid%genUniform(1d0,2d0,10)
    allocate(r(grid%nCell))
    allocate(v(grid%nCell))
    r(:)=1d0
    v(:)=0d0
    v(1:5)=1d0
    do l=1,1000
      v=v+0.4d-4/0.1d0*findDiffus(r,BIND_CELL,v,grid)
    end do
    write(*,*),v
  else
  end if
  call finalMPI()
end program
