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
  double precision,allocatable::u(:),gradU(:),v(:)
  
  call initMPI()
  if(pidMPI==0)then
    call grid%genUniform(1d0,2d0,10)
    allocate(u(grid%nCell))
    allocate(gradU(grid%nCell))
    allocate(v(grid%nCell))
    v(:)=1d0
    u(:)=0d0
    u(1)=1d0
    do l=1,10
      gradU(:)=0d0
      gradU(2:grid%nCell-1)=(u(3:grid%nCell)-u(1:grid%nCell-2))/0.2d0
      gradU(1)=(u(2)-u(1))/0.1d0
      gradU(10)=(u(10)-u(9))/0.1d0
      u=u+0.05d0/0.1d0*findConvect(u,v,BIND_CELL,grid,gradU)
      u(1)=1d0
    end do
    write(*,*),u
  else
  end if
  call finalMPI()
end program
