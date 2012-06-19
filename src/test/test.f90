!----------------------------------------------------------------------------- best with 100 columns

program tests
  call testReadGMSH()
end program

subroutine testReadGMSH()
  use moduleBasicDataStruct
  use moduleSimpleSetLogic
  use moduleFileIO
  use moduleGrid
  use moduleFVMGrad
  type(typeGrid)::grid
  double precision,allocatable::v(:),vv(:,:),vvv(:,:,:)
  
  open(12,file='bin/gridGMSH2.msh',status='old')
  call readGMSH(12,grid)
  close(12)
  allocate(v(grid%nBlock))
  allocate(vv(DIMS,grid%nBlock))
  allocate(vvv(DIMS,DIMS,grid%nBlock))
  call grid%updateBlockPos()
  do i=1,grid%nBlock
    v(i)=sin(3d0*grid%BlockPos(1,i))+sin(3d0*grid%BlockPos(2,i))+sin(3d0*grid%BlockPos(3,i))
  end do
  vv=findGrad(v,grid,bind=BIND_BLOCK)
  vvv=findGrad(vv,grid,bind=BIND_BLOCK)
  do i=1,grid%nBlock
    write(*,*),vvv(1,1,i),vvv(3,2,i)
  end do
  open(13,file='rst.msh',status='replace')
  call writeGMSH(13,grid)
  close(13)
end subroutine
