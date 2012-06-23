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
  use moduleFVMDiffus
  use moduleInterpolation
  type(typeGrid)::grid
  double precision,allocatable::v(:),vv(:,:)
  double precision t
  
  open(12,file='bin/gridGMSH3.msh',status='old')
  call readGMSH(12,grid)
  close(12)
  allocate(v(grid%nBlock))
  allocate(vv(DIMS,grid%nBlock))
  call grid%updateBlockPos()
  do i=1,grid%nBlock
    v(i)=merge(1d0,0d0,norm2(grid%BlockPos(:,i))>0.8d0)
  end do
  t=0d0
  open(13,file='rst.msh',status='replace')
  call writeGMSH(13,grid)
  call writeGMSH(13,v,grid,BIND_BLOCK,'name1',0,t)
  do i=1,100
    vv=findGrad(v,grid,BIND_BLOCK)
    v=v+findDiffus(v,grid,vv)*0.8d0
    t=t+0.8d0
    call writeGMSH(13,v,grid,BIND_BLOCK,'name1',i,t)
  end do
  close(13)
end subroutine
