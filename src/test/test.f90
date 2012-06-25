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
  use moduleFVMConvect
  use moduleInterpolation
  type(typeGrid)::grid
  double precision,allocatable::v(:),tempv(:),vv(:,:),u(:,:)
  double precision t
  
  open(12,file='bin/gridGMSH4.msh',status='old')
  call readGMSH(12,grid)
  close(12)
  call grid%updateIntf()
  allocate(v(grid%nBlock))
  allocate(tempv(grid%nBlock))
  allocate(vv(DIMS,grid%nBlock))
  call grid%updateBlockPos()
  do i=1,grid%nBlock
    v(i)=merge(100d0,000d0,grid%BlockPos(1,i)<0.1d0)
  end do
  allocate(u(DIMS,grid%nIntf))
  do i=1,grid%nIntf
    u(:,i)=[1d0,0d0,0d0]
  end do
  t=0d0
  open(13,file='rst.msh',status='replace')
  call writeGMSH(13,grid)
  call writeGMSH(13,v,grid,BIND_BLOCK,'name1',0,t)
  call grid%updateBlockVol()
  do i=1,100
    vv=findGrad(v,grid,BIND_BLOCK)
    tempv=v+5d-3*findConvect(v,u,BIND_INTF,grid,vv,limiter=minmod)/grid%BlockVol(:)
    do j=1,801,100
      tempv(j)=100d0
    end do
    do j=100,900,100
      tempv(j)=tempv(j)
    end do
    v(:)=tempv(:)
    t=t+5d-3
    if(mod(i,1)==0)then
      call writeGMSH(13,v,grid,BIND_BLOCK,'name1',i,t)
    end if
  end do
  close(13)
end subroutine
