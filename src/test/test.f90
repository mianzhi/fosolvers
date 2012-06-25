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
    v(i)=merge(100d0,100d0,grid%BlockPos(1,i)<0.5d0)
  end do
  allocate(u(DIMS,grid%nIntf))
  !do i=1,grid%nIntf
  !  u(:,i)=[1d0,0d0,0d0]
  !end do
  t=0d0
  open(13,file='rst.msh',status='replace')
  call writeGMSH(13,grid)
  call writeGMSH(13,v,grid,BIND_BLOCK,'name1',0,t)
  call grid%updateBlockVol()
  do i=1,10000
    vv=findGrad(v,grid,BIND_BLOCK)
    tempv=v+2d-5*findDiffus(v,grid,vv)/grid%BlockVol(:)
    do j=1,801,100
      tempv(j)=tempv(j)-2d0*v(j)/1d-4*2d-5
    end do
    do j=100,900,100
      tempv(j)=tempv(j)-(v(j)+dot_product(vv(:,j),[1d-2,0d0,0d0]))/1d-2*2d-5
    end do
    v(:)=tempv(:)
    t=t+2d-5
    if(mod(i,100)==0)then
      call writeGMSH(13,v,grid,BIND_BLOCK,'name1',i,t)
      write(*,*),i
    end if
  end do
  close(13)
end subroutine
