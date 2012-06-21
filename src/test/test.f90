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
  use moduleInterpolation
  type(typeGrid)::grid
  double precision,allocatable::v(:),vv(:,:),vvv(:,:,:)
  
  open(12,file='bin/gridGMSH1.msh',status='old')
  call readGMSH(12,grid)
  close(12)
  allocate(v(grid%nNode))
  allocate(vv(DIMS,grid%nNode))
  allocate(vvv(DIMS,DIMS,grid%nNode))
  do i=1,grid%nNode
    v(i)=sin(10d0*grid%NodePos(1,i))+sin(10d0*grid%NodePos(2,i))+sin(10d0*grid%NodePos(3,i))
  end do
  vv=findGrad(v,grid,BIND_NODE)
  vvv=findGrad(vv,grid,BIND_NODE)
  open(13,file='rst.msh',status='replace')
  call writeGMSH(13,grid)
  call writeGMSH(13,v,grid,BIND_NODE,'name1')
  call writeGMSH(13,vv,grid,BIND_NODE,'name2')
  call writeGMSH(13,vvv,grid,BIND_NODE,'name3')
  close(13)
end subroutine
