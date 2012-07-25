!----------------------------------------------------------------------------- best with 100 columns

program test
  use moduleBasicDataStruct
  use moduleSimpleSetLogic
  use moduleFileIO
  use moduleGrid
  use moduleGridOperation
  use moduleFVMGrad
  use moduleFVMDiffus
  use moduleFVMConvect
  use moduleInterpolation
  use moduleMPIComm
  type(typeGrid)::grid
  double precision,allocatable::v(:),gradv(:,:),disp(:,:),temp(:)
  double precision t
  
  call initMPI()
  if(pidMPI==0)then
    open(12,file='bin/gridGMSH5.msh',status='old')
    call readGMSH(12,grid)
    close(12)
    call grid%updateBlockPos()
    allocate(disp(DIMS,grid%nNode))
    allocate(v(grid%nBlock))
    allocate(temp(grid%nBlock))
    allocate(gradv(DIMS,grid%nBlock))
    disp(:,:)=0d0
    forall(i=1:grid%nNode)
      disp(1,i)=0.05d0-0.04d0*grid%NodePos(1,i)
    end forall
    do i=1,grid%nBlock
      v(i)=merge(100d0,000d0,grid%BlockPos(1,i)<0.9d0.and.grid%BlockPos(1,i)>0.7d0)
    end do
    t=0d0
    open(13,file='rst.msh',status='replace')
    call writeGMSH(13,grid)
    call writeGMSH(13,v,grid,BIND_BLOCK,'name1',0,t)
    call writeGMSH(13,0d0*disp,grid,BIND_NODE,'s',0,t)
    do i=1,10
      t=t+1d0
      gradv=findGrad(v,grid,BIND_BLOCK)
      temp=findDispConvect(v,BIND_BLOCK,disp,grid,gradv,limiter=vanLeer)
      call mvGrid(grid,disp)
      call grid%updateBlockVol()
      v=v+temp/grid%BlockVol
      call writeGMSH(13,v,grid,BIND_BLOCK,'name1',i,t)
      call writeGMSH(13,dble(i)*disp,grid,BIND_NODE,'s',i,t)
    end do
    close(13)
  else
  end if
  call finalMPI()
end program
