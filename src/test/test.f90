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
  double precision,allocatable::v(:),gradv(:,:),disp(:,:)
  double precision t
  
  call initMPI()
  if(pidMPI==0)then
    open(12,file='bin/gridGMSH4.msh',status='old')
    call readGMSH(12,grid)
    close(12)
    call grid%updateBlockPos()
    call grid%updateBlockVol()
    allocate(disp(DIMS,grid%nNode))
    allocate(v(grid%nBlock))
    allocate(gradv(DIMS,grid%nBlock))
    disp(:,:)=0d0
    disp(1,:)=0.0019d0
    do i=1,grid%nBlock
      v(i)=merge(100d0,000d0,grid%BlockPos(1,i)<0.9d0.and.grid%BlockPos(1,i)>0.6d0)
    end do
    t=0d0
    open(13,file='rst.msh',status='replace')
    call writeGMSH(13,grid)
    call writeGMSH(13,v,grid,BIND_BLOCK,'name1',0,t)
    do i=1,100
      t=t+1d0
      gradv=findGrad(v,grid,BIND_BLOCK)
      v=v+findDispConvect(v,disp,grid,gradv)/grid%BlockVol
      call writeGMSH(13,v,grid,BIND_BLOCK,'name1',i,t)
    end do
    close(13)
  else
  end if
  call finalMPI()
end program
