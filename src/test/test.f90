!----------------------------------------------------------------------------- best with 100 columns

program test
  use moduleBasicDataStruct
  use moduleSimpleSetLogic
  use moduleFileIO
  use moduleGrid
  use moduleFVMGrad
  use moduleFVMDiffus
  use moduleFVMConvect
  use moduleInterpolation
  use moduleMPIComm
  type(typeGrid)::grid
  
  call initMPI()
  if(pidMPI==0)then
    open(12,file='bin/gridGMSH1.msh',status='old')
    call readGMSH(12,grid)
    close(12)
    call sendDat(grid%Point(1),1)
  else
    allocate(grid%Point(20))
    call recvDat(grid%Point(1),0)
    write(*,*),'recv',grid%Point(1)%Ent,grid%Point(1)%Shp,grid%Point(1)%nNode
    write(*,*),grid%Point(1)%iNode,grid%Point(1)%Dmn,grid%Point(1)%Prt
    !write(*,*),v(:,1),v(:,2),v(:,3)
    !open(13,file='rst.msh',status='replace')
    !call writeGMSH(13,grid)
    !close(13)
  end if
  call finalMPI()
end program
