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
  
  call initMPI()
  if(pidMPI==0)then
    open(12,file='bin/gridGMSH4.msh',status='old')
    call readGMSH(12,grid)
    close(12)
    call grid%updateNodeVol()
    do i=1,grid%nNode
      write(*,*),i,grid%NodeVol(i)
    end do
    write(*,*),sum(grid%NodeVol(:)),sum(grid%BlockVol(:))
  else
  end if
  call finalMPI()
end program
