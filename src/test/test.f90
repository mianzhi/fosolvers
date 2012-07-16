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
    call grid%updateDualBlock()
    do i=1,grid%nNode
      write(*,*),i,':',grid%NodeVol(i)
      do j=1,size(grid%NodeNeibBlock(i)%dat)
        write(*,*),i,j,norm2(grid%NBAreaVect(i)%dat(:,j))
      end do
    end do
  else
  end if
  call finalMPI()
end program
