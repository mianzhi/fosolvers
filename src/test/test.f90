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
  type(typeGrid),allocatable::sgrid(:)
  type(typeHtr1DIArr),allocatable::r2o(:,:)
  
  call initMPI()
  if(pidMPI==0)then
    open(12,file='bin/gridGMSH1.msh',status='old')
    call readGMSH(12,grid)
    close(12)
    call splitGridPrt(grid,sgrid,r2o)
    call sendDat(sgrid(1),1)
  else
    call recvDat(grid,0)
    open(13,file='rst.msh',status='replace')
    call writeGMSH(13,grid)
    close(13)
  end if
  call finalMPI()
end program
