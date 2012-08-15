!----------------------------------------------------------------------------- best with 100 columns

program test
  use moduleBasicDataStruct
  use moduleSimpleSetLogic
  use moduleFileIO
  use moduleGrid
  use moduleGridInspection
  use moduleGridOperation
  use moduleFVMGrad
  use moduleFVMDiffus
  use moduleFVMConvect
  use moduleInterpolation
  use moduleMPIComm
  type(typeGrid)::grid
  double precision box(DIMS,2)
  
  call initMPI()
  if(pidMPI==0)then
    open(12,file='bin/gridGMSH5.msh',status='old')
    call readGMSH(12,grid)
    close(12)
    box=findBoundBox(grid)
    write(*,*),box(:,1)
    write(*,*),box(:,2)
  else
  end if
  call finalMPI()
end program
