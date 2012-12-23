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
  use moduleCondition
  use moduleMPIComm
  type(typeGrid)::grid
  double precision box(DIMS,2)
  type(typeCondition),allocatable::condition(:)
  
  call initMPI()
  if(pidMPI==0)then
    open(12,file='bin/gridGMSH5.msh',status='old')
    call readGMSH(12,grid)
    close(12)
    box=findBoundBox(grid)
    write(*,*),box(:,1)
    write(*,*),box(:,2)
    open(13,file='bin/condition1',status='old')
    call readCondition(13,condition)
    write(*,*),condition(1)%Ent,condition(1)%bind,condition(1)%dat%get('a_longer_name')
    write(*,*),condition(2)%Ent,condition(2)%bind,condition(2)%dat%get('str')
  else
  end if
  call finalMPI()
end program
