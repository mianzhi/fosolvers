!----------------------------------------------------------------------------- best with 100 columns

program tests
  call testReadGMSH()
end program

subroutine testReadGMSH()
  use moduleBasicDataStruct
  use moduleSimpleSetLogic
  use moduleFileIO
  use moduleGrid
  use moduleGridOperation
  type(typeGrid)::grid
  type(typeGrid),allocatable::sgrid(:)
  type(typeHtr1DIArr),allocatable::map(:,:)
  
  integer,allocatable::a(:)
  allocate(a(4))
  a(:)=[1,2,3,4]
  call applComplement(a,[2,4])
  write(*,*),a
  
  open(12,file='bin/gridGMSH1.msh',status='old')
  call readGMSH(12,grid)
  close(12)
  call splitGridPrt(grid,sgrid,map,overlap=.false.)
  open(13,file='rst.msh',status='replace')
  call writeGMSH(13,sgrid(2))
  close(13)
end subroutine
