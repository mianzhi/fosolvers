!----------------------------------------------------------------------------- best with 100 columns

program tests
  call testReadGMSH()
end program

subroutine testReadGMSH()
  use moduleFileIO
  use moduleGrid
  use moduleSimpleSetLogic
  type(typeGrid)::grid
  integer,allocatable::a(:),c(:)
  integer b(3)
  a=[1,2,3,4,5,6]
  b=[4,6,8]
  c=findUnion(a,b)
  write(*,*),c
  
  open(12,file='bin/gridGMSH1.msh',status='old')
  call readGMSH(12,grid)
  close(12)
  open(13,file='rst.msh',status='replace')
  call writeGMSH(13,grid)
  close(13)
end subroutine
