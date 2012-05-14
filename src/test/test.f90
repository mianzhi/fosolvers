!----------------------------------------------------------------------------- best with 100 columns

program tests
  call testReadGMSH()
end program

subroutine testReadGMSH()
  use moduleFileIO
  use moduleGrid
  type(typeGrid)::grid
  
  open(12,file='bin/gridGMSH2.msh',status='old')
  call readGMSH(12,grid)
  
  call grid%updateFacetArea()
  do i=1,grid%nFacet
    write(*,*),i,grid%FacetArea(i)
  end do
end subroutine
