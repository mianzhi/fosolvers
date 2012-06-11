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
  close(12)
  do i=1,grid%nFacet
    write(*,*),i,grid%Facet(i)%Ent,grid%Facet(i)%Prt(:)
  end do
  open(13,file='rst.msh',status='replace')
  call writeGMSH(13,grid)
  close(13)
end subroutine
