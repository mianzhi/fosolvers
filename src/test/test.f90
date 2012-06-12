!----------------------------------------------------------------------------- best with 100 columns

program tests
  call testReadGMSH()
end program

subroutine testReadGMSH()
  use moduleFileIO
  use moduleGrid
  type(typeGrid)::grid
  
  open(12,file='bin/gridGMSH1.msh',status='old')
  call readGMSH(12,grid)
  close(12)
  write(*,*),grid%nPrt
  open(13,file='rst.msh',status='replace')
  call writeGMSH(13,grid)
  close(13)
end subroutine
