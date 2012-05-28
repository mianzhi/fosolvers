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
  
  call grid%updateNodeNeibBlock()
  do i=1,grid%nNode
    write(*,*),i,grid%NodeNeibBlock(i)%dat(:)
  end do
end subroutine
