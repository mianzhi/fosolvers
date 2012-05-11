!----------------------------------------------------------------------------- best with 100 columns

program tests
  call testReadGMSH()
end program

subroutine testReadGMSH()
  use moduleFileIO
  use moduleGrid
  type(typeGrid)::grid
  
  open(12,file='bin/grid.msh',status='old')
  call readGMSH(12,grid)
  
  call grid%updateBlockPos()
  do i=1,grid%nBlock
    write(*,*),i,grid%BlockPos(:,i)
  end do
end subroutine
