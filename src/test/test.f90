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
  
  write(*,*),grid%nNode,grid%nPoint,grid%nLine,grid%nFacet,grid%nBlock
  do i=1,grid%nBlock
    write(*,*),i,grid%Block(i)%Shp,grid%Block(i)%Ent,grid%Block(i)%nNode,grid%Block(i)%Dmn(:),grid%Block(i)%Prt(:)
    write(*,*),grid%Block(i)%iNode(:)
  end do
end subroutine
