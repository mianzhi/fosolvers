!----------------------------------------------------------------------------- best with 100 columns

program tests
  call testReadGMSH()
end program

subroutine testReadGMSH()
  use moduleBasicDataStruct
  use moduleFileIO
  use moduleGrid
  use moduleGridOperation
  type(typeGrid)::grid
  type(typeGrid),allocatable::sgrid(:)
  type(typeHtr1DIArr),allocatable::map(:,:)
  
  open(12,file='bin/gridGMSH1.msh',status='old')
  call readGMSH(12,grid)
  close(12)
  call splitGridPrt(grid,sgrid,map,isOverlap=.false.)
  do i=1,grid%nPrt
    write(*,*),i,'::',size(map(MAP_FACET,i)%dat)
    do j=1,sgrid(i)%nFacet
      write(*,*),'    ',j,sgrid(i)%Facet(j)%Prt(:)
    end do
  end do
  open(13,file='rst.msh',status='replace')
  call writeGMSH(13,grid)
  close(13)
end subroutine
