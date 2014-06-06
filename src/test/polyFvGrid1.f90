!----------------------------------------------------------------------------- best with 100 columns

function polyFvGrid1() result(ierr)
  use modFileIO
  use modPolyFvGrid
  integer ierr
  type(polyFvGrid)::grid
  
  ierr=0
  open(10,file='/home/mianzhi/Desktop/grid_classical.vtk',action='read')
  call readVTK(10,grid)
  close(10)
  call grid%up()
  !do i=1,grid%nC
  !  write(*,*),'neib',i,grid%neib(:,i)
  !end do
  !do i=1,grid%nP
  !  write(*,*),'pair',i,grid%iEP(:,i)
  !end do
  call grid%clear()
end function
