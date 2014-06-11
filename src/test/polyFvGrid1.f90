!----------------------------------------------------------------------------- best with 100 columns

function polyFvGrid1() result(ierr)
  use modFileIO
  use modPolyFvGrid
  integer ierr
  type(polyFvGrid)::grid
  
  ierr=0
  call readCGNS('/home/mianzhi/Desktop/tut21.cgns',grid)
  call grid%up()
  do i=1,grid%nC
    !write(*,*),i,grid%gid(i)
  end do
  open(10,file='/home/mianzhi/Desktop/tut21.vtk',action='write')
  call writeVTK(10,grid)
  call writeVTK(10,grid,N_DATA)
  close(10)
  call grid%clear()
end function
