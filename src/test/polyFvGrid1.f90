!----------------------------------------------------------------------------- best with 100 columns

function polyFvGrid1() result(ierr)
  use modFileIO
  use modPolyFvGrid
  integer ierr
  type(polyFvGrid)::grid
  
  ierr=0
  call readCGNS('/home/mianzhi/Desktop/tut21.cgns',grid)
  call grid%up()
  open(10,file='/home/mianzhi/Desktop/tut21.vtk',action='write')
  call writeVTK(10,grid)
  close(10)
  call grid%clear()
end function
