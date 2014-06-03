!----------------------------------------------------------------------------- best with 100 columns

function polyGrid1() result(ierr)
  use modFileIO
  use modPolyGrid
  integer ierr
  type(polyGrid)::grid
  
  ierr=0
  open(10,file='data/grid_simple.vtk',action='read')
  call readVTK(10,grid)
  close(10)
  call grid%up()
  if(.not.abs(grid%v(1)-1d0)<1d-14)then
    ierr=ierr+1
  end if
  call grid%clear()
end function
