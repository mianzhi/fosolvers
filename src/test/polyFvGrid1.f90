!----------------------------------------------------------------------------- best with 100 columns

function polyFvGrid1() result(ierr)
  use modFileIO
  use modPolyFvGrid
  integer ierr
  type(polyFvGrid)::grid
  
  ierr=0
  call readCGNS('data/tut21.cgns',grid)
  call grid%up()
  open(10,file='polyFvGrid1.rst',action='write')
  do i=1,grid%nP
    write(10,'(i6,7g14.7)'),i,grid%aP(i),grid%pP(:,i),grid%normP(:,i)
  end do
  close(10)
end function
