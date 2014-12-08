!----------------------------------------------------------------------------- best with 100 columns

function gradient1() result(ierr)
  use modFileIO
  use modPolyFvGrid
  use modGradient
  integer ierr
  type(polyFvGrid)::grid
  double precision,allocatable::v(:),gradv(:,:)
  
  ierr=0
  call readCGNS('data/tut21.cgns',grid)
  call grid%up()
  allocate(v(grid%nC))
  do i=1,grid%nC
    v(i)=norm2(grid%p(:,i)+[0.1d0,0.1d0,0.1d0])
  end do
  call findGrad(grid,v,gradv)
  open(10,file='tut21.vtk',action='write')
  call writeVTK(10,grid)
  call writeVTK(10,grid,E_DATA)
  call writeVTK(10,'V',v)
  call writeVTK(10,'gradV',gradv)
  close(10)
  call grid%clear()
  deallocate(v,gradv)
end function
