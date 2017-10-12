!----------------------------------------------------------------------------- best with 100 columns

function gradient2() result(ierr)
  use modFileIO
  use modPolyFvGrid
  use modGradient
  integer ierr
  type(polyFvGrid)::grid
  double precision,allocatable::s(:),grads(:,:)
  double precision::p(3)
  
  ierr=0
  open(10,file='data/bar_tet.vtk',action='read')
  call readVTK(10,grid)
  close(10)
  call grid%up()
  allocate(s(grid%nE))
  allocate(grads(3,grid%nE))
  do i=1,grid%nC
    p=grid%p(:,i)
    s(i)=p(1)+2d0*p(2)+3d0*p(3)
  end do
  do i=1,grid%nP
    if(grid%iEP(2,i)>grid%nC)then
      p=grid%p(:,grid%iEP(1,i))
      s(grid%iEP(2,i))=p(1)+2d0*p(2)+3d0*p(3)
    end if
  end do
  forall(i=1:grid%nE)
    grads(:,i)=[1d0,2d0,3d0]
  end forall
  call findGrad(grid,s,grads)
  open(10,file='gradient2_rst.vtk',action='write')
  call writeVTK(10,grid)
  call writeVTK(10,grid,E_DATA)
  call writeVTK(10,'s',s)
  call writeVTK(10,'grads',grads)
  close(10)
  call grid%clear()
  deallocate(s,grads)
end function
