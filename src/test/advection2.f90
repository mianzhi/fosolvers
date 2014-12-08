!----------------------------------------------------------------------------- best with 100 columns

function advection2() result(ierr)
  use modFileIO
  use modPolyFvGrid
  use modAdvection
  integer ierr
  type(polyFvGrid)::grid
  double precision,allocatable::u(:,:),f(:,:),tmpu(:)
  double precision::p(3),dt
  
  ierr=0
  open(10,file='data/bar_tet.vtk',action='read')
  call readVTK(10,grid)
  close(10)
  call grid%up()
  allocate(u(3,grid%nC))
  allocate(f(3,grid%nC))
  u(:,:)=0d0
  do i=1,grid%nC
    p=grid%p(:,i)
    if(p(1)<0.5d0.and.p(1)>0.4d0)then
      u(1,i)=(p(1)-0.4d0)*10d0
    else if(p(1)<0.6d0.and.p(1)>0.5d0)then
      u(1,i)=(0.6d0-p(1))*10d0
    end if
  end do
  dt=0.0005d0
  do l=1,200
    forall(i=1:grid%nC)
      f(:,i)=0.5d0*u(1,i)*u(:,i)
    end forall
    call findAdv(grid,u(1,:),f,tmpu)
    u(1,:)=u(1,:)+dt*tmpu(:)/grid%v(:)
  end do
  open(10,file='burgers_rst_1.vtk',action='write')
  call writeVTK(10,grid)
  call writeVTK(10,grid,E_DATA)
  call writeVTK(10,'u',[u(1,:),[(0d0,i=1,grid%nE-grid%nC)]])
  close(10)
  do l=1,200
    forall(i=1:grid%nC)
      f(:,i)=0.5d0*u(1,i)*u(:,i)
    end forall
    call findAdv(grid,u(1,:),f,tmpu)
    u(1,:)=u(1,:)+dt*tmpu(:)/grid%v(:)
  end do
  open(10,file='burgers_rst_2.vtk',action='write')
  call writeVTK(10,grid)
  call writeVTK(10,grid,E_DATA)
  call writeVTK(10,'u',[u(1,:),[(0d0,i=1,grid%nE-grid%nC)]])
  close(10)
  call grid%clear()
  deallocate(u,f,tmpu)
end function
