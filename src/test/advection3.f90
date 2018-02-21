!----------------------------------------------------------------------------- best with 100 columns

function advection3() result(ierr)
  use modFileIO
  use modPolyFvGrid
  use modAdvection
  integer ierr
  type(polyFvGrid)::grid
  double precision,allocatable::rho(:),u(:,:),rhou(:,:),flow(:),adv(:)
  double precision::p(3),dt
  
  ierr=0
  open(10,file='data/bar_tet.vtk',action='read')
  call readVTK(10,grid)
  close(10)
  call grid%up()
  allocate(rho(grid%nE))
  allocate(u(3,grid%nE))
  allocate(rhou(3,grid%nE))
  forall(i=1:grid%nE)
    u(:,i)=[1d0,0d0,0d0]
  end forall
  do i=1,grid%nE
    p=grid%p(:,i)
    if(p(1)<0.6d0.and.p(1)>0.4d0)then
      rho(i)=1.1d0
    else if(p(1)<0.3d0.and.p(1)>0.1d0)then
      rho(i)=1d0+0.1d0*sin((p(1)-0.1d0)*40d0*atan(1d0))
    else
      rho(i)=1d0
    end if
  end do
  forall(i=1:grid%nE)
    rhou(:,i)=u(:,i)*rho(i)
  end forall
  dt=0.0005d0
  do l=1,200
    call findMassFlow(grid,rhou,flow)
    call findAdv(grid,flow,adv)
    rho(1:grid%nC)=rho(1:grid%nC)+dt*adv(1:grid%nC)/grid%v(:)
    forall(i=1:grid%nC)
      rhou(:,i)=u(:,i)*rho(i)
    end forall
  end do
  open(10,file='advection3_rst.vtk',action='write')
  call writeVTK(10,grid)
  call writeVTK(10,grid,E_DATA)
  call writeVTK(10,'rho',rho)
  close(10)
  call grid%clear()
  deallocate(rho,u,rhou,flow,adv)
end function
