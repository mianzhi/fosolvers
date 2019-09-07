!----------------------------------------------------------------------------- best with 100 columns

function euler() result(ierr)
  use modFileIO
  use modPolyFvGrid
  use modEuler
  integer ierr
  type(polyFvGrid)::grid
  double precision,allocatable::rho(:),rhou(:,:),rhoE(:),p(:),dRho(:),dRhou(:,:),dRhoE(:)
  double precision::pos(3),dt
  
  ierr=0
  open(10,file='data/bar_tet.vtk',action='read')
  call readVTK(10,grid)
  close(10)
  call grid%up()
  allocate(rho(grid%nE))
  allocate(rhou(3,grid%nE))
  allocate(rhoE(grid%nE))
  allocate(p(grid%nE))
  do i=1,grid%nE
    pos=grid%p(:,i)
    if(pos(1)<0.5d0)then
      p(i)=1d5
      rho(i)=1d0
    else
      p(i)=1d4
      rho(i)=0.125d0
    end if
    rhou(:,i)=0d0
    rhoE(i)=0d0+p(i)/0.4d0
  end do
  dt=1d-6
  do l=1,500
    call findEuler(grid,rho,rhou,rhoE,p,1.4d0,dRho,dRhou,dRhoE)
    forall(i=1:grid%nC)
      rho(i)=rho(i)+dt*dRho(i)/grid%v(i)
      rhou(:,i)=rhou(:,i)+dt*dRhou(:,i)/grid%v(i)
      rhoE(i)=rhoE(i)+dt*dRhoE(i)/grid%v(i)
      p(i)=0.4d0*(rhoE(i)-0.5d0*dot_product(rhou(:,i),rhou(:,i))/rho(i))
    end forall
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(n>grid%nC)then
        rho(n)=rho(m)
        rhou(:,n)=rhou(:,m)-2d0*dot_product(rhou(:,m),grid%normP(:,i))*grid%normP(:,i)
        rhoE(n)=rhoE(m)
        p(n)=p(m)
      end if
    end do
  end do
  open(10,file='euler_rst.vtk',action='write')
  call writeVTK(10,grid)
  call writeVTK(10,grid,E_DATA)
  call writeVTK(10,'rho',rho(:))
  close(10)
  call grid%clear()
  deallocate(rho,rhou,rhoE,p,dRho,dRhou,dRhoE)
end function
