!----------------------------------------------------------------------------- best with 100 columns

function eulerJac() result(ierr)
  use modFileIO
  use modPolyFvGrid
  use modAdvection
  integer ierr
  type(polyFvGrid)::grid
  double precision,allocatable::rho(:),rhou(:,:),rhoE(:),p(:),&
  &                             rho1(:),rhou1(:,:),rhoE1(:),p1(:),&
  &                             dRho(:),dRhou(:,:),dRhoE(:),&
  &                             JacP(:,:,:),JacC(:,:,:)
  double precision::pos(3),JacM(5,5),JacN(5,5)
  integer,parameter::M=8889,N=8896
  double precision,parameter::E1=1d-5,E2=1d0
  
  ierr=0
  open(10,file='data/bar_tet.vtk',action='read')
  call readVTK(10,grid)
  close(10)
  call grid%up()
  allocate(rho(grid%nE))
  allocate(rhou(3,grid%nE))
  allocate(rhoE(grid%nE))
  allocate(p(grid%nE))
  allocate(rho1(grid%nE))
  allocate(rhou1(3,grid%nE))
  allocate(rhoE1(grid%nE))
  allocate(p1(grid%nE))
  do i=1,grid%nE
    pos=grid%p(:,i)
    rho(i)=1d0+0.2d0*pos(1)
    rhou(:,i)=[200d0,0d0,0d0]*rho(i)
    rhoE(i)=(1d5+0.1d5*pos(1))/0.4d0
  end do
  ! numerical Jacobian
  ! no disturbance
  rho1(:)=rho(:)
  rhou1(:,:)=rhou(:,:)
  rhoE1(:)=rhoE(:)
  forall(i=1:grid%nC)
    p1(i)=(1.4d0-1d0)*(rhoE1(i)-0.5d0*dot_product(rhou1(:,i),rhou1(:,i))/rho1(i))
  end forall
  call findAdv(grid,rho1,rhou1,rhoE1,p1,1.4d0,dRho,dRhou,dRhoE)
  JacM(:,:)=spread([dRho(M),dRhou(:,M),dRhoE(M)],2,5)
  JacN(:,:)=spread([dRho(M),dRhou(:,M),dRhoE(M)],2,5)
  ! disturbance 1
  rho1(:)=rho(:)
  rho1(M)=rho(M)+E1
  rhou1(:,:)=rhou(:,:)
  rhoE1(:)=rhoE(:)
  forall(i=1:grid%nC)
    p1(i)=(1.4d0-1d0)*(rhoE1(i)-0.5d0*dot_product(rhou1(:,i),rhou1(:,i))/rho1(i))
  end forall
  call findAdv(grid,rho1,rhou1,rhoE1,p1,1.4d0,dRho,dRhou,dRhoE)
  JacM(:,1)=([dRho(M),dRhou(:,M),dRhoE(M)]-JacM(:,1))/E1
  ! disturbance 2
  rho1(:)=rho(:)
  rhou1(:,:)=rhou(:,:)
  rhou1(1,M)=rhou(1,M)+E1
  rhoE1(:)=rhoE(:)
  forall(i=1:grid%nC)
    p1(i)=(1.4d0-1d0)*(rhoE1(i)-0.5d0*dot_product(rhou1(:,i),rhou1(:,i))/rho1(i))
  end forall
  call findAdv(grid,rho1,rhou1,rhoE1,p1,1.4d0,dRho,dRhou,dRhoE)
  JacM(:,2)=([dRho(M),dRhou(:,M),dRhoE(M)]-JacM(:,2))/E1
  ! disturbance 3
  rho1(:)=rho(:)
  rhou1(:,:)=rhou(:,:)
  rhou1(2,M)=rhou(2,M)+E1
  rhoE1(:)=rhoE(:)
  forall(i=1:grid%nC)
    p1(i)=(1.4d0-1d0)*(rhoE1(i)-0.5d0*dot_product(rhou1(:,i),rhou1(:,i))/rho1(i))
  end forall
  call findAdv(grid,rho1,rhou1,rhoE1,p1,1.4d0,dRho,dRhou,dRhoE)
  JacM(:,3)=([dRho(M),dRhou(:,M),dRhoE(M)]-JacM(:,3))/E1
  ! disturbance 4
  rho1(:)=rho(:)
  rhou1(:,:)=rhou(:,:)
  rhou1(3,M)=rhou(3,M)+E1
  rhoE1(:)=rhoE(:)
  forall(i=1:grid%nC)
    p1(i)=(1.4d0-1d0)*(rhoE1(i)-0.5d0*dot_product(rhou1(:,i),rhou1(:,i))/rho1(i))
  end forall
  call findAdv(grid,rho1,rhou1,rhoE1,p1,1.4d0,dRho,dRhou,dRhoE)
  JacM(:,4)=([dRho(M),dRhou(:,M),dRhoE(M)]-JacM(:,4))/E1
  ! disturbance 5
  rho1(:)=rho(:)
  rhou1(:,:)=rhou(:,:)
  rhoE1(:)=rhoE(:)
  rhoE1(M)=rhoE(M)+E2
  forall(i=1:grid%nC)
    p1(i)=(1.4d0-1d0)*(rhoE1(i)-0.5d0*dot_product(rhou1(:,i),rhou1(:,i))/rho1(i))
  end forall
  call findAdv(grid,rho1,rhou1,rhoE1,p1,1.4d0,dRho,dRhou,dRhoE)
  JacM(:,5)=([dRho(M),dRhou(:,M),dRhoE(M)]-JacM(:,5))/E2
  ! disturbance 6
  rho1(:)=rho(:)
  rho1(N)=rho(N)+E1
  rhou1(:,:)=rhou(:,:)
  rhoE1(:)=rhoE(:)
  forall(i=1:grid%nC)
    p1(i)=(1.4d0-1d0)*(rhoE1(i)-0.5d0*dot_product(rhou1(:,i),rhou1(:,i))/rho1(i))
  end forall
  call findAdv(grid,rho1,rhou1,rhoE1,p1,1.4d0,dRho,dRhou,dRhoE)
  JacN(:,1)=([dRho(M),dRhou(:,M),dRhoE(M)]-JacN(:,1))/E1
  ! disturbance 7
  rho1(:)=rho(:)
  rhou1(:,:)=rhou(:,:)
  rhou1(1,N)=rhou(1,N)+E1
  rhoE1(:)=rhoE(:)
  forall(i=1:grid%nC)
    p1(i)=(1.4d0-1d0)*(rhoE1(i)-0.5d0*dot_product(rhou1(:,i),rhou1(:,i))/rho1(i))
  end forall
  call findAdv(grid,rho1,rhou1,rhoE1,p1,1.4d0,dRho,dRhou,dRhoE)
  JacN(:,2)=([dRho(M),dRhou(:,M),dRhoE(M)]-JacN(:,2))/E1
  ! disturbance 8
  rho1(:)=rho(:)
  rhou1(:,:)=rhou(:,:)
  rhou1(2,N)=rhou(2,N)+E1
  rhoE1(:)=rhoE(:)
  forall(i=1:grid%nC)
    p1(i)=(1.4d0-1d0)*(rhoE1(i)-0.5d0*dot_product(rhou1(:,i),rhou1(:,i))/rho1(i))
  end forall
  call findAdv(grid,rho1,rhou1,rhoE1,p1,1.4d0,dRho,dRhou,dRhoE)
  JacN(:,3)=([dRho(M),dRhou(:,M),dRhoE(M)]-JacN(:,3))/E1
  ! disturbance 9
  rho1(:)=rho(:)
  rhou1(:,:)=rhou(:,:)
  rhou1(3,N)=rhou(3,N)+E1
  rhoE1(:)=rhoE(:)
  forall(i=1:grid%nC)
    p1(i)=(1.4d0-1d0)*(rhoE1(i)-0.5d0*dot_product(rhou1(:,i),rhou1(:,i))/rho1(i))
  end forall
  call findAdv(grid,rho1,rhou1,rhoE1,p1,1.4d0,dRho,dRhou,dRhoE)
  JacN(:,4)=([dRho(M),dRhou(:,M),dRhoE(M)]-JacN(:,4))/E1
  ! disturbance 10
  rho1(:)=rho(:)
  rhou1(:,:)=rhou(:,:)
  rhoE1(:)=rhoE(:)
  rhoE1(N)=rhoE(N)+E2
  forall(i=1:grid%nC)
    p1(i)=(1.4d0-1d0)*(rhoE1(i)-0.5d0*dot_product(rhou1(:,i),rhou1(:,i))/rho1(i))
  end forall
  call findAdv(grid,rho1,rhou1,rhoE1,p1,1.4d0,dRho,dRhou,dRhoE)
  JacN(:,5)=([dRho(M),dRhou(:,M),dRhoE(M)]-JacN(:,5))/E2
  ! approximate analytical Jacobian
  call findAdvJac(grid,rho,rhou,rhoE,p,1.4d0,JacP,JacC)
  ! compare
  open(10,file='jac',action='write')
  write(10,*)'numerical'
  write(10,*)JacM(:,1)
  write(10,*)JacM(:,2)
  write(10,*)JacM(:,3)
  write(10,*)JacM(:,4)
  write(10,*)JacM(:,5)
  write(10,*)JacN(:,1)
  write(10,*)JacN(:,2)
  write(10,*)JacN(:,3)
  write(10,*)JacN(:,4)
  write(10,*)JacN(:,5)
  write(10,*)'approximate analytical'
  write(10,*)JacC(:,1,M)
  write(10,*)JacC(:,2,M)
  write(10,*)JacC(:,3,M)
  write(10,*)JacC(:,4,M)
  write(10,*)JacC(:,5,M)
  write(10,*)JacC(:,16,M)
  write(10,*)JacC(:,17,M)
  write(10,*)JacC(:,18,M)
  write(10,*)JacC(:,19,M)
  write(10,*)JacC(:,20,M)
  close(10)
  call grid%clear()
  deallocate(rho,rhou,rhoE,p,rho1,rhou1,rhoE1,p1,dRho,dRhou,dRhoE,JacP,JacC)
end function
