!----------------------------------------------------------------------------- best with 100 columns

!> environment for Euler problem
module modEuler
  use modPolyFvGrid
  use modCondition
  use iso_c_binding
  
  public
  
  integer,parameter::DIMS=3 !< three dimensions
  double precision,parameter::RTOL=1d-7 !< relative tolerance
  double precision,parameter::ATOL=1d-13 !< absolute tolerance
  
  type(polyFvGrid)::grid !< computational grid
  type(condTab)::bc !< boundary conditions
  
  double precision::t !< time
  double precision::tFinal !< final time
  double precision::tInt !< time interval of output
  double precision::tNext !< time for next output
  integer::iOut !< index of output
  integer,allocatable::iBC(:) !< indexes boundary conditions
  
  double precision::r !< specific gas constant
  double precision::gamm !< gamma of gas
  
  double precision,allocatable::rho(:) !< density
  double precision,allocatable::rhou(:,:) !< momentum
  double precision,allocatable::rhoE(:) !< total energy
  double precision,allocatable::u(:,:) !< velocity
  double precision,allocatable::p(:) !< pressure
  double precision,allocatable::temp(:) !< temperature
  double precision,allocatable::c(:) !< speed of sound
  
  double precision,allocatable::y(:) !< solution vector
  double precision,allocatable::dRho(:) !< time derivative of density
  double precision,allocatable::dRhou(:,:) !< time derivative of momentum
  double precision,allocatable::dRhoE(:) !< time derivative of total energy
  integer(C_LONG)::nEq !< number of equations
  integer(C_LONG)::iStat(100),iPar(1) !< integer solver outputs and parameters
  double precision::rStat(100),rPar(1) !< real solver outputs and parameters

contains
  
  !> initialize the simulation
  subroutine init()
    use modFileIO
    integer,parameter::FID=10
    double precision::p0,T0,u0(DIMS)
    
    ! read inputs
    open(FID,file='grid.msh',action='read')
    call readGMSH(FID,grid)
    call grid%up()
    close(FID)
    open(FID,file='bc',action='read')
    call readCondTab(FID,bc)
    close(FID)
    open(FID,file='ic',action='read')
    read(FID,*),p0
    read(FID,*),T0
    read(FID,*),u0(1)
    read(FID,*),u0(2)
    read(FID,*),u0(3)
    close(FID)
    open(FID,file='sim',action='read')
    read(FID,*),tFinal
    read(FID,*),tInt
    close(FID)
    open(FID,file='fl',action='read')
    read(FID,*),r
    read(FID,*),gamm
    close(FID)
    ! work space and initial state
    allocate(rho(grid%nE))
    allocate(rhou(DIMS,grid%nE))
    allocate(rhoE(grid%nE))
    allocate(u(DIMS,grid%nE))
    allocate(p(grid%nE))
    allocate(temp(grid%nE))
    allocate(c(grid%nE))
    forall(i=1:grid%nE)
      rho(i)=p0/r/T0
      rhou(:,i)=rho(i)*u0(:)
      rhoE(i)=rho(i)*(1d0/(gamm-1d0)*r*T0+0.5d0*dot_product(u0,u0))
      u(:,i)=u0(:)
      p(i)=p0
      temp(i)=T0
      c(i)=sqrt(gamm*r*T0)
    end forall
    t=0d0
    tNext=tInt
    iOut=0
    ! indexes of boundary conditions
    allocate(iBC(grid%nE))
    iBC(:)=0
    do i=1,grid%nE
      do j=1,size(bc%gid)
        if(grid%gid(i)==bc%gid(j))then
          iBC(i)=j
          exit
        end if
      end do
    end do
    ! initialize ODE solver
    nEq=5*grid%nC
    allocate(y(nEq))
    allocate(dRho(grid%nC))
    allocate(dRhou(DIMS,grid%nC))
    allocate(dRhoE(grid%nC))
    forall(i=1:grid%nC)
      y(i)=rho(i)
      y(grid%nC+i)=rhou(1,i)
      y(2*grid%nC+i)=rhou(2,i)
      y(3*grid%nC+i)=rhou(3,i)
      y(4*grid%nC+i)=rhoE(i)
      dRho(i)=0d0
      dRhou(:,i)=0d0
      dRhoE(i)=0d0
    end forall
    call fnvinits(1,nEq,ier)
    call fcvmalloc(t,y,2,2,1,RTOL,ATOL,iStat,rStat,iPar,rPar,ier)
    call fcvspgmr(0,1,0,0d0,ier)
    call fcvsetrin('MAX_STEP',tInt/50d0,ier)
  end subroutine
  
  !> clear the simulation environment
  subroutine clear()
    deallocate(iBC,rho,rhou,rhoE,u,p,temp,c)
    deallocate(y,dRho,dRhou,dRhoE)
    call grid%clear()
    call bc%clear()
  end subroutine
  
  !> synchronize the state with the solution vector
  subroutine syncState(x)
    double precision,intent(in)::x(*) !< solution vector
    
    forall(i=1:grid%nC)
      rho(i)=x(i)
      rhou(1,i)=x(grid%nC+i)
      rhou(2,i)=x(2*grid%nC+i)
      rhou(3,i)=x(3*grid%nC+i)
      rhoE(i)=x(4*grid%nC+i)
    end forall
    forall(i=1:grid%nE)
      u(:,i)=rhou(:,i)/rho(i)
      p(i)=(gamm-1d0)*(rhoE(i)-0.5d0*dot_product(rhou(:,i),rhou(:,i))/rho(i))
      temp(i)=p(i)/rho(i)/r
      c(i)=sqrt(gamm*r*temp(i))
    end forall
  end subroutine
  
  !> write the state
  subroutine writeState(fName)
    use modFileIO
    character(*),intent(in)::fName
    integer,parameter::FID=10
    
    open(FID,file=trim(fName),action='write')
    call writeVTK(FID,grid)
    call writeVTK(FID,grid,E_DATA)
    call writeVTK(FID,'density',rho)
    call writeVTK(FID,'velocity',u)
    call writeVTK(FID,'pressure',p)
    call writeVTK(FID,'temperature',temp)
    close(FID)
  end subroutine
  
end module

!> Euler solver
program foeuler
  use modEuler
  character(20)::tmpStr
  
  call init()
  write(tmpStr,*),iOut
  call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
  do while(t<tFinal)
    call fcvode(tFinal,t,y,2,ier)
    if(t>=tNext)then
      iOut=iOut+1
      write(tmpStr,*),iOut
      call syncState(y)
      call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
      tNext=tNext+tInt
    end if
  end do
  call clear()
end program

!> the ODE to be advanced
subroutine fcvfun(time,x,dxdt,iPara,rPara,ier)
  use modEuler
  use modAdvection
  use iso_c_binding
  double precision::time
  double precision::x(*)
  double precision::dxdt(*)
  integer(C_LONG)::iPara(*)
  double precision::rPara(*)
  integer::ier
  
  do i=1,grid%nP
    m=grid%iEP(1,i)
    n=grid%iEP(2,i)
    if(n>grid%nC)then
      ! TODO set boundary conditions
      rho(n)=rho(m)
      rhou(:,n)=rhou(:,m)-2d0*dot_product(rhou(:,m),grid%normP(:,i))*grid%normP(:,i)
      rhoE(n)=rhoE(m)
      p(n)=p(m)
    end if
  end do
  call syncState(x)
  call findAdv(grid,rho,rhou,rhoE,p,dRho,dRhou,dRhoE)
  forall(i=1:grid%nC)
    dxdt(i)=dRho(i)
    dxdt(grid%nC+i)=dRhou(1,i)
    dxdt(2*grid%nC+i)=dRhou(2,i)
    dxdt(3*grid%nC+i)=dRhou(3,i)
    dxdt(4*grid%nC+i)=dRhoE(i)
  end forall
  ier=0
end subroutine
