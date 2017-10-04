!----------------------------------------------------------------------------- best with 100 columns

!> environment for the compressible PISO solver
module modPiso
  use modPolyFvGrid
  use modCondition
  use modUDF
  use modNumerics
  use iso_c_binding
  
  public
  
  integer,parameter::DIMS=3 !< three dimensions
  integer,parameter::MAXIT_MOMENTUM=40 !< max number of momentum prediction equation iteration
  integer,parameter::MAXIT_PRESSURE=20 !< max number of pressure correction equation iteration
  integer,parameter::MAXIT_DENSITY=20 !< max number of density equation iteration
  integer,parameter::MAXIT_ENERGY=20 !< max number of energy equation iteration
  integer,parameter::MAXIT_PISO=20 !< max number of PISO iterations
  double precision,parameter::RTOL_MOMENTUM=1d-8 !< momentum prediction relative tolerance
  double precision,parameter::RTOL_PRESSURE=1d-8 !< pressure correction relative tolerance
  double precision,parameter::RTOL_DENSITY=1d-7 !< density equation relative tolerance
  double precision,parameter::RTOL_ENERGY=1d-8 !< energy equation relative tolerance
  
  type(polyFvGrid)::grid !< computational grid
  
  double precision::t !< time
  double precision::tFinal !< final time
  double precision::dt !< time step size [s]
  double precision::tInt !< time interval of output
  double precision::tNext !< time for next output
  integer::iOut !< index of output
  integer::nRetry !< number of retry for the current time step
  logical::needRetry !< if retry of the current time step is needed
  
  double precision,allocatable::rho(:) !< density [kg/m^3]
  double precision,allocatable::rhou(:,:) !< momentum [kg/s/m^2]
  double precision,allocatable::rhoE(:) !< total energy [J/m^3]
  double precision,allocatable::p(:) !< pressure [Pa]
  double precision,allocatable::u(:,:) !< velocity [m/s]
  double precision,allocatable::temp(:) !< temperature [K]
  double precision,allocatable::Y(:,:) !< mass fraction of species
  
  ! state at the beginning of a time step
  double precision,allocatable::rho0(:),rhou0(:,:),rhoE0(:),p0(:),u0(:,:),temp0(:),Y0(:,:)
  
  double precision,allocatable::visc(:) !< viscosity [Pa*s]
  double precision,allocatable::cond(:) !< thermal conductivity [W/m/K]
  double precision,allocatable::gradP(:,:) !< gradient of p [Pa/m]
  double precision,allocatable::laP(:) !< Laplacian of p [Pa/m^2] by diffusion scheme
  double precision,allocatable::viscF(:,:) !< viscous force applied on cell [N]
  double precision,allocatable::presF(:,:) !< pressure force applied on cell [N]
  double precision,allocatable::condQ(:) !< heat conduction into cell [W]
  double precision,allocatable::carrier(:,:) !< {rho,rhou,rhoH} of each cell [*]
  double precision,allocatable::flowRho(:) !< flow rate of rho [kg/s]
  double precision,allocatable::fluxRhou(:,:,:) !< flux of rhou [kg/m/s^2]
  double precision,allocatable::flowRhou(:,:) !< flow rate of rhou into cell [kg*m/s^2]
  double precision,allocatable::fluxRhoH(:,:) !< flux of rhoH [J/m^2/s]
  double precision,allocatable::flowRhoH(:) !< flow rate of rhoH into cell [J/s]
  
  ! auxiliary variables for the PISO iteration
  double precision,allocatable::rho1(:),p1(:),laP1(:),presF1(:,:),temp1(:),dRho(:)
  
  ! scales
  double precision,allocatable::rhoScale !< density scale [kg/m^3]
  double precision,allocatable::rhouScale !< momentum scale [kg/s/m^2]
  double precision,allocatable::rhoEScale !< total energy scale [J/m^3]
  double precision,allocatable::pScale !< pressure scale [Pa]
  double precision,allocatable::uScale !< velocity scale [m/s]
  double precision,allocatable::tempScale !< temperature scale [K]
  
  ! data for algebraic solver
  type(fixPt)::momentumEq !< momentum equation as a fix point problem
  type(NewtonKrylov)::pressureEq !< pressure correction problem to be solved by Newton-GMRES
  type(fixPt)::densityEq !< density equation as a fix point problem
  type(fixPt)::energyEq !< energy equation as a fix point problem
  integer::nItMomentum,nItPressure,nItDensity,nItEnergy,nItPISO !< number of iterations
  
contains
  
  !> initialize the simulation
  subroutine init()
    use modFileIO
    integer,parameter::FID=10
    double precision::pInit,uInit(DIMS),TInit,pE(DIMS)!,tmpD
    double precision,allocatable::YInit(:)
    
    ! read inputs
    open(FID,file='grid.msh',action='read')
    call readGMSH(FID,grid)
    call grid%up()
    close(FID)
    !open(FID,file='bc',action='read')
    !call readCondTab(FID,bc)
    !if(any(bc%t(:)==BC_IN_STATIC_UDF).or.&
    !&  any(bc%t(:)==BC_IN_TOTAL_UDF).or.&
    !&  any(bc%t(:)==BC_OUT_UDF).or.&
    !&  any(bc%t(:)==BC_FAR_UDF))then
    !  udfBc=1
    !else
    !  udfBc=0
    !end if
    !close(FID)
    !open(FID,file='ic',action='read')
    !read(FID,*),udfIC
    !if(udfIc==0)then
    !  read(FID,*),pInit
    !  read(FID,*),TInit
    !  read(FID,*),uInit(1)
    !  read(FID,*),uInit(2)
    !  read(FID,*),uInit(3)
    !else
    !  read(FID,*),tmpD
    !  iUdf(1)=int(tmpD)
    !  read(FID,*),tmpD
    !  iUdf(2)=int(tmpD)
    !  read(FID,*),tmpD
    !  iUdf(3)=int(tmpD)
    !  read(FID,*),tmpD
    !  iUdf(4)=int(tmpD)
    !  read(FID,*),tmpD
    !  iUdf(5)=int(tmpD)
    !end if
    !close(FID)
    open(FID,file='sim',action='read')
    read(FID,*)tFinal
    read(FID,*)tInt
    close(FID)
    !open(FID,file='fl',action='read')
    !read(FID,*),r
    !read(FID,*),gamm
    !close(FID)
    !if(udfIc/=0.or.udfBc/=0)then
    !  open(FID,file='udf',action='read')
    !  call readUDFTab(FID,udf)
    !  close(FID)
    !end if
    ! indexes of boundary conditions
    !allocate(iBC(grid%nE))
    !iBC(:)=0
    !do i=1,grid%nE
    !  do j=1,size(bc%gid)
    !    if(grid%gid(i)==bc%gid(j))then
    !      iBC(i)=j
    !      exit
    !    end if
    !  end do
    !end do
    ! initialize algebraic solver
    call momentumEq%init(grid%nC*DIMS,momentumRHS,maa=MAXIT_MOMENTUM)
    call momentumEq%setMaxIt(MAXIT_MOMENTUM)
    call pressureEq%init(grid%nC,pressureRHS)
    call pressureEq%setMaxIt(MAXIT_PRESSURE)
    call densityEq%init(grid%nC,densityRHS,maa=MAXIT_DENSITY)
    call densityEq%setMaxIt(MAXIT_DENSITY)
    call energyEq%init(grid%nC,energyRHS,maa=MAXIT_ENERGY)
    call energyEq%setMaxIt(MAXIT_ENERGY)
    ! work space and initial state
    allocate(rho(grid%nE))
    allocate(rho0(grid%nE))
    allocate(rho1(grid%nE))
    allocate(rhou(DIMS,grid%nE))
    allocate(rhou0(DIMS,grid%nE))
    allocate(rhoE(grid%nE))
    allocate(rhoE0(grid%nE))
    allocate(p(grid%nE))
    allocate(p0(grid%nE))
    allocate(p1(grid%nE))
    allocate(u(DIMS,grid%nE))
    allocate(u0(DIMS,grid%nE))
    allocate(temp(grid%nE))
    allocate(temp0(grid%nE))
    allocate(temp1(grid%nE))
    ! FIXME fix nSp
    allocate(Y(1,grid%nE))
    allocate(Y0(1,grid%nE))
    allocate(YInit(1))
    allocate(visc(grid%nE))
    allocate(cond(grid%nE))
    allocate(fluxRhou(DIMS,DIMS,grid%nE))
    allocate(fluxRhoH(DIMS,grid%nE))
    allocate(gradP(DIMS,grid%nC))
    allocate(laP(grid%nC))
    allocate(laP1(grid%nC))
    allocate(viscF(DIMS,grid%nC))
    allocate(presF(DIMS,grid%nC))
    allocate(presF1(DIMS,grid%nC))
    allocate(condQ(grid%nC))
    allocate(flowRho(grid%nC))
    allocate(flowRhou(DIMS,grid%nC))
    allocate(flowRhoH(grid%nC))
    allocate(dRho(grid%nC))
    do i=1,grid%nE
    !  if(udfIc/=0)then
    !    pE(:)=grid%p(:,i)
    !    pInit=udf%eval(iUdf(1),pE,0d0)
    !    TInit=udf%eval(iUdf(2),pE,0d0)
    !    uInit(1)=udf%eval(iUdf(3),pE,0d0)
    !    uInit(2)=udf%eval(iUdf(4),pE,0d0)
    !    uInit(3)=udf%eval(iUdf(5),pE,0d0)
    !  end if
      ! FIXME remove below
      pE(:)=grid%p(:,i)
      pInit=min(max(0.9d5,-0.5d5*(pE(3)-0.5d0)+0.95d5),1d5)
      uInit(:)=[0d0,0d0,0d0]
      TInit=298d0
      YInit=[1d0]
      ! FIXME remove above
      ! FIXME calculation of rho, rhoE based on cantera
      p(i)=pInit
      u(:,i)=uInit(:)
      temp(i)=TInit
      Y(:,i)=YInit(:)
    end do
    call recoverState(p,u,temp,Y,rho,rhou,rhoE)
    t=0d0
    dt=1d-7 ! TODO replace the good old 1d-7 with user input
    tNext=tInt
    iOut=0
    nRetry=0
    needRetry=.false.
    deallocate(YInit)
  end subroutine
  
  !> clear the simulation environment
  subroutine clear()
    call grid%clear()
    call momentumEq%clear()
    call pressureEq%clear()
    call densityEq%clear()
    call energyEq%clear()
  end subroutine
  
  !> record {rho,rhoU,rhoE,p,u,temp,Y} in {rho0,rhoU0,rhoE0,p0,u0,temp0,Y0}
  subroutine recordState0()
    rho0(:)=rho(:)
    rhou0(:,:)=rhou(:,:)
    rhoE0(:)=rhoE(:)
    p0(:)=p(:)
    u0(:,:)=u(:,:)
    temp0(:)=temp(:)
    Y0(:,:)=Y(:,:)
  end subroutine
  
  !> load {rho,rhoU,rhoE,p,u,temp,Y} from {rho0,rhoU0,rhoE0,p0,u0,temp0,Y0}
  subroutine loadState0()
    rho(:)=rho0(:)
    rhou(:,:)=rhou0(:,:)
    rhoE(:)=rhoE0(:)
    p(:)=p0(:)
    u(:,:)=u0(:,:)
    temp(:)=temp0(:)
    Y(:,:)=Y0(:,:)
  end subroutine
  
  !> derive primitive state {p,u,T} from conserved state {rho,rhou,rhoE}
  subroutine deriveState(rhoi,rhoui,rhoEi,Yi,po,uo,tempo)
    double precision,intent(in)::rhoi(:) !< density
    double precision,intent(in)::rhoui(:,:) !< momentum
    double precision,intent(in)::rhoEi(:) !< total energy
    double precision,intent(in)::Yi(:,:) !< mass fraction of species
    double precision,intent(inout)::po(:) !< pressure
    double precision,intent(inout)::uo(:,:) !< velocity
    double precision,intent(inout)::tempo(:) !< temperature
    
    !$omp parallel do default(shared)
    do i=1,grid%nE
      uo(:,i)=rhoui(:,i)/rhoi(i)
      ! FIXME calculation of p, T based on cantera
      po(i)=(1.4d0-1d0)*(rhoEi(i)-0.5d0*dot_product(rhoui(:,i),rhoui(:,i))/rhoi(i))
      tempo(i)=po(i)/rhoi(i)/287.058d0
    end do
    !$omp end parallel do
  end subroutine
  
  !> recover conserved state {rho,rhou,rhoE} from primitive state {p,u,T}
  subroutine recoverState(pi,ui,tempi,Yi,rhoo,rhouo,rhoEo)
    double precision,intent(in)::pi(:) !< pressure
    double precision,intent(in)::ui(:,:) !< velocity
    double precision,intent(in)::tempi(:) !< temperature
    double precision,intent(in)::Yi(:,:) !< mass fraction of species
    double precision,intent(inout)::rhoo(:) !< density
    double precision,intent(inout)::rhouo(:,:) !< momentum
    double precision,intent(inout)::rhoEo(:) !< total energy
    
    !$omp parallel do default(shared)
    do i=1,grid%nE
      ! FIXME calculation of rho, rhoE based on cantera
      rhoo(i)=pi(i)/287.058d0/tempi(i)
      rhouo(:,i)=rhoo(i)*ui(:,i)
      rhoEo(i)=rhoo(i)*(1d0/(1.4d0-1d0)*287.058d0*tempi(i)+0.5d0*dot_product(ui(:,i),ui(:,i)))
    end do
    !$omp end parallel do
  end subroutine
  
  !> extract primitive state {p,u,T} from variable vector
  subroutine extractVar(var,po,uo,tempo)
    double precision,intent(in)::var(:) !< variable vector of the nonlinear problem
    double precision,intent(inout)::po(:) !< pressure
    double precision,intent(inout)::uo(:,:) !< velocity
    double precision,intent(inout)::tempo(:) !< temperature
    
    !$omp parallel do default(shared)&
    !$omp& private(j)
    do i=1,grid%nC
      j=(i-1)*5
      po(i)=p0(i)+var(j+1)
      uo(:,i)=u0(:,i)+var(j+2:j+4)
      tempo(i)=temp0(i)+var(j+5)
    end do
    !$omp end parallel do
  end subroutine
  
  !> write the state to post-processing file
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
  
  !> calculate time step size, scales etc.
  subroutine preSolve()
    double precision,parameter::CFL_ACCOUSTIC=2d0
    double precision,parameter::CFL_FLOW=0.5d0
    double precision,parameter::CFL_DIFFUSION=0.5d0
    double precision,parameter::MINFRAC_OUT_ALIGN=0.05d0
    
    ! TODO calculate dt, cantera sound speed
    dt=min(dt,minval(CFL_ACCOUSTIC*grid%v(:)**(1d0/3d0)&
    &                     /sqrt(1.4d0*p(1:grid%nC)/rho(1:grid%nC))))
    dt=min(dt,minval(CFL_FLOW*grid%v(:)**(1d0/3d0)&
    &                /norm2(u(:,1:grid%nC),1)))
    dt=min(dt,minval(CFL_DIFFUSION*grid%v(:)**(2d0/3d0)&
    &                /(visc(1:grid%nC)/rho(1:grid%nC))))
    if(nRetry>0)then
      dt=dt/2d0**nRetry
      write(*,'(a,i2,a,g12.6)')'[i] starting retry No. ',nRetry,' at t: ',dt
    end if
    dt=max(min(dt,tNext-t),dt*MINFRAC_OUT_ALIGN)
    
    ! solution and residual scales
    pScale=max(maxval(p)-minval(p),maxval(0.5d0*rho*norm2(u,1)**2),1d0)
    uScale=max(maxval(norm2(u,1)),sqrt(2d0*pScale/minval(rho)))
    tempScale=max(maxval(temp)-minval(temp),pScale/minval(rho)/287.058d0)
    rhoScale=maxval(rho)
    rhouScale=rhoScale*uScale
    rhoEScale=max(rhoScale*287.058*tempScale,0.5d0*rhoScale*uScale**2)
    
    ! set tolerance and maximum number of iterations, and etc.
    call momentumEq%setTol(rhoScale*RTOL_MOMENTUM)
    call pressureEq%setTol(pScale*RTOL_PRESSURE)
    call densityEq%setTol(rhoScale*RTOL_DENSITY)
    call energyEq%setTol(rhoEScale*RTOL_ENERGY)
    nItMomentum=0
    nItPressure=0
    nItDensity=0
    nItEnergy=0
    
    write(*,'(a,g12.6,a,g12.6)')'[i] starting step, t: ',t,' dt: ',dt
  end subroutine
  
  !> guess next time step size, etc.
  subroutine postSolve()
    
    ! the time step size is adjusted to maintain half of the maximum number iterations
    dt=dt*minval(dble([MAXIT_MOMENTUM,MAXIT_PRESSURE,MAXIT_DENSITY,MAXIT_ENERGY,MAXIT_PISO])/&
    &            dble([nItMomentum,nItPressure,nItDensity,nItEnergy,nItPISO]))*0.5d0
    
    write(*,'(a,i2,a,i2,a,i2,a,i2,a,i2,a)')'[i] finished step, nIt[rhou,p,rho,rhoE,PISO]: [',&
    &    nItMomentum,',',nItPressure,',',nItDensity,',',nItEnergy,',',nItPISO,']'
  end subroutine
  
  !> set the boundary conditions
  subroutine setBC()
    
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(n>grid%nC)then
        if(.true.)then ! default wall boundary
          p(n)=p(m)
          u(:,n)=-u(:,m)
          temp(n)=temp(m)
          Y(:,n)=Y(:,m)
          rho(n)=rho(m)
          rhou(:,n)=-rhou(:,m)
          rhoE(n)=rhoE(m)
        end if
      end if
    end do
  end subroutine
  
  !> predict the rhou and u, during which the presF is updated
  subroutine predictMomentum()
    use modGradient
    use modPressure
    double precision,allocatable,save::tmpRhou(:)
    integer::info
    
    if(.not.allocated(tmpRhou))then
      allocate(tmpRhou(DIMS*grid%nC))
    else if(size(tmpRhou)/=DIMS*grid%nC)then
      deallocate(tmpRhou)
      allocate(tmpRhou(DIMS*grid%nC))
    end if
    tmpRhou=reshape(rhou(:,1:grid%nC),[DIMS*grid%nC])
    call setBC()
    call findPresForce(grid,p,presF)
    call findGrad(grid,p,gradP)
    call momentumEq%solve(tmpRhou,info=info)
    needRetry=info/=0
    call momentumEq%getNIt(n)
    nItMomentum=max(nItMomentum,n)
    rhou(:,1:grid%nC)=reshape(tmpRhou,[DIMS,grid%nC])
    forall(i=1:grid%nE)
      u(:,i)=rhou(:,i)/rho(i)
    end forall
  end subroutine
  
  !> solve the pressure, during which the laP is updated
  subroutine solvePressure()
    use modAdvection
    double precision,allocatable,save::tmpP(:)
    integer::info
    
    if(.not.allocated(tmpP))then
      allocate(tmpP(grid%nC))
    else if(size(tmpP)/=grid%nC)then
      deallocate(tmpP)
      allocate(tmpP(grid%nC))
    end if
    tmpP(:)=p(1:grid%nC)
    call findAdv(grid,rho,rhou,flowRho) ! vary only the Rhie-Chow part of flowRho in RHS
    call pressureEq%solve(tmpP,info=info)
    needRetry=info/=0
    call pressureEq%getNIt(n)
    nItPressure=max(nItPressure,n)
    p(1:grid%nC)=tmpP(:)
  end subroutine
  
  !> correct the momentum with updated pressure, during which presF and u is updated
  subroutine correctMomentum()
    use modPressure
    
    call setBC()
    call findPresForce(grid,p,presF)
    forall(i=1:grid%nC)
      rhou(:,i)=rhou(:,i)+dt/grid%v*(presF(:,i)-presF1(:,i))
    end forall
  end subroutine
  
  !> solve the density
  subroutine solveDensity()
    double precision,allocatable,save::tmpRho(:)
    integer::info
    
    if(.not.allocated(tmpRho))then
      allocate(tmpRho(grid%nC))
    else if(size(tmpRho)/=grid%nC)then
      deallocate(tmpRho)
      allocate(tmpRho(grid%nC))
    end if
    tmpRho(:)=rho(1:grid%nC)
    call densityEq%solve(tmpRho,info)
    needRetry=info/=0
    call densityEq%getNIt(n)
    nItDensity=max(nItDensity,n)
    rho(1:grid%nC)=tmpRho(:)
  end subroutine
  
  !> solve the (total) energy, during which T is updated
  subroutine solveEnergy()
    use modAdvection
    use modDiffusion
    double precision,allocatable,save::tmpRhoE(:)
    integer::info
    
    if(.not.allocated(tmpRhoE))then
      allocate(tmpRhoE(grid%nC))
    else if(size(tmpRhoE)/=grid%nC)then
      deallocate(tmpRhoE)
      allocate(tmpRhoE(grid%nC))
    end if
    tmpRhoE(:)=rhoE(1:grid%nC)
    call energyEq%solve(tmpRhoE,info)
    needRetry=info/=0
    call energyEq%getNIt(n)
    nItEnergy=max(nItEnergy,n)
    rhoE(1:grid%nC)=tmpRhoE(:)
  end subroutine
  
  !> RHS of the momentum equation in the form of a fix point problem
  function momentumRHS(oldVector,newVector,dat)
    use iso_c_binding
    use ieee_arithmetic
    use modNumerics
    use modAdvection
    use modNewtonian
    type(C_PTR),value::oldVector !< old N_Vector
    type(C_PTR),value::newVector !< new N_Vector
    type(C_PTR),value::dat !< optional user data object
    integer(C_INT)::momentumRHS !< error code
    double precision,pointer::x(:) !< fortran pointer associated with oldVector
    double precision,pointer::y(:) !< fortran pointer associated with newVector
    
    call associateVector(oldVector,x)
    call associateVector(newVector,y)
    
    rhou(:,1:grid%nC)=reshape(x,[DIMS,grid%nC])
    forall(i=1:grid%nE)
      u(:,i)=rhou(:,i)/rho(i)
    end forall
    call setBC()
    call findViscForce(grid,u,visc,viscF)
    viscF(:,:)=0d0 ! FIXME: enable viscous force
    forall(i=1:grid%nE,j=1:DIMS)
      fluxRhou(:,j,i)=rhou(j,i)*u(:,i)
    end forall
    call findAdv(grid,rhou,fluxRhou,flowRhou)
    forall(i=1:grid%nC)
      rhou(:,i)=rhou0(:,i)+dt/grid%v(i)*(flowRhou(:,i)+presF(:,i)+viscF(:,i))
    end forall
    y=reshape(rhou(:,1:grid%nC),[grid%nC*DIMS])
    momentumRHS=merge(1,0,any(ieee_is_nan(y).or.(.not.ieee_is_finite(y))))
  end function
  
  !> residual function of the pressure equation
  function pressureRHS(oldVector,newVector,dat)
    use iso_c_binding
    use ieee_arithmetic
    use modNumerics
    use modGradient
    use modDiffusion
    use modRhieChow
    type(C_PTR),value::oldVector !< old N_Vector
    type(C_PTR),value::newVector !< new N_Vector
    type(C_PTR),value::dat !< optional user data object
    integer(C_INT)::pressureRHS !< error code
    double precision,pointer::x(:) !< fortran pointer associated with oldVector
    double precision,pointer::y(:) !< fortran pointer associated with newVector
    double precision,allocatable,save::tmpFlowRho(:)
    double precision,parameter::R_AIR=287.058d0 ! TODO fluid data
    ! TODO R*T under-relaxation
    
    call associateVector(oldVector,x)
    call associateVector(newVector,y)
    
    p(1:grid%nC)=x(1:grid%nC)
    call setBC()
    call findGrad(grid,p,gradP)
    if(.not.allocated(tmpFlowRho))then
      allocate(tmpFlowRho(grid%nC))
    else if(size(tmpFlowRho)/=grid%nC)then
      deallocate(tmpFlowRho)
      allocate(tmpFlowRho(grid%nC))
    end if
    tmpFlowRho(:)=flowRho(:)
    call addRhieChow(grid,rho,p,gradP,rho,dt,tmpFlowRho)
    call findDiff(grid,p,[(1d0,i=1,grid%nC)],laP)
    forall(i=1:grid%nC)
      y(i)=p(i)-R_AIR*temp(i)*dt**2/grid%v(i)*(laP(i)-laP1(i))-R_AIR*temp(i)/R_AIR/temp1(i)*p1(i)&
      &    +R_AIR*temp(i)*(rho1(i)-rho0(i))-R_AIR*temp(i)*dt/grid%v(i)*tmpFlowRho(i)
    end forall
    pressureRHS=merge(1,0,any(ieee_is_nan(y).or.(.not.ieee_is_finite(y))))
  end function
  
  !> RHS of the density equation in the form of a fix point problem
  function densityRHS(oldVector,newVector,dat)
    use iso_c_binding
    use ieee_arithmetic
    use modNumerics
    use modAdvection
    use modRhieChow
    type(C_PTR),value::oldVector !< old N_Vector
    type(C_PTR),value::newVector !< new N_Vector
    type(C_PTR),value::dat !< optional user data object
    integer(C_INT)::densityRHS !< error code
    double precision,pointer::x(:) !< fortran pointer associated with oldVector
    double precision,pointer::y(:) !< fortran pointer associated with newVector
    
    call associateVector(oldVector,x)
    call associateVector(newVector,y)
    
    rho(1:grid%nC)=x(1:grid%nC)
    forall(i=1:grid%nE)
      rhou(:,i)=rho(i)*u(:,i)
    end forall
    call setBC()
    call findAdv(grid,rho,rhou,flowRho)
    call addRhieChow(grid,rho,p,gradP,rho,dt,flowRho)
    forall(i=1:grid%nC)
      y(i)=rho0(i)+dt/grid%v(i)*flowRho(i)
    end forall
    densityRHS=merge(1,0,any(ieee_is_nan(y).or.(.not.ieee_is_finite(y))))
  end function
  
  !> RHS of the energy equation in the form of a fix point problem
  function energyRHS(oldVector,newVector,dat)
    use iso_c_binding
    use ieee_arithmetic
    use modNumerics
    use modAdvection
    use modDiffusion
    use modRhieChow
    type(C_PTR),value::oldVector !< old N_Vector
    type(C_PTR),value::newVector !< new N_Vector
    type(C_PTR),value::dat !< optional user data object
    integer(C_INT)::energyRHS !< error code
    double precision,pointer::x(:) !< fortran pointer associated with oldVector
    double precision,pointer::y(:) !< fortran pointer associated with newVector
    double precision,parameter::R_AIR=287.058d0 ! TODO fluid data
    
    call associateVector(oldVector,x)
    call associateVector(newVector,y)
    
    rhoE(1:grid%nC)=x(1:grid%nC)
    forall(i=1:grid%nE)
      temp(i)=(rhoE(i)/rho(i)-0.5d0*dot_product(u(:,i),u(:,i)))/(1d0/(1.4d0-1d0))/R_AIR
    end forall
    call setBC()
    forall(i=1:grid%nE)
      fluxRhoH(:,i)=(rhoE(i)+p(i))*u(:,i)
    end forall
    call findAdv(grid,rhoE+p,fluxRhoH,flowRhoH)
    call addRhieChow(grid,rhoE+p,p,gradP,rho,dt,flowRhoH)
    call findDiff(grid,temp,cond,condQ)
    forall(i=1:grid%nC)
      y(i)=rhoE0(i)+dt/grid%v(i)*(flowRhoH(i)+condQ(i)) ! TODO add viscous heating
    end forall
    energyRHS=merge(1,0,any(ieee_is_nan(y).or.(.not.ieee_is_finite(y))))
  end function
  
end module