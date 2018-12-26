!----------------------------------------------------------------------------- best with 100 columns

!> environment for the pressure-based coupled solver
module modPbc
  use modPolyFvGrid
  use modCondition
  use modUDF
  use modNumerics
  use iso_c_binding
  
  public
  
  integer,parameter::DIMS=3 !< three dimensions
  
  integer,parameter::MAXIT_PBC=20 !< max number of momentum and pressure equation iterations
  integer,parameter::MAXIT_ENERGY=20 !< max number of energy equation iterations
  integer,parameter::MAXIT_OUTER=10 !< max number of outer iterations
  
  integer,parameter::BC_WALL_TEMP=0 !< wall boundary with prescribed temperature
  integer,parameter::BC_WALL_TEMP_UDF=5 !< wall boundary with prescribed temperature by UDF
  integer,parameter::BC_WALL_FLUX=1 !< wall boundary with prescribed heat flux
  integer,parameter::BC_WALL_FLUX_UDF=6 !< wall boundary with prescribed heat flux by UDF
  integer,parameter::BC_IN_STATIC=10 !< inflow boundary with static properties
  integer,parameter::BC_IN_STATIC_UDF=15 !< inflow boundary with static properties by UDF
  integer,parameter::BC_IN_TOTAL=11 !< inflow boundary with total properties
  integer,parameter::BC_IN_TOTAL_UDF=16 !< inflow boundary with total properties by UDF
  integer,parameter::BC_OUT=20 !< outflow boundary
  integer,parameter::BC_OUT_UDF=25 !< outflow boundary by UDF
  integer,parameter::BC_FAR=30 !< far-field boundary
  integer,parameter::BC_FAR_UDF=35 !< far-field boundary by UDF
  
  type(polyFvGrid)::grid !< computational grid
  type(condTab)::bc !< boundary conditions
  type(UDFTab)::udf !< UDF
  
  integer,allocatable::iBC(:) !< indexes boundary conditions
  
  double precision::t !< time
  double precision::tFinal !< final time
  double precision::dt !< time step size [s]
  double precision::tInt !< time interval of output
  double precision::tNext !< time for next output
  integer::iOut !< index of output
  
  double precision::Rgas !< specific gas constant [J/kg/K]
  double precision::gamm !< gamma of gas
  
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
  double precision,allocatable::viscF(:,:) !< viscous force applied on cell [N]
  double precision,allocatable::presF(:,:) !< pressure force applied on cell [N]
  double precision,allocatable::condQ(:) !< heat conduction into cell [W]
  
  ! scales
  double precision,allocatable::rhoScale !< density scale [kg/m^3]
  double precision,allocatable::rhouScale !< momentum scale [kg/s/m^2]
  double precision,allocatable::rhoEScale !< total energy scale [J/m^3]
  double precision,allocatable::pScale !< pressure scale [Pa]
  double precision,allocatable::uScale !< velocity scale [m/s]
  double precision,allocatable::tempScale !< temperature scale [K]
  
  ! data for algebraic solver
  type(NewtonKrylov)::pbcEq !< momentum and pressure equations (isothermal) solved by Newton-GMRES
  type(fixPt)::energyEq !< energy equation as a fix point problem
  integer::nItPBC,nItEnergy,nOuter !< number of iterations
  
contains
  
  !> initialize the simulation
  subroutine init()
    use modFileIO
    integer,parameter::FID=10
    double precision::pInit,uInit(DIMS),TInit,pE(DIMS),tmpD
    double precision,allocatable::YInit(:)
    integer::udfIc,udfBc,iUdf(5)
    
    ! read inputs
    open(FID,file='grid.msh',action='read')
      call readGMSH(FID,grid)
    close(FID)
    call grid%up()
    open(FID,file='bc',action='read')
      call readCondTab(FID,bc)
    close(FID)
    call mapCondTab(grid,bc,iBC)
    if(any(bc%t(:)==BC_WALL_TEMP_UDF).or.&
    &  any(bc%t(:)==BC_WALL_FLUX_UDF).or.&
    &  any(bc%t(:)==BC_IN_STATIC_UDF).or.&
    &  any(bc%t(:)==BC_IN_TOTAL_UDF).or.&
    &  any(bc%t(:)==BC_OUT_UDF).or.&
    &  any(bc%t(:)==BC_FAR_UDF))then
      udfBc=1
    else
      udfBc=0
    end if
    open(FID,file='ic',action='read')
      read(FID,*)udfIC
      if(udfIc==0)then
        read(FID,*)pInit
        read(FID,*)TInit
        read(FID,*)uInit(1)
        read(FID,*)uInit(2)
        read(FID,*)uInit(3)
      else
        read(FID,*)tmpD
        iUdf(1)=int(tmpD)
        read(FID,*)tmpD
        iUdf(2)=int(tmpD)
        read(FID,*)tmpD
        iUdf(3)=int(tmpD)
        read(FID,*)tmpD
        iUdf(4)=int(tmpD)
        read(FID,*)tmpD
        iUdf(5)=int(tmpD)
      end if
    close(FID)
    open(FID,file='sim',action='read')
      read(FID,*)tFinal
      read(FID,*)dt
      read(FID,*)tInt
    close(FID)
    open(FID,file='fl',action='read')
      read(FID,*)Rgas
      read(FID,*)gamm
    close(FID)
    if(udfIc/=0.or.udfBc/=0)then
      open(FID,file='udf',action='read')
        call readUDFTab(FID,udf)
      close(FID)
    end if
    ! initialize algebraic solver
    call pbcEq%init(grid%nC*(DIMS+1),pbcRHS)
    call pbcEq%setMaxIt(MAXIT_PBC)
    call energyEq%init(grid%nC,energyRHS,maa=MAXIT_ENERGY)
    call energyEq%setMaxIt(MAXIT_ENERGY)
    ! work space and initial state
    allocate(rho(grid%nE))
    allocate(rho0(grid%nE))
    allocate(rhou(DIMS,grid%nE))
    allocate(rhou0(DIMS,grid%nE))
    allocate(rhoE(grid%nE))
    allocate(rhoE0(grid%nE))
    allocate(p(grid%nE))
    allocate(p0(grid%nE))
    allocate(u(DIMS,grid%nE))
    allocate(u0(DIMS,grid%nE))
    allocate(temp(grid%nE))
    allocate(temp0(grid%nE))
    allocate(Y(1,grid%nE))! FIXME fix nSp
    allocate(Y0(1,grid%nE))! FIXME fix nSp
    allocate(YInit(1))! FIXME fix nSp
    allocate(visc(grid%nE))
    allocate(cond(grid%nE))
    allocate(gradP(DIMS,grid%nC))
    allocate(viscF(DIMS,grid%nC))
    allocate(presF(DIMS,grid%nC))
    allocate(condQ(grid%nC))
    do i=1,grid%nE
      if(udfIc/=0)then
        pE(:)=grid%p(:,i)
        pInit=udf%eval(iUdf(1),pE,0d0)
        TInit=udf%eval(iUdf(2),pE,0d0)
        uInit(1)=udf%eval(iUdf(3),pE,0d0)
        uInit(2)=udf%eval(iUdf(4),pE,0d0)
        uInit(3)=udf%eval(iUdf(5),pE,0d0)
      end if
      ! FIXME remove below
      YInit=[1d0]
      ! FIXME remove above
      p(i)=pInit
      u(:,i)=uInit(:)
      temp(i)=TInit
      Y(:,i)=YInit(:)
    end do
    call recoverState(p,u,temp,Y,rho,rhou,rhoE)
    t=0d0
    tNext=tInt
    iOut=0
    deallocate(YInit)
  end subroutine
  
  !> clear the simulation environment
  subroutine clear()
    call grid%clear()
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
      po(i)=(gamm-1d0)*(rhoEi(i)-0.5d0*dot_product(rhoui(:,i),rhoui(:,i))/rhoi(i))
      tempo(i)=po(i)/rhoi(i)/Rgas
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
      rhoo(i)=pi(i)/Rgas/tempi(i)
      rhouo(:,i)=rhoo(i)*ui(:,i)
      rhoEo(i)=rhoo(i)*(1d0/(gamm-1d0)*Rgas*tempi(i)+0.5d0*dot_product(ui(:,i),ui(:,i)))
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
  
  !> calculate time step size, scaling vectors and initial solution vector
  subroutine preSolve()
    double precision,parameter::CFL_ACCOUSTIC=10d0
    double precision,parameter::CFL_FLOW=0.5d0
    double precision,parameter::CFL_DIFFUSION=0.5d0
    
    ! TODO calculate dt, cantera sound speed
    dt=min(dt,minval(CFL_ACCOUSTIC*grid%v(:)**(1d0/3d0)&
    &                     /sqrt(gamm*p(1:grid%nC)/rho(1:grid%nC))))
    dt=min(dt,minval(CFL_FLOW*grid%v(:)**(1d0/3d0)&
    &                /norm2(u(:,1:grid%nC),1)))
    dt=min(dt,minval(CFL_DIFFUSION*grid%v(:)**(2d0/3d0)&
    &                /(visc(1:grid%nC)/rho(1:grid%nC))))
    dt=1d-4
    
    write(*,'(a,g12.6,a,g12.6)')'[i] starting step, t: ',t,' dt: ',dt
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
        end if
      end if
    end do
  end subroutine
  
  !> residual function of the momentum and pressure equations
  function pbcRHS(oldVector,newVector,dat)
    use iso_c_binding
    use ieee_arithmetic
    use modNumerics
    use modGradient
    use modAdvection
    use modDiffusion
    use modRhieChow
    type(C_PTR),value::oldVector !< old N_Vector
    type(C_PTR),value::newVector !< new N_Vector
    type(C_PTR),value::dat !< optional user data object
    integer(C_INT)::pbcRHS !< error code
    double precision,pointer::x(:) !< fortran pointer associated with oldVector
    double precision,pointer::y(:) !< fortran pointer associated with newVector
    
    call associateVector(oldVector,x)
    call associateVector(newVector,y)
    
    pbcRHS=merge(1,0,any(ieee_is_nan(y).or.(.not.ieee_is_finite(y))))
    if(c_associated(dat))then
    end if
  end function
  
  !> RHS of the energy equation in the form of a fix point problem
  function energyRHS(oldVector,newVector,dat)
    use iso_c_binding
    use ieee_arithmetic
    use modNumerics
    use modAdvection
    use modDiffusion
    type(C_PTR),value::oldVector !< old N_Vector
    type(C_PTR),value::newVector !< new N_Vector
    type(C_PTR),value::dat !< optional user data object
    integer(C_INT)::energyRHS !< error code
    double precision,pointer::x(:) !< fortran pointer associated with oldVector
    double precision,pointer::y(:) !< fortran pointer associated with newVector
    
    call associateVector(oldVector,x)
    call associateVector(newVector,y)
    
    !rhoE(1:grid%nC)=x(1:grid%nC)
    !forall(i=1:grid%nC)
    !  temp(i)=(rhoE(i)/rho(i)-0.5d0*dot_product(u(:,i),u(:,i)))/(1d0/(gamm-1d0))/Rgas
    !end forall
    !call setBC()
    !call findVarFlow(grid,(rhoE+p)/rho,flowRho,flowRhoH)
    !call findAdv(grid,flowRhoH,advRhoH)
    !call findDiff(grid,temp,cond,condQ)
    !forall(i=1:grid%nC)
    !  y(i)=rhoE0(i)+dt/grid%v(i)*(advRhoH(i)+condQ(i)) ! TODO add viscous heating
    !end forall
    energyRHS=merge(1,0,any(ieee_is_nan(y).or.(.not.ieee_is_finite(y))))
    if(c_associated(dat))then
    end if
  end function
  
end module
