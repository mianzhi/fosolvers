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
  
  integer,parameter::MAXIT_PBC=50 !< max number of momentum and pressure equation iterations
  integer,parameter::MAXIT_ENERGY=50 !< max number of energy equation iterations
  integer,parameter::MAXIT_OUTER=20 !< max number of outer iterations
  integer,parameter::MAXIT_FULL=50 !< max number of full NS equation iterations
  
  integer,parameter::BC_WALL_TEMP=1 !< wall boundary with prescribed temperature
  integer,parameter::BC_WALL_TEMP_UDF=101 !< wall boundary with prescribed temperature by UDF
  integer,parameter::BC_WALL_SLIP=2 !< adiabatic slip wall
  integer,parameter::BC_IN_STATIC=10 !< inflow boundary with static properties
  integer,parameter::BC_IN_STATIC_UDF=110 !< inflow boundary with static properties by UDF
  integer,parameter::BC_IN_TOTAL=11 !< inflow boundary with total properties
  integer,parameter::BC_IN_TOTAL_UDF=111 !< inflow boundary with total properties by UDF
  integer,parameter::BC_OUT=20 !< outflow boundary
  integer,parameter::BC_OUT_UDF=120 !< outflow boundary by UDF
  integer,parameter::BC_FAR=30 !< far-field boundary
  integer,parameter::BC_FAR_UDF=130 !< far-field boundary by UDF
  
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
  integer::nRetry !< number of retry for the current time step
  logical::needRetry !< if retry of the current time step is needed
  
  double precision::Rgas !< specific gas constant [J/kg/K]
  double precision::gamm !< gamma of gas
  
  double precision,allocatable::rho(:) !< density [kg/m^3]
  double precision,allocatable::rhou(:,:) !< momentum [kg/s/m^2]
  double precision,allocatable::rhoE(:) !< total energy [J/m^3]
  double precision,allocatable::p(:) !< pressure [Pa]
  double precision,allocatable::u(:,:) !< velocity [m/s]
  double precision,allocatable::temp(:) !< temperature [K]
  double precision,allocatable::massFrac(:,:) !< mass fraction of species
  
  ! state at the beginning of a time step
  double precision,allocatable::rho0(:),rhou0(:,:),rhoE0(:),p0(:),u0(:,:),temp0(:),massFrac0(:,:)
  
  ! state at the beginning of an outer iteration
  double precision,allocatable::rho1(:)
  double precision,allocatable::p1(:)
  
  double precision,allocatable::visc(:) !< viscosity [Pa*s]
  double precision,allocatable::cond(:) !< thermal conductivity [W/m/K]
  double precision,allocatable::flowRho(:) !< mass flow through pair [kg/s]
  double precision,allocatable::flowRhou(:,:) !< momentum flow through pair [kg*m/s^2]
  double precision,allocatable::flowRhoH(:) !< enthalpy flow through pair [W]
  double precision,allocatable::advRho(:) !< mass advection into cell [kg/s]
  double precision,allocatable::advRhou(:,:) !< momentum advection into cell [kg*m/s^2]
  double precision,allocatable::advRhoH(:) !< enthalpy advection into cell [W]
  double precision,allocatable::gradU(:,:,:) !< gradient of velocity [1/s]
  double precision,allocatable::gradP(:,:) !< gradient of pressure [Pa/m]
  double precision,allocatable::gradTemp(:,:) !< gradient of temp [K/m]
  double precision,allocatable::viscF(:,:) !< viscous force applied on cell [N]
  double precision,allocatable::presF(:,:) !< pressure force applied on cell [N]
  double precision,allocatable::condQ(:) !< heat conduction into cell [W]
  
  ! scales and scaling factors
  double precision,allocatable::rhoScale !< density scale [kg/m^3]
  double precision,allocatable::rhouScale !< momentum scale [kg/s/m^2]
  double precision,allocatable::rhoEScale !< total energy scale [J/m^3]
  double precision,allocatable::pScale !< pressure scale [Pa]
  double precision,allocatable::uScale !< velocity scale [m/s]
  double precision,allocatable::tempScale !< temperature scale [K]
  double precision,allocatable::pbcXFact(:) !< pbc solution scaling factor vector
  double precision,allocatable::pbcRFact(:) !< pbc residual scaling factor vector
  double precision,allocatable::fullXFact(:) !< full system solution scaling factor vector
  double precision,allocatable::fullRFact(:) !< full system residual scaling factor vector
  
  ! data for algebraic solver
  type(NewtonKrylov)::pbcEq !< momentum and pressure equations (isothermal) solved by Newton-GMRES
  type(fixPt)::energyEq !< energy equation as a fix point problem
  type(NewtonKrylov)::fullEq !< full NS equations solved by Newton-GMRES
  integer::nItPBC,nItEnergy,nItOuter,nItFull !< number of iterations
  
contains
  
  !> initialize the simulation
  subroutine init()
    use modFileIO
    integer,parameter::FID=10
    double precision::pInit,uInit(DIMS),TInit,pE(DIMS),tmpD
    double precision,allocatable::massFracInit(:)
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
    call pbcEq%init(grid%nC*(DIMS+1),pbcRHS,maxl=MAXIT_PBC)
    call pbcEq%setMaxIt(MAXIT_PBC)
    call energyEq%init(grid%nC,energyRHS,maa=MAXIT_ENERGY)
    call energyEq%setMaxIt(MAXIT_ENERGY)
    call fullEq%init(grid%nC*(DIMS+2),fullRHS,maxl=MAXIT_FULL)
    call fullEq%setMaxIt(MAXIT_FULL)
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
    allocate(massFrac(1,grid%nE))! FIXME fix nSp
    allocate(massFrac0(1,grid%nE))! FIXME fix nSp
    allocate(massFracInit(1))! FIXME fix nSp
    allocate(visc(grid%nE))
    allocate(cond(grid%nE))
    allocate(flowRho(grid%nP))
    allocate(flowRhou(DIMS,grid%nP))
    allocate(flowRhoH(grid%nP))
    allocate(advRho(grid%nC))
    allocate(advRhou(DIMS,grid%nC))
    allocate(advRhoH(grid%nC))
    allocate(gradU(DIMS,DIMS,grid%nC))
    allocate(gradP(DIMS,grid%nC))
    allocate(gradTemp(DIMS,grid%nC))
    allocate(viscF(DIMS,grid%nC))
    allocate(presF(DIMS,grid%nC))
    allocate(condQ(grid%nC))
    allocate(pbcXFact((DIMS+1)*grid%nC))
    allocate(pbcRFact((DIMS+1)*grid%nC))
    allocate(fullXFact((DIMS+2)*grid%nC))
    allocate(fullRFact((DIMS+2)*grid%nC))
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
      massFracInit=[1d0]
      ! FIXME remove above
      p(i)=pInit
      u(:,i)=uInit(:)
      temp(i)=TInit
      massFrac(:,i)=massFracInit(:)
    end do
    call recoverState(p,u,temp,massFrac,rho,rhou,rhoE)
    t=0d0
    tNext=tInt
    iOut=0
    deallocate(massFracInit)
  end subroutine
  
  !> clear the simulation environment
  subroutine clear()
    call grid%clear()
  end subroutine
  
  !> record {rho,rhoU,rhoE,p,u,temp,massFrac} in {rho0,rhoU0,rhoE0,p0,u0,temp0,massFrac0}
  subroutine recordState0()
    rho0(:)=rho(:)
    rhou0(:,:)=rhou(:,:)
    rhoE0(:)=rhoE(:)
    p0(:)=p(:)
    u0(:,:)=u(:,:)
    temp0(:)=temp(:)
    massFrac0(:,:)=massFrac(:,:)
  end subroutine
  
  !> load {rho,rhoU,rhoE,p,u,temp,massFrac} from {rho0,rhoU0,rhoE0,p0,u0,temp0,massFrac0}
  subroutine loadState0()
    rho(:)=rho0(:)
    rhou(:,:)=rhou0(:,:)
    rhoE(:)=rhoE0(:)
    p(:)=p0(:)
    u(:,:)=u0(:,:)
    temp(:)=temp0(:)
    massFrac(:,:)=massFrac0(:,:)
  end subroutine
  
  !> derive primitive state {p,u,T} from conserved state {rho,rhou,rhoE}
  subroutine deriveState(rhoi,rhoui,rhoEi,massFraci,po,uo,tempo)
    double precision,intent(in)::rhoi(:) !< density
    double precision,intent(in)::rhoui(:,:) !< momentum
    double precision,intent(in)::rhoEi(:) !< total energy
    double precision,intent(in)::massFraci(:,:) !< mass fraction of species
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
  subroutine recoverState(pi,ui,tempi,massFraci,rhoo,rhouo,rhoEo)
    double precision,intent(in)::pi(:) !< pressure
    double precision,intent(in)::ui(:,:) !< velocity
    double precision,intent(in)::tempi(:) !< temperature
    double precision,intent(in)::massFraci(:,:) !< mass fraction of species
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
    double precision::eps
    
    open(FID,file=trim(fName),action='write')
    call writeVTK(FID,grid)
    call writeVTK(FID,grid,E_DATA)
    call writeVTK(FID,'geoID',dble(grid%gid))
    call writeVTK(FID,'geoType',[(merge(0d0,1d0,i<=grid%nC),i=1,grid%nE)])
    if(iOut==0)then
      call writeVTK(FID,'density',rho)
      call writeVTK(FID,'velocity',u)
      call writeVTK(FID,'pressure',p)
      call writeVTK(FID,'temperature',temp)
    else
      eps=(t-tNext)/dt
      call writeVTK(FID,'density',eps*rho0+(1d0-eps)*rho)
      call writeVTK(FID,'velocity',eps*u0+(1d0-eps)*u)
      call writeVTK(FID,'pressure',eps*p0+(1d0-eps)*p)
      call writeVTK(FID,'temperature',eps*temp0+(1d0-eps)*temp)
    end if
    close(FID)
  end subroutine
  
  !> calculate time step size, scaling vectors and initial solution vector
  subroutine preSolve()
    double precision,parameter::CFL_ACCOUSTIC=50d0
    double precision,parameter::CFL_FLOW=0.5d0
    double precision,parameter::CFL_DIFFUSION=0.5d0
    
    ! TODO calculate dt, cantera sound speed
    dt=min(dt,minval(CFL_ACCOUSTIC*grid%v(:)**(1d0/3d0)&
    &                     /sqrt(gamm*p(1:grid%nC)/rho(1:grid%nC))))
    dt=min(dt,minval(CFL_FLOW*grid%v(:)**(1d0/3d0)&
    &                /norm2(u(:,1:grid%nC),1)))
    dt=min(dt,minval(CFL_DIFFUSION*grid%v(:)**(2d0/3d0)&
    &                /(visc(1:grid%nC)/rho(1:grid%nC))))
    if(nRetry>0)then
      dt=dt/2d0**nRetry
      write(*,'(a,i2,a,g12.6)')'[i] starting retry No. ',nRetry,' at t: ',t
    end if
    
    ! solution and residual scales
    pScale=max(maxval(p),maxval(0.5d0*rho*norm2(u,1)**2))
    uScale=max(maxval(norm2(u,1)),sqrt(2d0*(maxval(p)-minval(p))/minval(rho)),0.01d0)
    tempScale=maxval(temp)
    rhoScale=maxval(rho)
    rhouScale=rhoScale*uScale
    rhoEScale=max(rhoScale*Rgas*tempScale,0.5d0*rhoScale*uScale**2)
    forall(i=1:grid%nC)
      pbcXFact((DIMS+1)*(i-1)+1:(DIMS+1)*(i-1)+DIMS)=1d0/uScale
      pbcXFact((DIMS+1)*(i-1)+(DIMS+1))=1d0/pScale
      pbcRFact((DIMS+1)*(i-1)+1:(DIMS+1)*(i-1)+DIMS)=1d0/rhouScale
      pbcRFact((DIMS+1)*(i-1)+(DIMS+1))=1d0/rhoScale
      fullXFact((DIMS+2)*(i-1)+1:(DIMS+2)*(i-1)+DIMS)=1d0/uScale
      fullXFact((DIMS+2)*(i-1)+(DIMS+1))=1d0/pScale
      fullXFact((DIMS+2)*(i-1)+(DIMS+2))=1d0/tempScale
      fullRFact((DIMS+2)*(i-1)+1:(DIMS+2)*(i-1)+DIMS)=1d0/rhouScale
      fullRFact((DIMS+2)*(i-1)+(DIMS+1))=1d0/rhoScale
      fullRFact((DIMS+2)*(i-1)+(DIMS+2))=1d0/rhoEScale
    end forall
    call pbcEq%setScale(pbcXFact,pbcRFact)
    call fullEq%setScale(fullXFact,fullRfact)
    
    nItPBC=0
    nItEnergy=0
    nItFull=0
    
    write(*,'(a,g12.6,a,g12.6)')'[i] starting step, t: ',t,' dt: ',dt
  end subroutine
  
  !> guess next time step size, etc.
  subroutine postSolve()
    double precision,parameter::MIN_EFFORT=0.3d0 !< increase dt if "effort" lower than lower limit
    double precision,parameter::MAX_EFFORT=0.5d0 !< reduce dt if "effort" higher than higher limit
    double precision,parameter::MAX_GROWTH=1.2d0 !< new dt is no more than 120% of the current dt
    double precision::effort !< the "effort"
  
    ! adjust time step size to maintain effort
    effort=maxval(dble([nItPBC,nItEnergy,nItOuter])/dble([MAXIT_PBC,MAXIT_ENERGY,MAXIT_OUTER]))
    if(effort<MIN_EFFORT)then
      dt=dt*min(MIN_EFFORT/effort,MAX_GROWTH)
    else if(effort>MAX_EFFORT)then
      dt=dt*MAX_EFFORT/effort
    end if
    write(*,'(a,i2,a,i2,a,i2,a,i2,a)')'[i] finished step, nIt[PBC,Energy,Outer]: [',&
    &                                 nItPBC,',',nItEnergy,',',nItOuter,']'
  end subroutine
  
  !> set the boundary conditions
  subroutine setBC()
    double precision::pGst,TGst,uGst(DIMS),rhoGst,rhouGst(DIMS),rhoEGst,Mach,pP(DIMS),Tt,pt,&
    &                 un(DIMS),ut(DIMS),ui(DIMS),charI,charO
    
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(n>grid%nC)then
        ! default adiabatic wall boundary
        pGst=p(m)
        uGst(:)=-u(:,m)
        TGst=temp(m)
        rhoGst=rho(m)
        rhouGst=-rhou(:,m)
        rhoEGst=rhoE(m)
        if(iBC(n)>0)then
          select case(bc%t(iBC(n)))
          case(BC_WALL_TEMP,BC_WALL_TEMP_UDF) ! wall temperature boundary
            ! TODO implement this
          case(BC_WALL_SLIP) ! adiabatic slip wall boundary
            pGst=p(m)
            TGst=temp(m)
            uGst(:)=u(:,m)-2d0*grid%normP(:,i)*dot_product(grid%normP(:,i),u(:,m))
          case(BC_IN_STATIC,BC_IN_STATIC_UDF) ! inflow boundary with static properties
            if(bc%t(iBC(n))==BC_IN_STATIC)then
              pGst=bc%p(1,iBC(n))
              TGst=bc%p(2,iBC(n))
              uGst(:)=bc%p(3:5,iBC(n))
            else
              pP(:)=grid%pP(:,i)
              pGst=udf%eval(int(bc%p(1,iBC(n))),pP,t)
              TGst=udf%eval(int(bc%p(2,iBC(n))),pP,t)
              uGst(1)=udf%eval(int(bc%p(3,iBC(n))),pP,t)
              uGst(2)=udf%eval(int(bc%p(4,iBC(n))),pP,t)
              uGst(3)=udf%eval(int(bc%p(5,iBC(n))),pP,t)
            end if
            Mach=norm2(uGst)/sqrt(gamm*Rgas*TGst)
            if(Mach<1d0)then
              if(dot_product(uGst,uGst)>0d0)then
                uGst(:)=dot_product(u(:,m),uGst(:))*uGst/dot_product(uGst,uGst)
              else
                uGst(:)=u(:,m)
              end if
            end if
          case(BC_IN_TOTAL,BC_IN_TOTAL_UDF) ! inflow boundary with total properties
            if(bc%t(iBC(n))==BC_IN_TOTAL)then
              pt=bc%p(1,iBC(n))
              Tt=bc%p(2,iBC(n))
              uGst=bc%p(3:5,iBC(n))
            else
              pP(:)=grid%pP(:,i)
              pt=udf%eval(int(bc%p(1,iBC(n))),pP,t)
              Tt=udf%eval(int(bc%p(2,iBC(n))),pP,t)
              uGst(1)=udf%eval(int(bc%p(3,iBC(n))),pP,t)
              uGst(2)=udf%eval(int(bc%p(4,iBC(n))),pP,t)
              uGst(3)=udf%eval(int(bc%p(5,iBC(n))),pP,t)
            end if
            TGst=Tt-0.5d0*(gamm-1d0)/gamm/Rgas*dot_product(uGst,uGst)
            Mach=norm2(uGst)/sqrt(gamm*Rgas*TGst)
            if(Mach<1d0)then
              if(dot_product(uGst,uGst)>0d0)then
                uGst(:)=dot_product(u(:,m),uGst(:))*uGst/dot_product(uGst,uGst)
              else
                uGst(:)=u(:,m)
              end if
              TGst=Tt-0.5d0*(gamm-1d0)/gamm/Rgas*dot_product(uGst,uGst)
              Mach=norm2(uGst)/sqrt(gamm*Rgas*TGst)
            end if
            pGst=pt*(1d0+0.5d0*(gamm-1d0)*Mach**2)**(-gamm/(gamm-1d0))
          case(BC_OUT,BC_OUT_UDF) ! outflow boundary
            if(bc%t(iBC(n))==BC_OUT)then
              pGst=bc%p(1,iBC(n))
            else
              pP(:)=grid%pP(:,i)
              pGst=udf%eval(int(bc%p(1,iBC(n))),pP,t)
            end if
            uGst(:)=u(:,m)
            TGst=temp(m)
          case(BC_FAR,BC_FAR_UDF) ! far-field boundary
            if(bc%t(iBC(n))==BC_FAR)then
              pGst=bc%p(1,iBC(n))
              TGst=bc%p(2,iBC(n))
              uGst(:)=bc%p(3:5,iBC(n))
            else
              pP(:)=grid%pP(:,i)
              pGst=udf%eval(int(bc%p(1,iBC(n))),pP,t)
              TGst=udf%eval(int(bc%p(2,iBC(n))),pP,t)
              uGst(1)=udf%eval(int(bc%p(3,iBC(n))),pP,t)
              uGst(2)=udf%eval(int(bc%p(4,iBC(n))),pP,t)
              uGst(3)=udf%eval(int(bc%p(5,iBC(n))),pP,t)
            end if
            charO=dot_product(uGst,grid%normP(:,i))-2d0*sqrt(gamm*Rgas*TGst)/(gamm-1d0)
            ui(:)=u(:,m)
            charI=dot_product(ui,grid%normP(:,i))+2d0*sqrt(gamm*Rgas*temp(m))/(gamm-1d0)
            un(:)=0.5d0*grid%normP(:,i)*(charI+charO)
            if(dot_product(un(:),grid%normP(:,i))>0d0)then
              ut(:)=ui(:)-grid%normP(:,i)*dot_product(grid%normP(:,i),ui(:))
            else
              ut(:)=uGst(:)-un
            end if
            uGst(:)=un(:)+ut(:)
          case default
          end select
          rhoGst=pGst/Rgas/TGst
          rhouGst(:)=rhoGst*uGst(:)
          rhoEGst=rhoGst*(1d0/(gamm-1d0)*Rgas*TGst+0.5d0*dot_product(uGst,uGst))
        end if
        ! apply ghost values
        p(n)=pGst
        u(:,n)=uGst(:)
        temp(n)=TGst
        rho(n)=rhoGst
        rhou(:,n)=rhouGst(:)
        rhoE(n)=rhoEGst
      end if
    end do
  end subroutine
  
  !> solve the coupled momentum and pressure equations while keeping isentropic
  subroutine solvePBC()
    use modGradient
    double precision,allocatable,save::x(:)
    integer::info
    
    if(.not.allocated(x))then
      allocate(x((DIMS+1)*grid%nC))
    else if(size(x)/=(DIMS+1)*grid%nC)then
      deallocate(x)
      allocate(x((DIMS+1)*grid%nC))
    end if
    forall(i=1:grid%nC)
      x((DIMS+1)*(i-1)+1:(DIMS+1)*(i-1)+DIMS)=u(:,i)
      x((DIMS+1)*(i-1)+(DIMS+1))=p(i)
    end forall
    call setBC()
    call findGrad(grid,p,gradP)
    call findGrad(grid,u,gradU)
    call pbcEq%solve(x,info=info)
    needRetry=info<0
    call pbcEq%getNIt(n)
    nItPBC=max(nItPBC,n)
    forall(i=1:grid%nC)
      u(:,i)=x((DIMS+1)*(i-1)+1:(DIMS+1)*(i-1)+DIMS)
      p(i)=x((DIMS+1)*(i-1)+(DIMS+1))
    end forall
  end subroutine
  
  !> residual function of the momentum and pressure equations
  function pbcRHS(oldVector,newVector,dat)
    use iso_c_binding
    use ieee_arithmetic
    use modNumerics
    use modAdvection
    use modDiffusion
    use modPressure
    use modNewtonian
    type(C_PTR),value::oldVector !< old N_Vector
    type(C_PTR),value::newVector !< new N_Vector
    type(C_PTR),value::dat !< optional user data object
    integer(C_INT)::pbcRHS !< error code
    double precision,pointer::x(:) !< fortran pointer associated with oldVector
    double precision,pointer::y(:) !< fortran pointer associated with newVector
    
    call associateVector(oldVector,x)
    call associateVector(newVector,y)
    
    !$omp parallel workshare
    forall(i=1:grid%nC)
      u(:,i)=x((DIMS+1)*(i-1)+1:(DIMS+1)*(i-1)+DIMS)
      p(i)=x((DIMS+1)*(i-1)+(DIMS+1))
      rho(i)=p1(i)/Rgas/temp(i)*(p(i)/p1(i))**(1d0/gamm)
      rhou(:,i)=rho(i)*u(:,i)
    end forall
    !$omp end parallel workshare
    call setBC()
    call findPresForce(grid,p,gradP,presF)
    call findViscForce(grid,u,gradU,visc,viscF)
    call findMassFlow(grid,rho,u,p,presF,dt,flowRho)
    call findVarFlow(grid,u,flowRho,flowRhou)
    call findAdv(grid,flowRho,advRho)
    call findAdv(grid,flowRhou,advRhou)
    !$omp parallel workshare
    forall(i=1:grid%nC)
      y((DIMS+1)*(i-1)+1:(DIMS+1)*(i-1)+DIMS)=& ! momentum equation residual
      &  rhou(:,i)-rhou0(:,i)-dt/grid%v(i)*(advRhou(:,i)+presF(:,i)+viscF(:,i))
      y((DIMS+1)*(i-1)+(DIMS+1))=& ! pressure equation residual
      &  rho(i)-rho0(i)-dt/grid%v(i)*advRho(i)
    end forall
    !$omp end parallel workshare
    pbcRHS=merge(1,0,any(ieee_is_nan(y).or.(.not.ieee_is_finite(y))))
    if(c_associated(dat))then
    end if
  end function
  
  !> solve the energy equation
  subroutine solveEnergy()
    use modGradient
    double precision,allocatable,save::x(:)
    integer::info
    
    if(.not.allocated(x))then
      allocate(x(grid%nC))
    else if(size(x)/=grid%nC)then
      deallocate(x)
      allocate(x(grid%nC))
    end if
    x(1:grid%nC)=rhoE(1:grid%nC)
    call setBC()
    call findGrad(grid,temp,gradTemp)
    call energyEq%solve(x,info=info)
    needRetry=info<0
    call energyEq%getNIt(n)
    nItEnergy=max(nItEnergy,n)
    forall(i=1:grid%nC)
      rhoE(i)=x(i)
      temp(i)=(rhoE(i)/rho(i)-0.5d0*dot_product(u(:,i),u(:,i)))/(1d0/(gamm-1d0))/Rgas
    end forall
  end subroutine
  
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
    
    rhoE(1:grid%nC)=x(1:grid%nC)
    !$omp parallel workshare
    forall(i=1:grid%nC)
      temp(i)=(rhoE(i)/rho(i)-0.5d0*dot_product(u(:,i),u(:,i)))/(1d0/(gamm-1d0))/Rgas
    end forall
    !$omp end parallel workshare
    call setBC()
    call findVarFlow(grid,(rhoE+p)/rho,flowRho,flowRhoH)
    call findAdv(grid,flowRhoH,advRhoH)
    call findDiff(grid,temp,gradTemp,cond,condQ)
    !$omp parallel workshare
    forall(i=1:grid%nC)
      y(i)=rhoE0(i)+dt/grid%v(i)*(advRhoH(i)+condQ(i)) ! TODO add viscous heating
    end forall
    !$omp end parallel workshare
    energyRHS=merge(1,0,any(ieee_is_nan(y).or.(.not.ieee_is_finite(y))))
    if(c_associated(dat))then
    end if
  end function
  
  !> solve the full NS equations
  subroutine solveFull()
    double precision,allocatable,save::x(:)
    integer::info
    
    if(.not.allocated(x))then
      allocate(x((DIMS+2)*grid%nC))
    else if(size(x)/=(DIMS+2)*grid%nC)then
      deallocate(x)
      allocate(x((DIMS+2)*grid%nC))
    end if
    forall(i=1:grid%nC)
      x((DIMS+2)*(i-1)+1:(DIMS+2)*(i-1)+DIMS)=u(:,i)
      x((DIMS+2)*(i-1)+(DIMS+1))=p(i)
      x((DIMS+2)*(i-1)+(DIMS+2))=temp(i)
    end forall
    call setBC()
    call fullEq%solve(x,info=info)
    needRetry=info<0
    call fullEq%getNIt(nItFull)
    forall(i=1:grid%nC)
      u(:,i)=x((DIMS+2)*(i-1)+1:(DIMS+2)*(i-1)+DIMS)
      p(i)=x((DIMS+2)*(i-1)+(DIMS+1))
      temp(i)=x((DIMS+2)*(i-1)+(DIMS+2))
    end forall
  end subroutine
  
  !> residual function of the full NS equations
  function fullRHS(oldVector,newVector,dat)
    use iso_c_binding
    use ieee_arithmetic
    use modNumerics
    use modGradient
    use modAdvection
    use modDiffusion
    use modPressure
    use modNewtonian
    type(C_PTR),value::oldVector !< old N_Vector
    type(C_PTR),value::newVector !< new N_Vector
    type(C_PTR),value::dat !< optional user data object
    integer(C_INT)::fullRHS !< error code
    double precision,pointer::x(:) !< fortran pointer associated with oldVector
    double precision,pointer::y(:) !< fortran pointer associated with newVector
    
    call associateVector(oldVector,x)
    call associateVector(newVector,y)
    
    !$omp parallel workshare
    forall(i=1:grid%nC)
      u(:,i)=x((DIMS+2)*(i-1)+1:(DIMS+2)*(i-1)+DIMS)
      p(i)=x((DIMS+2)*(i-1)+(DIMS+1))
      temp(i)=x((DIMS+2)*(i-1)+(DIMS+2))
    end forall
    !$omp end parallel workshare
    call recoverState(p,u,temp,massFrac,rho,rhou,rhoE)
    call setBC()
    call findGrad(grid,p,gradP)
    call findGrad(grid,u,gradU)
    call findPresForce(grid,p,gradP,presF)
    call findViscForce(grid,u,gradU,visc,viscF)
    call findMassFlow(grid,rho,u,p,presF,dt,flowRho)
    call findVarFlow(grid,u,flowRho,flowRhou)
    call findVarFlow(grid,(rhoE+p)/rho,flowRho,flowRhoH)
    call findAdv(grid,flowRho,advRho)
    call findAdv(grid,flowRhou,advRhou)
    call findAdv(grid,flowRhoH,advRhoH)
    call findDiff(grid,temp,gradTemp,cond,condQ)
    !$omp parallel workshare
    forall(i=1:grid%nC)
      y((DIMS+2)*(i-1)+1:(DIMS+2)*(i-1)+DIMS)=& ! momentum equation residual
      &  rhou(:,i)-rhou0(:,i)-dt/grid%v(i)*(advRhou(:,i)+presF(:,i)+viscF(:,i))
      y((DIMS+2)*(i-1)+(DIMS+1))=& ! pressure equation residual
      &  rho(i)-rho0(i)-dt/grid%v(i)*advRho(i)
      y((DIMS+2)*(i-1)+(DIMS+2))=& ! energy equation residual
      &  rhoE(i)-rhoE0(i)-dt/grid%v(i)*(advRhoH(i)+condQ(i))
    end forall
    !$omp end parallel workshare
    fullRHS=merge(1,0,any(ieee_is_nan(y).or.(.not.ieee_is_finite(y))))
    if(c_associated(dat))then
    end if
  end function
  
end module
