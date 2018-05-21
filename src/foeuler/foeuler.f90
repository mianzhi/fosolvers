!----------------------------------------------------------------------------- best with 100 columns

!> environment for Euler problem
module modEuler
  use modPolyFvGrid
  use modCondition
  use modUDF
  use modNumerics
  use iso_c_binding
  
  public
  
  integer,parameter::DIMS=3 !< three dimensions
  double precision,parameter::RTOL=1d-6 !< relative tolerance
  double precision,parameter::ATOL=1d-12 !< absolute tolerance
  
  integer,parameter::BC_WALL=0 !< wall boundary
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
  
  double precision::t !< time
  double precision::dt !< time step size determined by ODE solver
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
  double precision,allocatable::JacP(:,:,:) !< Jacobian matrix on pairs
  double precision,allocatable::JacC(:,:,:) !< Jacobian matrix in cells
  
  double precision,allocatable::y(:) !< solution vector
  double precision,allocatable::dRho(:) !< time derivative of density
  double precision,allocatable::dRhou(:,:) !< time derivative of momentum
  double precision,allocatable::dRhoE(:) !< time derivative of total energy
  integer::nEq !< number of equations
  real,allocatable::JacS(:,:,:) !< single precision version of JacC
  integer,allocatable::precPiv(:,:) !< pivoting arrays for Lapack
  real,allocatable::precRhs(:,:) !< RHS arrays for Lapack
  type(BDFNewtonKrylov)::ode !< the ODE problem

contains
  
  !> initialize the simulation
  subroutine init()
    use modFileIO
    integer,parameter::FID=10
    double precision::p0,T0,u0(DIMS),pE(DIMS),tmpD
    integer::udfIc,udfBc,iUdf(5)
    
    ! read inputs
    open(FID,file='grid.msh',action='read')
    call readGMSH(FID,grid)
    call grid%up()
    close(FID)
    open(FID,file='bc',action='read')
    call readCondTab(FID,bc)
    if(any(bc%t(:)==BC_IN_STATIC_UDF).or.&
    &  any(bc%t(:)==BC_IN_TOTAL_UDF).or.&
    &  any(bc%t(:)==BC_OUT_UDF).or.&
    &  any(bc%t(:)==BC_FAR_UDF))then
      udfBc=1
    else
      udfBc=0
    end if
    close(FID)
    open(FID,file='ic',action='read')
    read(FID,*)udfIC
    if(udfIc==0)then
      read(FID,*)p0
      read(FID,*)T0
      read(FID,*)u0(1)
      read(FID,*)u0(2)
      read(FID,*)u0(3)
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
    read(FID,*)tInt
    close(FID)
    open(FID,file='fl',action='read')
    read(FID,*)r
    read(FID,*)gamm
    close(FID)
    if(udfIc/=0.or.udfBc/=0)then
      open(FID,file='udf',action='read')
      call readUDFTab(FID,udf)
      close(FID)
    end if
    ! work space and initial state
    allocate(rho(grid%nE))
    allocate(rhou(DIMS,grid%nE))
    allocate(rhoE(grid%nE))
    allocate(u(DIMS,grid%nE))
    allocate(p(grid%nE))
    allocate(temp(grid%nE))
    allocate(c(grid%nE))
    do i=1,grid%nE
      if(udfIc/=0)then
        pE(:)=grid%p(:,i)
        p0=udf%eval(iUdf(1),pE,0d0)
        T0=udf%eval(iUdf(2),pE,0d0)
        u0(1)=udf%eval(iUdf(3),pE,0d0)
        u0(2)=udf%eval(iUdf(4),pE,0d0)
        u0(3)=udf%eval(iUdf(5),pE,0d0)
      end if
      rho(i)=p0/r/T0
      rhou(:,i)=rho(i)*u0(:)
      rhoE(i)=rho(i)*(1d0/(gamm-1d0)*r*T0+0.5d0*dot_product(u0,u0))
      u(:,i)=u0(:)
      p(i)=p0
      temp(i)=T0
      c(i)=sqrt(gamm*r*T0)
    end do
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
    call ode%init(nEq,rhs)
    call ode%setTol(RTOL,ATOL)
    call ode%setMaxSteps(-1)
    call ode%setIV(t,y)
    ! initialize preconditioning solver
    allocate(precPiv(5,grid%nC))
    allocate(precRhs(5,grid%nC))
  end subroutine
  
  !> clear the simulation environment
  subroutine clear()
    deallocate(iBC,rho,rhou,rhoE,u,p,temp,c)
    deallocate(y,dRho,dRhou,dRhoE)
    deallocate(JacP,JacC,JacS)
    deallocate(precPiv,precRhs)
    call grid%clear()
    call bc%clear()
  end subroutine
  
  !> synchronize the state with the solution vector
  subroutine syncState(x)
    double precision,intent(in)::x(*) !< solution vector
    
    !$omp parallel do default(shared)
    do i=1,grid%nC
      rho(i)=x(i)
      rhou(1,i)=x(grid%nC+i)
      rhou(2,i)=x(2*grid%nC+i)
      rhou(3,i)=x(3*grid%nC+i)
      rhoE(i)=x(4*grid%nC+i)
    end do
    !$omp end parallel do
    !$omp parallel do default(shared)
    do i=1,grid%nE
      u(:,i)=rhou(:,i)/rho(i)
      p(i)=(gamm-1d0)*(rhoE(i)-0.5d0*dot_product(rhou(:,i),rhou(:,i))/rho(i))
      temp(i)=p(i)/rho(i)/r
      c(i)=sqrt(gamm*r*temp(i))
    end do
    !$omp end parallel do
  end subroutine
  
  !> set the boundary conditions
  subroutine setBC(x)
    double precision,intent(in)::x(*) !< solution vector
    double precision::pGst,TGst,ptGst,TtGst,uGst(DIMS),rhoGst,Mach,&
    &                 un(DIMS),ut(DIMS),ui(DIMS),charI,charO,a,pP(DIMS)
    
    do i=1,grid%nP
      m=grid%iEP(1,i)
      n=grid%iEP(2,i)
      if(n>grid%nC)then
        if(iBC(n)==0)then ! default wall boundary
          rho(n)=x(m)
          rhou(:,n)=x(grid%nC*[1,2,3]+m)&
          &         -2d0*dot_product(x(grid%nC*[1,2,3]+m),grid%normP(:,i))*grid%normP(:,i)
          rhoE(n)=x(4*grid%nC+m)
        else
          select case(bc%t(iBC(n)))
          case(BC_WALL) ! wall boundary
            rho(n)=x(m)
            rhou(:,n)=x(grid%nC*[1,2,3]+m)&
            &         -2d0*dot_product(x(grid%nC*[1,2,3]+m),grid%normP(:,i))*grid%normP(:,i)
            rhoE(n)=x(4*grid%nC+m)
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
            Mach=norm2(uGst)/sqrt(gamm*r*TGst)
            rhoGst=pGst/r/TGst
            if(Mach<1d0)then
              uGst(:)=x(grid%nC*[1,2,3]+m)/rhoGst
            end if
            rho(n)=rhoGst
            rhou(:,n)=rhoGst*uGst(:)
            rhoE(n)=rhoGst*(1d0/(gamm-1d0)*r*TGst+0.5d0*dot_product(uGst,uGst))
          case(BC_IN_TOTAL,BC_IN_TOTAL_UDF) ! inflow boundary with total properties
            if(bc%t(iBC(n))==BC_IN_TOTAL)then
              ptGst=bc%p(1,iBC(n))
              TtGst=bc%p(2,iBC(n))
              uGst=bc%p(3:5,iBC(n))
            else
              pP(:)=grid%pP(:,i)
              ptGst=udf%eval(int(bc%p(1,iBC(n))),pP,t)
              TtGst=udf%eval(int(bc%p(2,iBC(n))),pP,t)
              uGst(1)=udf%eval(int(bc%p(3,iBC(n))),pP,t)
              uGst(2)=udf%eval(int(bc%p(4,iBC(n))),pP,t)
              uGst(3)=udf%eval(int(bc%p(5,iBC(n))),pP,t)
            end if
            TGst=TtGst-0.5d0*(gamm-1d0)/gamm/r*dot_product(uGst,uGst)
            Mach=norm2(uGst)/sqrt(gamm*r*TGst)
            pGst=ptGst*(1d0+0.5d0*(gamm-1d0)*Mach**2)**(-gamm/(gamm-1d0))
            rhoGst=pGst/r/TGst
            if(Mach<1d0)then
              uGst(:)=dot_product(x(grid%nC*[1,2,3]+m)/x(m),uGst(:))*uGst/dot_product(uGst,uGst)
              TGst=TtGst-0.5d0*(gamm-1d0)/gamm/r*dot_product(uGst,uGst)
              Mach=norm2(uGst)/sqrt(gamm*r*TGst)
              pGst=ptGst*(1d0+0.5d0*(gamm-1d0)*Mach**2)**(-gamm/(gamm-1d0))
              rhoGst=pGst/r/TGst
            end if
            rho(n)=rhoGst
            rhou(:,n)=rhoGst*uGst(:)
            rhoE(n)=rhoGst*(1d0/(gamm-1d0)*r*TGst+0.5d0*dot_product(uGst,uGst))
          case(BC_OUT,BC_OUT_UDF) ! outflow boundary
            if(bc%t(iBC(n))==BC_OUT)then
              pGst=bc%p(1,iBC(n))
            else
              pP(:)=grid%pP(:,i)
              pGst=udf%eval(int(bc%p(1,iBC(n))),pP,t)
            end if
            rho(n)=x(m)
            rhou(:,n)=x(grid%nC*[1,2,3]+m)
            rhoE(n)=1d0/(gamm-1d0)*pGst&
            &       +0.5d0*dot_product(x(grid%nC*[1,2,3]+m),x(grid%nC*[1,2,3]+m))/x(m)
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
            charO=dot_product(uGst,grid%normP(:,i))-2d0*sqrt(gamm*r*TGst)/(gamm-1d0)
            ui(:)=x(grid%nC*[1,2,3]+m)/x(m)
            charI=dot_product(ui,grid%normP(:,i))&
            &     +2d0*sqrt(gamm/x(m)/(gamm-1d0)*(x(4*grid%nC+m)-0.5d0*dot_product(ui,ui)*x(m)))
            un(:)=0.5d0*grid%normP(:,i)*(charI+charO)
            a=0.25d0*(gamm-1d0)*(charI-charO)
            if(dot_product(un(:),grid%normP(:,i))>0d0)then
              ut(:)=ui(:)-grid%normP(:,i)*dot_product(grid%normP(:,i),ui(:))
              Mach=norm2(ui)/&
              &    sqrt(gamm/x(m)*(gamm-1d0)*(x(4*grid%nC+m)-0.5d0*dot_product(ui,ui)*x(m)))
              ptGst=(gamm-1d0)*(x(4*grid%nC+m)-0.5d0*dot_product(ui,ui)*x(m))*&
              &     (1d0+0.5d0*(gamm-1d0)*Mach**2)**(gamm/(gamm-1d0))
            else
              ut(:)=uGst(:)-un
              Mach=norm2(uGst)/sqrt(gamm*r*TGst)
              ptGst=pGst*(1d0+0.5d0*(gamm-1d0)*Mach**2)**(gamm/(gamm-1d0))
            end if
            uGst(:)=un(:)+ut(:)
            TGst=a**2/gamm/r
            Mach=norm2(uGst)/a
            pGst=ptGst*(1d0+0.5d0*(gamm-1d0)*Mach**2)**(-gamm/(gamm-1d0))
            rhoGst=pGst/r/TGst
            rho(n)=rhoGst
            rhou(:,n)=rhoGst*uGst(:)
            rhoE(n)=rhoGst*(1d0/(gamm-1d0)*r*TGst+0.5d0*dot_product(uGst,uGst))
          case default
          end select
        end if
      end if
    end do
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
  
  !> the ODE to be advanced
  function rhs(c_time,c_x,c_dxdt,dat)
    use modAdvection
    real(C_DOUBLE),value::c_time
    type(C_PTR),value::c_x,c_dxdt,dat
    integer(C_INT)::rhs
    double precision::time
    double precision,pointer::x(:)
    double precision,pointer::dxdt(:)
    type(C_PTR)::foo
    
    time=c_time
    call associateVector(c_x,x)
    call associateVector(c_dxdt,dxdt)
    
    call setBC(x)
    call syncState(x)
    call findAdv(grid,rho,rhou,rhoE,p,gamm,dRho,dRhou,dRhoE)
    forall(i=1:grid%nC)
      dxdt(i)=dRho(i)/grid%v(i)
      dxdt(grid%nC+i)=dRhou(1,i)/grid%v(i)
      dxdt(2*grid%nC+i)=dRhou(2,i)/grid%v(i)
      dxdt(3*grid%nC+i)=dRhou(3,i)/grid%v(i)
      dxdt(4*grid%nC+i)=dRhoE(i)/grid%v(i)
    end forall
    rhs=0
    foo=dat
  end function
  
end module

!> Euler solver
program foeuler
  use modEuler
  character(20)::tmpStr
  integer::ier
  
  call init()
  write(tmpStr,*)iOut
  call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
  do while(t<tFinal)
    call ode%step(tNext,t,y)
    call ode%getDt(dt)
    write(*,'(a,g12.5,a,g12.5)')'[i] finished time step: t=',t,' dt=',dt
    if(t+tiny(1d0)>=tNext)then
      iOut=iOut+1
      write(tmpStr,*)iOut
      call syncState(y)
      write(*,'(a)')'[i] writing: rst_'//trim(adjustl(tmpStr))//'.vtk'
      call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
      tNext=tNext+tInt
    end if
  end do
  call clear()
end program

!> setup/factor preconditioning matrix
subroutine fcvpset(time,x,fx,Jok,Jcur,pGamm,h,iPara,rPara,work1,work2,work3,ier)
  use modEuler
  use modAdvection
  use iso_c_binding
  double precision::time
  double precision::x(*)
  double precision::fx(*)
  integer::Jok
  integer::Jcur
  double precision::pGamm
  double precision::h
  integer(C_LONG)::iPara(*)
  double precision::rPara(*)
  double precision::work1(*)
  double precision::work2(*)
  double precision::work3(*)
  integer::ier
  
  if(Jok==1)then
    Jcur=0
  else
    call setBC(x)
    call syncState(x)
    call findAdvJac(grid,rho,rhou,rhoE,p,gamm,JacP,JacC)
    if(.not.allocated(JacS))then
      allocate(JacS(size(JacC,1),size(JacC,2),size(JacC,3)))
    end if
    JacS(:,:,:)=real(JacC(:,:,:))
    Jcur=1
  end if
  JacS(:,:,:)=-real(pGamm)*JacS(:,:,:)
  forall(i=1:grid%nC)
    forall(j=1:5)
      JacS(j,j,i)=JacS(j,j,i)+1e0
    end forall
  end forall
  !$omp parallel do default(shared)&
  !$omp& private(ier)
  do i=1,grid%nC
    call SGETRF(5,5,JacS(:,1:5,i),5,precPiv(:,i),ier)
  end do
  !$end omp parallel do
  ier=0
end subroutine

!> solve the preconditioning problem
subroutine fcvpsol(time,x,fx,res,z,pGamm,delta,lr,iPara,rPara,work,ier)
  use modEuler
  use iso_c_binding
  double precision::time
  double precision::x(*)
  double precision::fx(*)
  double precision::res(*)
  double precision::z(*)
  double precision::pGamm
  double precision::delta
  integer::lr
  integer(C_LONG)::iPara(*)
  double precision::rPara(*)
  double precision::work(*)
  integer::ier
  
  n=max(floor(-log10(pGamm)),0)
  ! adaptive number of Jacobi iteration
  if(n>=6)then
    z(1:5*grid%nC)=res(1:5*grid%nC)
    ier=0
    return
  else if(n>=5)then
    n=0
  else
    n=1
  end if
  !$omp parallel do default(shared)&
  !$omp& private(ier)
  do i=1,grid%nC
    precRhs(:,i)=real(res([0,1,2,3,4]*grid%nC+i))
    call SGETRS('N',5,1,JacS(:,1:5,i),5,precPiv(:,i),precRhs(:,i),5,ier)
    z([0,1,2,3,4]*grid%nC+i)=dble(precRhs(:,i))
  end do
  !$end omp parallel do
  do l=1,n ! additional Jacobi loops
    !$omp parallel do default(shared)&
    !$omp& private(j,m,ier)
    do i=1,grid%nC
      precRhs(:,i)=real(res([0,1,2,3,4]*grid%nC+i))
      do j=1,size(grid%neib,1)
        m=grid%neib(j,i)
        if(1<=m.and.m<=grid%nC)then
          precRhs(:,i)=precRhs(:,i)-matmul(JacS(:,j*5+1:j*5+5,i),real(z([0,1,2,3,4]*grid%nC+m)))
        end if
      end do
    end do
    !$end omp parallel do
    !$omp parallel do default(shared)&
    !$omp& private(ier)
    do i=1,grid%nC
      call SGETRS('N',5,1,JacS(:,1:5,i),5,precPiv(:,i),precRhs(:,i),5,ier)
      z([0,1,2,3,4]*grid%nC+i)=dble(precRhs(:,i))
    end do
    !$end omp parallel do
  end do
  ier=0
end subroutine
