!----------------------------------------------------------------------------- best with 100 columns

!> environment for the pressure-based coupled solver
module modPbc
  use modPolyFvGrid
  use modCondition
  use modUDF
  use iso_c_binding
  
  public
  
  integer,parameter::DIMS=3 !< three dimensions
  
  type(polyFvGrid)::grid !< computational grid
  
  double precision::t !< time
  double precision::tFinal !< final time
  double precision::tInt !< time interval of output
  double precision::tNext !< time for next output
  integer::iOut !< index of output
  
  double precision,allocatable::rho(:) !< density [kg/m^3]
  double precision,allocatable::rhou(:,:) !< momentum [kg/s/m^2]
  double precision,allocatable::rhoE(:) !< total energy [J/m^3]
  double precision,allocatable::p(:) !< pressure [Pa]
  double precision,allocatable::u(:,:) !< velocity [m/s]
  double precision,allocatable::temp(:) !< temperature [K]
  double precision,allocatable::Y(:,:) !< mass fraction of species
  
  integer(C_LONG)::nEq !< number of equations
  integer(C_LONG)::ioutFKIN(100) !< integer output of FKINSOL
  double precision::routFKIN(100) !< real output of FKINSOL
  
contains
  
  !> initialize the simulation
  subroutine init()
    use modFileIO
    integer,parameter::FID=10
    double precision::p0,u0(DIMS),T0,pE(DIMS)!,tmpD
    double precision,allocatable::Y0(:)
    integer::ier!,udfIc,udfBc,iUdf(5)
    
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
    !  read(FID,*),p0
    !  read(FID,*),T0
    !  read(FID,*),u0(1)
    !  read(FID,*),u0(2)
    !  read(FID,*),u0(3)
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
    read(FID,*),tFinal
    read(FID,*),tInt
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
    ! initialize nonlinear solver
    nEq=5*grid%nC
    call fnvinits(3,nEq,ier)
    call fkinmalloc(ioutFKIN,routFKIN,ier)
    call fkinspgmr(50,20,ier)
    call fkinsetrin('FNORM_TOL',1d-5,ier)
    call fkinsetrin('SSTEP_TOL',1d-9,ier)
    ! work space and initial state
    allocate(rho(grid%nE))
    allocate(rhou(DIMS,grid%nE))
    allocate(rhoE(grid%nE))
    allocate(p(grid%nE))
    allocate(u(DIMS,grid%nE))
    allocate(temp(grid%nE))
    ! FIXME fix nSp
    allocate(Y(1,grid%nE))
    allocate(Y0(1))
    do i=1,grid%nE
    !  if(udfIc/=0)then
    !    pE(:)=grid%p(:,i)
    !    p0=udf%eval(iUdf(1),pE,0d0)
    !    T0=udf%eval(iUdf(2),pE,0d0)
    !    u0(1)=udf%eval(iUdf(3),pE,0d0)
    !    u0(2)=udf%eval(iUdf(4),pE,0d0)
    !    u0(3)=udf%eval(iUdf(5),pE,0d0)
    !  end if
      ! FIXME remove below
      pE(:)=grid%p(:,i)
      p0=merge(1d5,0.9d5,pE(3)<0.5d0)
      u0(:)=[0d0,0d0,0d0]
      T0=298d0
      Y0=[1d0]
      ! FIXME remove above
      ! FIXME calculation of rho, rhoE based on cantera
      p(i)=p0
      u(:,i)=u0(:)
      temp(i)=T0
      Y(:,i)=Y0(:)
    end do
    call recoverState(p,u,temp,Y,rho,rhou,rhoE)
    t=0d0
    tNext=tInt
    iOut=0
    deallocate(Y0)
  end subroutine
  
  !> clear the simulation environment
  subroutine clear()
    call grid%clear()
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
  
  !> scale and fill primitive state {p,u,T} into variable vector
  subroutine fillVar(pi,ui,tempi,var)
    double precision,intent(in)::pi(:) !< pressure
    double precision,intent(in)::ui(:,:) !< velocity
    double precision,intent(in)::tempi(:) !< temperature
    double precision,intent(inout)::var(:) !< variable vector of the nonlinear problem
    
    !$omp parallel do default(shared)&
    !$omp& private(j)
    do i=1,grid%nE
      j=(i-1)*5
      ! FIXME scale the variable
      var(j+1)=pi(i)*1d0
      var(j+2:j+4)=ui(:,i)*1d0
      var(j+5)=tempi(i)*1d0
    end do
    !$omp end parallel do
  end subroutine
  
  !> extract and de-scale primitive state {p,u,T} from variable vector
  subroutine extractVar(var,po,uo,tempo)
    double precision,intent(in)::var(:) !< variable vector of the nonlinear problem
    double precision,intent(inout)::po(:) !< pressure
    double precision,intent(inout)::uo(:,:) !< velocity
    double precision,intent(inout)::tempo(:) !< temperature
    
    !$omp parallel do default(shared)&
    !$omp& private(j)
    do i=1,grid%nE
      j=(i-1)*5
      ! FIXME de-scale the variable
      po(i)=var(j+1)/1d0
      uo(:,i)=var(j+2:j+4)/1d0
      tempo(i)=var(j+5)/1d0
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
  
end module
