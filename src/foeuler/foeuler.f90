!----------------------------------------------------------------------------- best with 100 columns

!> environment for Euler problem
module modEuler
  use modPolyFvGrid
  use modCondition
  use iso_c_binding
  
  public
  
  integer,parameter::DIMS=3 !< three dimensions
  
  type(polyFvGrid)::grid !< computational grid
  type(condTab)::bc !< boundary conditions
  double precision::t !< time
  double precision::tOut !< time interval of output
  double precision::r !< specific gas constant
  double precision::gamm !< gamma of gas
  
  double precision,allocatable::rho(:) !< density
  double precision,allocatable::rhou(:,:) !< momentum
  double precision,allocatable::rhoE(:) !< total energy
  double precision,allocatable::p(:) !< pressure
  double precision,allocatable::temp(:) !< temperature
  double precision,allocatable::c(:) !< speed of sound
  
  integer,allocatable::iBC(:) !< indexes boundary conditions

contains
  
  !> initialize the simulation
  subroutine init()
    use modFileIO
    integer,parameter::FID=10
    double precision::p0,T0,u0(DIMS)
    
    ! read inputs
    open(FID,file='grid.msh',action='read')
    call readGMSH(FID,grid)
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
    open(FID,file='fl',action='read')
    read(FID,*),r
    read(FID,*),gamm
    close(FID)
    ! work space and initial state
    allocate(rho(grid%nE))
    allocate(rhou(DIMS,grid%nE))
    allocate(rhoE(grid%nE))
    allocate(p(grid%nE))
    allocate(temp(grid%nE))
    allocate(c(grid%nE))
    forall(i=1:grid%nE)
      rho(i)=p0/r/T0
      rhou(:,i)=rho(i)*u0(:)
      rhoE(i)=rho(i)*(1d0/(gamm-1d0)*r*T0+0.5d0*dot_product(u0,u0))
      p(i)=p0
      temp(i)=T0
      c(i)=sqrt(gamm*r*T0)
    end forall
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
  end subroutine
  
end module

!> Euler solver
program foeuler
  use modEuler
  
  call init()
end program
