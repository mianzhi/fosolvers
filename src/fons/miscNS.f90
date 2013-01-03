!----------------------------------------------------------------------------- best with 100 columns

!> miscellaneous variables and utilities for fons
module miscNS
  use moduleGrid
  use moduleCondition
  public
  
  type(typeGrid)::grid !< the grid
  type(typeCondition),allocatable::condition(:) !< the simulation conditions
  double precision,allocatable::rho(:) !< density
  double precision,allocatable::u(:,:) !< velocity
  double precision,allocatable::p(:) !< pressure
  double precision,allocatable::Temp(:) !< temperature
  double precision,allocatable::E(:) !< internal_energy and kinetic_energy
  double precision,allocatable::IE(:) !< internal_energy
  double precision,allocatable::rhou(:,:) !< momentum per unit volume
  double precision,allocatable::rhoE(:) !< rho*E
  double precision,allocatable::Mass(:) !< block mass
  double precision,allocatable::Mom(:,:) !< extensive node momentum
  double precision,allocatable::Energy(:) !< extensive block energy
  double precision,allocatable::IEnergy(:) !< extensive block internal energy
  double precision,allocatable::rhoNode(:) !< density at node
  double precision,allocatable::uBlock(:,:) !< velocity at block
  double precision,allocatable::uIntf(:,:) !< velocity at interface
  double precision,allocatable::pIntf(:) !< pressure at interface
  double precision,allocatable::visc(:) !< dynamic viscosity
  double precision,allocatable::viscRate(:) !< viscosity rate lambda/mu
  double precision,allocatable::thermK(:) !< thermal conductivity
  double precision,allocatable::tao(:,:,:) !< viscous stress tensor
  double precision,allocatable::taoIntf(:,:,:) !< viscous stress tensor at interface
  double precision,allocatable::oldU(:,:) !< old velocity
  double precision,allocatable::oldP(:) !< old pressure
  double precision,allocatable::oldMom(:,:) !< old extensive momentum
  double precision,allocatable::oldEnergy(:) !< old extensive energy
  double precision,allocatable::preU(:,:) !< velocity of previous pressure coupling iteration
  double precision,allocatable::gradRho(:,:) !< gradient of rho
  double precision,allocatable::gradRhou(:,:,:) !< gradient of rhou
  double precision,allocatable::gradRhoE(:,:) !< gradient of rhoE
  double precision,allocatable::gradT(:,:) !< gradient of temperature
  double precision gamm !< gamma=c_p/c_v
  double precision t !< current time
  double precision dt !< time step size
  double precision tFinal !< final time
  double precision tWrite !< interval of result output
  integer iWrite !< index of result output snapshot
  
  contains
  
  !> read opened simulation control file id
  subroutine readSim(id)
    integer,intent(in)::id !< the file id
    integer ierr
    integer,parameter::DFLT_STR_LEN=400
    character(DFLT_STR_LEN)::tempStr
    
    ierr=0
    rewind(id,iostat=ierr)
    do while(ierr==0)
      read(id,*,iostat=ierr),tempStr
      if(tempStr(1:15)=='$SimulationTime')then
        m=index(tempStr,'(')
        n=index(tempStr,')')
        read(tempStr(m+1:n-1),*),tFinal
      end if
      if(tempStr(1:16)=='$RstTimeInterval')then
        m=index(tempStr,'(')
        n=index(tempStr,')')
        read(tempStr(m+1:n-1),*),tWrite
      end if
    end do
  end subroutine
  
  !> find the stress tensor field from given velocity field Uin
  function findTao(Uin)
    use moduleGrid
    use moduleFVMGrad
    use moduleInterpolation
    double precision,intent(in)::Uin(:,:) !< the input velocity field (binding with nodes)
    double precision,allocatable::findTao(:,:,:) !< the stress tensor
    double precision,allocatable::gradU(:,:,:),gradUBlock(:,:,:)
    
    allocate(findTao(DIMS,DIMS,grid%nBlock))
    allocate(gradU(DIMS,DIMS,grid%nNode))
    allocate(gradUBlock(DIMS,DIMS,grid%nBlock))
    gradU=findGrad(Uin,grid,BIND_NODE)
    gradUBlock=itplNode2Block(gradU,grid)
    forall(l=1:grid%nBlock)
      findTao(:,:,l)=visc(l)*(gradUBlock(:,:,l)+transpose(gradUBlock(:,:,l)))
      forall(i=1:DIMS)
        findTao(i,i,l)=findTao(i,i,l)+viscRate(l)*visc(l)*sum([(gradUBlock(j,j,l),j=1,DIMS)])
      end forall
    end forall
    deallocate(gradU)
    deallocate(gradUBlock)
  end function
  
end module
