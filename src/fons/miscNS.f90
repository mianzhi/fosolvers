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
  double precision,allocatable::E(:) !< internal_energy+kinetic_energy
  double precision,allocatable::rhou(:,:) !< momentum per unit volume
  double precision,allocatable::rhoE(:) !< rho*E
  double precision,allocatable::Mass(:) !< block mass
  double precision,allocatable::Mom(:,:) !< extensive node momentum
  double precision,allocatable::Energy(:) !< extensive block energy
  double precision,allocatable::rhoNode(:) !< density at node
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
  
end module
