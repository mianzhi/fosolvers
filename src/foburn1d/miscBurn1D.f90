!----------------------------------------------------------------------------- best with 100 columns

!> miscellaneous variables and utilities for foburn1d
module miscBurn1D
  use moduleGrid1D
  public
  
  type(typeGrid1D)::grid !< the grid
  double precision t !< current time
  double precision dt !< time step size
  double precision tFinal !< simulation time
  double precision,allocatable::rho(:) !< density
  double precision,allocatable::Y(:) !< mass fraction of reactant
  double precision,allocatable::u(:) !< velocity
  double precision,allocatable::p(:) !< pressure
  double precision,allocatable::Temp(:) !< temperature
  double precision,allocatable::E(:) !< internal_energy and kinetic_energy
  double precision,allocatable::IE(:) !< internal_energy
  double precision,allocatable::rhou(:) !< momentum per unit volume
  double precision,allocatable::rhoE(:) !< rho*E
  double precision,allocatable::Mass(:) !< cell mass
  double precision,allocatable::Mom(:) !< extensive cell momentum
  double precision,allocatable::Energy(:) !< extensive cell energy
  double precision,allocatable::IEnergy(:) !< extensive cell internal energy
  double precision,allocatable::rhox(:) !< derivative of rho along x
  double precision gamm
  
contains
  
  !> setup environment
  subroutine setEnv()
    allocate(rho(grid%nCell))
    allocate(Y(grid%nCell))
    allocate(u(grid%nCell))
    allocate(p(grid%nCell))
    allocate(Temp(grid%nCell))
    allocate(E(grid%nCell))
    allocate(IE(grid%nCell))
    allocate(rhou(grid%nCell))
    allocate(rhoE(grid%nCell))
    allocate(Mass(grid%nCell))
    allocate(Mom(grid%nCell))
    allocate(Energy(grid%nCell))
    allocate(IEnergy(grid%nCell))
    allocate(rhox(grid%nCell))
  end subroutine
  
  !> clear environment
  subroutine clearEnv()
    deallocate(rho)
    deallocate(Y)
    deallocate(u)
    deallocate(p)
    deallocate(Temp)
    deallocate(E)
    deallocate(IE)
    deallocate(rhou)
    deallocate(rhoE)
    deallocate(Mass)
    deallocate(Mom)
    deallocate(Energy)
    deallocate(IEnergy)
    deallocate(rhox)
  end subroutine
  
end module
