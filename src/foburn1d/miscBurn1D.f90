!----------------------------------------------------------------------------- best with 100 columns

!> miscellaneous variables and utilities for foburn1d
module miscBurn1D
  use moduleGrid1D
  public
  
  type(typeGrid1D)::grid !< the grid
  double precision t !< current time
  double precision dt !< time step size
  double precision tFinal !< simulation time
  double precision width !< width of the domain
  double precision,allocatable::rho(:) !< density
  double precision,allocatable::Y(:) !< mass fraction of reactant
  double precision,allocatable::u(:) !< velocity
  double precision,allocatable::p !< pressure
  double precision,allocatable::Temp(:) !< temperature
  double precision,allocatable::burnR(:) !< burn rate
  double precision,allocatable::Mass(:) !< mass
  double precision gamm !< ratio of heat capacities
  double precision Dm !< mass diffusivity
  double precision alpha !< thermal diffusivity
  double precision mw !< molecular weight
  double precision R !< gas constant for product and reactant
  double precision Cv,Cp !< heat capacities
  double precision Q !< heating value
  
  double precision,parameter::RU=8.3144621d0 !< universal gas constant
  
contains
  
  !> setup environment
  subroutine setEnv()
    allocate(rho(grid%nCell))
    allocate(Y(grid%nCell))
    allocate(u(grid%nNode))
    allocate(Temp(grid%nCell))
    allocate(burnR(grid%nCell))
    allocate(Mass(grid%nCell))
  end subroutine
  
  !> clear environment
  subroutine clearEnv()
    deallocate(rho)
    deallocate(Y)
    deallocate(u)
    deallocate(Temp)
    deallocate(burnR)
    deallocate(Mass)
  end subroutine
  
end module
