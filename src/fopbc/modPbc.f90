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
  
contains
  
  !> initialize the simulation
  subroutine init()
    use modFileIO
    integer,parameter::FID=10
  end subroutine
  
  !> clear the simulation environment
  subroutine clear()
    call grid%clear()
  end subroutine
  
end module
