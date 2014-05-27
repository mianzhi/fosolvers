!----------------------------------------------------------------------------- best with 100 columns

!> polyhedron grid module
module modPolyGrid
  use modPolyX
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  
  !> polyhedron grid type
  type,extends(polyX),public::polyGrid
  contains
  end type
  
contains
end module
