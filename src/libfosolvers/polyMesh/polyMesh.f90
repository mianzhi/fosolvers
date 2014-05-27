!----------------------------------------------------------------------------- best with 100 columns

!> polygon surface mesh module
module modPolyMesh
  use modPolyX
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  
  !> polygon surface mesh type
  type,extends(polyX),public::polyMesh
  contains
  end type
  
contains
end module
