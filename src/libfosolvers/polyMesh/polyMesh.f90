!----------------------------------------------------------------------------- best with 100 columns

!> polygon surface mesh module
module modPolyMesh
  use modPolyX
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  integer,public,parameter::TRI=10 !< triangle
  integer,public,parameter::TRI_N=3 !< 3 nodes per triangle
  integer,public,parameter::QUAD=11 !< quadrilateral
  integer,public,parameter::QUAD_N=3 !< 3 nodes per quadrilateral
  
  !> polygon surface mesh type
  type,extends(polyX),public::polyMesh
    double precision,allocatable::a(:) !< surface area
  contains
    procedure,public::init=>initPolyMesh
    procedure,public::clear=>clearPolyMesh
    !FIXME:final::purgePolyMesh
  end type
  
contains
  
  !> initialize this polyMesh
  elemental subroutine initPolyMesh(this,nN,nE,m)
    class(polyMesh),intent(inout)::this !< this polyMesh
    integer,intent(in)::nN !< number of nodes
    integer,intent(in)::nE !< number of elements
    integer,intent(in)::m !< maximum number of nodes per elements
    
    call this%polyX%init(nN,nE,m)
    allocate(this%a(nE))
  end subroutine
  
  !> clear this polyMesh
  elemental subroutine clearPolyMesh(this)
    class(polyMesh),intent(inout)::this !< this polyMesh
    
    call this%polyX%clear()
    if(allocated(this%a)) deallocate(this%a)
  end subroutine
  
  !> destructor of polyMesh
  elemental subroutine purgePolyMesh(this)
    type(polyMesh),intent(inout)::this !< this polyMesh
    
    call this%clear()
  end subroutine
  
end module
