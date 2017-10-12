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
  integer,public,parameter::QUAD_N=4 !< 4 nodes per quadrilateral
  
  !> polygon surface mesh type
  type,extends(polyX),public::polyMesh
    double precision,allocatable::a(:) !< surface area
    double precision,allocatable::norm(:,:) !< normal vector
  contains
    procedure,public::init=>initPolyMesh
    procedure,public::clear=>clearPolyMesh
    procedure,public::up=>upPolyMesh
    final::purgePolyMesh
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
    allocate(this%norm(DIMS,nE))
  end subroutine
  
  !> clear this polyMesh
  elemental subroutine clearPolyMesh(this)
    class(polyMesh),intent(inout)::this !< this polyMesh
    
    call this%polyX%clear()
    if(allocated(this%a)) deallocate(this%a)
    if(allocated(this%norm)) deallocate(this%norm)
  end subroutine
  
  !> update this polyMesh
  subroutine upPolyMesh(this)
    use modGeometry
    class(polyMesh),intent(inout)::this !< this polyMesh
    
    if(.not.this%isUp)then
      call this%polyX%up()
      do i=1,this%nE
        select case(this%sE(i))
        case(TRI)
          this%a(i)=a3p(this%pN(:,this%iNE(:,i)))
          this%norm(:,i)=n3p(this%pN(:,this%iNE(:,i)))
        case(QUAD)
          this%a(i)=a4p(this%pN(:,this%iNE(:,i)))
          this%norm(:,i)=n4p(this%pN(:,this%iNE(:,i)))
        case default
          this%a(i)=0d0
          this%norm(:,i)=0d0
        end select
      end do
      this%isUp=.true.
    end if
  end subroutine
  
  !> destructor of polyMesh
  elemental subroutine purgePolyMesh(this)
    type(polyMesh),intent(inout)::this !< this polyMesh
    
    call this%clear()
  end subroutine
  
end module
