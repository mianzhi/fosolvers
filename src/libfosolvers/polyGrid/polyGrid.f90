!----------------------------------------------------------------------------- best with 100 columns

!> polyhedron grid module
module modPolyGrid
  use modPolyX
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  integer,public,parameter::TET=20 !< tetrahedron
  integer,public,parameter::TET_N=4 !< 4 nodes per tetrahedron
  integer,public,parameter::HEX=21 !< hexahedron
  integer,public,parameter::HEX_N=8 !< 8 nodes per hexahedron
  
  !> polyhedron grid type
  type,extends(polyX),public::polyGrid
    logical::isUp !< if auxiliary data is updated
    double precision,allocatable::v(:) !< volume
  contains
    procedure,public::init=>initPolyGrid
    procedure,public::clear=>clearPolyGrid
    procedure,public::up=>upPolyGrid
    !FIXME:final::purgePolyGrid
  end type
  
contains
  
  !> initialize this polyGrid
  elemental subroutine initPolyGrid(this,nN,nE,m)
    class(polyGrid),intent(inout)::this !< this polyGrid
    integer,intent(in)::nN !< number of nodes
    integer,intent(in)::nE !< number of elements
    integer,intent(in)::m !< maximum number of nodes per elements
    
    call this%polyX%init(nN,nE,m)
    allocate(this%v(nE))
    this%isUp=.false.
  end subroutine
  
  !> clear this polyGrid
  elemental subroutine clearPolyGrid(this)
    class(polyGrid),intent(inout)::this !< this polyGrid
    
    call this%polyX%clear()
    if(allocated(this%v)) deallocate(this%v)
  end subroutine
  
  !> update this polyGrid
  elemental subroutine upPolyGrid(this)
    use modGeometry
    class(polyGrid),intent(inout)::this !< this polyGrid
    
    if(.not.this%isUp)then
      do i=1,this%nE
        select case(this%sE(i))
        case(TET)
          this%v(i)=v4p(this%pN(:,this%iNE(:,i)))
        case(HEX)
          this%v(i)=v8p(this%pN(:,this%iNE(:,i)))
        case default
          this%v(i)=0d0
        end select
      end do
      this%isUp=.true.
    end if
  end subroutine
  
  !> destructor of polyGrid
  elemental subroutine purgePolyGrid(this)
    type(polyGrid),intent(inout)::this !< this polyGrid
    
    call this%clear()
  end subroutine
  
end module
