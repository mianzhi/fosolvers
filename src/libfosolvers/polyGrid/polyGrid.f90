!----------------------------------------------------------------------------- best with 100 columns

!> polyhedron and polygon grid module
module modPolyGrid
  use modPolyX
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  integer,public,parameter::TRI=10 !< triangle
  integer,public,parameter::TRI_N=3 !< 3 nodes per triangle
  integer,public,parameter::QUAD=11 !< quadrilateral
  integer,public,parameter::QUAD_N=4 !< 4 nodes per quadrilateral
  integer,public,parameter::TET=20 !< tetrahedron
  integer,public,parameter::TET_N=4 !< 4 nodes per tetrahedron
  integer,public,parameter::HEX=21 !< hexahedron
  integer,public,parameter::HEX_N=8 !< 8 nodes per hexahedron
  
  !> polyhedron and polygon grid type
  type,extends(polyX),public::polyGrid
    integer,allocatable::gid(:) !< geometric group identifier
    logical::isUp !< if auxiliary data is updated
    double precision,allocatable::v(:) !< volume
  contains
    procedure,public::init=>initPolyGrid
    procedure,public::clear=>clearPolyGrid
    procedure,public::up=>upPolyGrid
    !FIXME:final::purgePolyGrid
  end type
  
  ! public procedures
  public::nF
  public::nNF
  public::getINF
  
contains
  
  !> initialize this polyGrid
  elemental subroutine initPolyGrid(this,nN,nE,m)
    class(polyGrid),intent(inout)::this !< this polyGrid
    integer,intent(in)::nN !< number of nodes
    integer,intent(in)::nE !< number of elements
    integer,intent(in)::m !< maximum number of nodes per elements
    
    call this%polyX%init(nN,nE,m)
    allocate(this%gid(nE))
    allocate(this%v(nE))
    this%isUp=.false.
  end subroutine
  
  !> clear this polyGrid
  elemental subroutine clearPolyGrid(this)
    class(polyGrid),intent(inout)::this !< this polyGrid
    
    call this%polyX%clear()
    if(allocated(this%gid)) deallocate(this%gid)
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
  
  !> number of faces per given shape
  elemental function nF(s)
    integer,intent(in)::s !< shape
    integer::nF !< result
    
    select case(s)
    case(TRI)
      nF=1
    case(QUAD)
      nF=1
    case(TET)
      nF=4
    case(HEX)
      nF=6
    case default
      nF=0
    end select
  end function
  
  !> number of nodes per given shape's i_th face
  elemental function nNF(s,i)
    integer,intent(in)::s !< shape
    integer,intent(in)::i !< face index
    integer::nNF !< result
    
    nNF=0
    select case(s)
    case(TRI)
      select case(i)
      case(1)
        nNF=3
      case default
      end select
    case(QUAD)
      select case(i)
      case(1)
        nNF=4
      case default
      end select
    case(TET)
      select case(i)
      case(1:4)
        nNF=3
      case default
      end select
    case(HEX)
      select case(i)
      case(1:6)
        nNF=4
      case default
      end select
    case default
    end select
  end function
  
  !> get the indices of nodes of given shape's i_th face
  pure subroutine getINF(s,i,ind)
    integer,intent(in)::s !< shape
    integer,intent(in)::i !< face index
    integer,intent(inout)::ind(:) !< result
    
    select case(s)
    case(TRI)
      select case(i)
      case(1)
        ind(1:3)=[1,2,3]
      case default
      end select
    case(QUAD)
      select case(i)
      case(1)
        ind(1:4)=[1,2,3,4]
      case default
      end select
    case(TET)
      select case(i)
      case(1)
        ind(1:3)=[1,3,2]
      case(2)
        ind(1:3)=[1,2,4]
      case(3)
        ind(1:3)=[1,4,3]
      case(4)
        ind(1:3)=[2,3,4]
      case default
      end select
    case(HEX)
      select case(i)
      case(1)
        ind(1:4)=[2,3,7,6]
      case(2)
        ind(1:4)=[1,5,8,4]
      case(3)
        ind(1:4)=[3,4,8,7]
      case(4)
        ind(1:4)=[1,2,6,5]
      case(5)
        ind(1:4)=[5,6,7,8]
      case(6)
        ind(1:4)=[1,4,3,2]
      case default
      end select
    case default
    end select
  end subroutine
  
end module
