!----------------------------------------------------------------------------- best with 100 columns

!> grid
module moduleGrid
  use moduleBasicDataStruct
  use moduleSimpleGeometry
  private
  
  ! constants
  integer,parameter,public::DIMS=3 !< dimensions
  
  integer,parameter,public::POINT_TYPE=15 !< point type
  integer,parameter,public::LINE_TYPE=1 !< line type
  integer,parameter,public::TRI_TYPE=2 !< tri type
  integer,parameter,public::QUAD_TYPE=3 !< quad type
  integer,parameter,public::TET_TYPE=4 !< tet type
  integer,parameter,public::HEX_TYPE=5 !< hex type
  
  integer,parameter,public::POINT_NODE_NUM=1 !< number of nodes per point
  integer,parameter,public::LINE_NODE_NUM=2 !< number of nodes per line
  integer,parameter,public::TRI_NODE_NUM=3 !< number of nodes per tri
  integer,parameter,public::QUAD_NODE_NUM=4 !< number of nodes per quad
  integer,parameter,public::TET_NODE_NUM=4 !< number of nodes per tet
  integer,parameter,public::HEX_NODE_NUM=8 !< number of nodes per hex
  
  !> basic information of elements
  type,public::typeEle
    integer Ent !< geometric entity number
    integer Shp !< shape type
    integer nNode !< number of nodes per element
    integer,allocatable::iNode(:) !< node index
    integer,allocatable::Dmn(:) !< physical domain numbers
    integer,allocatable::Prt(:) !< partition numbers
  contains
    procedure,public::init=>initEle
    procedure,public::clear=>clearEle
    !FIXME:final::purgeEle
  end type
  
  !> grid data and procedures
  type,public::typeGrid
    ! basic grid data
    integer nNode !< number of nodes
    double precision,allocatable::NodePos(:,:) !< node position
    integer nPoint !< number of points
    type(typeEle),allocatable::Point(:) !< basic information of points
    integer nLine !< number of lines
    type(typeEle),allocatable::Line(:) !< basic information of lines
    integer nFacet !< number of facets
    type(typeEle),allocatable::Facet(:) !< basic information of facets
    integer nBlock !< number of blocks
    type(typeEle),allocatable::Block(:) !< basic information of facets
    integer nDmn !< number of physical domains
    integer,allocatable::lDmn(:) !< list of physical domains
    integer nPrt !< number of partitions
    integer,allocatable::lPrt(:) !< list of partitions
    ! auxiliary grid data
  contains
    ! basic grid procedures
    procedure,public::init=>initGrid
    procedure,public::clear=>clearGrid
    !FIXME:final::purgeGrid
    ! auxiliary grid procedures
  end type
  
contains
  
  !> initialize this Ele
  elemental subroutine initEle(this)
    class(typeEle),intent(inout)::this !< this element
    
    this%Ent=0
    this%Shp=0
    this%nNode=0
    call this%clear()
  end subroutine
  
  !> clear this Ele
  elemental subroutine clearEle(this)
    class(typeEle),intent(inout)::this !< this element
    
    if(allocated(this%iNode))then
      deallocate(this%iNode)
    end if
    if(allocated(this%Dmn))then
      deallocate(this%Dmn)
    end if
    if(allocated(this%Prt))then
      deallocate(this%Prt)
    end if
  end subroutine
  
  !> destructor of typeEle
  elemental subroutine purgeEle(this)
    type(typeEle),intent(inout)::this !< this element
    
    call this%clear()
  end subroutine
  
  !> initialize this Grid
  elemental subroutine initGrid(this)
    class(typeGrid),intent(inout)::this !< this grid
    
    this%nNode=0
    this%nPoint=0
    this%nLine=0
    this%nFacet=0
    this%nBlock=0
    this%nDmn=0
    this%nPrt=0
    call this%clear()
  end subroutine
  
  !> clear this Grid
  elemental subroutine clearGrid(this)
    class(typeGrid),intent(inout)::this !< this grid
    
    if(allocated(this%NodePos))then
      deallocate(this%NodePos)
    end if
    if(allocated(this%Point))then
      deallocate(this%Point)
    end if
    if(allocated(this%Line))then
      deallocate(this%Line)
    end if
    if(allocated(this%Facet))then
      deallocate(this%Facet)
    end if
    if(allocated(this%Block))then
      deallocate(this%Block)
    end if
    if(allocated(this%lDmn))then
      deallocate(this%lDmn)
    end if
    if(allocated(this%lPrt))then
      deallocate(this%lPrt)
    end if
  end subroutine
  
  !> destructor of typeGrid
  elemental subroutine purgeGrid(this)
    type(typeGrid),intent(inout)::this !< this grid
    
    call this%clear()
  end subroutine
  
end module
