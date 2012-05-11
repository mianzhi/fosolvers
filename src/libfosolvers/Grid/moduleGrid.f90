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
    logical isUpPointPos !< if the point position is updated
    double precision,allocatable::PointPos(:,:) !< point position
    logical isUpLinePos !< if the line position is updated
    double precision,allocatable::LinePos(:,:) !< line position
    logical isUpFacetPos !< if the facet position is updated
    double precision,allocatable::FacetPos(:,:) !< facet position
    logical isUpBlockPos !< if the block position is updated
    double precision,allocatable::BlockPos(:,:) !< block position
  contains
    ! basic grid procedures
    procedure,public::init=>initGrid
    procedure,public::clear=>clearGrid
    !FIXME:final::purgeGrid
    ! auxiliary grid procedures
    procedure,public::updatePointPos
    procedure,public::updateLinePos
    procedure,public::updateFacetPos
    procedure,public::updateBlockPos
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
    
    if(allocated(this%iNode)) deallocate(this%iNode)
    if(allocated(this%Dmn)) deallocate(this%Dmn)
    if(allocated(this%Prt)) deallocate(this%Prt)
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
    this%isUpPointPos=.false.
    this%isUpLinePos=.false.
    this%isUpFacetPos=.false.
    this%isUpBlockPos=.false.
    call this%clear()
  end subroutine
  
  !> clear this Grid
  elemental subroutine clearGrid(this)
    class(typeGrid),intent(inout)::this !< this grid
    
    if(allocated(this%NodePos)) deallocate(this%NodePos)
    if(allocated(this%Point)) deallocate(this%Point)
    if(allocated(this%Line)) deallocate(this%Line)
    if(allocated(this%Facet)) deallocate(this%Facet)
    if(allocated(this%Block)) deallocate(this%Block)
    if(allocated(this%lDmn)) deallocate(this%lDmn)
    if(allocated(this%lPrt)) deallocate(this%lPrt)
    if(allocated(this%PointPos)) deallocate(this%PointPos)
    if(allocated(this%LinePos)) deallocate(this%LinePos)
    if(allocated(this%FacetPos)) deallocate(this%FacetPos)
    if(allocated(this%BlockPos)) deallocate(this%BlockPos)
  end subroutine
  
  !> destructor of typeGrid
  elemental subroutine purgeGrid(this)
    type(typeGrid),intent(inout)::this !< this grid
    
    call this%clear()
  end subroutine
  
  !> update the point position
  elemental subroutine updatePointPos(this)
    class(typeGrid),intent(inout)::this !< this grid
    
    if(.not.this%isUpPointPos)then
      if(allocated(this%PointPos))then
        if(size(this%PointPos,2)/=this%nPoint)then
          deallocate(this%PointPos)
          allocate(this%PointPos(DIMS,this%nPoint))
        end if
      else
        allocate(this%PointPos(DIMS,this%nPoint))
      end if
      forall(i=1:this%nPoint)
        this%PointPos(:,i)=this%NodePos(:,this%Point(i)%iNode(1))
      end forall
      this%isUpPointPos=.true.
    end if
  end subroutine
  
  !> update the line position
  elemental subroutine updateLinePos(this)
    class(typeGrid),intent(inout)::this !< this grid
    
    if(.not.this%isUpLinePos)then
      if(allocated(this%LinePos))then
        if(size(this%LinePos,2)/=this%nLine)then
          deallocate(this%LinePos)
          allocate(this%LinePos(DIMS,this%nLine))
        end if
      else
        allocate(this%LinePos(DIMS,this%nLine))
      end if
      forall(i=1:this%nLine)
        this%LinePos(:,i)=sum(this%NodePos(:,this%Line(i)%iNode(:)),2)/dble(this%Line(i)%nNode)
      end forall
      this%isUpLinePos=.true.
    end if
  end subroutine
  
  !> update the facet position
  elemental subroutine updateFacetPos(this)
    class(typeGrid),intent(inout)::this !< this grid
    
    if(.not.this%isUpFacetPos)then
      if(allocated(this%FacetPos))then
        if(size(this%FacetPos,2)/=this%nFacet)then
          deallocate(this%FacetPos)
          allocate(this%FacetPos(DIMS,this%nFacet))
        end if
      else
        allocate(this%FacetPos(DIMS,this%nFacet))
      end if
      forall(i=1:this%nFacet)
        this%FacetPos(:,i)=sum(this%NodePos(:,this%Facet(i)%iNode(:)),2)/dble(this%Facet(i)%nNode)
      end forall
      this%isUpFacetPos=.true.
    end if
  end subroutine
  
  !> update the block position
  elemental subroutine updateBlockPos(this)
    class(typeGrid),intent(inout)::this !< this grid
    
    if(.not.this%isUpBlockPos)then
      if(allocated(this%BlockPos))then
        if(size(this%BlockPos,2)/=this%nBlock)then
          deallocate(this%BlockPos)
          allocate(this%BlockPos(DIMS,this%nBlock))
        end if
      else
        allocate(this%BlockPos(DIMS,this%nBlock))
      end if
      forall(i=1:this%nBlock)
        this%BlockPos(:,i)=sum(this%NodePos(:,this%Block(i)%iNode(:)),2)/dble(this%Block(i)%nNode)
      end forall
      this%isUpBlockPos=.true.
    end if
  end subroutine
  
end module
