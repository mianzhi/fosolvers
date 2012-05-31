!----------------------------------------------------------------------------- best with 100 columns

!> grid
module moduleGrid
  use moduleBasicDataStruct
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
  
  integer,parameter,public::TET_SURF_NUM=4 !< number of surfaces per tet
  integer,parameter,public::HEX_SURF_NUM=6 !< number of surfaces per hex
  
  integer,parameter,public::INTF_NEIB_BLOCK_NUM=2 !< number of neighbour blocks per interface
  
  integer,public::TET_SURF_TAB(TRI_NODE_NUM,TET_SURF_NUM) !< table of surface nodes for tet
  parameter(TET_SURF_TAB=reshape([1,3,2,1,2,4,1,4,3,2,3,4],[TRI_NODE_NUM,TET_SURF_NUM]))
  integer,public::HEX_SURF_TAB(QUAD_NODE_NUM,HEX_SURF_NUM) !< table of surface nodes for hex
  parameter(HEX_SURF_TAB=reshape([2,3,7,6,1,5,8,4,3,4,8,7,1,2,6,5,5,6,7,8,1,4,3,2],&
  &                              [QUAD_NODE_NUM,HEX_SURF_NUM]))
  
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
    logical isUpNodeNeibBlock !< if the node neighbour block is updated
    type(typeHtr1DIArr),allocatable::NodeNeibBlock(:) !< node neighbour block
    
    logical isUpIntf !< if the interface between blocks is updated
    integer nIntf !< number of interfaces between blocks
    type(typeEle),allocatable::Intf(:) !< interface between blocks
    integer,allocatable::IntfNeibBlock(:,:) !< neighbour blocks of interface
    
    logical isUpPointPos !< if the point position is updated
    double precision,allocatable::PointPos(:,:) !< point position
    
    logical isUpLinePos !< if the line position is updated
    double precision,allocatable::LinePos(:,:) !< line position
    
    logical isUpFacetPos !< if the facet position is updated
    double precision,allocatable::FacetPos(:,:) !< facet position
    
    logical isUpBlockPos !< if the block position is updated
    double precision,allocatable::BlockPos(:,:) !< block position
    
    logical isUpLineLen !< if the line length is updated
    double precision,allocatable::LineLen(:) !< line length
    
    logical isUpFacetArea !< if the facet area is updated
    double precision,allocatable::FacetArea(:) !< facet area
    
    logical isUpBlockVol !< if the block volume is updated
    double precision,allocatable::BlockVol(:) !< block volume
    
    logical isUpFacetNorm !< if the facet normal vector is updated
    double precision,allocatable::FacetNorm(:,:) !< facet normal vector
  contains
    ! basic grid procedures
    procedure,public::init=>initGrid
    procedure,public::clear=>clearGrid
    !FIXME:final::purgeGrid
    ! auxiliary grid procedures
    procedure,public::updateNodeNeibBlock
    procedure,public::updateIntf
    procedure,public::updatePointPos
    procedure,public::updateLinePos
    procedure,public::updateFacetPos
    procedure,public::updateBlockPos
    procedure,public::updateLineLen
    procedure,public::updateFacetArea
    procedure,public::updateBlockVol
    procedure,public::updateFacetNorm
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
    this%isUpNodeNeibBlock=.false.
    this%isUpIntf=.false.
    this%isUpPointPos=.false.
    this%isUpLinePos=.false.
    this%isUpFacetPos=.false.
    this%isUpBlockPos=.false.
    this%isUpLineLen=.false.
    this%isUpFacetArea=.false.
    this%isUpBlockVol=.false.
    this%isUpFacetNorm=.false.
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
    if(allocated(this%NodeNeibBlock)) deallocate(this%NodeNeibBlock)
    if(allocated(this%Intf)) deallocate(this%Intf)
    if(allocated(this%IntfNeibBlock)) deallocate(this%IntfNeibBlock)
    if(allocated(this%PointPos)) deallocate(this%PointPos)
    if(allocated(this%LinePos)) deallocate(this%LinePos)
    if(allocated(this%FacetPos)) deallocate(this%FacetPos)
    if(allocated(this%BlockPos)) deallocate(this%BlockPos)
    if(allocated(this%LineLen)) deallocate(this%LineLen)
    if(allocated(this%FacetArea)) deallocate(this%FacetArea)
    if(allocated(this%BlockVol)) deallocate(this%BlockVol)
    if(allocated(this%FacetNorm)) deallocate(this%FacetNorm)
  end subroutine
  
  !> destructor of typeGrid
  elemental subroutine purgeGrid(this)
    type(typeGrid),intent(inout)::this !< this grid
    
    call this%clear()
  end subroutine
  
  !> update the node neighbour block
  elemental subroutine updateNodeNeibBlock(this)
    class(typeGrid),intent(inout)::this !< this grid
    
    if(.not.this%isUpNodeNeibBlock)then
      call reallocArr(this%NodeNeibBlock,this%nNode)
      do i=1,this%nBlock
        do j=1,this%Block(i)%nNode
          call pushArr(this%NodeNeibBlock(this%Block(i)%iNode(j))%dat,i)
        end do
      end do
      this%isUpNodeNeibBlock=.true.
    end if
  end subroutine
  
  !> update the interface between blocks
  elemental subroutine updateIntf(this)
    class(typeGrid),intent(inout)::this !< this grid
    
    if(.not.this%isUpIntf)then
      !TODO:count interface num
      this%nIntf=3
      if(allocated(this%Intf))then
        if(size(this%Intf)/=this%nIntf)then
          deallocate(this%Intf)
          allocate(this%Intf(this%nIntf))
        end if
      else
        allocate(this%Intf(this%nIntf))
      end if
      call reallocArr(this%IntfNeibBlock,INTF_NEIB_BLOCK_NUM,this%nIntf)
      !TODO:save interfaces
      this%isUpIntf=.true.
    end if
  end subroutine
  
  !> update the point position
  elemental subroutine updatePointPos(this)
    class(typeGrid),intent(inout)::this !< this grid
    
    if(.not.this%isUpPointPos)then
      call reallocArr(this%PointPos,DIMS,this%nPoint)
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
      call reallocArr(this%LinePos,DIMS,this%nLine)
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
      call reallocArr(this%FacetPos,DIMS,this%nFacet)
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
      call reallocArr(this%BlockPos,DIMS,this%nBlock)
      forall(i=1:this%nBlock)
        this%BlockPos(:,i)=sum(this%NodePos(:,this%Block(i)%iNode(:)),2)/dble(this%Block(i)%nNode)
      end forall
      this%isUpBlockPos=.true.
    end if
  end subroutine
  
  !> update the line length
  elemental subroutine updateLineLen(this)
    class(typeGrid),intent(inout)::this !< this grid
    
    if(.not.this%isUpLineLen)then
      call reallocArr(this%LineLen,this%nLine)
      do i=1,this%nLine
        select case(this%Line(i)%Shp)
        case(LINE_TYPE)
          this%LineLen(i)=norm2(this%NodePos(:,this%Line(i)%iNode(1))&
          &                    -this%NodePos(:,this%Line(i)%iNode(2)))
        case default
        end select
      end do
      this%isUpLineLen=.true.
    end if
  end subroutine
  
  !> update the facet area
  elemental subroutine updateFacetArea(this)
    use moduleSimpleGeometry
    class(typeGrid),intent(inout)::this !< this grid
    
    if(.not.this%isUpFacetArea)then
      call reallocArr(this%FacetArea,this%nFacet)
      do i=1,this%nFacet
        select case(this%Facet(i)%Shp)
        case(TRI_TYPE)
          this%FacetArea(i)=find3PArea(this%NodePos(:,this%Facet(i)%iNode(:)))
        case(QUAD_TYPE)
          this%FacetArea(i)=find4PArea(this%NodePos(:,this%Facet(i)%iNode(:)))
        case default
        end select
      end do
      this%isUpFacetArea=.true.
    end if
  end subroutine
  
  !> update the block volume
  elemental subroutine updateBlockVol(this)
    use moduleSimpleGeometry
    class(typeGrid),intent(inout)::this !< this grid
    
    if(.not.this%isUpBlockVol)then
      call reallocArr(this%BlockVol,this%nBlock)
      do i=1,this%nBlock
        select case(this%Block(i)%Shp)
        case(TET_TYPE)
          this%BlockVol(i)=find4PVol(this%NodePos(:,this%Block(i)%iNode(:)))
        case(HEX_TYPE)
          this%BlockVol(i)=find8PVol(this%NodePos(:,this%Block(i)%iNode(:)))
        case default
        end select
      end do
      this%isUpBlockVol=.true.
    end if
  end subroutine
  
  !> update the facet normal vector
  elemental subroutine updateFacetNorm(this)
    use moduleSimpleGeometry
    class(typeGrid),intent(inout)::this !< this grid
    
    if(.not.this%isUpFacetNorm)then
      call reallocArr(this%FacetNorm,DIMS,this%nFacet)
      do i=1,this%nFacet
        select case(this%Facet(i)%Shp)
        case(TRI_TYPE)
          this%FacetNorm(:,i)=find3PNorm(this%NodePos(:,this%Facet(i)%iNode(:)))
        case(QUAD_TYPE)
          this%FacetNorm(:,i)=find3PNorm(this%NodePos(:,this%Facet(i)%iNode([1,2,3])))
        case default
        end select
      end do
      this%isUpFacetNorm=.true.
    end if
  end subroutine
  
end module
