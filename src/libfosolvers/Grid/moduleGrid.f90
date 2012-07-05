!----------------------------------------------------------------------------- best with 100 columns

!> grid
module moduleGrid
  use moduleBasicDataStruct
  private
  
  ! constants
  integer,parameter,public::DIMS=3 !< dimensions
  
  integer,parameter,public::BIND_NODE=1 !< bind with node
  integer,parameter,public::BIND_POINT=2 !< bind with point
  integer,parameter,public::BIND_LINE=3 !< bind with line
  integer,parameter,public::BIND_FACET=4 !< bind with facet
  integer,parameter,public::BIND_BLOCK=5 !< bind with block
  integer,parameter,public::BIND_INTF=6 !< bind with interface
  
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
    integer nPrt !< number of partitions
    
    ! auxiliary grid data
    logical isUpNodeNeib !< if the node neighbour facet and block are updated
    type(typeHtr1DIArr),allocatable::NodeNeibFacet(:) !< node neighbour block
    type(typeHtr1DIArr),allocatable::NodeNeibBlock(:) !< node neighbour block
    
    logical isUpBlockNeib !< if the block neighbour facet and neighbour block is updated
    type(typeHtr1DIArr),allocatable::BlockNeibFacet(:) !< block neighbour facet
    type(typeHtr1DIArr),allocatable::BlockNeibBlock(:) !< block neighbour block
    
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
    
    logical isUpIntfPos !< if the interface position is updated
    double precision,allocatable::IntfPos(:,:) !< interface position
    
    logical isUpLineLen !< if the line length is updated
    double precision,allocatable::LineLen(:) !< line length
    
    logical isUpFacetArea !< if the facet area is updated
    double precision,allocatable::FacetArea(:) !< facet area
    
    logical isUpIntfArea !< if the interface area is updated
    double precision,allocatable::IntfArea(:) !< interface area
    
    logical isUpBlockVol !< if the block volume is updated
    double precision,allocatable::BlockVol(:) !< block volume
    
    logical isUpNodeVol !< if the node/median-dual block volume is updated
    double precision,allocatable::NodeVol(:) !< node/median-dual block volume
    
    logical isUpFacetNorm !< if the facet normal vector is updated
    double precision,allocatable::FacetNorm(:,:) !< facet normal vector
    
    logical isUpIntfNorm !< if the interface normal vector is updated
    double precision,allocatable::IntfNorm(:,:) !< interface normal vector
  contains
    ! basic grid procedures
    procedure,public::init=>initGrid
    procedure,public::clear=>clearGrid
    !FIXME:final::purgeGrid
    ! auxiliary grid procedures
    procedure,public::updateNodeNeib
    procedure,public::updateBlockNeib
    procedure,public::updateIntf
    procedure,public::updatePointPos
    procedure,public::updateLinePos
    procedure,public::updateFacetPos
    procedure,public::updateBlockPos
    procedure,public::updateIntfPos
    procedure,public::updateLineLen
    procedure,public::updateFacetArea
    procedure,public::updateIntfArea
    procedure,public::updateBlockVol
    procedure,public::updateNodeVol
    procedure,public::updateFacetNorm
    procedure,public::updateIntfNorm
  end type
  
  ! individual procedures
  public::getBlockSurfNum
  public::getBlockSurfNode
  
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
  
  !> get the number of surfaces of a block
  elemental function getBlockSurfNum(block)
    class(typeEle),intent(in)::block !< block to be get the number of surfaces of
    integer getBlockSurfNum !< number of surfaces
    
    select case(block%Shp)
    case(TET_TYPE)
      getBlockSurfNum=TET_SURF_NUM
    case(HEX_TYPE)
      getBlockSurfNum=HEX_SURF_NUM
    case default
      getBlockSurfNum=0
    end select
  end function
  
  !> get the nodes of the k_th surface of a block
  pure function getBlockSurfNode(block,k)
    class(typeEle),intent(in)::block !< block to be get the nodes of k_th surface of
    integer,intent(in)::k !< index of surface
    integer,allocatable::getBlockSurfNode(:) !< surface nodes
    
    select case(block%Shp)
    case(TET_TYPE)
      allocate(getBlockSurfNode(TRI_NODE_NUM))
      getBlockSurfNode(:)=block%iNode(TET_SURF_TAB(:,k))
    case(HEX_TYPE)
      allocate(getBlockSurfNode(QUAD_NODE_NUM))
      getBlockSurfNode(:)=block%iNode(HEX_SURF_TAB(:,k))
    case default
    end select
  end function
  
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
    this%isUpNodeNeib=.false.
    this%isUpBlockNeib=.false.
    this%isUpIntf=.false.
    this%isUpPointPos=.false.
    this%isUpLinePos=.false.
    this%isUpFacetPos=.false.
    this%isUpBlockPos=.false.
    this%isUpIntfPos=.false.
    this%isUpLineLen=.false.
    this%isUpFacetArea=.false.
    this%isUpIntfArea=.false.
    this%isUpBlockVol=.false.
    this%isUpNodeVol=.false.
    this%isUpFacetNorm=.false.
    this%isUpIntfNorm=.false.
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
    if(allocated(this%NodeNeibFacet)) deallocate(this%NodeNeibFacet)
    if(allocated(this%NodeNeibBlock)) deallocate(this%NodeNeibBlock)
    if(allocated(this%BlockNeibFacet)) deallocate(this%BlockNeibFacet)
    if(allocated(this%BlockNeibBlock)) deallocate(this%BlockNeibBlock)
    if(allocated(this%Intf)) deallocate(this%Intf)
    if(allocated(this%IntfNeibBlock)) deallocate(this%IntfNeibBlock)
    if(allocated(this%PointPos)) deallocate(this%PointPos)
    if(allocated(this%LinePos)) deallocate(this%LinePos)
    if(allocated(this%FacetPos)) deallocate(this%FacetPos)
    if(allocated(this%BlockPos)) deallocate(this%BlockPos)
    if(allocated(this%IntfPos)) deallocate(this%IntfPos)
    if(allocated(this%LineLen)) deallocate(this%LineLen)
    if(allocated(this%FacetArea)) deallocate(this%FacetArea)
    if(allocated(this%IntfArea)) deallocate(this%IntfArea)
    if(allocated(this%BlockVol)) deallocate(this%BlockVol)
    if(allocated(this%NodeVol)) deallocate(this%NodeVol)
    if(allocated(this%FacetNorm)) deallocate(this%FacetNorm)
    if(allocated(this%IntfNorm)) deallocate(this%IntfNorm)
  end subroutine
  
  !> destructor of typeGrid
  elemental subroutine purgeGrid(this)
    type(typeGrid),intent(inout)::this !< this grid
    
    call this%clear()
  end subroutine
  
  !> update the node neighbour facet and block
  elemental subroutine updateNodeNeib(this)
    class(typeGrid),intent(inout)::this !< this grid
    
    if(.not.this%isUpNodeNeib)then
      call reallocArr(this%NodeNeibFacet,this%nNode)
      do i=1,this%nFacet
        do j=1,this%Facet(i)%nNode
          call pushArr(this%NodeNeibFacet(this%Facet(i)%iNode(j))%dat,i)
        end do
      end do
      call reallocArr(this%NodeNeibBlock,this%nNode)
      do i=1,this%nBlock
        do j=1,this%Block(i)%nNode
          call pushArr(this%NodeNeibBlock(this%Block(i)%iNode(j))%dat,i)
        end do
      end do
      this%isUpNodeNeib=.true.
    end if
  end subroutine
  
  !> update the block neighbour facet and neighbour block
  elemental subroutine updateBlockNeib(this)
    class(typeGrid),intent(inout)::this !< this grid
    integer,allocatable::SurfNode(:),ScanList(:)
    logical,allocatable::mask(:)
    
    if(.not.this%isUpBlockNeib)then
      call reallocArr(this%BlockNeibFacet,this%nBlock)
      call reallocArr(this%BlockNeibBlock,this%nBlock)
      call this%updateNodeNeib()
      do i=1,this%nBlock
        m=getBlockSurfNum(this%Block(i))
        call reallocArr(this%BlockNeibFacet(i)%dat,m)
        this%BlockNeibFacet(i)%dat(:)=0
        call reallocArr(this%BlockNeibBlock(i)%dat,m)
        this%BlockNeibBlock(i)%dat(:)=0
        do j=1,m
          SurfNode=getBlockSurfNode(this%Block(i),j)
          allocate(mask(size(SurfNode)))
          if(allocated(this%NodeNeibBlock(SurfNode(1))%dat))then
            ScanList=this%NodeNeibBlock(SurfNode(1))%dat
            do k=1,size(ScanList)
              if(ScanList(k)/=i)then
                mask(:)=.false.
                forall(l=1:size(SurfNode))
                  mask(l)=any(this%Block(ScanList(k))%iNode(:)==SurfNode(l))
                end forall
                if(all(mask(:)))then
                  this%BlockNeibBlock(i)%dat(j)=ScanList(k)
                  exit
                end if
              end if
            end do
            deallocate(ScanList)
          end if
          if(allocated(this%NodeNeibFacet(SurfNode(1))%dat))then
            ScanList=this%NodeNeibFacet(SurfNode(1))%dat
            do k=1,size(ScanList)
              mask(:)=.false.
              forall(l=1:size(SurfNode))
                mask(l)=any(this%Facet(ScanList(k))%iNode(:)==SurfNode(l))
              end forall
              if(all(mask(:)))then
                this%BlockNeibFacet(i)%dat(j)=ScanList(k)
                if(.not.(this%BlockNeibBlock(i)%dat(j)/=0.and.&
                &        this%Facet(this%BlockNeibFacet(i)%dat(j))%Ent>=0))then
                  exit
                end if
              end if
            end do
            deallocate(ScanList)
          end if
          deallocate(SurfNode)
          deallocate(mask)
        end do
      end do
      this%isUpBlockNeib=.true.
    end if
  end subroutine
  
  !> update the interface between blocks
  elemental subroutine updateIntf(this)
    class(typeGrid),intent(inout)::this !< this grid
    
    if(.not.this%isUpIntf)then
      call this%updateBlockNeib()
      this%nIntf=0
      do i=1,this%nBlock
        do j=1,getBlockSurfNum(this%Block(i))
          if(this%BlockNeibBlock(i)%dat(j)>i.and.this%BlockNeibFacet(i)%dat(j)==0)then
            this%nIntf=this%nIntf+1
          end if
        end do
      end do
      if(allocated(this%Intf))then
        if(size(this%Intf)/=this%nIntf)then
          deallocate(this%Intf)
          allocate(this%Intf(this%nIntf))
        end if
      else
        allocate(this%Intf(this%nIntf))
      end if
      call reallocArr(this%IntfNeibBlock,INTF_NEIB_BLOCK_NUM,this%nIntf)
      k=0
      do i=1,this%nBlock
        do j=1,getBlockSurfNum(this%Block(i))
          if(this%BlockNeibBlock(i)%dat(j)>i.and.this%BlockNeibFacet(i)%dat(j)==0)then
            k=k+1
            select case(this%Block(i)%Shp)
            case(TET_TYPE)
              this%Intf(k)%Shp=TRI_TYPE
              this%Intf(k)%nNode=TRI_NODE_NUM
            case(HEX_TYPE)
              this%Intf(k)%Shp=QUAD_TYPE
              this%Intf(k)%nNode=QUAD_NODE_NUM
            case default
            end select
            call reallocArr(this%Intf(k)%iNode,this%Intf(k)%nNode)
            this%Intf(k)%iNode=getBlockSurfNode(this%Block(i),j)
            call reallocArr(this%Intf(k)%Dmn,size(this%Block(i)%Dmn))
            this%Intf(k)%Dmn(:)=this%Block(i)%Dmn(:)
            call reallocArr(this%Intf(k)%Prt,1)
            this%Intf(k)%Prt(1)=maxval(this%Block(i)%Prt(:))
            this%IntfNeibBlock(:,k)=[i,this%BlockNeibBlock(i)%dat(j)]
          end if
        end do
      end do
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
  
  !> update the interface position
  elemental subroutine updateIntfPos(this)
    class(typeGrid),intent(inout)::this !< this grid
    
    if(.not.this%isUpIntfPos)then
      call this%updateIntf()
      call reallocArr(this%IntfPos,DIMS,this%nIntf)
      forall(i=1:this%nIntf)
        this%IntfPos(:,i)=sum(this%NodePos(:,this%Intf(i)%iNode(:)),2)/dble(this%Intf(i)%nNode)
      end forall
      this%isUpIntfPos=.true.
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
  
  !> update the interface area
  elemental subroutine updateIntfArea(this)
    use moduleSimpleGeometry
    class(typeGrid),intent(inout)::this !< this grid
    
    if(.not.this%isUpIntfArea)then
      call this%updateIntf()
      call reallocArr(this%IntfArea,this%nIntf)
      do i=1,this%nIntf
        select case(this%Intf(i)%Shp)
        case(TRI_TYPE)
          this%IntfArea(i)=find3PArea(this%NodePos(:,this%Intf(i)%iNode(:)))
        case(QUAD_TYPE)
          this%IntfArea(i)=find4PArea(this%NodePos(:,this%Intf(i)%iNode(:)))
        case default
        end select
      end do
      this%isUpIntfArea=.true.
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
  
  !> update the node/median-dual block volume
  elemental subroutine updateNodeVol(this)
    use moduleSimpleGeometry
    class(typeGrid),intent(inout)::this !< this grid
    
    if(.not.this%isUpNodeVol)then
      call this%updateBlockVol()
      call reallocArr(this%NodeVol,this%nNode)
      this%NodeVol(:)=0d0
      do i=1,this%nBlock
        this%NodeVol(this%Block(i)%iNode(:))=this%NodeVol(this%Block(i)%iNode(:))&
        &                                    +this%BlockVol(i)/dble(this%Block(i)%nNode)
      end do
      this%isUpNodeVol=.true.
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
  
  !> update the interface normal vector
  elemental subroutine updateIntfNorm(this)
    use moduleSimpleGeometry
    class(typeGrid),intent(inout)::this !< this grid
    
    if(.not.this%isUpIntfNorm)then
      call this%updateIntf()
      call reallocArr(this%IntfNorm,DIMS,this%nIntf)
      do i=1,this%nIntf
        select case(this%Intf(i)%Shp)
        case(TRI_TYPE)
          this%IntfNorm(:,i)=find3PNorm(this%NodePos(:,this%Intf(i)%iNode(:)))
        case(QUAD_TYPE)
          this%IntfNorm(:,i)=find3PNorm(this%NodePos(:,this%Intf(i)%iNode([1,2,3])))
        case default
        end select
      end do
      this%isUpIntfNorm=.true.
    end if
  end subroutine
  
end module
