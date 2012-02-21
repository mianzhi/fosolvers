!----------------------------------------------------------------------------- best with 100 columns

!*****************************************
! data and elementary procedures for grid
!*****************************************
module moduleGrid
  use moduleMiscDataStruct
  use moduleSimpleGeometry
  private
  
  ! some constants
  integer,parameter,public::DIMS=3
  
  integer,parameter,public::POINT_TYPE=15
  integer,parameter,public::LINE_TYPE=1
  integer,parameter,public::TRI_TYPE=2
  integer,parameter,public::QUAD_TYPE=3
  integer,parameter,public::TET_TYPE=4
  integer,parameter,public::HEX_TYPE=5
  
  integer,parameter,public::POINT_NODE_NUM=1
  integer,parameter,public::LINE_NODE_NUM=2
  integer,parameter,public::TRI_NODE_NUM=3
  integer,parameter,public::QUAD_NODE_NUM=4
  integer,parameter,public::TET_NODE_NUM=4
  integer,parameter,public::HEX_NODE_NUM=8
  
  integer,parameter,public::TET_SURF_NUM=4
  integer,parameter,public::HEX_SURF_NUM=6
  
  integer,parameter,public::TET_EDGE_NUM=6
  integer,parameter,public::HEX_EDGE_NUM=12
  
  integer,parameter,public::FACET_NEIB_BLOCK_NUM=2
  
  integer,parameter,public::BIND_NODE=1
  integer,parameter,public::BIND_FACET=2
  integer,parameter,public::BIND_BLOCK=3
  
  ! table of surface nodes for all kinds of block
  integer,public::TET_SURF_TAB(TET_SURF_NUM,TRI_NODE_NUM),HEX_SURF_TAB(HEX_SURF_NUM,QUAD_NODE_NUM)
  parameter(TET_SURF_TAB=reshape([1,1,1,2,3,2,4,3,2,4,3,4],[TET_SURF_NUM,TRI_NODE_NUM]))
  parameter(HEX_SURF_TAB=reshape([2,1,3,1,5,1,3,5,4,2,6,4,7,8,8,6,7,3,6,4,7,5,8,2],&
  &                              [HEX_SURF_NUM,QUAD_NODE_NUM]))
  
  ! table of edge nodes for all kinds of block
  integer,public::TET_EDGE_TAB(TET_EDGE_NUM,LINE_NODE_NUM),HEX_EDGE_TAB(HEX_EDGE_NUM,LINE_NODE_NUM)
  parameter(TET_EDGE_TAB=reshape([1,1,1,2,3,4,2,3,4,3,4,2],[TET_EDGE_NUM,LINE_NODE_NUM]))
  parameter(HEX_EDGE_TAB=reshape([1,2,3,4,5,6,7,8,1,2,3,4,2,3,4,1,6,7,8,5,5,6,7,8],&
  &                              [HEX_EDGE_NUM,LINE_NODE_NUM]))
  
  ! bounding box of the geometry
  double precision,public,save::BoundBox(2,DIMS)
  
  ! stand-along procedures
  public updateAllNodeFacetInd
  public updateAllNodeBlockInd
  public updateGrid
  
  !----------
  ! typeNode
  !----------
  type,public::typeNode
    integer Ind
    double precision Pos(DIMS)
    integer,allocatable::FacetInd(:)
    integer,allocatable::BlockInd(:)
    !integer,allocatable::NeibNode(:)
  contains
    procedure,public::updateFacetInd=>updateNodeFacetInd
    procedure,public::updateBlockInd=>updateNodeBlockInd
    !procedure,public::updateNeibNode=>updateNodeNeibNode
    !TODO: wait for gcc to implement
    !final::cleanNode
  end type
  type(typeNode),public,allocatable,save::Node(:)
  integer,public,save::nNode
  
  !-----------
  ! typePoint
  !-----------
  type,public::typePoint
    integer Ind
    integer NodeInd
    integer GeoEnti
    double precision Pos(DIMS)
  contains
    procedure,public::updatePos=>updatePointPos
  end type
  type(typePoint),public,allocatable,save::Point(:)
  integer,public,save::nPoint
  
  !----------
  ! typeLine
  !----------
  type,public::typeLine
    integer Ind
    integer ShapeType
    integer NodeNum
    integer,allocatable::NodeInd(:)
    integer GeoEnti
    double precision PC(DIMS)
    double precision Length
  contains
    procedure,public::specify=>specifyLine
    procedure,public::updatePC=>updateLinePC
    procedure,public::updateLength=>updateLineLength
    !TODO:final::cleanLine
  end type
  type(typeLine),public,allocatable,save::Line(:)
  integer,public,save::nLine
  
  !-----------
  ! typeFacet
  !-----------
  type,public::typeFacet
    integer Ind
    integer ShapeType
    integer NodeNum
    integer,allocatable::NodeInd(:)
    integer GeoEnti
    integer NeibBlock(FACET_NEIB_BLOCK_NUM)
    double precision PC(DIMS)
    double precision Area
    double precision Norm(DIMS)
  contains
    procedure,public::specify=>specifyFacet
    procedure,public::updatePC=>updateFacetPC
    procedure,public::updateArea=>updateFacetArea
    procedure,public::updateNorm=>updateFacetNorm
    procedure,public::updateNeibBlock=>updateFacetNeibBlock
    !TODO:final::cleanFacet
  end type
  type(typeFacet),public,allocatable,save::Facet(:)
  integer,public,save::nFacet
  
  !-----------
  ! typeBlock
  !-----------
  type,public::typeBlock
    integer Ind
    integer ShapeType
    integer NodeNum
    integer EdgeNum
    integer SurfNum
    integer,allocatable::NodeInd(:)
    integer GeoEnti
    integer,allocatable::Prt(:)
    integer,allocatable::Neib(:)
    double precision PC(DIMS)
    double precision Vol
    integer,allocatable::EdgeNodeInd(:,:)
    integer,allocatable::SurfNodeInd(:,:)
    double precision,allocatable::SurfPC(:,:)
    double precision,allocatable::SurfArea(:)
    double precision,allocatable::SurfNorm(:,:)
  contains
    procedure,public::specify=>specifyBlock
    procedure,public::updatePC=>updateBlockPC
    procedure,public::updateVol=>updateBlockVol
    procedure,public::updateNeib=>updateBlockNeib
    procedure,public::updateEdgeNodeInd=>updateBlockEdgeNodeInd
    procedure,public::updateSurfNodeInd=>updateBlockSurfNodeInd
    procedure,public::updateSurfPC=>updateBlockSurfPC
    procedure,public::updateSurfArea=>updateBlockSurfArea
    procedure,public::updateSurfNorm=>updateBlockSurfNorm
    !TODO:final::cleanBlock
  end type
  type(typeBlock),public,allocatable,save::Block(:)
  integer,public,save::nBlock,nPrt
  
contains
  
  !-------------------------------------------------------
  ! update the index of facets in which this node is used
  !-------------------------------------------------------
  ! Note: not recommend to use in common cases. should use updateAllNodeFacetInd()
  elemental subroutine updateNodeFacetInd(this)
    class(typeNode),intent(inout)::this
    
    if(allocated(this%FacetInd))then
      deallocate(this%FacetInd)
    end if
    do i=1,nFacet
      if(any(Facet(i)%NodeInd==this%Ind))then
        call extendArray(this%FacetInd,1)
        this%FacetInd(size(this%FacetInd))=i
      end if
    end do
  end subroutine
  
  !--------------------------------------------------------------------
  ! update the index of facets in which one node is used for all nodes
  !--------------------------------------------------------------------
  subroutine updateAllNodeFacetInd()
    !$omp parallel do
    do i=1,nNode
      if(allocated(Node(i)%FacetInd))then
        deallocate(Node(i)%FacetInd)
      end if
    end do
    !$omp end parallel do
    do i=1,nFacet
      do j=1,Facet(i)%NodeNum
        call extendArray(Node(Facet(i)%NodeInd(j))%FacetInd,1)
        Node(Facet(i)%NodeInd(j))%FacetInd(size(Node(Facet(i)%NodeInd(j))%FacetInd))=i
      end do
    end do
  end subroutine
  
  !-------------------------------------------------------
  ! update the index of blocks in which this node is used
  !-------------------------------------------------------
  ! Note: not recommend to use in common cases. should use updateAllNodeBlockInd()
  elemental subroutine updateNodeBlockInd(this)
    class(typeNode),intent(inout)::this
    
    if(allocated(this%BlockInd))then
      deallocate(this%BlockInd)
    end if
    do i=1,nBlock
      if(any(Block(i)%NodeInd==this%Ind))then
        call extendArray(this%BlockInd,1)
        this%BlockInd(size(this%BlockInd))=i
      end if
    end do
  end subroutine
  
  !--------------------------------------------------------------------
  ! update the index of blocks in which one node is used for all nodes
  !--------------------------------------------------------------------
  subroutine updateAllNodeBlockInd()
    !$omp parallel do
    do i=1,nNode
      if(allocated(Node(i)%BlockInd))then
        deallocate(Node(i)%BlockInd)
      end if
    end do
    !$omp end parallel do
    do i=1,nBlock
      do j=1,Block(i)%NodeNum
        call extendArray(Node(Block(i)%NodeInd(j))%BlockInd,1)
        Node(Block(i)%NodeInd(j))%BlockInd(size(Node(Block(i)%NodeInd(j))%BlockInd))=i
      end do
    end do
  end subroutine
  
  !------------------------------------------------------------------
  ! update the index of neighbour nodes (the other end of the edges)
  !------------------------------------------------------------------
  elemental subroutine updateNodeNeibNode(this)
    class(typeNode),intent(inout)::this
    
    if(allocated(this%BlockInd))then
      do i=1,size(this%BlockInd)
        !TODO:implement
      end do
    end if
  end subroutine
  
  !------------------------
  ! destructor of typeNode
  !------------------------
  elemental subroutine cleanNode(this)
    type(typeNode),intent(inout)::this
    
    if(allocated(this%FacetInd))then
      deallocate(this%FacetInd)
    end if
    if(allocated(this%BlockInd))then
      deallocate(this%BlockInd)
    end if
  end subroutine
  
  !-----------------------------------
  ! update the position of this point
  !-----------------------------------
  elemental subroutine updatePointPos(this)
    class(typePoint),intent(inout)::this
    
    this%Pos(:)=Node(this%NodeInd)%Pos(:)
  end subroutine
  
  !------------------------------
  ! specify a shape to this line
  !------------------------------
  elemental subroutine specifyLine(this,shapetype)
    use moduleUtility
    class(typeLine),intent(inout)::this
    integer,intent(in)::shapetype
    integer::nodenum
    
    nodenum=0
    select case(shapetype)
      case(LINE_TYPE)
        nodenum=LINE_NODE_NUM
      case default
    end select
    this%ShapeType=shapetype
    this%NodeNum=nodenum
    if(allocated(this%NodeInd))then
      deallocate(this%NodeInd)
    end if
    allocate(this%NodeInd(nodenum))
  end subroutine
  
  !-----------------------------------------
  ! update the centre position of this line
  !-----------------------------------------
  elemental subroutine updateLinePC(this)
    class(typeLine),intent(inout)::this
    
    this%PC(:)=[(sum(Node(this%NodeInd(:))%Pos(i)),i=1,DIMS)]/dble(this%NodeNum)
  end subroutine
  
  !--------------------------------
  ! update the length of this line
  !--------------------------------
  elemental subroutine updateLineLength(this)
    class(typeLine),intent(inout)::this
    
    select case(this%ShapeType)
      case(LINE_TYPE)
        this%Length=norm2(Node(this%NodeInd(1))%Pos(:)-Node(this%NodeInd(2))%Pos(:))
      case default
    end select
  end subroutine
  
  !------------------------
  ! destructor of typeLine
  !------------------------
  elemental subroutine cleanLine(this)
    type(typeLine),intent(inout)::this
    
    if(allocated(this%NodeInd))then
      deallocate(this%NodeInd)
    end if
  end subroutine
  
  !-------------------------------
  ! specify a shape to this facet
  !-------------------------------
  elemental subroutine specifyFacet(this,shapetype)
    use moduleUtility
    class(typeFacet),intent(inout)::this
    integer,intent(in)::shapetype
    integer::nodenum
    
    nodenum=0
    select case(shapetype)
      case(TRI_TYPE)
        nodenum=TRI_NODE_NUM
      case(QUAD_TYPE)
        nodenum=QUAD_NODE_NUM
      case default
    end select
    this%ShapeType=shapetype
    this%NodeNum=nodenum
    if(allocated(this%NodeInd))then
      deallocate(this%NodeInd)
    end if
    allocate(this%NodeInd(nodenum))
  end subroutine
  
  !------------------------------------------
  ! update the centre position of this facet
  !------------------------------------------
  elemental subroutine updateFacetPC(this)
    class(typeFacet),intent(inout)::this
    
    this%PC(:)=[(sum(Node(this%NodeInd(:))%Pos(i)),i=1,DIMS)]/dble(this%NodeNum)
  end subroutine
  
  !-------------------------------
  ! update the area of this facet
  !-------------------------------
  elemental subroutine updateFacetArea(this)
    class(typeFacet),intent(inout)::this
    
    select case(this%ShapeType)
      case(TRI_TYPE)
        this%Area=find3PArea(Node(this%NodeInd(1))%Pos(:),Node(this%NodeInd(2))%Pos(:),&
        &                    Node(this%NodeInd(3))%Pos(:))
      case(QUAD_TYPE)
        this%Area=find3PArea(Node(this%NodeInd(1))%Pos(:),Node(this%NodeInd(2))%Pos(:),&
        &                    Node(this%NodeInd(3))%Pos(:))&
        &        +find3PArea(Node(this%NodeInd(3))%Pos(:),Node(this%NodeInd(4))%Pos(:),&
        &                    Node(this%NodeInd(1))%Pos(:))
      case default
    end select
  end subroutine
  
  !----------------------------------------
  ! update the normal vector of this facet
  !----------------------------------------
  elemental subroutine updateFacetNorm(this)
    class(typeFacet),intent(inout)::this
    
    this%Norm(:)=find3PNorm(Node(this%NodeInd(1))%Pos(:),Node(this%NodeInd(2))%Pos(:),&
    &                       Node(this%NodeInd(3))%Pos(:))
  end subroutine
  
  !------------------------------------------
  ! update the neighbour block of this facet
  !------------------------------------------
  ! Note: one facet may have 2 neighbour blocks
  elemental subroutine updateFacetNeibBlock(this)
    class(typeFacet),intent(inout)::this
    logical,allocatable::mask(:)
    
    allocate(mask(this%NodeNum))
    this%NeibBlock(:)=0
    l=1
    if(allocated(Node(this%NodeInd(1))%BlockInd))then
      do i=1,size(Node(this%NodeInd(1))%BlockInd)
        k=Node(this%NodeInd(1))%BlockInd(i)
        mask(:)=.false.
        forall(j=1:this%NodeNum)
          mask(j)=any(Block(k)%NodeInd(:)==this%NodeInd(j))
        end forall
        if(all(mask(:)))then
          this%NeibBlock(l)=k
          if(l==FACET_NEIB_BLOCK_NUM)then
            exit
          end if
          l=l+1
        end if
      end do
    end if
    deallocate(mask)
  end subroutine
  
  !-------------------------
  ! destructor of typeFacet
  !-------------------------
  elemental subroutine cleanFacet(this)
    type(typeFacet),intent(inout)::this
    
    if(allocated(this%NodeInd))then
      deallocate(this%NodeInd)
    end if
  end subroutine
  
  !-------------------------------
  ! specify a shape to this block
  !-------------------------------
  elemental subroutine specifyBlock(this,shapetype)
    use moduleUtility
    class(typeBlock),intent(inout)::this
    integer,intent(in)::shapetype
    integer::nodenum,surfnum,surfnodenum,edgenum,edgenodenum
    
    nodenum=0
    surfnum=0
    surfnodenum=0
    edgenum=0
    edgenodenum=0
    select case(shapetype)
      case(TET_TYPE)
        nodenum=TET_NODE_NUM
        surfnum=TET_SURF_NUM
        surfnodenum=TRI_NODE_NUM
        edgenum=TET_EDGE_NUM
        edgenodenum=LINE_NODE_NUM
      case(HEX_TYPE)
        nodenum=HEX_NODE_NUM
        surfnum=HEX_SURF_NUM
        surfnodenum=QUAD_NODE_NUM
        edgenum=HEX_EDGE_NUM
        edgenodenum=LINE_NODE_NUM
      case default
    end select
    this%ShapeType=shapetype
    this%NodeNum=nodenum
    this%EdgeNum=edgenum
    this%SurfNum=surfnum
    if(allocated(this%NodeInd))then
      deallocate(this%NodeInd)
    end if
    allocate(this%NodeInd(nodenum))
    if(allocated(this%Neib))then
      deallocate(this%Neib)
    end if
    allocate(this%Neib(surfnum))
    if(allocated(this%EdgeNodeInd))then
      deallocate(this%EdgeNodeInd)
    end if
    allocate(this%EdgeNodeInd(edgenum,edgenodenum))
    if(allocated(this%SurfNodeInd))then
      deallocate(this%SurfNodeInd)
    end if
    allocate(this%SurfNodeInd(surfnum,surfnodenum))
    if(allocated(this%SurfPC))then
      deallocate(this%SurfPC)
    end if
    allocate(this%SurfPC(surfnum,DIMS))
    if(allocated(this%SurfArea))then
      deallocate(this%SurfArea)
    end if
    allocate(this%SurfArea(surfnum))
    if(allocated(this%SurfNorm))then
      deallocate(this%SurfNorm)
    end if
    allocate(this%SurfNorm(surfnum,DIMS))
  end subroutine
  
  !------------------------------------------
  ! update the centre position of this block
  !------------------------------------------
  elemental subroutine updateBlockPC(this)
    class(typeBlock),intent(inout)::this
    
    this%PC(:)=[(sum(Node(this%NodeInd(:))%Pos(i)),i=1,DIMS)]/dble(this%NodeNum)
  end subroutine
  
  !---------------------------------
  ! update the volume of this block
  !---------------------------------
  elemental subroutine updateBlockVol(this)
    class(typeBlock),intent(inout)::this
    
    select case(this%ShapeType)
      case(TET_TYPE)
        this%Vol=find4PVol(Node(this%NodeInd(1))%Pos(:),Node(this%NodeInd(2))%Pos(:),&
        &                  Node(this%NodeInd(3))%Pos(:),Node(this%NodeInd(4))%Pos(:))
      case(HEX_TYPE)
        this%Vol=find4PVol(Node(this%NodeInd(1))%Pos(:),Node(this%NodeInd(2))%Pos(:),&
        &                  Node(this%NodeInd(3))%Pos(:),Node(this%NodeInd(6))%Pos(:))&
        &       +find4PVol(Node(this%NodeInd(1))%Pos(:),Node(this%NodeInd(3))%Pos(:),&
        &                  Node(this%NodeInd(4))%Pos(:),Node(this%NodeInd(6))%Pos(:))&
        &       +find4PVol(Node(this%NodeInd(1))%Pos(:),Node(this%NodeInd(4))%Pos(:),&
        &                  Node(this%NodeInd(5))%Pos(:),Node(this%NodeInd(6))%Pos(:))&
        &       +find4PVol(Node(this%NodeInd(4))%Pos(:),Node(this%NodeInd(5))%Pos(:),&
        &                  Node(this%NodeInd(6))%Pos(:),Node(this%NodeInd(7))%Pos(:))&
        &       +find4PVol(Node(this%NodeInd(3))%Pos(:),Node(this%NodeInd(4))%Pos(:),&
        &                  Node(this%NodeInd(6))%Pos(:),Node(this%NodeInd(7))%Pos(:))&
        &       +find4PVol(Node(this%NodeInd(4))%Pos(:),Node(this%NodeInd(5))%Pos(:),&
        &                  Node(this%NodeInd(7))%Pos(:),Node(this%NodeInd(8))%Pos(:))
      case default
    end select
  end subroutine
  
  !-----------------------------------------------------
  ! update the neighbour blocks or facets of this block
  !-----------------------------------------------------
  ! Note: if facet are found, the negative value of the facet index would be stored.
  ! Note: if a facet between two partitions is found, then this facet is the neighbour,
  !       otherwise if a block is found, then this block is the neighbour,
  !       otherwise if a boundary facet is found, then this facet is the neighbour.
  elemental subroutine updateBlockNeib(this)
    class(typeBlock),intent(inout)::this
    logical,allocatable::mask(:)
    integer,allocatable::scanlist(:)
    
    allocate(mask(this%NodeNum))
    this%Neib(:)=0
    
    allocate(scanlist(0))
    do i=1,this%NodeNum
      if(allocated(Node(this%NodeInd(i))%BlockInd))then
        do j=1,size(Node(this%NodeInd(i))%BlockInd)
          if(all(scanlist(:)/=Node(this%NodeInd(i))%BlockInd(j))&
          &  .and.Node(this%NodeInd(i))%BlockInd(j)/=this%Ind)then
            call extendArray(scanlist,1)
            scanlist(size(scanlist))=Node(this%NodeInd(i))%BlockInd(j)
          end if
        end do
      end if
    end do
    do i=1,size(scanlist)
      k=scanlist(i)
      mask(:)=.false.
      forall(j=1:this%NodeNum)
        mask(j)=any(Block(k)%NodeInd(:)==this%NodeInd(j))
      end forall
      do j=1,this%SurfNum
        select case(this%ShapeType)
          case(TET_TYPE)
            if(all(mask(TET_SURF_TAB(j,:))))then
              this%Neib(j)=k
              exit
            end if
          case(HEX_TYPE)
            if(all(mask(HEX_SURF_TAB(j,:))))then
              this%Neib(j)=k
              exit
            end if
          case default
        end select
      end do
      if(all(this%Neib(:)/=0))then
        exit
      end if
    end do
    deallocate(scanlist)
    
    allocate(scanlist(0))
    do i=1,this%NodeNum
      if(allocated(Node(this%NodeInd(i))%FacetInd))then
        do j=1,size(Node(this%NodeInd(i))%FacetInd)
          if(all(scanlist(:)/=Node(this%NodeInd(i))%FacetInd(j)))then
            call extendArray(scanlist,1)
            scanlist(size(scanlist))=Node(this%NodeInd(i))%FacetInd(j)
          end if
        end do
      end if
    end do
    do i=1,size(scanlist)
      k=scanlist(i)
      mask(:)=.false.
      forall(j=1:this%NodeNum)
        mask(j)=any(Facet(k)%NodeInd(:)==this%NodeInd(j))
      end forall
      do j=1,this%SurfNum
        select case(this%ShapeType)
          case(TET_TYPE)
            if(all(mask(TET_SURF_TAB(j,:))))then
              if(this%Neib(j)==0.or.Facet(k)%GeoEnti<0)then ! partition interface has GeoEnti<0
                this%Neib(j)=-k
                exit
              end if
            end if
          case(HEX_TYPE)
            if(all(mask(HEX_SURF_TAB(j,:))))then
              if(this%Neib(j)==0.or.Facet(k)%GeoEnti<0)then ! partition interface has GeoEnti<0
                this%Neib(j)=-k
                exit
              end if
            end if
          case default
        end select
      end do
    end do
    deallocate(scanlist)
    
    deallocate(mask)
  end subroutine
  
  !--------------------------------------------------
  ! update the node index of the edges of this block
  !--------------------------------------------------
  elemental subroutine updateBlockEdgeNodeInd(this)
    class(typeBlock),intent(inout)::this
    
    select case(this%ShapeType)
      case(TET_TYPE)
        forall(i=1:this%EdgeNum)
          this%EdgeNodeInd(i,:)=this%NodeInd(TET_EDGE_TAB(i,:))
        end forall
      case(HEX_TYPE)
        forall(i=1:this%EdgeNum)
          this%EdgeNodeInd(i,:)=this%NodeInd(HEX_EDGE_TAB(i,:))
        end forall
      case default
    end select
  end subroutine
  
  !-----------------------------------------------------
  ! update the node index of the surfaces of this block
  !-----------------------------------------------------
  elemental subroutine updateBlockSurfNodeInd(this)
    class(typeBlock),intent(inout)::this
    
    select case(this%ShapeType)
      case(TET_TYPE)
        forall(i=1:this%SurfNum)
          this%SurfNodeInd(i,:)=this%NodeInd(TET_SURF_TAB(i,:))
        end forall
      case(HEX_TYPE)
        forall(i=1:this%SurfNum)
          this%SurfNodeInd(i,:)=this%NodeInd(HEX_SURF_TAB(i,:))
        end forall
      case default
    end select
  end subroutine
  
  !----------------------------------------------------------
  ! update the centre position of the surfaces of this block
  !----------------------------------------------------------
  elemental subroutine updateBlockSurfPC(this)
    class(typeBlock),intent(inout)::this
    type(typeFacet)::tempFacet
    
    select case(this%ShapeType)
      case(TET_TYPE)
        call tempFacet%specify(TRI_TYPE)
        do i=1,this%SurfNum
          tempFacet%NodeInd(:)=this%SurfNodeInd(i,:)
          call tempFacet%updatePC
          this%SurfPC(i,:)=tempFacet%PC(:)
        end do
      case(HEX_TYPE)
        call tempFacet%specify(QUAD_TYPE)
        do i=1,this%SurfNum
          tempFacet%NodeInd(:)=this%SurfNodeInd(i,:)
          call tempFacet%updatePC
          this%SurfPC(i,:)=tempFacet%PC(:)
        end do
      case default
    end select
  end subroutine
  
  !----------------------------------------
  ! update the surface areas of this block
  !----------------------------------------
  elemental subroutine updateBlockSurfArea(this)
    class(typeBlock),intent(inout)::this
    type(typeFacet)::tempFacet
    
    select case(this%ShapeType)
      case(TET_TYPE)
        call tempFacet%specify(TRI_TYPE)
        do i=1,this%SurfNum
          tempFacet%NodeInd(:)=this%SurfNodeInd(i,:)
          call tempFacet%updateArea
          this%SurfArea(i)=tempFacet%Area
        end do
      case(HEX_TYPE)
        call tempFacet%specify(QUAD_TYPE)
        do i=1,this%SurfNum
          tempFacet%NodeInd(:)=this%SurfNodeInd(i,:)
          call tempFacet%updateArea
          this%SurfArea(i)=tempFacet%Area
        end do
      case default
    end select
  end subroutine
  
  !--------------------------------------------------------
  ! update the normal vector of the surfaces of this block
  !--------------------------------------------------------
  elemental subroutine updateBlockSurfNorm(this)
    class(typeBlock),intent(inout)::this
    type(typeFacet)::tempFacet
    
    select case(this%ShapeType)
      case(TET_TYPE)
        call tempFacet%specify(TRI_TYPE)
        do i=1,this%SurfNum
          tempFacet%NodeInd(:)=this%SurfNodeInd(i,:)
          call tempFacet%updateNorm
          this%SurfNorm(i,:)=tempFacet%Norm(:)
        end do
      case(HEX_TYPE)
        call tempFacet%specify(QUAD_TYPE)
        do i=1,this%SurfNum
          tempFacet%NodeInd(:)=this%SurfNodeInd(i,:)
          call tempFacet%updateNorm
          this%SurfNorm(i,:)=tempFacet%Norm(:)
        end do
      case default
    end select
  end subroutine
  
  !-------------------------
  ! destructor of typeBlock
  !-------------------------
  elemental subroutine cleanBlock(this)
    type(typeBlock),intent(inout)::this
    
    if(allocated(this%NodeInd))then
      deallocate(this%NodeInd)
    end if
    if(allocated(this%Prt))then
      deallocate(this%Prt)
    end if
    if(allocated(this%Neib))then
      deallocate(this%Neib)
    end if
    if(allocated(this%EdgeNodeInd))then
      deallocate(this%EdgeNodeInd)
    end if
    if(allocated(this%SurfNodeInd))then
      deallocate(this%SurfNodeInd)
    end if
    if(allocated(this%SurfPC))then
      deallocate(this%SurfPC)
    end if
    if(allocated(this%SurfArea))then
      deallocate(this%SurfArea)
    end if
    if(allocated(this%SurfNorm))then
      deallocate(this%SurfNorm)
    end if
  end subroutine
  
  !-------------------------------------------
  ! update all components of grid information
  !-------------------------------------------
  subroutine updateGrid()
    Node(:)%Ind=[(i,i=1,nNode)]
    call updateAllNodeFacetInd()
    call updateAllNodeBlockInd()
    !$omp parallel do
    do i=1,nPoint
      Point(i)%Ind=i
      call Point(i)%updatePos()
    end do
    !$omp end parallel do
    !$omp parallel do
    do i=1,nLine
      Line(i)%Ind=i
      call Line(i)%updatePC()
      call Line(i)%updateLength()
    end do
    !$omp end parallel do
    !$omp parallel do
    do i=1,nFacet
      Facet(i)%Ind=i
      call Facet(i)%updatePC()
      call Facet(i)%updateArea()
      call Facet(i)%updateNorm()
      call Facet(i)%updateNeibBlock()
    end do
    !$omp end parallel do
    !$omp parallel do
    do i=1,nBlock
      Block(i)%Ind=i
      call Block(i)%updatePC()
      call Block(i)%updateVol()
      call Block(i)%updateNeib()
      call Block(i)%updateEdgeNodeInd()
      call Block(i)%updateSurfNodeInd()
      call Block(i)%updateSurfPC()
      call Block(i)%updateSurfArea()
      call Block(i)%updateSurfNorm()
    end do
    !$omp end parallel do
    nPrt=0
    if(nBlock>0)then
      nPrt=maxval([(maxval(Block(i)%Prt(:)),i=1,nBlock)])
    end if
    forall(i=1:DIMS)
      BoundBox(1,i)=minval([(Node(j)%Pos(i),j=1,nNode)])
      BoundBox(2,i)=maxval([(Node(j)%Pos(i),j=1,nNode)])
    end forall
  end subroutine
  
end module
