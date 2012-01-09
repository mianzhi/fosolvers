!----------------------------------------------------------------------------- best with 100 columns

!*****************************************
! multi-processing related data (level 2)
!*****************************************
module moduleMP2
  use moduleMiscDataStruct
  private
  
  ! data list to be communicated
  type(typePtrScalArray),allocatable,public,save::commNodeScal(:),commFacetScal(:),commBlockScal(:)
  type(typePtrVectArray),allocatable,public,save::commNodeVect(:),commFacetVect(:),commBlockVect(:)
  type(typePtrTensArray),allocatable,public,save::commNodeTens(:),commFacetTens(:),commBlockTens(:)
  
  ! index map
  integer,allocatable,public,save::mapNode(:),mapPoint(:),mapLine(:),mapFacet(:),mapBlock(:)
  
  !--------------------------------
  ! generic add communication list
  !--------------------------------
  interface addComm
    module procedure::addCommScal
    module procedure::addCommVect
    module procedure::addCommTens
  end interface
  public addComm
  
  ! other procedures
  public distriPrt
  public recvPrt
  
contains
  
  !---------------------------------------
  ! add scaler data to communication list
  !---------------------------------------
  subroutine addCommScal(v,binding)
    use moduleGrid
    double precision,target,intent(in)::v(:)
    integer,intent(in),optional::binding
    logical bindnode,bindfacet,bindblock
    type(typePtrScalArray),allocatable::temp(:)
    
    bindnode=.false.
    bindfacet=.false.
    bindblock=.false.
    if(present(binding))then
      select case(binding)
        case(BIND_NODE)
          bindnode=.true.
        case(BIND_FACET)
          bindfacet=.true.
        case(BIND_BLOCK)
          bindblock=.true.
        case default
      end select
    else
      bindnode=(size(v,1)==nNode)
      bindfacet=(size(v,1)==nFacet)
      bindblock=(size(v,1)==nBlock)
    end if
    
    if(bindnode)then
      if(allocated(commNodeScal))then
        allocate(temp(size(commNodeScal)+1))
        temp(1:size(commNodeScal))=commNodeScal(:)
        temp(size(temp))%ptr=>v
        call move_alloc(temp,commNodeScal)
      else
        allocate(commNodeScal(1))
        commNodeScal(1)%ptr=>v
      end if
    end if
    if(bindfacet)then
      if(allocated(commFacetScal))then
        allocate(temp(size(commFacetScal)+1))
        temp(1:size(commFacetScal))=commFacetScal(:)
        temp(size(temp))%ptr=>v
        call move_alloc(temp,commFacetScal)
      else
        allocate(commFacetScal(1))
        commFacetScal(1)%ptr=>v
      end if
    end if
    if(bindblock)then
      if(allocated(commBlockScal))then
        allocate(temp(size(commBlockScal)+1))
        temp(1:size(commBlockScal))=commBlockScal(:)
        temp(size(temp))%ptr=>v
        call move_alloc(temp,commBlockScal)
      else
        allocate(commBlockScal(1))
        commBlockScal(1)%ptr=>v
      end if
    end if
  end subroutine
  
  !---------------------------------------
  ! add vector data to communication list
  !---------------------------------------
  subroutine addCommVect(v,binding)
    use moduleGrid
    double precision,target,intent(in)::v(:,:)
    integer,intent(in),optional::binding
    logical bindnode,bindfacet,bindblock
    type(typePtrVectArray),allocatable::temp(:)
    
    bindnode=.false.
    bindfacet=.false.
    bindblock=.false.
    if(present(binding))then
      select case(binding)
        case(BIND_NODE)
          bindnode=.true.
        case(BIND_FACET)
          bindfacet=.true.
        case(BIND_BLOCK)
          bindblock=.true.
        case default
      end select
    else
      bindnode=(size(v,1)==nNode)
      bindfacet=(size(v,1)==nFacet)
      bindblock=(size(v,1)==nBlock)
    end if
    
    if(bindnode)then
      if(allocated(commNodeVect))then
        allocate(temp(size(commNodeVect)+1))
        temp(1:size(commNodeVect))=commNodeVect(:)
        temp(size(temp))%ptr=>v
        call move_alloc(temp,commNodeVect)
      else
        allocate(commNodeVect(1))
        commNodeVect(1)%ptr=>v
      end if
    end if
    if(bindfacet)then
      if(allocated(commFacetVect))then
        allocate(temp(size(commFacetVect)+1))
        temp(1:size(commFacetVect))=commFacetVect(:)
        temp(size(temp))%ptr=>v
        call move_alloc(temp,commFacetVect)
      else
        allocate(commFacetVect(1))
        commFacetVect(1)%ptr=>v
      end if
    end if
    if(bindblock)then
      if(allocated(commBlockVect))then
        allocate(temp(size(commBlockVect)+1))
        temp(1:size(commBlockVect))=commBlockVect(:)
        temp(size(temp))%ptr=>v
        call move_alloc(temp,commBlockVect)
      else
        allocate(commBlockVect(1))
        commBlockVect(1)%ptr=>v
      end if
    end if
  end subroutine
  
  !---------------------------------------
  ! add tensor data to communication list
  !---------------------------------------
  subroutine addCommTens(v,binding)
    use moduleGrid
    double precision,target,intent(in)::v(:,:,:)
    integer,intent(in),optional::binding
    logical bindnode,bindfacet,bindblock
    type(typePtrTensArray),allocatable::temp(:)
    
    bindnode=.false.
    bindfacet=.false.
    bindblock=.false.
    if(present(binding))then
      select case(binding)
        case(BIND_NODE)
          bindnode=.true.
        case(BIND_FACET)
          bindfacet=.true.
        case(BIND_BLOCK)
          bindblock=.true.
        case default
      end select
    else
      bindnode=(size(v,1)==nNode)
      bindfacet=(size(v,1)==nFacet)
      bindblock=(size(v,1)==nBlock)
    end if
    
    if(bindnode)then
      if(allocated(commNodeTens))then
        allocate(temp(size(commNodeTens)+1))
        temp(1:size(commNodeTens))=commNodeTens(:)
        temp(size(temp))%ptr=>v
        call move_alloc(temp,commNodeTens)
      else
        allocate(commNodeTens(1))
        commNodeTens(1)%ptr=>v
      end if
    end if
    if(bindfacet)then
      if(allocated(commFacetTens))then
        allocate(temp(size(commFacetTens)+1))
        temp(1:size(commFacetTens))=commFacetTens(:)
        temp(size(temp))%ptr=>v
        call move_alloc(temp,commFacetTens)
      else
        allocate(commFacetTens(1))
        commFacetTens(1)%ptr=>v
      end if
    end if
    if(bindblock)then
      if(allocated(commBlockTens))then
        allocate(temp(size(commBlockTens)+1))
        temp(1:size(commBlockTens))=commBlockTens(:)
        temp(size(temp))%ptr=>v
        call move_alloc(temp,commBlockTens)
      else
        allocate(commBlockTens(1))
        commBlockTens(1)%ptr=>v
      end if
    end if
  end subroutine
  
  !-------------------------------------
  ! distribute partition k to process p
  !-------------------------------------
  subroutine distriPrt(k,p)
    use moduleGrid
    use moduleCond
    use moduleTime
    use moduleMP1
    use moduleMiscDataStruct
    integer,intent(in)::k,p
    integer nNodePrt,nPointPrt,nLinePrt,nFacetPrt,nBlockPrt
    type(typeNode),allocatable::buffNode(:)
    type(typePoint),allocatable::buffPoint(:)
    type(typeLine),allocatable::buffLine(:)
    type(typeFacet),allocatable::buffFacet(:)
    type(typeBlock),allocatable::buffBlock(:)
    type(typeDataSet),allocatable::buffCondNode(:)
    type(typeDataSet),allocatable::buffCondFacet(:)
    type(typeDataSet),allocatable::buffCondBlock(:)
    integer,allocatable::buffMapNode(:),buffMapPoint(:),buffMapLine(:),buffMapFacet(:),&
    &                    buffMapBlock(:)
    
    ! copy blocks first (only blocks carry partition information)
    if(.not.allocated(mapBlock))then
      allocate(mapBlock(nBlock))
    end if
    mapBlock(:)=0
    nBlockPrt=count([(any(Block(i)%Prt(:)==k),i=1,nBlock)])
    allocate(buffBlock(nBlockPrt))
    allocate(buffMapBlock(nBlockPrt))
    nBlockPrt=0
    do i=1,nBlock
      if(any(Block(i)%Prt(:)==k))then
        nBlockPrt=nBlockPrt+1
        mapBlock(i)=nBlockPrt
        buffMapBlock(nBlockPrt)=i
      end if
    end do
    buffBlock(:)=Block(buffMapBlock(:))
    ! copy nodes then (decides if other elements belongs to partition k)
    if(.not.allocated(mapNode))then
      allocate(mapNode(nNode))
    end if
    mapNode(:)=0
    nNodePrt=0
    do i=1,nBlockPrt
      do j=1,buffBlock(i)%NodeNum
        if(mapNode(buffBlock(i)%NodeInd(j))==0)then
          nNodePrt=nNodePrt+1
          mapNode(buffBlock(i)%NodeInd(j))=nNodePrt
        end if
      end do
    end do
    allocate(buffNode(nNodePrt))
    allocate(buffMapNode(nNodePrt))
    forall(i=1:nNode,mapNode(i)>0)
      buffMapNode(mapNode(i))=i
    end forall
    buffNode(:)=Node(buffMapNode(:))
    ! copy points
    if(.not.allocated(mapPoint))then
      allocate(mapPoint(nPoint))
    end if
    mapPoint(:)=0
    nPointPrt=0
    do i=1,nPoint
      if(any(buffMapNode(:)==Point(i)%NodeInd))then
        nPointPrt=nPointPrt+1
        mapPoint(i)=nPointPrt
      end if
    end do
    allocate(buffPoint(nPointPrt))
    allocate(buffMapPoint(nPointPrt))
    forall(i=1:nPoint,mapPoint(i)>0)
      buffMapPoint(mapPoint(i))=i
    end forall
    buffPoint(:)=Point(buffMapPoint(:))
    ! copy lines
    if(.not.allocated(mapLine))then
      allocate(mapLine(nLine))
    end if
    mapLine(:)=0
    nLinePrt=0
    do i=1,nLine
      if(all([(any(buffMapNode(:)==Line(i)%NodeInd(j)),j=1,LINE_NODE_NUM)]))then
        nLinePrt=nLinePrt+1
        mapLine(i)=nLinePrt
      end if
    end do
    allocate(buffLine(nLinePrt))
    allocate(buffMapLine(nLinePrt))
    forall(i=1:nLine,mapLine(i)>0)
      buffMapLine(mapLine(i))=i
    end forall
    buffLine(:)=Line(buffMapLine(:))
    ! copy facets
    if(.not.allocated(mapFacet))then
      allocate(mapFacet(nFacet))
    end if
    mapFacet(:)=0
    nFacetPrt=0
    do i=1,nFacet
      if(all([(any(buffMapNode(:)==Facet(i)%NodeInd(j)),j=1,Facet(i)%NodeNum)]))then
        nFacetPrt=nFacetPrt+1
        mapFacet(i)=nFacetPrt
      end if
    end do
    allocate(buffFacet(nFacetPrt))
    allocate(buffMapFacet(nFacetPrt))
    forall(i=1:nFacet,mapFacet(i)>0)
      buffMapFacet(mapFacet(i))=i
    end forall
    buffFacet(:)=Facet(buffMapFacet(:))
    
    ! copy node conditions
    allocate(buffCondNode(nNodePrt))
    buffCondNode(:)=condNode(buffMapNode(:))
    ! copy facet conditions
    allocate(buffCondFacet(nFacetPrt))
    buffCondFacet(:)=condFacet(buffMapFacet(:))
    ! copy block conditions
    allocate(buffCondBlock(nBlockPrt))
    buffCondBlock(:)=condBlock(buffMapBlock(:))
    
    ! correct indices
    
    ! send date
    call sendData(nNodePrt,p,p)
    call sendData(buffNode,p,p)
    
    ! clean ups
    deallocate(buffNode,buffMapNode)
    deallocate(buffPoint,buffMapPoint)
    deallocate(buffLine,buffMapLine)
    deallocate(buffFacet,buffMapFacet)
    deallocate(buffBlock,buffMapBlock)
    deallocate(buffCondNode,buffCondFacet,buffCondBlock)
  end subroutine
  
  !-------------------
  ! receive partition
  !-------------------
  subroutine recvPrt()
    use moduleGrid
    use moduleCond
    use moduleTime
    use moduleMP1
    
    ! receive data
    call recvData(nNode,ROOT_PID,pidMPI)
    call recvData(Node,ROOT_PID,pidMPI,realloc=.true.)
  end subroutine
  
end module
