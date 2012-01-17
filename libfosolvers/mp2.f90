!----------------------------------------------------------------------------- best with 100 columns

! buffer of grid maping
module moduleMP2buffMap
  integer,allocatable,save::buffMapNode(:),buffMapPoint(:),buffMapLine(:),buffMapFacet(:),&
  &                         buffMapBlock(:)
end module

!*****************************************
! multi-processing related data (level 2)
!*****************************************
module moduleMP2
  use moduleMiscDataStruct
  private
  
  ! index map
  integer,allocatable,public,save::mapNode(:),mapPoint(:),mapLine(:),mapFacet(:),mapBlock(:)
  
  ! data list to be communicated
  type(typePtrScalArray),allocatable,public,save::commNodeScal(:),commFacetScal(:),commBlockScal(:)
  type(typePtrVectArray),allocatable,public,save::commNodeVect(:),commFacetVect(:),commBlockVect(:)
  type(typePtrTensArray),allocatable,public,save::commNodeTens(:),commFacetTens(:),commBlockTens(:)
  
  !------------------------
  ! generic send partition
  !------------------------
  interface sendPrt
    module procedure::sendPrtGridCondOnly
    module procedure::sendPrtWithScal
    module procedure::sendPrtWithVect
    module procedure::sendPrtWithTens
  end interface
  public sendPrt
  
  !---------------------------
  ! generic receive partition
  !---------------------------
  interface recvPrt
    module procedure::recvPrtGridCondOnly
    module procedure::recvPrtWithScal
    module procedure::recvPrtWithVect
    module procedure::recvPrtWithTens
  end interface
  public recvPrt
  
  ! other procedures
  public retnPrtData
  public gathPrtData
  
contains
  
  !---------------------------------------
  ! add scaler data to communication list
  !---------------------------------------
  subroutine addCommScal(v,binding)
    use moduleGrid
    double precision,target,intent(in)::v(:)
    integer,intent(in)::binding
    type(typePtrScalArray),allocatable::temp(:)
    
    select case(binding)
      case(BIND_NODE)
        if(allocated(commNodeScal))then
          if(.not.any([(associated(commNodeScal(i)%ptr,v),i=1,size(commNodeScal))]))then
            allocate(temp(size(commNodeScal)+1))
            temp(1:size(commNodeScal))=commNodeScal(:)
            temp(size(temp))%ptr=>v
            call move_alloc(temp,commNodeScal)
          end if
        else
          allocate(commNodeScal(1))
          commNodeScal(1)%ptr=>v
        end if
      case(BIND_FACET)
        if(allocated(commFacetScal))then
          if(.not.any([(associated(commFacetScal(i)%ptr,v),i=1,size(commFacetScal))]))then
            allocate(temp(size(commFacetScal)+1))
            temp(1:size(commFacetScal))=commFacetScal(:)
            temp(size(temp))%ptr=>v
            call move_alloc(temp,commFacetScal)
          end if
        else
          allocate(commFacetScal(1))
          commFacetScal(1)%ptr=>v
        end if
      case(BIND_BLOCK)
        if(allocated(commBlockScal))then
          if(.not.any([(associated(commBlockScal(i)%ptr,v),i=1,size(commBlockScal))]))then
            allocate(temp(size(commBlockScal)+1))
            temp(1:size(commBlockScal))=commBlockScal(:)
            temp(size(temp))%ptr=>v
            call move_alloc(temp,commBlockScal)
          end if
        else
          allocate(commBlockScal(1))
          commBlockScal(1)%ptr=>v
        end if
      case default
    end select
  end subroutine
  
  !---------------------------------------
  ! add vector data to communication list
  !---------------------------------------
  subroutine addCommVect(v,binding)
    use moduleGrid
    double precision,target,intent(in)::v(:,:)
    integer,intent(in)::binding
    type(typePtrVectArray),allocatable::temp(:)
    
    select case(binding)
      case(BIND_NODE)
        if(allocated(commNodeVect))then
          if(.not.any([(associated(commNodeVect(i)%ptr,v),i=1,size(commNodeVect))]))then
            allocate(temp(size(commNodeVect)+1))
            temp(1:size(commNodeVect))=commNodeVect(:)
            temp(size(temp))%ptr=>v
            call move_alloc(temp,commNodeVect)
          end if
        else
          allocate(commNodeVect(1))
          commNodeVect(1)%ptr=>v
        end if
      case(BIND_FACET)
        if(allocated(commFacetVect))then
          if(.not.any([(associated(commFacetVect(i)%ptr,v),i=1,size(commFacetVect))]))then
            allocate(temp(size(commFacetVect)+1))
            temp(1:size(commFacetVect))=commFacetVect(:)
            temp(size(temp))%ptr=>v
            call move_alloc(temp,commFacetVect)
          end if
        else
          allocate(commFacetVect(1))
          commFacetVect(1)%ptr=>v
        end if
      case(BIND_BLOCK)
        if(allocated(commBlockVect))then
          if(.not.any([(associated(commBlockVect(i)%ptr,v),i=1,size(commBlockVect))]))then
            allocate(temp(size(commBlockVect)+1))
            temp(1:size(commBlockVect))=commBlockVect(:)
            temp(size(temp))%ptr=>v
            call move_alloc(temp,commBlockVect)
          end if
        else
          allocate(commBlockVect(1))
          commBlockVect(1)%ptr=>v
        end if
      case default
    end select
  end subroutine
  
  !---------------------------------------
  ! add tensor data to communication list
  !---------------------------------------
  subroutine addCommTens(v,binding)
    use moduleGrid
    double precision,target,intent(in)::v(:,:,:)
    integer,intent(in)::binding
    type(typePtrTensArray),allocatable::temp(:)
    
    select case(binding)
      case(BIND_NODE)
        if(allocated(commNodeTens))then
          if(.not.any([(associated(commNodeTens(i)%ptr,v),i=1,size(commNodeTens))]))then
            allocate(temp(size(commNodeTens)+1))
            temp(1:size(commNodeTens))=commNodeTens(:)
            temp(size(temp))%ptr=>v
            call move_alloc(temp,commNodeTens)
          end if
        else
          allocate(commNodeTens(1))
          commNodeTens(1)%ptr=>v
        end if
      case(BIND_FACET)
        if(allocated(commFacetTens))then
          if(.not.any([(associated(commFacetTens(i)%ptr,v),i=1,size(commFacetTens))]))then
            allocate(temp(size(commFacetTens)+1))
            temp(1:size(commFacetTens))=commFacetTens(:)
            temp(size(temp))%ptr=>v
            call move_alloc(temp,commFacetTens)
          end if
        else
          allocate(commFacetTens(1))
          commFacetTens(1)%ptr=>v
        end if
      case(BIND_BLOCK)
        if(allocated(commBlockTens))then
          if(.not.any([(associated(commBlockTens(i)%ptr,v),i=1,size(commBlockTens))]))then
            allocate(temp(size(commBlockTens)+1))
            temp(1:size(commBlockTens))=commBlockTens(:)
            temp(size(temp))%ptr=>v
            call move_alloc(temp,commBlockTens)
          end if
        else
          allocate(commBlockTens(1))
          commBlockTens(1)%ptr=>v
        end if
      case default
    end select
  end subroutine
  
  !------------------------------------------------------
  ! send grid and conditions of partition k to process p
  !------------------------------------------------------
  subroutine sendPrtGridCondOnly(k,p)
    use moduleGrid
    use moduleCond
    use moduleTime
    use moduleMP1
    use moduleMiscDataStruct
    use moduleMP2buffMap
    use mpi
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
    
    integer,save::lastk=-1,lastp=-1
    
    if(k/=lastk.and.p/=lastp)then ! need to send the grid and conditions
      call MPI_send(.true.,1,MPI_logical,p,p,MPI_comm_world,errMPI)
      ! copy blocks first (only blocks carry partition information)
      if(.not.allocated(mapBlock))then
        allocate(mapBlock(nBlock))
      end if
      mapBlock(:)=0
      nBlockPrt=count([(any(Block(i)%Prt(:)==k),i=1,nBlock)])
      allocate(buffBlock(nBlockPrt))
      if(allocated(buffMapBlock))then
        deallocate(buffMapBlock)
      end if
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
      if(allocated(buffMapNode))then
        deallocate(buffMapNode)
      end if
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
      if(allocated(buffMapPoint))then
        deallocate(buffMapPoint)
      end if
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
      if(allocated(buffMapLine))then
        deallocate(buffMapLine)
      end if
      allocate(buffMapLine(nLinePrt))
      forall(i=1:nLine,mapLine(i)>0)
        buffMapLine(mapLine(i))=i
      end forall
      buffLine(:)=Line(buffMapLine(:))
      ! copy facets (depend on blocks instead of nodes)
      if(.not.allocated(mapFacet))then
        allocate(mapFacet(nFacet))
      end if
      mapFacet(:)=0
      nFacetPrt=0
      do i=1,nFacet
        do j=1,FACET_NEIB_BLOCK_NUM
          if(Facet(i)%NeibBlock(j)>0)then
            if(mapBlock(Facet(i)%NeibBlock(j))>0)then
              nFacetPrt=nFacetPrt+1
              mapFacet(i)=nFacetPrt
              exit
            end if
          end if
        end do
      end do
      allocate(buffFacet(nFacetPrt))
      if(allocated(buffMapFacet))then
        deallocate(buffMapFacet)
      end if
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
      do i=1,nNodePrt
        buffNode(i)%Ind=i
        if(allocated(buffNode(i)%FacetInd))then
          buffNode(i)%FacetInd(:)=mapFacet(buffNode(i)%FacetInd)
        end if
        if(allocated(buffNode(i)%BlockInd))then
          buffNode(i)%BlockInd(:)=mapBlock(buffNode(i)%BlockInd)
        end if
      end do
      forall(i=1:nPointPrt)
        buffPoint(i)%Ind=i
        buffPoint(i)%NodeInd=mapNode(buffPoint(i)%NodeInd)
      end forall
      forall(i=1:nLinePrt)
        buffLine(i)%Ind=i
        buffLine(i)%NodeInd(:)=mapNode(buffLine(i)%NodeInd(:))
      end forall
      forall(i=1:nFacetPrt)
        buffFacet(i)%Ind=i
        buffFacet(i)%NodeInd(:)=mapNode(buffFacet(i)%NodeInd(:))
      end forall
      do i=1,nFacetPrt
        forall(j=1:FACET_NEIB_BLOCK_NUM,buffFacet(i)%NeibBlock(j)/=0)
          buffFacet(i)%NeibBlock(j)=mapBlock(buffFacet(i)%NeibBlock(j))
        end forall
      end do
      forall(i=1:nBlockPrt)
        buffBlock(i)%Ind=i
        buffBlock(i)%NodeInd(:)=mapNode(buffBlock(i)%NodeInd(:))
      end forall
      do i=1,nBlockPrt
        if(allocated(buffBlock(i)%Neib))then
          do j=1,buffBlock(i)%SurfNum
            if(buffBlock(i)%Neib(j)>0)then
              buffBlock(i)%Neib(j)=mapBlock(buffBlock(i)%Neib(j))
            else if(buffBlock(i)%Neib(j)<0)then
              buffBlock(i)%Neib(j)=-mapFacet(-buffBlock(i)%Neib(j))
            end if
          end do
        end if
      end do
      
      ! send grid and conditions
      call sendData(nNodePrt,p,p)
      call sendData(buffNode,p,p)
      call sendData(buffMapNode,p,p)
      call sendData(nPointPrt,p,p)
      call sendData(buffPoint,p,p)
      call sendData(buffMapPoint,p,p)
      call sendData(nLinePrt,p,p)
      call sendData(buffLine,p,p)
      call sendData(buffMapLine,p,p)
      call sendData(nFacetPrt,p,p)
      call sendData(buffFacet,p,p)
      call sendData(buffMapFacet,p,p)
      call sendData(nBlockPrt,p,p)
      call sendData(buffBlock,p,p)
      call sendData(buffMapBlock,p,p)
      call sendData(buffCondNode,p,p)
      call sendData(buffCondFacet,p,p)
      call sendData(buffCondBlock,p,p)
      
      ! clean ups
      deallocate(buffNode)
      deallocate(buffPoint)
      deallocate(buffLine)
      deallocate(buffFacet)
      deallocate(buffBlock)
      deallocate(buffCondNode,buffCondFacet,buffCondBlock)
    else
      call MPI_send(.false.,1,MPI_logical,p,p,MPI_comm_world,errMPI)
    end if
    
    ! record k and p
    lastk=k
    lastp=p
  end subroutine
  
  !-------------------------------------------------------------------------
  ! send grid, conditions and other scalar data of partition k to process p
  !-------------------------------------------------------------------------
  subroutine sendPrtWithScal(k,p,scalV,binding)
    use moduleMP1
    use moduleMP2buffMap
    use moduleGrid
    integer,intent(in)::k,p
    double precision,intent(in)::scalV(:)
    integer,intent(in)::binding
    
    call sendPrtGridCondOnly(k,p)
    select case(binding)
      case(BIND_NODE)
        call sendData(scalV(buffMapNode(:)),p,p)
      case(BIND_FACET)
        call sendData(scalV(buffMapFacet(:)),p,p)
      case(BIND_BLOCK)
        call sendData(scalV(buffMapBlock(:)),p,p)
      case default
    end select
    call sendData(binding,p,p)
    call addCommScal(scalV,binding=binding)
  end subroutine
  
  !-------------------------------------------------------------------------
  ! send grid, conditions and other vector data of partition k to process p
  !-------------------------------------------------------------------------
  subroutine sendPrtWithVect(k,p,vectV,binding)
    use moduleMP1
    use moduleMP2buffMap
    use moduleGrid
    integer,intent(in)::k,p
    double precision,intent(in)::vectV(:,:)
    integer,intent(in)::binding
    
    call sendPrtGridCondOnly(k,p)
    select case(binding)
      case(BIND_NODE)
        call sendData(vectV(buffMapNode(:),:),p,p)
      case(BIND_FACET)
        call sendData(vectV(buffMapFacet(:),:),p,p)
      case(BIND_BLOCK)
        call sendData(vectV(buffMapBlock(:),:),p,p)
      case default
    end select
    call sendData(binding,p,p)
    call addCommVect(vectV,binding=binding)
  end subroutine
  
  !-------------------------------------------------------------------------
  ! send grid, conditions and other tensor data of partition k to process p
  !-------------------------------------------------------------------------
  subroutine sendPrtWithTens(k,p,tensV,binding)
    use moduleMP1
    use moduleMP2buffMap
    use moduleGrid
    integer,intent(in)::k,p
    double precision,intent(in)::tensV(:,:,:)
    integer,intent(in)::binding
    
    call sendPrtGridCondOnly(k,p)
    select case(binding)
      case(BIND_NODE)
        call sendData(size(buffMapNode),p,p)
        do i=1,DIMS
          call sendData(tensV(buffMapNode(:),:,i),p,p)
        end do
      case(BIND_FACET)
        call sendData(size(buffMapFacet),p,p)
        do i=1,DIMS
          call sendData(tensV(buffMapFacet(:),:,i),p,p)
        end do
      case(BIND_BLOCK)
        call sendData(size(buffMapBlock),p,p)
        do i=1,DIMS
          call sendData(tensV(buffMapBlock(:),:,i),p,p)
        end do
      case default
    end select
    call sendData(binding,p,p)
    call addCommTens(tensV,binding=binding)
  end subroutine
  
  !--------------------------------------------
  ! receive grid and conditions of a partition
  !--------------------------------------------
  subroutine recvPrtGridCondOnly()
    use moduleGrid
    use moduleCond
    use moduleTime
    use moduleMP1
    use mpi
    logical::needRecv
    
    call MPI_recv(needRecv,1,MPI_logical,ROOT_PID,pidMPI,MPI_comm_world,statMPI,errMPI)
    if(needRecv)then ! need to receive grid and conditions
      call recvData(nNode,ROOT_PID,pidMPI)
      call recvData(Node,ROOT_PID,pidMPI,realloc=.true.)
      call recvData(mapNode,ROOT_PID,pidMPI,realloc=.true.)
      call recvData(nPoint,ROOT_PID,pidMPI)
      call recvData(Point,ROOT_PID,pidMPI,realloc=.true.)
      call recvData(mapPoint,ROOT_PID,pidMPI,realloc=.true.)
      call recvData(nLine,ROOT_PID,pidMPI)
      call recvData(Line,ROOT_PID,pidMPI,realloc=.true.)
      call recvData(mapLine,ROOT_PID,pidMPI,realloc=.true.)
      call recvData(nFacet,ROOT_PID,pidMPI)
      call recvData(Facet,ROOT_PID,pidMPI,realloc=.true.)
      call recvData(mapFacet,ROOT_PID,pidMPI,realloc=.true.)
      call recvData(nBlock,ROOT_PID,pidMPI)
      call recvData(Block,ROOT_PID,pidMPI,realloc=.true.)
      call recvData(mapBlock,ROOT_PID,pidMPI,realloc=.true.)
      call recvData(condNode,ROOT_PID,pidMPI,realloc=.true.)
      call recvData(condFacet,ROOT_PID,pidMPI,realloc=.true.)
      call recvData(condBlock,ROOT_PID,pidMPI,realloc=.true.)
    end if
  end subroutine
  
  !---------------------------------------------------------
  ! receive grid, conditions and scalar data of a partition
  !---------------------------------------------------------
  subroutine recvPrtWithScal(scalV)
    use moduleMP1
    double precision,allocatable,intent(inout)::scalV(:)
    
    call recvPrtGridCondOnly()
    call recvData(scalV,ROOT_PID,pidMPI,realloc=.true.)
    call recvData(k,ROOT_PID,pidMPI)
    call addCommScal(scalV,binding=k)
  end subroutine
  
  !---------------------------------------------------------
  ! receive grid, conditions and vector data of a partition
  !---------------------------------------------------------
  subroutine recvPrtWithVect(vectV)
    use moduleMP1
    double precision,allocatable,intent(inout)::vectV(:,:)
    
    call recvPrtGridCondOnly()
    call recvData(vectV,ROOT_PID,pidMPI,realloc=.true.)
    call recvData(k,ROOT_PID,pidMPI)
    call addCommVect(vectV,binding=k)
  end subroutine
  
  !---------------------------------------------------------
  ! receive grid, conditions and tensor data of a partition
  !---------------------------------------------------------
  subroutine recvPrtWithTens(tensV)
    use moduleMP1
    use moduleGrid
    double precision,allocatable,intent(inout)::tensV(:,:,:)
    
    call recvPrtGridCondOnly()
    if(allocated(tensV))then
      deallocate(tensV)
    end if
    call recvData(n,ROOT_PID,pidMPI)
    allocate(tensV(n,DIMS,DIMS))
    do i=1,DIMS
      call recvData(tensV(:,:,i),ROOT_PID,pidMPI)
    end do
    call recvData(k,ROOT_PID,pidMPI)
    call addCommTens(tensV,binding=k)
  end subroutine
  
  !-------------------------------------
  ! return data bind with the partition
  !-------------------------------------
  subroutine retnPrtData()
    use moduleMP1
    use moduleMP2buffMap
    use moduleGrid
    use mpi
    
    call sendData(mapNode,ROOT_PID,pidMPI)
    call sendData(mapFacet,ROOT_PID,pidMPI)
    call sendData(mapBlock,ROOT_PID,pidMPI)
    
    if(allocated(commNodeScal))then
      call MPI_send(.true.,1,MPI_logical,ROOT_PID,pidMPI,MPI_comm_world,errMPI)
      call sendData(size(commNodeScal),ROOT_PID,pidMPI)
      do i=1,size(commNodeScal)
        call sendData(commNodeScal(i)%ptr,ROOT_PID,pidMPI)
      end do
    else
      call MPI_send(.false.,1,MPI_logical,ROOT_PID,pidMPI,MPI_comm_world,errMPI)
    end if
  end subroutine
  
  !----------------------------------------------------------------------------
  ! try to gather data from any process, received partition prt from process p
  !----------------------------------------------------------------------------
  subroutine gathPrtData(prt,p)
    use moduleMP1
    use moduleGrid
    use mpi
    
    integer,intent(out)::prt,p
    integer,allocatable::tempMapNode(:),tempMapFacet(:),tempMapBlock(:)
    double precision,allocatable::buffVect(:),buffMat(:,:)
    logical isAllocated
    
    call recvData(tempMapNode,MPI_any_source,MPI_any_tag,realloc=.true.)
    p=statMPI(MPI_source)
    call recvData(tempMapFacet,p,p,realloc=.true.)
    call recvData(tempMapBlock,p,p,realloc=.true.)
    prt=maxval(Block(tempMapBlock(1))%Prt(:))
    
    call MPI_recv(isAllocated,1,MPI_logical,p,p,MPI_comm_world,statMPI,errMPI)
    if(isAllocated)then
      call recvData(n,p,p)
      do i=1,n
        call recvData(buffVect,p,p,realloc=.true.)
        commNodeScal(i)%ptr(tempMapNode(:))=buffVect(:)
      end do
    end if
    
    deallocate(tempMapNode,tempMapFacet,tempMapBlock)
    if(allocated(buffVect))then
      deallocate(buffVect)
    end if
    if(allocated(buffMat))then
      deallocate(buffMat)
    end if
  end subroutine
  
end module
