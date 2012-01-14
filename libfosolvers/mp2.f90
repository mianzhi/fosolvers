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
  
  ! procedures
  !public sendPrt
  !public recvPrt
  
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
  
contains
  
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
      call sendData(nPointPrt,p,p)
      call sendData(buffPoint,p,p)
      call sendData(nLinePrt,p,p)
      call sendData(buffLine,p,p)
      call sendData(nFacetPrt,p,p)
      call sendData(buffFacet,p,p)
      call sendData(nBlockPrt,p,p)
      call sendData(buffBlock,p,p)
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
      call recvData(nPoint,ROOT_PID,pidMPI)
      call recvData(Point,ROOT_PID,pidMPI,realloc=.true.)
      call recvData(nLine,ROOT_PID,pidMPI)
      call recvData(Line,ROOT_PID,pidMPI,realloc=.true.)
      call recvData(nFacet,ROOT_PID,pidMPI)
      call recvData(Facet,ROOT_PID,pidMPI,realloc=.true.)
      call recvData(nBlock,ROOT_PID,pidMPI)
      call recvData(Block,ROOT_PID,pidMPI,realloc=.true.)
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
  end subroutine
  
  !---------------------------------------------------------
  ! receive grid, conditions and vector data of a partition
  !---------------------------------------------------------
  subroutine recvPrtWithVect(vectV)
    use moduleMP1
    double precision,allocatable,intent(inout)::vectV(:,:)
    
    call recvPrtGridCondOnly()
    call recvData(vectV,ROOT_PID,pidMPI,realloc=.true.)
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
  end subroutine
  
end module
