!----------------------------------------------------------------------------- best with 100 columns

!> grid operations
module moduleGridOperation
  use moduleGrid
  private
  
  ! constants
  integer,parameter::NPLFB=5 !< node, point, line, facet & block, 5 things need to be mapped
  integer,parameter,public::MAP_NODE=1 !< map index of node
  integer,parameter,public::MAP_POINT=2 !< map index of Point
  integer,parameter,public::MAP_LINE=3 !< map index of Line
  integer,parameter,public::MAP_FACET=4 !< map index of Facet
  integer,parameter,public::MAP_BLOCK=5 !< map index of Block
  
  ! public procedures
  public::splitGridPrt
  
contains
  
  !> split the original grid into resulting grids with respect to partition
  subroutine splitGridPrt(org,rst,r2o,overlap)
    use moduleBasicDataStruct
    use moduleSimpleSetLogic
    type(typeGrid),intent(in)::org !< the original grid
    type(typeGrid),intent(inout),allocatable::rst(:) !< the resulting grids
    type(typeHtr1DIArr),intent(inout),allocatable::r2o(:,:) !< mapping from rst to org
    logical,intent(in),optional::overlap !< whether the partitions are overlapping
    integer o2rNodeMap(org%nNode,org%nPrt)
    logical isOverlap
    
    if(present(overlap))then
      isOverlap=overlap
    else
      isOverlap=.false.
    end if
    ! initiate output data
    if(allocated(rst)) deallocate(rst)
    if(allocated(r2o)) deallocate(r2o)
    allocate(rst(org%nPrt))
    allocate(r2o(NPLFB,org%nPrt))
    do i=1,org%nPrt
      call rst(i)%init()
      rst(i)%nPrt=1
    end do
    o2rNodeMap(:,:)=0
    ! generate mapping
    do i=1,org%nPoint
      do j=1,size(org%Point(i)%Prt)
        if(org%Point(i)%Prt(j)/=0)then
          call pushArr(r2o(MAP_POINT,org%Point(i)%Prt(j))%dat,i)
          rst(org%Point(i)%Prt(j))%nPoint=rst(org%Point(i)%Prt(j))%nPoint+1
        end if
      end do
    end do
    do i=1,org%nLine
      do j=1,size(org%Line(i)%Prt)
        if(org%Line(i)%Prt(j)/=0)then
          call pushArr(r2o(MAP_LINE,org%Line(i)%Prt(j))%dat,i)
          rst(org%Line(i)%Prt(j))%nLine=rst(org%Line(i)%Prt(j))%nLine+1
        end if
      end do
    end do
    do i=1,org%nFacet
      do j=1,size(org%Facet(i)%Prt)
        if(org%Facet(i)%Prt(j)/=0)then
          call pushArr(r2o(MAP_FACET,org%Facet(i)%Prt(j))%dat,i)
          rst(org%Facet(i)%Prt(j))%nFacet=rst(org%Facet(i)%Prt(j))%nFacet+1
        end if
      end do
    end do
    do i=1,org%nBlock
      do j=1,size(org%Block(i)%Prt)
        if(isOverlap)then
          if(org%Block(i)%Prt(j)/=0)then
            call pushArr(r2o(MAP_BLOCK,abs(org%Block(i)%Prt(j)))%dat,i)
            rst(abs(org%Block(i)%Prt(j)))%nBlock=rst(abs(org%Block(i)%Prt(j)))%nBlock+1
            do k=1,org%Block(i)%nNode
              call applUnion(r2o(MAP_NODE,abs(org%Block(i)%Prt(j)))%dat,[org%Block(i)%iNode(k)])
              if(rst(abs(org%Block(i)%Prt(j)))%nNode/=&
              &  size(r2o(MAP_NODE,abs(org%Block(i)%Prt(j)))%dat))then
                rst(abs(org%Block(i)%Prt(j)))%nNode=size(r2o(MAP_NODE,abs(org%Block(i)%Prt(j)))%dat)
                o2rNodeMap(org%Block(i)%iNode(k),abs(org%Block(i)%Prt(j)))=&
                & rst(abs(org%Block(i)%Prt(j)))%nNode
              end if
            end do
          end if
        else
          if(org%Block(i)%Prt(j)>0)then
            call pushArr(r2o(MAP_BLOCK,org%Block(i)%Prt(j))%dat,i)
            rst(org%Block(i)%Prt(j))%nBlock=rst(org%Block(i)%Prt(j))%nBlock+1
            do k=1,org%Block(i)%nNode
              call applUnion(r2o(MAP_NODE,org%Block(i)%Prt(j))%dat,[org%Block(i)%iNode(k)])
              if(rst(org%Block(i)%Prt(j))%nNode/=size(r2o(MAP_NODE,org%Block(i)%Prt(j))%dat))then
                rst(org%Block(i)%Prt(j))%nNode=size(r2o(MAP_NODE,org%Block(i)%Prt(j))%dat)
                o2rNodeMap(org%Block(i)%iNode(k),org%Block(i)%Prt(j))=rst(org%Block(i)%Prt(j))%nNode
              end if
            end do
          end if
        end if
      end do
    end do
    do i=1,org%nPrt
      rst(i)%nNode=size(r2o(MAP_NODE,i)%dat)
    end do
    ! copy items
    do i=1,org%nPrt
      allocate(rst(i)%NodePos(DIMS,rst(i)%nNode))
      rst(i)%NodePos(:,:)=org%NodePos(:,r2o(MAP_NODE,i)%dat(:))
      allocate(rst(i)%Point(rst(i)%nPoint))
      rst(i)%Point(:)=org%Point(r2o(MAP_POINT,i)%dat(:))
      allocate(rst(i)%Line(rst(i)%nLine))
      rst(i)%Line(:)=org%Line(r2o(MAP_LINE,i)%dat(:))
      allocate(rst(i)%Facet(rst(i)%nFacet))
      rst(i)%Facet(:)=org%Facet(r2o(MAP_FACET,i)%dat(:))
      allocate(rst(i)%Block(rst(i)%nBlock))
      rst(i)%Block(:)=org%Block(r2o(MAP_BLOCK,i)%dat(:))
    end do
    ! correct node reference
    do i=1,org%nPrt
      forall(j=1:rst(i)%nPoint)
        rst(i)%Point(j)%iNode(:)=o2rNodeMap(rst(i)%Point(j)%iNode(:),i)
      end forall
      forall(j=1:rst(i)%nLine)
        rst(i)%Line(j)%iNode(:)=o2rNodeMap(rst(i)%Line(j)%iNode(:),i)
      end forall
      forall(j=1:rst(i)%nFacet)
        rst(i)%Facet(j)%iNode(:)=o2rNodeMap(rst(i)%Facet(j)%iNode(:),i)
      end forall
      forall(j=1:rst(i)%nBlock)
        rst(i)%Block(j)%iNode(:)=o2rNodeMap(rst(i)%Block(j)%iNode(:),i)
      end forall
    end do
  end subroutine
  
end module
