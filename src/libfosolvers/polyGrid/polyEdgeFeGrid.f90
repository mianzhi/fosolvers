!----------------------------------------------------------------------------- best with 100 columns

!> edge-based finite element grid module
module modPolyEdgeFeGrid
  use modPolyFeGrid
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  integer,parameter,public::TRI_N_EDGE=3 !< 3 edges per TRI
  integer,parameter,public::TRI_NI_EDGE(2,TRI_N_EDGE)=reshape(& !< node index of each edge in TRI
  &                                                   [1,2,1,3,2,3],[2,TRI_N_EDGE])
  integer,parameter,public::TET_N_EDGE=6 !< 6 edges per TET
  integer,parameter,public::TET_NI_EDGE(2,TET_N_EDGE)=reshape(& !< node index of each edge in TET
  &                                                   [1,2,1,3,1,4,2,3,2,4,3,4],[2,TET_N_EDGE])
  
  !> edge-based finite volume grid type
  type,extends(polyFeGrid),public::polyEdgeFeGrid
    integer::nEdge !< number of edges
    integer,allocatable::iNEdge(:,:) !< node index of each edge
    integer,allocatable::iEdgeE(:,:) !< edge index of each element
    integer,allocatable::nEdgeE(:) !< number of edges in each element
    logical,allocatable::sameDir(:,:) !< if the element edge and grid edge has the same direction
  contains
    procedure,public::clear=>clearPolyEdgeFeGrid
    procedure,public::up=>upPolyEdgeFeGrid
    final::purgePolyEdgeFeGrid
  end type
  
contains
  
  !> clear this polyEdgeFeGrid
  elemental subroutine clearPolyEdgeFeGrid(this)
    class(polyEdgeFeGrid),intent(inout)::this !< this polyEdgeFeGrid
    
    call this%polyFeGrid%clear()
    if(allocated(this%iNEdge)) deallocate(this%iNEdge)
    if(allocated(this%iEdgeE)) deallocate(this%iEdgeE)
    if(allocated(this%nEdgeE)) deallocate(this%nEdgeE)
    if(allocated(this%sameDir)) deallocate(this%sameDir)
  end subroutine
  
  !> update this polyEdgeFeGrid
  subroutine upPolyEdgeFeGrid(this)
    use modPolyGrid
    class(polyEdgeFeGrid),intent(inout)::this !< this polyEdgeFeGrid
    integer,allocatable::iNEdge(:,:)
    integer,allocatable::iEdgeE(:,:)
    integer,allocatable::nEdgeE(:)
    logical,allocatable::sameDir(:,:)
    integer,parameter::MAX_EDGE_PER_E=20 !< maximum 20 edges per element
    integer::nEdge,existingEdge
    
    if(.not.this%isUp)then
      call this%polyFeGrid%up()
      allocate(iNEdge(2,this%nE*MAX_EDGE_PER_E),source=0)
      allocate(iEdgeE(MAX_EDGE_PER_E,this%nE),source=0)
      allocate(nEdgeE(this%nE),source=0)
      allocate(sameDir(MAX_EDGE_PER_E,this%nE),source=.true.)
      nEdge=0
      do i=1,this%nE
        select case(this%sE(i))
        case(TET)
          nEdgeE(i)=TET_N_EDGE
          do j=1,TET_N_EDGE
            m=this%iNE(TET_NI_EDGE(1,j),i)
            n=this%iNE(TET_NI_EDGE(2,j),i)
            existingEdge=0
            do k=1,nEdge ! check if edge already exists
              if(all(iNEdge(:,k)==[min(m,n),max(m,n)]))then
                existingEdge=k
                exit
              end if
            end do
            if(existingEdge>0)then
              iEdgeE(j,i)=existingEdge
              sameDir(j,i)=m<n
            else
              nEdge=nEdge+1
              iNEdge(:,nEdge)=[min(m,n),max(m,n)]
              iEdgeE(j,i)=nEdge
              sameDir(j,i)=m<n
            end if
          end do
        case(TRI)
          nEdgeE(i)=TRI_N_EDGE
          do j=1,TRI_N_EDGE
            m=this%iNE(TRI_NI_EDGE(1,j),i)
            n=this%iNE(TRI_NI_EDGE(2,j),i)
            existingEdge=0
            do k=1,nEdge ! check if edge already exists
              if(all(iNEdge(:,k)==[min(m,n),max(m,n)]))then
                existingEdge=k
                exit
              end if
            end do
            if(existingEdge>0)then
              iEdgeE(j,i)=existingEdge
              sameDir(j,i)=m<n
            else
              nEdge=nEdge+1
              iNEdge(:,nEdge)=[min(m,n),max(m,n)]
              iEdgeE(j,i)=nEdge
              sameDir(j,i)=m<n
            end if
          end do
        case default
        end select
      end do
      if(allocated(this%iNEdge)) deallocate(this%iNEdge)
      if(allocated(this%iEdgeE)) deallocate(this%iEdgeE)
      if(allocated(this%nEdgeE)) deallocate(this%nEdgeE)
      if(allocated(this%sameDir)) deallocate(this%sameDir)
      this%nEdge=nEdge
      allocate(this%iNEdge,source=iNEdge(:,1:nEdge))
      allocate(this%iEdgeE,source=iEdgeE(1:maxval(nEdgeE(:)),:))
      allocate(this%nEdgeE,source=nEdgeE(:))
      allocate(this%sameDir,source=sameDir(1:maxval(nEdgeE(:)),:))
      deallocate(iNEdge,iEdgeE,nEdgeE,sameDir)
    end if
  end subroutine
  
  !> destructor of polyEdgeFeGrid
  elemental subroutine purgePolyEdgeFeGrid(this)
    type(polyEdgeFeGrid),intent(inout)::this !< this polyEdgeFeGrid
    
    call this%clear()
  end subroutine
end module
