!----------------------------------------------------------------------------- best with 100 columns

!> polyhedron and polygon finite volume grid module
module modPolyFvGrid
  use modPolyGrid
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  
  !> polyhedron and polygon finite volume grid type
  type,extends(polyGrid),public::polyFvGrid
    integer::nP !< number of pairs of elements
    integer,allocatable::iEP(:,:) !< indices of elements of each pair
    integer,allocatable::neib(:,:) !< neighbor list
  contains
    procedure,public::clear=>clearPolyFvGrid
    procedure,public::up=>upPolyFvGrid
    !FIXME:final::purgePolyFvGrid
  end type
  
contains
  
  !> clear this polyFvGrid
  elemental subroutine clearPolyFvGrid(this)
    class(polyFvGrid),intent(inout)::this !< this polyFvGrid
    
    call this%polyGrid%clear()
    if(allocated(this%iEP)) deallocate(this%iEP)
    if(allocated(this%neib)) deallocate(this%neib)
  end subroutine
  
  !> update this polyFvGrid
  elemental subroutine upPolyFvGrid(this)
    class(polyFvGrid),intent(inout)::this !< this polyFvGrid
    
    if(.not.this%isUp)then
      call this%polyGrid%up()
      call getNeibPolyFvGrid(this)
    end if
  end subroutine
  
  !> get neighbor list and pairs
  elemental subroutine getNeibPolyFvGrid(grid)
    class(polyFvGrid),intent(inout)::grid !< the polyFvGrid
    integer,allocatable::iEN(:,:),nEN(:),nFace(:),iNF(:),iNF1(:),iEP(:,:)
    
    ! indices of elements containing each node
    allocate(nEN(grid%nN))
    nEN(:)=0
    do i=1,grid%nE
      do j=1,grid%nNE(i)
        nEN(grid%iNE(j,i))=nEN(grid%iNE(j,i))+1
      end do
    end do
    allocate(iEN(maxval(nEN),grid%nN))
    iEN(:,:)=0
    nEN(:)=0
    do i=1,grid%nE
      do j=1,grid%nNE(i)
        nEN(grid%iNE(j,i))=nEN(grid%iNE(j,i))+1
        iEN(nEN(grid%iNE(j,i)),grid%iNE(j,i))=i
      end do
    end do
    ! prepare storage for the face with most nodes
    allocate(nFace(grid%nE))
    nFace(:)=nF(grid%sE(:))
    allocate(iNF(maxval([(maxval([(nNF(grid%sE(i),j),j=1,nFace(1))]),i=1,grid%nE)])))
    allocate(iNF1(size(iNF)))
    ! find neighbor through each face
    if(allocated(grid%neib)) deallocate(grid%neib)
    if(allocated(grid%iEP)) deallocate(grid%iEP)
    allocate(grid%neib(maxval(nFace),grid%nC))
    allocate(iEP(2,sum(nFace)))
    grid%neib(:,:)=0
    grid%nP=0
    do i=1,grid%nC
      do j=1,nFace(i)
        if(grid%neib(j,i)==0)then
          ! find list of nodes of interface
          call getINF(grid%sE(i),j,iNF)
          forall(k=1:size(iNF),iNF(k)>0)
            iNF(k)=grid%iNE(iNF(k),i)
          end forall
          ! 1_st round find facets only
          do k=1,nEN(iNF(1))
            m=iEN(k,iNF(1)) ! m is an element with the 1_st node in iNF
            if(m/=i.and.m>grid%nC)then ! facets only
              do l=1,nFace(m)
                ! find list of nodes on face l of candidate element m
                call getINF(grid%sE(m),l,iNF1)
                forall(ii=1:size(iNF1),iNF1(ii)>0)
                  iNF1(ii)=grid%iNE(iNF1(ii),m)
                end forall
                if(all([(any(iNF1(:)==iNF(ii)),ii=1,size(iNF))]))then ! matching pair
                  grid%neib(j,i)=m
                  grid%nP=grid%nP+1
                  iEP(:,grid%nP)=[i,m]
                end if
              end do
            end if
          end do
          ! move on and omit neighbor cells if found neighbor facet
          if(grid%neib(j,i)/=0)then
            cycle
          end if
          ! 2_nd round find cells only
          do k=1,nEN(iNF(1))
            m=iEN(k,iNF(1)) ! m is an element with the 1_st node in iNF
            if(m/=i.and.m<=grid%nC)then ! cells only
              do l=1,nFace(m)
                ! find list of nodes on face l of candidate element m
                call getINF(grid%sE(m),l,iNF1)
                forall(ii=1:size(iNF1),iNF1(ii)>0)
                  iNF1(ii)=grid%iNE(iNF1(ii),m)
                end forall
                if(all([(any(iNF1(:)==iNF(ii)),ii=1,size(iNF))]))then ! matching pair
                  grid%neib(j,i)=m
                  grid%neib(l,m)=i
                  grid%nP=grid%nP+1
                  iEP(:,grid%nP)=[i,m]
                  exit
                end if
              end do
            end if
          end do
        end if
      end do
    end do
    allocate(grid%iEP(2,grid%nP),source=iEP(:,1:grid%nP))!FIXME:remove work-around
    deallocate(iEN)
    deallocate(nEN)
    deallocate(nFace)
    deallocate(iNF)
    deallocate(iNF1)
    deallocate(iEP)
  end subroutine
  
  !> destructor of polyGrid
  elemental subroutine purgePolyFvGrid(this)
    type(polyFvGrid),intent(inout)::this !< this polyFvGrid
    
    call this%clear()
  end subroutine
  
end module
