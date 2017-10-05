!----------------------------------------------------------------------------- best with 100 columns

!> polyhedron and polygon finite volume grid module
module modPolyFvGrid
  use modPolyGrid
  private
  
  ! constants
  integer,parameter::DIMS=3 !< dimensions
  integer,parameter::MIN_N_NEAR=8 !< minimum nearby elements for each cell
  integer,parameter::MAX_N_NEAR=25 !< maximum nearby elements for each cell
  
  !> polyhedron and polygon finite volume grid type
  type,extends(polyGrid),public::polyFvGrid
    integer::nP !< number of pairs of elements
    integer,allocatable::iEP(:,:) !< indices of elements of each pair
    integer,allocatable::neib(:,:) !< list of neighbor (adjacent via pair) elements to cells
    integer,allocatable::near(:,:) !< list of nearby elements to cells
    double precision,allocatable::aP(:) !< area of pair interfaces
    double precision,allocatable::normP(:,:) !< normal vector of pair interfaces
    double precision,allocatable::pP(:,:) !< center position of pair interfaces
  contains
    procedure,public::clear=>clearPolyFvGrid
    procedure,public::up=>upPolyFvGrid
    final::purgePolyFvGrid
  end type
  
contains
  
  !> clear this polyFvGrid
  elemental subroutine clearPolyFvGrid(this)
    class(polyFvGrid),intent(inout)::this !< this polyFvGrid
    
    call this%polyGrid%clear()
    if(allocated(this%iEP)) deallocate(this%iEP)
    if(allocated(this%neib)) deallocate(this%neib)
    if(allocated(this%near)) deallocate(this%near)
    if(allocated(this%aP)) deallocate(this%aP)
    if(allocated(this%normP)) deallocate(this%normP)
    if(allocated(this%pP)) deallocate(this%pP)
  end subroutine
  
  !> update this polyFvGrid
  elemental subroutine upPolyFvGrid(this)
    class(polyFvGrid),intent(inout)::this !< this polyFvGrid
    
    if(.not.this%isUp)then
      call this%polyGrid%up()
      call getNeibPolyFvGrid(this)
      call getNearPolyFvGrid(this)
    end if
  end subroutine
  
  !> get neighbor list and pairs
  elemental subroutine getNeibPolyFvGrid(grid)
    use modGeometry
    class(polyFvGrid),intent(inout)::grid !< the polyFvGrid
    integer,allocatable::iEN(:,:),nEN(:),nFace(:),iNF(:),iNF1(:),iEP(:,:)
    double precision,allocatable::aP(:),normP(:,:),pP(:,:)
    
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
    allocate(iNF(maxval([(maxval([(nNF(grid%sE(i),j),j=1,nFace(i))]),i=1,grid%nE)])))
    allocate(iNF1(size(iNF)))
    ! find neighbor through each face
    if(allocated(grid%neib)) deallocate(grid%neib)
    if(allocated(grid%iEP)) deallocate(grid%iEP)
    if(allocated(grid%aP)) deallocate(grid%aP)
    if(allocated(grid%normP)) deallocate(grid%normP)
    if(allocated(grid%pP)) deallocate(grid%pP)
    allocate(grid%neib(maxval(nFace),grid%nC))
    allocate(iEP(2,sum(nFace)))
    allocate(aP(sum(nFace)))
    allocate(normP(DIMS,sum(nFace)))
    allocate(pP(DIMS,sum(nFace)))
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
                  select case(count(iNF(:)>0))
                  case(3)
                    aP(grid%nP)=a3p(grid%pN(:,iNF(1:3)))
                    normP(:,grid%nP)=n3p(grid%pN(:,iNF(1:3)))
                    pP(:,grid%nP)=sum(grid%pN(:,iNF(1:3)),2)/3d0
                  case(4)
                    aP(grid%nP)=a4p(grid%pN(:,iNF(1:4)))
                    normP(:,grid%nP)=n4p(grid%pN(:,iNF(1:4)))
                    pP(:,grid%nP)=sum(grid%pN(:,iNF(1:4)),2)/4d0
                  case default
                    aP(grid%nP)=0d0
                    normP(:,grid%nP)=0d0
                    pP(:,grid%nP)=0d0
                  end select
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
                  select case(count(iNF(:)>0))
                  case(3)
                    aP(grid%nP)=a3p(grid%pN(:,iNF(1:3)))
                    normP(:,grid%nP)=n3p(grid%pN(:,iNF(1:3)))
                    pP(:,grid%nP)=sum(grid%pN(:,iNF(1:3)),2)/3d0
                  case(4)
                    aP(grid%nP)=a4p(grid%pN(:,iNF(1:4)))
                    normP(:,grid%nP)=n4p(grid%pN(:,iNF(1:4)))
                    pP(:,grid%nP)=sum(grid%pN(:,iNF(1:4)),2)/4d0
                  case default
                    aP(grid%nP)=0d0
                    normP(:,grid%nP)=0d0
                    pP(:,grid%nP)=0d0
                  end select
                  exit
                end if
              end do
            end if
          end do
        end if
      end do
    end do
    allocate(grid%iEP(2,grid%nP),source=iEP(:,1:grid%nP))!FIXME:remove work-around
    allocate(grid%aP(grid%nP),source=aP(1:grid%nP))!FIXME:remove work-around
    allocate(grid%normP(DIMS,grid%nP),source=normP(:,1:grid%nP))!FIXME:remove work-around
    allocate(grid%pP(DIMS,grid%nP),source=pP(:,1:grid%nP))!FIXME:remove work-around
    deallocate(iEN)
    deallocate(nEN)
    deallocate(nFace)
    deallocate(iNF)
    deallocate(iNF1)
    deallocate(iEP)
    deallocate(aP)
    deallocate(normP)
  end subroutine
  
  !> get the list of nearby elements for each cell
  elemental subroutine getNearPolyFvGrid(grid)
    use modGeometry
    class(polyFvGrid),intent(inout)::grid !< the polyFvGrid
    integer::nNear
    
    if(allocated(grid%near)) deallocate(grid%near)
    allocate(grid%near(MAX_N_NEAR,grid%nC))
    grid%near(:,:)=0
    do i=1,grid%nC
      nNear=0
      ! neighbors are nearby
      do j=1,nF(grid%sE(i))
        if(nNear<MAX_N_NEAR.and.grid%neib(j,i)>0)then
          nNear=nNear+1
          grid%near(nNear,i)=grid%neib(j,i)
        end if
      end do
      ! find more nearby elements if needed
      do while(nNear<MIN_N_NEAR)
        n=nNear
        do k=1,n
          l=grid%near(k,i)
          if(l<=grid%nC)then
            do j=1,nF(grid%sE(l))
              if(nNear<MAX_N_NEAR.and.grid%neib(j,l)>0&
              &  .and.all(grid%near(1:nNear,i)/=grid%neib(j,l))&
              &  .and.grid%neib(j,l)/=i)then
                nNear=nNear+1
                grid%near(nNear,i)=grid%neib(j,l)
              end if
            end do
          end if
        end do
      end do
    end do
  end subroutine
  
  !> destructor of polyGrid
  elemental subroutine purgePolyFvGrid(this)
    type(polyFvGrid),intent(inout)::this !< this polyFvGrid
    
    call this%clear()
  end subroutine
  
end module
